// The wax family's layers, printed so each can be ported and checked.
//
// **The fluid is the first layer, and it is the one a port cannot build.** A wax system is
// made from a TBP fraction and a plus fraction: `addPlusFraction` plus
// `characterisePlusFraction()` splits the plus into pseudo-components with a Pedersen plus
// model, and `getWaxModel().addTBPWax()` then splits *those* again into parallel `wax<name>`
// components and sets each one's `waxFormer`, `heatOfFusion` and `triplePointTemperature`.
// So the components below - their molar masses, normal liquid densities, critical constants
// and the two wax temperatures - are what a port has to be *given*, because the subsystem
// that invents them is 16,000 lines of plus-fraction characterisation.
//
// The second layer is the wax phase: `ComponentWax.fugcoef2` builds the solid's fugacity from
// a **reference liquid phase** of the same component, the heat of fusion, the triple point and
// a heat-capacity difference, and the phase's composition follows. Both are printed.
//
//   javac -proc:none -cp neqsim-3.20.0.jar WaxProbe.java
//   java -cp .:neqsim-3.20.0.jar WaxProbe > captures/wax_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.util.database.NeqSimDataBase;

public class WaxProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    // NeqSim's own wax test fluid at 261 K and 5 bara, where it reports a wax phase.
    double[] temperatures = {285.0, 275.0, 261.0};
    for (double temperatureK : temperatures) {
      report(temperatureK, 5.0);
    }
  }

  private static SystemInterface build() {
    NeqSimDataBase.setCreateTemporaryTables(true);
    SystemInterface system = new SystemSrkEos(298.0, 10.0);
    system.addComponent("methane", 6.78);
    system.addTBPfraction("C19", 10.13, 170.0 / 1000.0, 0.7814);
    system.addPlusFraction("C20", 10.62, 381.0 / 1000.0, 0.850871882888);
    system.getCharacterization().characterisePlusFraction();
    system.getWaxModel().addTBPWax();
    system.createDatabase(true);
    system.setMixingRule(2);
    // **The flag has to be on again here, and the turn-off at the end of this method is
    // what made the instrument irreproducible for an hour.** `addSolidComplexPhase` builds
    // a hydrate phase and then a wax phase from the system's component *names*, and every
    // one of those is resolved against the database - the wax pseudo-components exist only
    // in the temporary table, so a run that reaches this line with the flag off throws
    // "waxPC1_PC not found in database" and the capture comes out with none of the wax rows
    // in it. NeqSim's own test never notices because it sets the flag once and never clears
    // it; this probe clears it, so it has to put it back.
    NeqSimDataBase.setCreateTemporaryTables(true);
    system.addSolidComplexPhase("wax");
    system.setMultiphaseWaxCheck(true);
    system.setMultiPhaseCheck(true);
    NeqSimDataBase.setCreateTemporaryTables(false);
    system.init(0);
    system.init(1);
    return system;
  }

  private static void report(double temperatureK, double pressureBara) {
    System.out.printf("# wax fluid at T = %.15g K, P = %.15g bara%n", temperatureK, pressureBara);
    try {
      SystemInterface fluid = build();
      // **The fluid a port has to be given.** Printed before the flash so it is the same at
      // every state, and once in enough detail to rebuild: the characterisation is what azoth
      // does not have, and these rows are its output.
      for (int i = 0; i < fluid.getPhase(0).getNumberOfComponents(); i++) {
        ComponentInterface component = fluid.getPhase(0).getComponent(i);
        String name = component.getName();
        row("component[" + i + "].mw", component.getMolarMass());
        row("component[" + i + "].normal_liquid_density", component.getNormalLiquidDensity());
        row("component[" + i + "].tc", component.getTC());
        row("component[" + i + "].pc", component.getPC());
        row("component[" + i + "].acentric", component.getAcentricFactor());
        row("component[" + i + "].moles", component.getNumberOfmoles());
        System.out.printf("component[%d] = %s%n", i, name);
        row("component[" + i + "].wax_former", component.isWaxFormer() ? 1.0 : 0.0);
        // **The Peneloux shift, which a wax calculation turns out to need.** `addTBPfraction`
        // fits `racketZ` from a reference flash and `ComponentSrk.getVolumeCorrection` turns
        // it into the volume translation the cubic and the fugacity coefficient carry; a port
        // without it is a per cent out, and the per cent grows with the cut's molar mass.
        row("component[" + i + "].racket_z", component.getRacketZ());
        row("component[" + i + "].volume_correction", component.getVolumeCorrection());
        row("component[" + i + "].heat_of_fusion", component.getHeatOfFusion());
        row("component[" + i + "].triple_point_temperature", component.getTriplePointTemperature());
        if (component.getAttractiveTerm() != null) {
          row("component[" + i + "].attractive_m", component.getAttractiveTerm().getm());
        }
      }

      fluid.setTemperature(temperatureK);
      fluid.setPressure(pressureBara);
      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      operations.TPflash();

      row("phases", fluid.getNumberOfPhases());
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        System.out.printf("phase[%d] = %s%n", p, phase.getType());
        row("phase[" + p + "].beta", phase.getBeta());
        for (int i = 0; i < phase.getNumberOfComponents(); i++) {
          String name = phase.getComponent(i).getName();
          row("phase[" + p + "].x[" + name + "]", phase.getComponent(i).getx());
          // **The wax model's own coefficient, and only in the wax phase.** `ComponentWax`
          // is what the wax phase's components are; the same component in the gas is an
          // ordinary cubic one, and its coefficient there is the cubic's. Printing both for
          // one name would make the model look like it agreed with itself.
          // **Every phase's coefficient for every component**, because the question the
          // capture has to answer is whether the wax is in equilibrium with the phases
          // beside it: `x_i phi_i` must be one number across them, and a model that solved
          // only its own equality would look the same in its own rows.
          row(
              "fugcoef[" + p + "][" + name + "]",
              phase.getComponent(i).getFugacityCoefficient());
        }
      }

      // **A material balance, because a wax fraction is one.** Counting each component
      // across the phases against the feed, as `HydrateFractionProbe` does for its own.
      int n = fluid.getPhase(0).getNumberOfComponents();
      for (int i = 0; i < n; i++) {
        double inPhases = 0.0;
        for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
          inPhases += fluid.getPhase(p).getBeta() * fluid.getPhase(p).getComponent(i).getx();
        }
        String name = fluid.getPhase(0).getComponent(i).getName();
        row("balance_error[" + name + "]", inPhases - fluid.getPhase(0).getComponent(i).getz());
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
