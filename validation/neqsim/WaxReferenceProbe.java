package neqsim.thermo.component;

// The reference liquid `ComponentWax.fugcoef2` builds, read directly out of the component.
//
// **This one is in NeqSim's own package, and that is the whole point.** `fugcoef2` builds a
// one-component phase of the host's class, takes its fugacity coefficient and its molar
// volume, and multiplies by an exponential; the two halves are not separable from outside,
// because `ComponentSolid.refPhase` is package-private. Every captured coefficient so far has
// been the *product*, so a divergence in it could be either half. Read from here they are two
// numbers, and a port can be wrong about one of them and right about the other.
//
// The other half is pure arithmetic from numbers `WaxProbe` already prints - heat of fusion,
// triple point, molar mass - so with this the two are separated completely.
//
//   javac -proc:none -cp neqsim-3.20.0.jar -d . WaxReferenceProbe.java
//   java -cp .:neqsim-3.20.0.jar neqsim.thermo.component.WaxReferenceProbe

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.util.database.NeqSimDataBase;

public class WaxReferenceProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    for (double temperatureK : new double[] {285.0, 275.0, 261.0}) {
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
    System.out.printf("# wax reference at T = %.15g K, P = %.15g bara%n", temperatureK, pressureBara);
    try {
      SystemInterface fluid = build();
      fluid.setTemperature(temperatureK);
      fluid.setPressure(pressureBara);
      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      operations.TPflash();

      PhaseInterface wax = fluid.getPhaseOfType("wax");
      row("wax_beta", wax.getBeta());
      for (int i = 0; i < wax.getNumberOfComponents(); i++) {
        ComponentInterface component = wax.getComponent(i);
        if (!component.isWaxFormer()) {
          continue;
        }
        String name = component.getName();
        // **The two halves.** `refPhase` is the one-component liquid the coefficient is built
        // from; what `fugcoef2` reports is this fugacity coefficient times its exponential.
        PhaseInterface reference = ((ComponentSolid) component).refPhase;
        if (reference == null) {
          System.out.printf("note = %s carries no reference phase%n", name);
          continue;
        }
        row("ref_phi[" + name + "]", reference.getComponent(0).getFugacityCoefficient());
        row("ref_molar_volume[" + name + "]", reference.getMolarVolume());
        row("ref_z[" + name + "]", reference.getZ());
        row("ref_pressure[" + name + "]", reference.getPressure());
        row("ref_temperature[" + name + "]", reference.getTemperature());
        // And the reference's own component, because a port that read the *parent's* critical
        // constants would be 21-40 per cent out on the lightest cut.
        row("ref_tc[" + name + "]", reference.getComponent(0).getTC());
        row("ref_pc[" + name + "]", reference.getComponent(0).getPC());
        row("ref_acentric[" + name + "]", reference.getComponent(0).getAcentricFactor());
        row("ref_m[" + name + "]", reference.getComponent(0).getAttractiveTerm().getm());
        // The phase's own coefficient, so the two rows can be checked against each other.
        row("wax_phi[" + name + "]", component.getFugacityCoefficient());
        // `SolidFug` and `x` are the two factors the coefficient is built from, and they are
        // what make the exponential recoverable exactly: `SolidFug / (x * f_liq)` is
        // `exp(exponent)` with nothing else in it.
        row("wax_solid_fug[" + name + "]", ((ComponentSolid) component).SolidFug);
        row("wax_x[" + name + "]", component.getx());
        row("wax_heat_of_fusion[" + name + "]", component.getHeatOfFusion());
        row("wax_triple_point[" + name + "]", component.getTriplePointTemperature());
        row("wax_molar_mass_g[" + name + "]", component.getMolarMass() * 1000.0);
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
