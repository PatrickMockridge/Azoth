// The scale family's layers, printed so each can be ported and checked.
//
// `CheckScalePotential` is a **saturation ratio per salt**, computed by walking NeqSim's own
// `compsalt` table: a solubility product `Ksp(T)` from five fitted coefficients, a pressure
// correction from the molar volume change, and an ion activity product from the aqueous
// phase's molalities and activity coefficients.
//
// **The ions are supplied directly, and that is the layer this probe settles.** NeqSim's own
// example drives a `SystemElectrolyteCPAstatoil` with `chemicalReactionInit()`, which
// speciates the carbonate system through a reactive equilibrium this library has none of. But
// `CheckScalePotential` reads the phase's components and nothing else, so a brine whose ions
// are *given* needs no reaction solve - which is what the example's commented-out lines do.
// The rows below are that: CO3-- and HCO3- stated, not solved for.
//
//   javac -proc:none -cp neqsim-f0c7436.jar ScaleProbe.java
//   java -cp .:neqsim-f0c7436.jar ScaleProbe > captures/scale_probe.tsv

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class ScaleProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    for (double temperatureC : new double[] {25.0, 50.0, 80.0}) {
      report(temperatureC, 10.0);
    }
  }

  private static void report(double temperatureC, double pressureBara) {
    System.out.printf("# brine at %.15g C, %.15g bara%n", temperatureC, pressureBara);
    try {
      SystemInterface fluid = new SystemSrkEos(273.15 + temperatureC, pressureBara);
      fluid.addComponent("water", 1.0);
      fluid.addComponent("Na+", 4e-5);
      fluid.addComponent("Cl-", 4e-5);
      fluid.addComponent("OH-", 220e-5);
      fluid.addComponent("Fe++", 110.1e-5);
      fluid.addComponent("CO3--", 0.95e-6);
      fluid.addComponent("HCO3-", 1.0e-4);
      fluid.addComponent("Ca++", 2.0e-3);
      fluid.setMixingRule(2);
      fluid.setMultiPhaseCheck(true);
      fluid.init(0);
      fluid.init(1);

      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      operations.checkScalePotential(1);

      row("phases", fluid.getNumberOfPhases());
      PhaseInterface aqueous = fluid.getPhase(1);
      System.out.printf("phase[1] = %s%n", aqueous.getType());
      for (int i = 0; i < aqueous.getNumberOfComponents(); i++) {
        String name = aqueous.getComponent(i).getName();
        row("x[" + name + "]", aqueous.getComponent(i).getx());
        // **The activity coefficients the saturation ratio is built from**, and they are the
        // phase model's own: a port that read a different gamma would produce a plausible
        // ratio for a different brine.
        row("gamma[" + name + "]", aqueous.getActivityCoefficient(i, 0));
      }

      // The operation's own table: salt name against relative solubility, as it reports it.
      String[][] table = operations.getResultTable();
      if (table == null) {
        System.out.println("note = no result table");
      } else {
        for (String[] entries : table) {
          System.out.printf("salt = %s%n", entries[0]);
          try {
            row("sr[" + entries[0] + "]", Double.parseDouble(entries[1]));
          } catch (NumberFormatException error) {
            System.out.printf("sr[%s] = %s%n", entries[0], entries[1]);
          }
        }
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
