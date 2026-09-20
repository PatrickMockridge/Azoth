// What a built `SystemSoreideWhitson` computes, before any of it is ported.
//
//     javac -proc:none -cp neqsim-3.20.0.jar SoreideWhitsonArithmetic.java
//     java -cp .:neqsim-3.20.0.jar SoreideWhitsonArithmetic
//
// The second of F's two instruments. `SoreideWhitsonProbe` answers what happens to the six
// interaction rows written with a comma; this prints the model those rows feed - the PR
// phase with a salinity-modified water alpha and the Whitson-Soreide kij matrix - so the
// port has an oracle at a state rather than at a table.
//
// Three things the source does not say and this prints:
//
//   1. **The alpha function is water-only and everything else is PR78.** The tranche already
//      ports the water form as `eos.soreide_whitson_alpha`, so the phase's content is which
//      component takes it, and the answer is the component's *name*.
//   2. **The aqueous and non-aqueous paths are different kij sources.** A water-gas pair
//      goes through the parameterisation - legacy, Chabab 2019 or Burgoyne-Nielsen 2026 -
//      and everything else through the interaction table's `KIJWhitsonSoriede`.
//   3. **What the salinity actually is.** `PhaseSoreideWhitson.getSalinityConcentration` is
//      the number the alpha and the aqueous correlations both read, and it is a *molality*
//      over the water in the phase, so it moves with the flash.

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSoreideWhitson;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class SoreideWhitsonArithmetic {

  private static SystemInterface build(String[] names, double[] z, double salinity, double tC,
      double pBara) {
    SystemSoreideWhitson system = new SystemSoreideWhitson(298.0, 20.0);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i], "mole/sec");
    }
    system.addSalinity(salinity, "mole/sec");
    system.setTotalFlowRate(15, "mole/sec");
    system.setMixingRule(11);
    system.setTemperature(tC, "C");
    system.setPressure(pBara, "bara");
    return system;
  }

  private static void alphaOf(SystemInterface system, String label) {
    PhaseInterface phase = system.getPhase(1);
    System.out.printf("%n  alpha against temperature, %s:%n", label);
    System.out.printf("    %10s", "T / K");
    for (int i = 0; i < phase.getNumberOfComponents(); i++) {
      System.out.printf(" %14s", phase.getComponent(i).getComponentName());
    }
    System.out.println();
    for (double t : new double[] {273.15, 298.15, 323.15, 373.15, 473.15}) {
      System.out.printf("    %10.2f", t);
      for (int i = 0; i < phase.getNumberOfComponents(); i++) {
        System.out.printf(" %14.10g",
            phase.getComponent(i).getAttractiveTerm().alpha(t));
      }
      System.out.println();
    }
  }

  private static void flash(SystemInterface system, String label) {
    System.out.printf("%n=== %s ===%n", label);
    try {
      // `calcSalinity` turns the salt moles into the molality both the alpha and the
      // aqueous correlations read, and nothing calls it for us.
      ((SystemSoreideWhitson) system).calcSalinity();
      System.out.printf("  salt input = %.12g mol/sec, phase concentration = %.12g%n",
          ((SystemSoreideWhitson) system).getSalinity(),
          ((neqsim.thermo.phase.PhaseSoreideWhitson) system.getPhase(1))
              .getSalinityConcentration());
      ThermodynamicOperations ops = new ThermodynamicOperations(system);
      ops.TPflash();
      system.initProperties();
      System.out.printf("  after the flash: salt = %.12g, phase concentration = %.12g%n",
          ((SystemSoreideWhitson) system).getSalinity(),
          ((neqsim.thermo.phase.PhaseSoreideWhitson) system.getPhase(1))
              .getSalinityConcentration());
      System.out.printf("  phases = %d   T = %.2f K   P = %.4f bar%n",
          system.getNumberOfPhases(), system.getTemperature(), system.getPressure());
      for (int p = 0; p < system.getNumberOfPhases(); p++) {
        PhaseInterface phase = system.getPhase(p);
        System.out.printf("  [%d] %-12s beta = %.12g  z = %.15g  salinity = %.12g%n", p,
            phase.getClass().getSimpleName(), phase.getBeta(), phase.getZ(),
            phase instanceof neqsim.thermo.phase.PhaseSoreideWhitson
                ? ((neqsim.thermo.phase.PhaseSoreideWhitson) phase).getSalinityConcentration()
                : Double.NaN);
        for (int i = 0; i < phase.getNumberOfComponents(); i++) {
          // **`ln phi` and not a gamma**: this is an equation of state, and
          // `Phase.getActivityCoefficient` is not its surface - it re-runs `PhaseEos.init`
          // on a phase whose mixing-rule arrays are not built for it and throws.
          System.out.printf("        %-10s x = %.12g  ln phi = %.12g%n",
              phase.getComponent(i).getComponentName(), phase.getComponent(i).getx(),
              Math.log(phase.getComponent(i).getFugacityCoefficient()));
        }
      }
    } catch (Throwable error) {
      System.out.println("  " + error.getClass().getSimpleName() + ": " + error.getMessage());
      for (StackTraceElement frame : error.getStackTrace()) {
        if (frame.getClassName().startsWith("neqsim.")) {
          System.out.println("        at " + frame);
        }
      }
    }
  }

  public static void main(String[] args) {
    // The shipped test's own mixture, at its own state, and its own expected compositions.
    String[] shipped = {"nitrogen", "CO2", "methane", "ethane", "water"};
    double[] moles = {0.1, 0.2, 0.3, 0.3, 0.1};
    flash(build(shipped, moles, 0.0, 45.0, 40.0), "the shipped test: N2 CO2 C1 C2 water, 45 C, 40 bara");

    // With salinity, which is what the model is for.
    flash(build(shipped, moles, 0.5, 45.0, 40.0),
        "the same brine at 0.5 mol/sec salinity");

    // A mixture whose pairs are the interaction table's rather than water's.
    String[] condensate = {"CO2", "methane", "n-butane", "n-heptane", "water"};
    double[] condensateMoles = {0.15, 0.4, 0.15, 0.1, 0.2};
    flash(build(condensate, condensateMoles, 0.0, 60.0, 100.0),
        "a condensate: CO2 C1 nC4 nC7 water, 60 C, 100 bara");

    SystemInterface system = build(shipped, moles, 0.0, 45.0, 40.0);
    system.init(0);
    system.init(1);
    alphaOf(system, "the shipped mixture");
  }
}
