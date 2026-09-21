// What a Furst electrolyte phase actually uses for a TEG pair, against what its own fit says.
//
//     javac -proc:none -cp neqsim-f0c7436.jar FurstTegBranch.java
//     java -cp .:neqsim-f0c7436.jar FurstTegBranch
//
// `EosMixingRuleHandler.ElectrolyteMixRule.calcWij` builds the short-range `Wij` table from
// a chain of `if/else if` tests on the *other* component's name. `TEG` is tested twice:
//
//     } else if (solventName.equals("TEG") || solventName.equals("triethylene glycol")) {
//       wij[0][i][j] = FurstElectrolyteConstants.getFurstParamMEG(2) * stokesDiam
//           + FurstElectrolyteConstants.getFurstParamMEG(3);      // <- MEG, not TEG
//     ...
//     } else if (solventName.equals("TEG")) {
//       wij[0][i][j] = FurstElectrolyteConstants.getFurstParamTEG(2) * stokesDiam
//           + FurstElectrolyteConstants.getFurstParamTEG(3);      // <- never reached
//
// The first match wins, so a `TEG` pair is assigned **the MEG parameter set** and the branch
// written for TEG below it cannot run. The two sets differ where it matters: `[2]` is
// `8.0e-5` against `4.98e-5` and `[3]` is `-1.15e-4` against `-1.22e-4`, and those are the
// slope and intercept of `Wij = p2 d + p3`.
//
// This prints, for a phase carrying both solvents in turn, the `Wij(Na+, solvent)` the model
// holds against the one each solvent's own parameters give - so the size of the substitution
// is a measured number rather than an argument about which branch is reachable.

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseModifiedFurstElectrolyteEos;
import neqsim.thermo.system.SystemFurstElectrolyteEos;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class FurstTegBranch {

  /**
   * `FurstElectrolyteConstants.getFurstParamX(i)`, by the set's name.
   *
   * `getWijParameter` is what the phase holds; this is what each solvent's own fit says it
   * should hold. The two are only comparable through `Wij = p2 d + p3`, which is the form
   * the branch uses.
   */
  private static double param(String set, int i) {
    try {
      java.lang.reflect.Method m = Class
          .forName("neqsim.thermo.util.constants.FurstElectrolyteConstants")
          .getMethod("getFurstParam" + set, int.class);
      return (double) m.invoke(null, i);
    } catch (Throwable e) {
      return Double.NaN;
    }
  }

  /**
   * The `Wij` a phase holds for a pair, from the rule's own public accessor.
   *
   * `ElectrolyteMixRule` is an *inner* class of the handler, so its `wij` array belongs to
   * the outer instance and a `getDeclaredField` on the rule finds nothing - which is how the
   * first attempt at this probe returned `NaN` for every pair.
   */
  private static double wij(PhaseInterface phase, int i, int j) {
    try {
      return ((PhaseModifiedFurstElectrolyteEos) phase).getElectrolyteMixingRule()
          .getWijParameter(i, j);
    } catch (Throwable e) {
      return Double.NaN;
    }
  }

  private static void report(String solvent) {
    SystemInterface system = new SystemFurstElectrolyteEos(298.15, 10.01325);
    system.addComponent("methane", 0.1);
    system.addComponent("water", 1.0);
    system.addComponent(solvent, 0.5);
    system.addComponent("Na+", 0.001);
    system.addComponent("Cl-", 0.001);
    system.setMixingRule(4);
    try {
      ThermodynamicOperations ops = new ThermodynamicOperations(system);
      ops.TPflash();
      system.initProperties();
      PhaseInterface phase = system.getPhase(1);
      int cation = phase.getComponent("Na+").getComponentNumber();
      int index = phase.getComponent(solvent).getComponentNumber();
      double d = phase.getComponent(cation).getStokesCationicDiameter();

      double held = wij(phase, cation, index);
      double own = param(solvent, 2) * d + param(solvent, 3);
      double meg = param("MEG", 2) * d + param("MEG", 3);
      System.out.printf("%n=== %s ===%n", solvent);
      System.out.printf("  stokes diameter of Na+       = %.12g%n", d);
      System.out.printf("  Wij(Na+, %-4s) the phase holds = %.12g%n", solvent, held);
      System.out.printf("  from the %s parameters        = %.12g%n", solvent, own);
      System.out.printf("  from the MEG parameters        = %.12g%n", meg);
      System.out.printf("  ratio held / own               = %.6f%n",
          (own != 0.0) ? held / own : Double.NaN);
    } catch (Throwable error) {
      System.out.printf("%n=== %s failed ===%n", solvent);
      System.out.println("  " + error.getClass().getSimpleName() + ": " + error.getMessage());
    }
  }

  public static void main(String[] args) {
    report("TEG");
    report("MEG");
  }
}
