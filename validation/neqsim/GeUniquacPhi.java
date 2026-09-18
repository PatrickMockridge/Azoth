// What NeqSim 3.20.0's UNIQUAC liquid reports - and, since #3774, it refuses to exist.
//
//     javac -proc:none -cp neqsim-3.20.0.jar GeUniquacPhi.java
//     java -cp .:neqsim-3.20.0.jar GeUniquacPhi
//
// `eos.ge_uniquac_phase` is `gamma_i * P0_i / P`, so an oracle needs NeqSim to produce a
// UNIQUAC `gamma_i`. It cannot, and it no longer pretends to.
//
// **Before commit `c5ec5fb` (PR #3774, closing upstream #3770)** a bare `PhaseGEUniquac`
// constructed, `getExcessGibbsEnergy` called the eight-argument `getGamma` which was
// `return 0.0`, and `fugcoef` read the inherited `gamma` field that nothing ever wrote.
// The driver printed those zeros.
//
// **After it**, both standalone constructors throw
// `UnsupportedOperationException` - "activity-coefficient implementation and parameter
// data are incomplete" - and `ComponentGEUniquac.getGamma` throws where it returned 0.0.
// The first section below now records that, rather than a coefficient: the constructor is
// what stops a bare UNIQUAC phase, so the second throw is not reachable from here.
//
// `ComponentGEUniquac`'s constructor reads `rUNIQUAQ`/`qUNIQUAQ` from the `unifaccomp`
// table, which is 0.0 for 109 of its 112 rows - one of the incompletenesses the refusal
// cites. The group sums `eos.uniquac_activity_coefficients` uses are the only usable r and
// q, and the second section prints those, from `ComponentGEUnifac.getR`/`getQ`, which is
// the half of this phase that *does* have an oracle.
//
// The third section is the contrast that explains why UNIFAC works and UNIQUAC does not:
// `PhaseGEUnifac extends PhaseGEUniquac` and overrides `getExcessGibbsEnergy` to call the
// *five*-argument `getGamma`, which virtual dispatch resolves to `ComponentGEUnifac`'s
// real implementation.

import neqsim.thermo.component.ComponentGEUnifac;
import neqsim.thermo.component.ComponentGEUniquac;
import neqsim.thermo.phase.PhaseGEInterface;
import neqsim.thermo.phase.PhaseGEUnifac;
import neqsim.thermo.phase.PhaseGEUniquac;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemUNIFAC;

public class GeUniquacPhi {

  /** A bare UNIQUAC phase, which upstream now refuses to build. */
  static void bare(String[] names, double[] x, double temperatureK, double pressureBar) {
    System.out.println("bare PhaseGEUniquac  T=" + temperatureK + " P=" + pressureBar);
    try {
      PhaseGEUniquac phase = new PhaseGEUniquac();
      phase.setTemperature(temperatureK);
      phase.setPressure(pressureBar);
      for (int i = 0; i < names.length; i++) {
        phase.addComponent(names[i], x[i], x[i], i);
      }
      System.out.println("  constructed - upstream no longer refuses this");
      System.out.println("  getExcessGibbsEnergy = "
          + ((PhaseGEInterface) phase).getExcessGibbsEnergy(phase, names.length, temperatureK,
              pressureBar, PhaseType.LIQUID));
    } catch (UnsupportedOperationException refused) {
      System.out.println("  PhaseGEUniquac()          -> " + refused.getClass().getName());
      System.out.println("    " + refused.getMessage());
    }
    try {
      new ComponentGEUniquac(names[0], x[0], x[0], 0);
      System.out.println("  constructed a ComponentGEUniquac - upstream no longer refuses this");
    } catch (UnsupportedOperationException refused) {
      System.out.println("  new ComponentGEUniquac()  -> " + refused.getClass().getName());
      System.out.println("    " + refused.getMessage());
    }
  }

  /** The subclass that works, for contrast, and the r/q oracle. */
  static void unifac(String[] names, double[] x, double temperatureK, double pressureBar) {
    SystemUNIFAC system = new SystemUNIFAC(temperatureK, pressureBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], x[i]);
    }
    // The array is synchronised before *every* read of it - the first `init(0)` included,
    // which is the order `GeUnifacPhi.java` measured.
    system.createDatabase(true);
    synchronise(system);
    system.init(0);
    synchronise(system);
    system.setMixingRule("classic");
    synchronise(system);
    system.init(0);
    synchronise(system);

    PhaseInterface ge = null;
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      if (system.getPhase(p) instanceof PhaseGEInterface) {
        ge = system.getPhase(p);
      }
    }
    if (ge == null) {
      System.out.println("no GE phase");
      return;
    }
    ((PhaseGEInterface) ge).getExcessGibbsEnergy(ge, ge.getNumberOfComponents(), temperatureK,
        pressureBar, PhaseType.LIQUID);
    System.out.println("SystemUNIFAC (" + ge.getClass().getSuperclass().getSimpleName()
        + " subclass)  T=" + temperatureK + " P=" + pressureBar);
    for (int i = 0; i < ge.getNumberOfComponents(); i++) {
      ge.getComponent(i).fugcoef(ge);
      ComponentGEUnifac component = (ComponentGEUnifac) ge.getComponent(i);
      System.out.println("  " + component.getName()
          + "  x=" + component.getx()
          + "  getR()=" + component.getR()
          + "  getQ()=" + component.getQ()
          + "  gamma=" + component.getGamma()
          + "  phi=" + component.getFugacityCoefficient());
    }
  }

  /** `ComponentGEUnifac`'s group array is empty until it is synchronised from the list. */
  static void synchronise(neqsim.thermo.system.SystemInterface system) {
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      if (!(system.getPhase(p) instanceof PhaseGEUnifac)) {
        continue;
      }
      PhaseInterface ge = system.getPhase(p);
      for (int i = 0; i < ge.getNumberOfComponents(); i++) {
        neqsim.thermo.component.ComponentGEUnifac component =
            (neqsim.thermo.component.ComponentGEUnifac) ge.getComponent(i);
        component.setUnifacGroups(component.getUnifacGroups2());
      }
    }
  }

  public static void main(String[] args) {
    // The withdrawal, at the two states `eos.ge_uniquac_phase`'s cases use.
    bare(new String[] {"methanol", "water"}, new double[] {0.5, 0.5}, 298.15, 1.0);
    bare(new String[] {"water", "nc10"}, new double[] {0.5, 0.5}, 350.0, 1.0);
    // The r/q oracle: `getR()`/`getQ()` are the UNIFAC group sums, which is what
    // `eos.uniquac_activity_coefficients` resolves `r` and `q` to. Printed for the
    // substances the phase's cases use, so the pair can be compared independently of the
    // composition that has no oracle.
    System.out.println("--- r and q, from ComponentGEUnifac.getR/getQ ---");
    unifac(new String[] {"methanol", "water"}, new double[] {0.5, 0.5}, 298.15, 1.0);
    unifac(new String[] {"water", "nc10"}, new double[] {0.5, 0.5}, 350.0, 1.0);
    unifac(new String[] {"benzene", "n-hexane"}, new double[] {0.4, 0.6}, 320.0, 1.0);
  }
}
