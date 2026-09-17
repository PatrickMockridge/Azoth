// What NeqSim 3.20.0's UNIQUAC liquid reports - and what it reports is nothing.
//
//     javac -proc:none -cp neqsim-3.20.0.jar GeUniquacPhi.java
//     java -cp .:neqsim-3.20.0.jar GeUniquacPhi
//
// `eos.ge_uniquac_phase` is `gamma_i * P0_i / P`, so an oracle needs NeqSim to produce a
// UNIQUAC `gamma_i`. It cannot. This drives the two paths that look like they should and
// prints what each gives:
//
//   1. A bare `PhaseGEUniquac` - the phase itself, with `ComponentGEUniquac` components.
//      `getExcessGibbsEnergy` calls the eight-argument `getGamma`, which is `return 0.0`.
//   2. `SystemUNIFAC`, whose `PhaseGEUnifac extends PhaseGEUniquac` - to show the contrast:
//      `PhaseGEUnifac` overrides `getExcessGibbsEnergy` and calls the *five*-argument
//      `getGamma`, which virtual dispatch resolves to `ComponentGEUnifac`'s real
//      implementation. That is why UNIFAC works and UNIQUAC does not.
//
// `fugcoef` is the five-argument path, which returns the inherited `gamma` field that
// nothing in `ComponentGEUniquac` ever writes - so it is 0.0 with or without the phase.
//
// Separately, `ComponentGEUniquac`'s constructor reads `rUNIQUAQ`/`qUNIQUAQ` from the
// `unifaccomp` table, which is the pair of columns measured to be 0.0 for 109 of its 112
// rows. So even a working `getGamma` would divide by zero for almost every mixture; the
// group sums `eos.uniquac_activity_coefficients` uses are the only usable r and q, and
// this prints those too, from `ComponentGEUnifac.getR`/`getQ`.

import neqsim.thermo.component.ComponentGEUniquac;
import neqsim.thermo.component.ComponentGEUnifac;
import neqsim.thermo.phase.PhaseGEInterface;
import neqsim.thermo.phase.PhaseGEUnifac;
import neqsim.thermo.phase.PhaseGEUniquac;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemUNIFAC;

public class GeUniquacPhi {

  /** A bare UNIQUAC phase, at a stated composition. */
  static void bare(String[] names, double[] x, double temperatureK, double pressureBar) {
    PhaseGEUniquac phase = new PhaseGEUniquac();
    phase.setTemperature(temperatureK);
    phase.setPressure(pressureBar);
    for (int i = 0; i < names.length; i++) {
      phase.addComponent(names[i], x[i], x[i], i);
    }
    for (int i = 0; i < names.length; i++) {
      phase.getComponent(i).setx(x[i]);
    }

    System.out.println("bare PhaseGEUniquac  T=" + temperatureK + " P=" + pressureBar);
    System.out.println("  getExcessGibbsEnergy = "
        + ((PhaseGEInterface) phase).getExcessGibbsEnergy(phase, names.length, temperatureK,
            pressureBar, PhaseType.LIQUID));
    for (int i = 0; i < names.length; i++) {
      ComponentGEUniquac component = (ComponentGEUniquac) phase.getComponent(i);
      double gamma8 = component.getGamma(phase, names.length, temperatureK, pressureBar,
          PhaseType.LIQUID, null, null, null, null);
      component.fugcoef(phase);
      System.out.println("  " + component.getName()
          + "  x=" + component.getx()
          + "  rUNIQUAQ=" + component.getr()
          + "  qUNIQUAQ=" + component.getq()
          + "  getGamma(8-arg)=" + gamma8
          + "  getGamma()=" + component.getGamma()
          + "  getFugacityCoefficient=" + component.getFugacityCoefficient());
    }
  }

  /** The subclass that works, for contrast. */
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
    bare(new String[] {"methanol", "water"}, new double[] {0.5, 0.5}, 298.15, 1.0);
    bare(new String[] {"water", "nc10"}, new double[] {0.5, 0.5}, 350.0, 1.0);
    unifac(new String[] {"methanol", "water"}, new double[] {0.5, 0.5}, 298.15, 1.0);
    // The r/q oracle: `getR()`/`getQ()` are the UNIFAC group sums, which is what
    // `eos.uniquac_activity_coefficients` resolves `r` and `q` to. Printed for the
    // substances the phase's cases use, so the pair can be compared independently of
    // the composition that has no oracle.
    System.out.println("--- r and q, from ComponentGEUnifac.getR/getQ ---");
    unifac(new String[] {"water", "nc10"}, new double[] {0.5, 0.5}, 350.0, 1.0);
    unifac(new String[] {"benzene", "n-hexane"}, new double[] {0.4, 0.6}, 320.0, 1.0);
  }
}
