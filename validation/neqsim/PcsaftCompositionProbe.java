package neqsim.thermo.component;

// NeqSim's composition derivatives against its own `F`.
//
// `ln phi_i = d(nF)/dn_i - ln Z`, and `ComponentPCSAFT.dFdN` is the analytic first half.
// A central difference of `getF()` at fixed temperature and volume is that same number
// computed with none of that code, so the two rows together say whether a disagreement is
// in the derivative or in the state the derivative was taken at. `getF()` is *extensive*
// here - `F_HC_SAFT` and both dispersion terms carry `getNumberOfMolesInPhase()` - which is
// the one thing a reader of `dFdN` has to know to difference it correctly.
//
// The second block forces the phase to a volume given here rather than the one the EOS
// chose, which is how a port that converged a slightly different root is compared without
// the volume difference standing in for a derivative difference.
//
// **`initSAFTDerivatives` is package-private**, which is why this probe declares NeqSim's
// package; `WaxReferenceProbe` is here for the same reason.
//
// `v` is m^3/mol, converted from the phase's own `molarVolume` field by the `1.0e-5` the
// PC-SAFT phase itself uses. `Z_from_v` is `p v / R T` with `p` in bara and NeqSim's own
// `R = 8.3144621`, and it reproduces the `Z` field in the first block exactly - which is
// what licenses it in the second, where `Z` is a stale field.
//
//   javac -proc:none -cp neqsim-f0c7436.jar -d . PcsaftCompositionProbe.java
//   java -cp .:neqsim-f0c7436.jar neqsim.thermo.component.PcsaftCompositionProbe

import neqsim.thermo.phase.PhasePCSAFT;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPCSAFT;

public class PcsaftCompositionProbe {

  /** The molar volume azoth's own root holds at this state, in m^3/mol. */
  private static final double AZOTH_V = 8.141097344315069e-4;

  /** NeqSim's `ThermodynamicConstantsInterface.R`. */
  private static final double R = 8.3144621;

  private static final double[] STEPS = {1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6};

  private static double t;
  private static double p;

  public static void main(String[] args) {
    System.out.println("# azoth PcsaftCompositionProbe - NeqSim master's analytic `dFdN` and");
    System.out.println("# `dFdNdN` against central differences of its own `getF()` and `dFdN` at");
    System.out.println("# fixed T and V. `fd_*` rows use (plus - minus) / (2h), with `h` in moles.");
    System.out.println("# `getF()` is already extensive, so `fd_dFdN[i]` meets `dFdN[i]` and");
    System.out.println("# `fd_dFdNdN[i][j]` meets `dFdNdN[i][j]`, with no scale factor either way.");
    t = 350.0;
    p = 30.0;
    SystemInterface system = new SystemPCSAFT(t, p);
    system.addComponent("methane", 0.6);
    system.addComponent("n-butane", 0.4);
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);
    system.init(2);

    PhasePCSAFT phase = (PhasePCSAFT) system.getPhase(0);
    System.out.printf("%n# methane/n-butane at T=%g K, P=%g bara, as the EOS chose it%n", t, p);
    report(phase);
    System.out.printf("%n# the same phase forced to azoth's molar volume %.15g%n", AZOTH_V);
    PhasePCSAFT forced = (PhasePCSAFT) phase.clone();
    forced.setMolarVolume(AZOTH_V * 1.0e5);
    forced.volInit();
    reinit(forced);
    report(forced);
  }

  private static void report(PhasePCSAFT phase) {
    int n = phase.getNumberOfComponents();
    row("n", phase.getNumberOfMolesInPhase());
    row("v", phase.getMolarVolume() * 1.0e-5 / phase.getNumberOfMolesInPhase());
    row("Z", phase.getZ());
    row("Z_from_v", p * phase.getMolarVolume() / (R * t));
    row("eta", phase.getNSAFT());
    row("F", phase.getF());
    row("F_hc", phase.F_HC_SAFT());
    row("F_disp1", phase.F_DISP1_SAFT());
    row("F_disp2", phase.F_DISP2_SAFT());
    for (int i = 0; i < n; i++) {
      ComponentPCSAFT c = (ComponentPCSAFT) phase.getComponent(i);
      row("dFdN[" + i + "]", c.dFdN(phase, n, t, p));
      row("dFdN_hc[" + i + "]", c.dF_HC_SAFTdN(phase, n, t, p));
      row("dFdN_disp1[" + i + "]", c.dF_DISP1_SAFTdN(phase, n, t, p));
      row("dFdN_disp2[" + i + "]", c.dF_DISP2_SAFTdN(phase, n, t, p));
      row("lnPhi[" + i + "]", Math.log(c.getFugacityCoefficient()));
      for (int j = 0; j < n; j++) {
        row("dFdNdN[" + i + "][" + j + "]", c.dFdNdN(j, phase, n, t, p));
      }
    }
    for (double h : STEPS) {
      for (int i = 0; i < n; i++) {
        row("fd_dFdN[" + i + "] h=" + h, (extensive(phase, i, h) - extensive(phase, i, -h)) / (2.0 * h));
      }
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          row("fd_dFdNdN[" + i + "][" + j + "] h=" + h,
              (dFdN(perturbed(phase, j, h), i) - dFdN(perturbed(phase, j, -h), i)) / (2.0 * h));
        }
      }
    }
  }

  /** `getF()` on a clone whose `k`th mole number moved by `h` at the original's volume. */
  private static double extensive(PhasePCSAFT phase, int k, double h) {
    return perturbed(phase, k, h).getF();
  }

  private static double dFdN(PhasePCSAFT phase, int i) {
    return ((ComponentPCSAFT) phase.getComponent(i)).dFdN(phase, phase.getNumberOfComponents(), t, p);
  }

  /**
   * A clone whose `k`th mole number moved by `h` at the original's total volume.
   *
   * **The mole fractions are written directly, and that is not optional.** `Component`'s
   * `setNumberOfMolesInPhase` takes the phase total and stores `total * x_i`, not a
   * component's own moles, and `x` is otherwise a set-once field - so a clone perturbed
   * through that setter computes `F` at the *old* composition, and the difference that
   * comes out of it is not the composition derivative at all.
   */
  private static PhasePCSAFT perturbed(PhasePCSAFT phase, int k, double h) {
    int n = phase.getNumberOfComponents();
    double[] moles = new double[n];
    double total = 0.0;
    for (int i = 0; i < n; i++) {
      moles[i] = phase.getComponent(i).getNumberOfMolesInPhase();
      total += moles[i];
    }
    moles[k] += h;
    double moved = total + h;
    double volume = phase.getMolarVolume() * total;
    PhasePCSAFT clone = (PhasePCSAFT) phase.clone();
    for (int i = 0; i < n; i++) {
      ComponentPCSAFT c = (ComponentPCSAFT) clone.getComponent(i);
      c.numberOfMolesInPhase = moles[i];
      c.x = moles[i] / moved;
    }
    clone.numberOfMolesInPhase = moved;
    clone.setMolarVolume(volume / moved);
    clone.volInit();
    reinit(clone);
    return clone;
  }

  private static void reinit(PhasePCSAFT phase) {
    for (int i = 0; i < phase.getNumberOfComponents(); i++) {
      ((ComponentPCSAFT) phase.getComponent(i)).initSAFTDerivatives(phase, phase.getNumberOfComponents(), t, p);
    }
  }

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }
}
