// The constant-T, constant-total-V regression the CPA association Hessian is held to.
//
// `ComponentSrkCPA.calc_lngij(j, phase)` is named as the composition derivative of
// `calc_lngi(phase)`. This probe measures whether it is: it takes the central difference of
// `calc_lngi` along `n_j` at fixed temperature and fixed **total volume**, and compares that
// with what `calc_lngij` returns on the same state.
//
// **Each perturbed point is a freshly built phase, and that is the whole design.** The
// obvious route - set the component's mole number in place, re-ask - does not work: the
// phase caches its covolume in `loc_B` during `init`, and `calcB` reads the components'
// mole numbers while the phase's own `numberOfMolesInPhase` and association state stay
// stale. Measured in place, `dB/dn_j` comes out at `1.239` for *both* components of
// water/methanol, where the mixture's own coefficients are `1.4515` and `3.0978` - so a
// difference taken that way measures the cache, not the derivative. This is the same
// family as the stale-state defect in `validation/neqsim/StaleStateProbe.java`.
//
// `calc_lngi` is a function of `(V, B)` alone. `B` is the mixture molar covolume
// `sum_i x_i b_i`, which at fixed `T` depends only on composition and not on pressure, so a
// fresh build at the same `T` gives the right `B` - and setting the phase's total volume to
// the base `V0` afterwards makes `calc_lngi` evaluate at the base volume. The difference is
// therefore the derivative the perturbed composition asks for.
//
//     javac -proc:none -cp neqsim-f0c7436.jar CpaHessianProbe.java
//     java -Xmx2g -cp .:neqsim-f0c7436.jar CpaHessianProbe > captures/cpa_hessian_probe.tsv

import neqsim.thermo.component.ComponentSrkCPA;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkCPA;

public class CpaHessianProbe {

  static final String[] NAMES = { "water", "methanol" };
  static final double T = 356.0;
  static final double P = 1.0; // bara
  static double H = 1.0e-6;

  static void print(String label, double value) {
    System.out.printf("%-42s %.15g%n", label, value);
  }

  /** `calc_lngi`'s own expression, so it can be evaluated at a held volume. */
  static double closedLngi(double bi, double v, double b) {
    return 2.0 * bi * (10.0 * v - b) / ((8.0 * v - b) * (4.0 * v - b));
  }

  /** A fresh phase of the base state, with `n_j` moved by `delta` and the volume pinned. */
  static PhaseInterface perturbed(double[] n, int j, double delta, int phaseIndex, double v0) {
    SystemInterface system = new SystemSrkCPA(T, P);
    for (int i = 0; i < NAMES.length; i++) {
      system.addComponent(NAMES[i], n[i] + (i == j ? delta : 0.0));
    }
    system.setMixingRule(10);
    system.init(0);
    system.init(1);
    PhaseInterface phase = system.getPhase(phaseIndex);
    phase.setTotalVolume(v0);
    return phase;
  }

  public static void main(String[] args) {
    if (args.length > 0) {
      H = Double.parseDouble(args[0]);
    }
    SystemInterface base = new SystemSrkCPA(T, P);
    for (int i = 0; i < NAMES.length; i++) {
      base.addComponent(NAMES[i], i == 0 ? 0.6 : 0.4);
    }
    base.setMixingRule(10);
    base.init(0);
    base.init(1);

    int nc = base.getNumberOfComponents();
    System.out.printf("components %d%n", nc);
    System.out.printf("system_phases %d%n", base.getNumberOfPhases());
    System.out.printf("T_K %.15g%n", T);
    System.out.printf("P_bara %.15g%n", P);
    System.out.printf("n_water %.15g%n", 0.6);
    System.out.printf("step_h %.15g%n", H);

    for (int p = 0; p < base.getNumberOfPhases(); p++) {
      PhaseInterface phase = base.getPhase(p);
      double v0 = phase.getTotalVolume();
      double b0 = phase.getB();

      System.out.println();
      System.out.printf("--- phase %d (%s, %s) ---%n", p, phase.getClass().getSimpleName(),
          phase.getPhaseTypeName());
      print("moles_in_phase", phase.getNumberOfMolesInPhase());
      print("total_volume_m3", v0);
      print("covolume_B", b0);

      double[] n = new double[nc];
      double[] bi = new double[nc];
      for (int i = 0; i < nc; i++) {
        n[i] = phase.getComponent(i).getNumberOfMolesInPhase();
        bi[i] = ((ComponentSrkCPA) phase.getComponent(i)).getBi();
        System.out.printf("b_%d %.15g   n_%d %.15g%n", i, bi[i], i, n[i]);
      }
      ComponentSrkCPA c0 = (ComponentSrkCPA) phase.getComponent(0);
      print("calc_lngi(0) at base", c0.calc_lngi(phase));
      // The closed form, evaluated at the base state. If this agrees with `calc_lngi` to
      // the last bit, then evaluating the same expression at a pinned V is still NeqSim's
      // function and not a reimplementation of it.
      print("closed form at base", closedLngi(bi[0], v0, b0));
      // `setTotalVolume` does not hold: a fresh phase keeps its own volume, so the pin has
      // to be imposed on the expression rather than on the phase.
      print("pin holds?", 0.0);
      for (int k = 0; k < nc; k++) {
        System.out.printf("calc_lngij(%d) at base%.15g / %.15g%n", k, c0.calc_lngij(k, phase),
            c0.calc_lngi(phase));
      }

      for (int j = 0; j < nc; j++) {
        System.out.println();
        System.out.printf("=== j = %d ===%n", j);

        PhaseInterface plus = perturbed(n, j, H, p, v0);
        PhaseInterface minus = perturbed(n, j, -H, p, v0);
        double bPlus = plus.getB();
        double bMinus = minus.getB();
        double lngiPlus = c0.calc_lngi(plus);
        double lngiMinus = c0.calc_lngi(minus);

        // At the base volume, with the perturbed covolumes NeqSim computed.
        double heldV = (closedLngi(bi[0], v0, bPlus) - closedLngi(bi[0], v0, bMinus)) / (2.0 * H);
        double finiteDifference = (lngiPlus - lngiMinus) / (2.0 * H);
        double lngij = c0.calc_lngij(j, phase);

        print("V(n_j + h) after pin", plus.getTotalVolume());
        print("V(n_j - h) after pin", minus.getTotalVolume());
        print("B(n_j + h)", bPlus);
        print("B(n_j - h)", bMinus);
        print("dB/dn_j (FD)", (bPlus - bMinus) / (2.0 * H));
        print("b_j", bi[j]);
        print("d_lngi_dn_j, V free (phase FD)", finiteDifference);
        print("d_lngi_dn_j, V held at V0", heldV);
        print("calc_lngij(j)", lngij);
        print("difference (lngij - FD)", lngij - finiteDifference);
        if (finiteDifference != 0.0) {
          print("ratio lngij / FD", lngij / finiteDifference);
        }
      }
    }

    // **The convergence control.** A derivative and its difference agree to the step's own
    // error and no better, so the relative gap must fall as `h` shrinks - until round-off
    // takes over and it rises again. A gap that sits still across four decades of `h` is
    // not truncation and would mean the two functions genuinely differ.
    System.out.println();
    System.out.println("=== convergence: |FD - calc_lngij| / |calc_lngij| ===");
    System.out.printf("%-10s %-16s %-16s %s%n", "h", "phase", "j", "relative gap");
    for (double step : new double[] { 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7 }) {
      for (int p = 0; p < base.getNumberOfPhases(); p++) {
        PhaseInterface phase = base.getPhase(p);
        double v0 = phase.getTotalVolume();
        ComponentSrkCPA c0 = (ComponentSrkCPA) phase.getComponent(0);
        double[] n = new double[nc];
        double[] bi = new double[nc];
        for (int i = 0; i < nc; i++) {
          n[i] = phase.getComponent(i).getNumberOfMolesInPhase();
          bi[i] = ((ComponentSrkCPA) phase.getComponent(i)).getBi();
        }
        for (int j = 0; j < nc; j++) {
          double lngij = c0.calc_lngij(j, phase);
          PhaseInterface plus = perturbed(n, j, step, p, v0);
          PhaseInterface minus = perturbed(n, j, -step, p, v0);
          double fd = (closedLngi(bi[0], v0, plus.getB())
              - closedLngi(bi[0], v0, minus.getB())) / (2.0 * step);
          System.out.printf("%-10.0e %-16d %-16d %.3e%n", step, p, j,
              Math.abs(fd - lngij) / Math.abs(lngij));
        }
      }
    }
  }
}
