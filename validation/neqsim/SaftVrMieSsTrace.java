import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSAFTVRMie;

/**
 * NeqSim's successive-substitution loop for a SAFT-VR-Mie flash, run in the open.
 *
 * <p>
 * {@code TPflashSAFT} reports a single phase for methane/n-butane at every state tried,
 * including states its own Wilson seed calls two-phase - so this reproduces its loop step by
 * step, from the same seed and with the same nested per-phase solve, and prints what the
 * convergence does. It is the port's specification: the scheme is short enough to write out,
 * and writing it out is what says whether the class's answer is the algorithm's or a failure
 * inside it.
 *
 * <pre>
 * javac -proc:none -cp neqsim-3.20.0.jar SaftVrMieSsTrace.java
 * java -cp .:neqsim-3.20.0.jar SaftVrMieSsTrace
 * </pre>
 */
public final class SaftVrMieSsTrace {

  private SaftVrMieSsTrace() {}

  /** NeqSim's own nested solve: a fresh system at the trial composition, one phase. */
  private static double[] phi(String[] names, double t, double pBar, double[] comp, PhaseType pt) {
    SystemInterface sys = new SystemSAFTVRMie(t, pBar);
    for (int i = 0; i < names.length; i++) {
      sys.addComponent(names[i], comp[i]);
    }
    sys.setMixingRule("classic");
    sys.init(0);
    sys.setPhaseType(0, pt);
    sys.init(1);
    double[] out = new double[names.length];
    for (int i = 0; i < names.length; i++) {
      out[i] = sys.getPhase(0).getComponent(i).getFugacityCoefficient();
    }
    return out;
  }

  private static double rachfordRice(double[] z, double[] k) {
    double lo = 0.0;
    double hi = 1.0;
    for (int iter = 0; iter < 200; iter++) {
      double mid = 0.5 * (lo + hi);
      double f = 0.0;
      for (int i = 0; i < z.length; i++) {
        f += z[i] * (k[i] - 1.0) / (1.0 + mid * (k[i] - 1.0));
      }
      if (f > 0.0) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    return 0.5 * (lo + hi);
  }

  private static void trace(double t, double pBar, String[] names, double[] z) {
    SystemInterface system = new SystemSAFTVRMie(t, pBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.setMixingRule("classic");
    system.init(0);

    int nc = names.length;
    double[] k = new double[nc];
    for (int i = 0; i < nc; i++) {
      ComponentInterface c = system.getPhase(0).getComponent(i);
      k[i] = c.getK() > 0 && !Double.isNaN(c.getK()) ? c.getK() : 1.0;
    }

    System.out.printf("%n# %s at T=%g K, P=%g bara%n", String.join("/", names), t, pBar);
    System.out.printf("seed K = [%.15g, %.15g]%n", k[0], k[1]);
    double[] x = new double[nc];
    double[] y = new double[nc];
    double beta = 0.5;
    for (int iter = 0; iter < 50; iter++) {
      beta = rachfordRice(z, k);
      if (beta < 1.0e-10 || beta > 1.0 - 1.0e-10) {
        System.out.printf("iter %d: beta = %.15g is outside the bounds, loop breaks%n", iter, beta);
        return;
      }
      for (int i = 0; i < nc; i++) {
        x[i] = z[i] / (1.0 + beta * (k[i] - 1.0));
        y[i] = k[i] * x[i];
      }
      double sumx = x[0] + x[1];
      double sumy = y[0] + y[1];
      for (int i = 0; i < nc; i++) {
        x[i] /= sumx;
        y[i] /= sumy;
      }
      double[] phiL = phi(names, t, pBar, x, PhaseType.LIQUID);
      double[] phiG = phi(names, t, pBar, y, PhaseType.GAS);
      double maxDelta = 0.0;
      for (int i = 0; i < nc; i++) {
        double knew = phiL[i] / phiG[i];
        if (Double.isNaN(knew) || knew <= 0) {
          knew = k[i];
        }
        maxDelta = Math.max(maxDelta, Math.abs(knew - k[i]) / Math.max(Math.abs(k[i]), 1.0e-10));
        k[i] = knew;
      }
      System.out.printf(
          "iter %d: beta=%.15g x=[%.10g, %.10g] y=[%.10g, %.10g] phiL=[%.10g, %.10g] "
              + "phiG=[%.10g, %.10g] K=[%.10g, %.10g] maxDeltaK=%.6g%n",
          iter, beta, x[0], x[1], y[0], y[1], phiL[0], phiL[1], phiG[0], phiG[1], k[0], k[1],
          maxDelta);
      if (maxDelta < 1.0e-6) {
        System.out.println("converged");
        return;
      }
    }
    System.out.println("fifty iterations without converging");
  }

  public static void main(String[] args) {
    System.out.println("# azoth SaftVrMieSsTrace - NeqSim's TPflashSAFT loop, in the open.");
    for (double[] state : new double[][] {{250.0, 30.0}, {200.0, 30.0}}) {
      trace(state[0], state[1], new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4});
    }
  }
}
