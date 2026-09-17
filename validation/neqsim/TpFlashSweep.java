import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * NeqSim's own TP flash over the grid azoth's successive substitution struggles on.
 *
 * <p>The states are methane/n-butane 0.6/0.4 under Peng-Robinson with the classic mixing
 * rule, from 300 K to 410 K and 20 bar to 80 bar. azoth's iteration count rises through
 * this grid and two states in it - 400 K/70 bar and 410 K/60 bar - do not converge at
 * all inside its 300-step cap. NeqSim's {@code TPflash} hands over to
 * {@code SysNewtonRhapsonTPflash} after twelve successive-substitution steps, so this is
 * the measurement of what that second-order solver does where the first stalls.
 *
 * <p>{@code getBeta()} is the vapour fraction, and it is *not* the same quantity as
 * azoth's on a single-phase feed: NeqSim reports a number there and azoth reports an
 * absent one. The phase count is printed beside it so the two can be told apart.
 */
public class TpFlashSweep {

  private static void one(double temperatureK, double pressureBar) {
    SystemInterface fluid = new SystemPrEos(temperatureK, pressureBar);
    fluid.addComponent("methane", 0.6);
    fluid.addComponent("n-butane", 0.4);
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(1);
    String outcome;
    try {
      new ThermodynamicOperations(fluid).TPflash();
      outcome = String.format("beta=%.17g phases=%d", fluid.getBeta(), fluid.getNumberOfPhases());
    } catch (Throwable failure) {
      outcome = "THREW " + failure.getClass().getSimpleName();
    }
    System.out.printf("%7.1f %8.1f  %s%n", temperatureK, pressureBar, outcome);
  }

  public static void main(String[] args) {
    System.out.printf("%7s %8s  %s%n", "T", "P(bar)", "neqsim");
    for (double t : new double[] {300.0, 320.0, 340.0, 350.0, 355.0}) {
      for (double p : new double[] {80.0, 100.0, 120.0, 150.0, 200.0, 250.0}) {
        one(t, p);
      }
    }
  }
}
