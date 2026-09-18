import neqsim.thermo.component.ComponentEos;
import neqsim.thermo.phase.PhaseEos;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkCPA;

/**
 * For {@code SystemSrkCPA}, {@code ComponentEos.dFdN[i]} is not {@code d getF()/d n_i}.
 *
 * <p>
 * {@code ComponentEos.fugcoef} is
 *
 * <pre>
 * logFugacityCoefficient = dFdN(phase, n, T, P) - log(P V / (R T))
 * </pre>
 *
 * so a component's fugacity coefficient is the derivative of its Helmholtz energy. That makes
 * {@code dFdN} checkable against {@code getF()} directly, with no other code involved: {@code F} is
 * a function of {@code (T, V, n)}, so along any path
 *
 * <pre>
 * dF = sum_i (dF/dn_i)|_{T,V} dn_i + (dF/dV)|_{T,n} dV
 * </pre>
 *
 * and perturbing one mole number at constant {@code T} and {@code P} isolates {@code (dF/dn_i)}
 * once the volume's share is subtracted. No prescribed-volume state is needed.
 *
 * <p>
 * This prints the measurement for an associating mixture and for a plain cubic, which is the
 * control: a plain cubic is self-consistent to ten digits, and an associating one is not.
 *
 * <p>
 * Usage:
 *
 * <pre>
 * javac -proc:none -cp neqsim-3.20.0.jar CpaFdMismatch.java
 * java -cp .:neqsim-3.20.0.jar CpaFdMismatch
 * </pre>
 */
public final class CpaFdMismatch {

  private CpaFdMismatch() {}

  private static final double D = 1.0e-7;

  /** One state: the Helmholtz energy, its volume derivative, the volume, and `dFdN`. */
  private record State(double f, double fv, double v, double[] dFdN) {}

  private static State at(
      boolean cpa, double temperature, double pressure, double first, double second) {
    // Built fresh every time: `init(1)` solves the volume from `T` and `P`, so a perturbed
    // composition has to be a new system rather than a mutated one.
    SystemInterface system =
        cpa ? new SystemSrkCPA(temperature, pressure) : new SystemPrEos(temperature, pressure);
    if (cpa) {
      system.addComponent("water", first);
      system.addComponent("methanol", second);
      system.setMixingRule(10);
    } else {
      system.addComponent("methane", first);
      system.addComponent("n-butane", second);
      system.setMixingRule("classic");
      system.setAttractiveTerm(1);
    }
    system.init(0);
    system.init(1);
    // `init(3)` is the flag that fills the composition and temperature derivatives. Run it
    // too, because a mismatch that vanishes once the system is fully initialised would be a
    // defect in *this probe* rather than in NeqSim - which is the first thing a reader would
    // and should suspect.
    system.init(3);
    PhaseInterface p = system.getPhase(0);
    int n = p.getNumberOfComponents();
    double[] dfdn = new double[n];
    for (int i = 0; i < n; i++) {
      dfdn[i] = ((ComponentEos) p.getComponent(i)).dFdN(p, n, temperature, pressure);
    }
    return new State(((PhaseEos) p).getF(), p.FV(), p.getTotalVolume(), dfdn);
  }

  /** `(dF/dn_i)|_{T,V}` at constant `T` and `P`, once the volume's share is removed. */
  private static double[] isolated(
      boolean cpa, double temperature, double pressure, double first, double second) {
    State base = at(cpa, temperature, pressure, first, second);
    double[] out = new double[base.dFdN.length];
    for (int i = 0; i < out.length; i++) {
      double dFirst = i == 0 ? D : 0.0;
      double dSecond = i == 1 ? D : 0.0;
      State up = at(cpa, temperature, pressure, first + dFirst, second + dSecond);
      State down = at(cpa, temperature, pressure, first - dFirst, second - dSecond);
      out[i] = (up.f - down.f) / (2.0 * D) - base.fv * (up.v - down.v) / (2.0 * D);
    }
    return out;
  }

  private static void report(
      String label, boolean cpa, double temperature, double pressure, double first, double second) {
    State base = at(cpa, temperature, pressure, first, second);
    State doubled = at(cpa, temperature, pressure, 2.0 * first, 2.0 * second);
    double[] isolated = isolated(cpa, temperature, pressure, first, second);

    System.out.printf("%n%s%n", label);
    System.out.printf("  getF()        = %.15g%n", base.f);
    System.out.printf("  F(2n)/F(n)    = %.15g   (2 for an extensive F)%n", doubled.f / base.f);
    System.out.printf("  FV()          = %.15g%n", base.fv);
    for (int i = 0; i < base.dFdN.length; i++) {
      double rel = Math.abs(isolated[i] / base.dFdN[i] - 1.0);
      System.out.printf(
          "  component %d: dFdN = %-22.15g dF/dn_i = %-22.15g  disagree by %5.1f%%%n",
          i, base.dFdN[i], isolated[i], 100.0 * rel);
    }
  }

  public static void main(String[] args) {
    // The associating mixture: `SystemSrkCPA` over water and methanol, the fluid the tranche
    // this was found in works with.
    report(
        "SystemSrkCPA  water/methanol 0.6/0.4 at 356 K, 1 bara",
        true,
        356.0,
        1.0,
        0.6,
        0.4);

    // The control: an ordinary cubic, where the same test passes to ten digits.
    report(
        "SystemPrEos  methane/n-butane 0.6/0.4 at 330 K, 25 bara  (control)",
        false,
        330.0,
        25.0,
        0.6,
        0.4);
  }
}
