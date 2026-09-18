import neqsim.thermo.component.ComponentEos;
import neqsim.thermo.phase.PhaseEos;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkCPA;

/**
 * For a {@code SystemSrkCPA}, {@code ComponentEos.dFdN[i]} is {@code d getF()/d n_i} at
 * constant {@code T} and {@code V} - measured, not read, and to ten digits.
 *
 * <p>
 * {@code ComponentEos.fugcoef} is
 *
 * <pre>
 * logFugacityCoefficient = dFdN(phase, n, T, P) - log(P V / (R T))
 * </pre>
 *
 * so a component's fugacity coefficient is the derivative of the Helmholtz energy. That makes
 * {@code dFdN} checkable against {@code getF()} directly, with no other code involved: {@code F} is
 * a function of {@code (T, V, n)}, so along any path
 *
 * <pre>
 * dF = sum_i (dF/dn_i)|_{T,V} dn_i + (dF/dV)|_{T,n} dV
 * </pre>
 *
 * and perturbing one mole number at constant {@code T} and {@code P} isolates
 * {@code (dF/dn_i)|_{T,V}} once the volume's share is subtracted. No prescribed-volume state is
 * needed.
 *
 * <p>
 * <b>The subtlety that makes this probe worth keeping is which volume derivative to subtract.</b>
 * {@code PhaseEos.dFdV()} is {@code return FV();}, and {@code PhaseSrkCPA.dFdV()} overrides it as
 * {@code super.dFdV() + cpaon * dFCPAdV()}. For a CPA phase the two are not close: at 356 K and
 * 1 bara the association's {@code dFCPAdV} is {@code 1.694e-5} against a cubic {@code FV()} of
 * {@code 9.75e-7}, seventeen times larger. Subtracting {@code FV()} leaves the association's share
 * of {@code dV} in the residual, and it reports {@code dFdN} as 38-54% away from the derivative of
 * {@code F} - which is what this probe printed before it used {@code dFdV()}, and the reading was
 * taken as a defect in NeqSim rather than in the probe.
 *
 * <p>
 * Both derivatives are therefore printed on every line: they are the two candidate corrections, and
 * on the plain-cubic control they are the same number, which is exactly why a control that passes
 * cannot find this.
 *
 * <p>
 * Usage:
 *
 * <pre>
 * javac -proc:none -cp neqsim-3.20.0.jar CpaFdProbe.java
 * java -cp .:neqsim-3.20.0.jar CpaFdProbe
 * </pre>
 */
public final class CpaFdProbe {

  private CpaFdProbe() {}

  private static final double D = 1.0e-7;

  /** One state: the Helmholtz energy, both volume derivatives, the volume, and `dFdN`. */
  private record State(double f, double fv, double dfdv, double v, double[] dFdN) {}

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
    // `init(3)` fills the composition and temperature derivative caches. It is run because
    // agreement that appeared only once the system was fully initialised would say `dFdN`
    // and `getF` are fed from different caches - and it is *measured* to move the state:
    // `getF()` goes from -0.0521800519231376 to -0.0521800519442854 and `FCPA` from
    // -0.0494373456962294 to -0.0494373457173083, a relative 4e-10. So these readings are
    // this program's, at `init(3)`, and `CpaSweep`'s are at `init(1)`; comparing the two
    // across the tenth digit compares the two initialisations rather than the two quantities.
    system.init(3);
    PhaseInterface p = system.getPhase(0);
    int n = p.getNumberOfComponents();
    double[] dfdn = new double[n];
    for (int i = 0; i < n; i++) {
      dfdn[i] = ((ComponentEos) p.getComponent(i)).dFdN(p, n, temperature, pressure);
    }
    return new State(((PhaseEos) p).getF(), p.FV(), p.dFdV(), p.getTotalVolume(), dfdn);
  }

  /**
   * `(dF/dn_i)|_{T,V}` at constant `T` and `P`, once the volume's share is removed - by
   * `dFdV()` in row 0 and by `FV()` in row 1, so the two candidate corrections are both on
   * the record. On a plain cubic the rows are identical.
   */
  private static double[][] isolated(
      boolean cpa, double temperature, double pressure, double first, double second) {
    State base = at(cpa, temperature, pressure, first, second);
    double[][] out = new double[2][base.dFdN.length];
    for (int i = 0; i < base.dFdN.length; i++) {
      double dFirst = i == 0 ? D : 0.0;
      double dSecond = i == 1 ? D : 0.0;
      State up = at(cpa, temperature, pressure, first + dFirst, second + dSecond);
      State down = at(cpa, temperature, pressure, first - dFirst, second - dSecond);
      double dF = (up.f - down.f) / (2.0 * D);
      double dV = (up.v - down.v) / (2.0 * D);
      out[0][i] = dF - base.dfdv * dV;
      out[1][i] = dF - base.fv * dV;
    }
    return out;
  }

  private static void report(
      String label, boolean cpa, double temperature, double pressure, double first, double second) {
    State base = at(cpa, temperature, pressure, first, second);
    State doubled = at(cpa, temperature, pressure, 2.0 * first, 2.0 * second);
    double[][] isolated = isolated(cpa, temperature, pressure, first, second);

    System.out.printf("%n%s%n", label);
    System.out.printf("  getF()        = %.15g%n", base.f);
    System.out.printf("  F(2n)/F(n)    = %.15g   (2 for an extensive F)%n", doubled.f / base.f);
    System.out.printf("  FV()          = %.15g   (the cubic part alone)%n", base.fv);
    System.out.printf(
        "  dFdV()        = %.15g   (= FV() + dFCPAdV() for a CPA phase)%n", base.dfdv);
    for (int i = 0; i < base.dFdN.length; i++) {
      System.out.printf("  component %d: dFdN = %-22.15g%n", i, base.dFdN[i]);
      System.out.printf("               dF/dn_i less dFdV()*dV = %-22.15g  agree to %5.1e%n",
          isolated[0][i], Math.abs(isolated[0][i] / base.dFdN[i] - 1.0));
      System.out.printf("               dF/dn_i less FV()*dV   = %-22.15g  agree to %5.1e%n",
          isolated[1][i], Math.abs(isolated[1][i] / base.dFdN[i] - 1.0));
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

    // The control: an ordinary cubic, where `FV()` and `dFdV()` are the same method, so this
    // test cannot tell the two apart and passes either way.
    report(
        "SystemPrEos  methane/n-butane 0.6/0.4 at 330 K, 25 bara  (control)",
        false,
        330.0,
        25.0,
        0.6,
        0.4);
  }
}
