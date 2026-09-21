import neqsim.thermo.component.ComponentEos;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhasePrCPA;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrCPA;
import neqsim.thermo.system.SystemPrEos;

/**
 * Which fitted set NeqSim's {@code SystemPrCPA} puts into the Peng-Robinson cubic.
 *
 * <p>
 * {@code PhasePrCPA} never sets its association-site count, so the class associates over
 * nothing and its flash is a Peng-Robinson - which is what {@code validation/neqsim/
 * PrCpaFlash.java} records and what {@code specs/cases/eos/pr_cpa_phase.toml} says. That
 * much was known. <b>What was not is that its attraction and covolume are not the
 * Peng-Robinson correlation's either.</b>
 *
 * <p>
 * This prints, for a pair neither of whose components associates, the plain
 * {@code SystemPrEos} component beside the {@code SystemPrCPA} one: {@code a}, {@code b},
 * the acentric factor and the attractive term number. The two differ - on n-butane by
 * `4.3%` in `a` - and neither the composition nor the association can be the cause, since
 * both are zero-sited here.
 *
 * <p>
 * <b>It is not a Peng-Robinson component at all.</b> Reading {@code a} and {@code b} back
 * off the two states below gives, for methane,
 *
 * <pre>
 *   Omega_a = a 1e-5 Pc /(R^2 Tc^2) = 0.427480233540342
 *   Omega_b = b 1e-5 Pc /(R Tc)     = 0.08664034996495783
 * </pre>
 *
 * which are the **Soave-Redlich-Kwong** constants, not Peng-Robinson's `0.45724`/`0.07780`
 * - to fifteen digits, and T-independent across 250 and 350 K. On top of that, where a
 * component has a fitted CPA set the SRK-CPA one is used: n-butane's `a` and `b` are
 * `aCPA_SRK`/`bCPA_SRK` verbatim (`131427.4`/`7.2081`), and its acentric factor is
 * `mCPA_SRK` inverted through Peng-Robinson's own `m(w) = 0.37464 + 1.54226 w - 0.26992 w^2`
 * (`0.224807298178772` to sixteen digits).
 *
 * <p>
 * azoth's {@code eos.pr_cpa_phase} is an actual Peng-Robinson reading the
 * {@code acpa_pr}/{@code bcpa_pr} columns - zero throughout the vendored table - so on this
 * non-associating pair it applies no substitution at all and the two disagree by 4% in
 * `Z`. NeqSim's class is a third thing: SRK constants, an SRK-CPA fit and a PR cubic.
 *
 * <pre>
 * javac -proc:none -cp neqsim-f0c7436.jar PrCpaSubstitution.java
 * java -cp .:neqsim-f0c7436.jar PrCpaSubstitution
 * </pre>
 */
public final class PrCpaSubstitution {

  private PrCpaSubstitution() {}

  private static void show(String label, SystemInterface system) {
    system.setMixingRule(10);
    system.init(0);
    system.init(1);
    PhaseInterface phase = system.getPhase(0);
    System.out.printf("%s%n", label);
    for (int i = 0; i < system.getNumberOfComponents(); i++) {
      ComponentInterface c = phase.getComponent(i);
      ComponentEos eos = (ComponentEos) c;
      System.out.printf(
          "   %-9s term=%d Tc=%.15g Pc=%.15g omega=%.15g a=%.15g b=%.15g sites=%d%n",
          c.getComponentName(), eos.getAttractiveTermNumber(), c.getTC(), c.getPC(),
          c.getAcentricFactor(), eos.geta(), eos.getb(), c.getNumberOfAssociationSites());
    }
    String sites = "-";
    String fcpa = "-";
    if (phase instanceof PhasePrCPA) {
      sites = Integer.toString(((PhasePrCPA) phase).getTotalNumberOfAccociationSites());
      fcpa = String.format("%.15g", ((PhasePrCPA) phase).FCPA());
    }
    System.out.printf("   %s sites=%s FCPA=%s Z=%.15g%n", phase.getClass().getSimpleName(), sites,
        fcpa, phase.getZ());
  }

  public static void main(String[] args) {
    System.out.println("# azoth PrCpaSubstitution - which fitted set PhasePrCPA carries.");
    System.out.println("# NeqSim master. `a`/`b` are in the class's own internal scale.");
    for (double[] state : new double[][] {{350.0, 30.0}, {250.0, 30.0}}) {
      SystemInterface pr = new SystemPrEos(state[0], state[1]);
      pr.addComponent("methane", 0.6);
      pr.addComponent("n-butane", 0.4);
      show(String.format("# SystemPrEos  at T=%g K, P=%g bara", state[0], state[1]), pr);

      SystemInterface cpa = new SystemPrCPA(state[0], state[1]);
      cpa.addComponent("methane", 0.6);
      cpa.addComponent("n-butane", 0.4);
      show(String.format("# SystemPrCPA at T=%g K, P=%g bara", state[0], state[1]), cpa);
    }
  }
}
