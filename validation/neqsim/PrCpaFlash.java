import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrCPA;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * NeqSim's PR-CPA flash over the state ladder, one row per state.
 *
 * <p>
 * The oracle for {@code eos.pr_cpa_phase}. It is deliberately small: the numbers a case
 * records are the phase's own - {@code beta}, {@code Z}, the composition and the
 * fugacity coefficient - and every one of those is a method of {@code PhaseInterface},
 * so this needs no CPA-specific accessor and no cast. Reading the fitted parameters and
 * the association's terms at each layer is what {@link CpaSweep} is for, and it is
 * hard-coded to {@code SystemSrkCPA} because generalising its fifteen casts is only
 * worth doing when a divergence needs localising. A row that none of these keys can
 * explain is that job.
 *
 * <p>
 * Usage, the same shape {@code CpaSweep} takes so the two are read the same way:
 *
 * <pre>
 * javac -proc:none -cp neqsim-3.20.0.jar PrCpaFlash.java
 * java -cp .:neqsim-3.20.0.jar PrCpaFlash [T_K P_bara n_water]...
 * </pre>
 */
public final class PrCpaFlash {

  private PrCpaFlash() {}

  /** The states the ladder covers, the same ones the SRK case does. */
  private static final double[][] DEFAULT_CASES = {
    {350.0, 1.0, 0.6},
    {354.0, 1.0, 0.6},
    {356.0, 1.0, 0.6},
    {358.0, 1.0, 0.6},
    {360.0, 1.0, 0.6},
    {364.0, 1.0, 0.6},
    {368.0, 1.0, 0.6},
  };

  private static void row(double temperature, double pressure, double nWater) {
    SystemInterface system = new SystemPrCPA(temperature, pressure);
    system.addComponent("water", nWater);
    system.addComponent("methanol", 1.0 - nWater);
    system.setMixingRule(10);
    system.init(0);
    system.init(1);

    ThermodynamicOperations operations = new ThermodynamicOperations(system);
    operations.TPflash();

    int phases = system.getNumberOfPhases();
    int n = system.getNumberOfComponents();
    StringBuilder out = new StringBuilder();
    out.append(String.format("T=%.15g\tP=%.15g\tn_water=%.15g\tphases=%d", temperature, pressure, nWater,
        phases));
    // **Whether the association is on at all.** The two counts are not the same number and
    // the difference is the whole finding: `PhaseSrkCPA` sums its components' sites in its
    // init path, and `PhasePrCPA` has the field and the setter and no block that ever sets
    // it. So water carries four sites on its `ComponentSrkCPA` while the phase's own total
    // is zero, and every association term it computes is computed over nothing. A phase
    // count of zero with a component count of four is a PR run wearing a CPA name.
    PhaseInterface first = system.getPhase(0);
    out.append(String.format("\tphase_sites=%d",
        ((neqsim.thermo.phase.PhaseCPAInterface) first).getTotalNumberOfAccociationSites()));
    StringBuilder perComponent = new StringBuilder();
    for (int c = 0; c < n; c++) {
      perComponent.append(c == 0 ? "" : ",");
      perComponent.append(first.getComponent(c).getNumberOfAssociationSites());
    }
    out.append(String.format("\tcomponent_sites=%s", perComponent));
    for (int i = 0; i < phases; i++) {
      PhaseInterface phase = system.getPhase(i);
      out.append(String.format("\tbeta[%d]=%.15g", i, phase.getBeta()));
      out.append(String.format("\tZ[%d]=%.15g", i, phase.getZ()));
      for (int c = 0; c < n; c++) {
        out.append(String.format("\tx[%d][%d]=%.15g", i, c, phase.getComponent(c).getx()));
        out.append(String.format("\tlnPhi[%d][%d]=%.15g", i, c,
            Math.log(phase.getComponent(c).getFugacityCoefficient())));
      }
    }
    System.out.println(out);
  }

  public static void main(String[] args) {
    double[][] cases;
    if (args.length >= 3) {
      cases = new double[args.length / 3][];
      for (int i = 0; i < cases.length; i++) {
        cases[i] =
            new double[] {
              Double.parseDouble(args[3 * i]),
              Double.parseDouble(args[3 * i + 1]),
              Double.parseDouble(args[3 * i + 2]),
            };
      }
    } else {
      cases = DEFAULT_CASES;
    }
    System.out.println("# azoth PrCpaFlash - NeqSim 3.20.0's SystemPrCPA with setMixingRule(10).");
    System.out.println("# Phase 0 is whatever the flash ordered first; `PhaseInterface.getBeta`");
    System.out.println("# is the phase's fraction, so the vapour is the smaller one at these");
    System.out.println("# states. `lnPhi` is `log(ComponentInterface.getFugacityCoefficient())`.");
    for (double[] c : cases) {
      try {
        row(c[0], c[1], c[2]);
      } catch (Throwable thrown) {
        System.out.printf("T=%.15g\tP=%.15g\tn_water=%.15g\terror=%s%n", c[0], c[1], c[2],
            thrown.getClass().getSimpleName());
      }
    }
  }
}
