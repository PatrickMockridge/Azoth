import java.util.LinkedHashMap;
import java.util.Map;

import neqsim.thermo.component.ComponentEos;
import neqsim.thermo.component.ComponentSrkCPA;
import neqsim.thermo.phase.PhaseCPAInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseSrkCPA;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkCPA;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * NeqSim's CPA internals over a grid of states, one row per state, in the order the model
 * is built.
 *
 * <p>
 * <b>Every key is declared in the header with the class and method it comes from and the
 * azoth field it corresponds to.</b> That is not documentation, it is the instrument. An
 * earlier probe printed {@code getAi} beside {@code calca} and it was read as the pure
 * component's attraction; {@code EosMixingRuleHandler.calcAi} shows it is in fact the
 * interaction <em>row sum</em>, {@code sum_j n_j sqrt(a_i a_j)}, so a whole session was
 * spent on an alpha "divergence" that never existed. A quantity whose provenance is not on
 * its own line is a quantity that will be misread.
 *
 * <p>
 * <b>Scales are converted to SI and said so.</b> NeqSim carries {@code a} in
 * {@code Pa*m^6/mol^2 * 1e5} and {@code b} in {@code m^3/mol * 1e5}
 * ({@code Component.java:526-531}); the header states the division, so a comparison with
 * azoth never has to guess a factor.
 *
 * <p>
 * The rows are ordered by layer - the fluid's pure constants, then the mixture's, then the
 * root, then the distribution function, then the site fractions, then the Helmholtz energy,
 * then the fugacity, then the flash - so that a disagreement can be attributed to the first
 * layer that shows one rather than to the number that happened to be compared.
 *
 * <p>
 * The fluid is water/methanol, the one this tranche works in: {@code SystemSrkCPA} builds
 * {@code ComponentSrkCPA} for every component it carries, and a general fluid list would
 * make "which component" the first question rather than the association.
 *
 * <p>
 * {@code CpaProbe} remains the single-state kernel probe the Rust tests cite; this is the
 * grid driver, and it reproduces every {@code CpaProbe} reading as a one-case run.
 *
 * <p>
 * Usage:
 *
 * <pre>
 * javac -proc:none -cp neqsim-3.20.0.jar CpaSweep.java
 * java -cp .:neqsim-3.20.0.jar CpaSweep [T_K P_bara n_water]...
 * </pre>
 *
 * With no arguments it sweeps a default grid: the states this tranche has measured, plus a
 * pressure and temperature ladder around the two-phase region.
 */
public final class CpaSweep {

  /** The default grid, as {@code {T_K, P_bara, n_water}} triples. */
  private static final double[][] DEFAULT_CASES = {
    // The state the kernel was first verified at, and the state the root was verified at.
    {300.0, 100.0, 0.6},
    {300.0, 1.0, 0.6},
    // The low-pressure region where the root finder was wrong: NeqSim splits 0.6/0.4.
    {356.0, 1.0, 0.6},
    // Single phase in NeqSim as well - the control that says the instrument can say "no".
    {340.0, 1.0, 0.6},
    {356.0, 1.0, 0.8340697447584025},
    // A ladder through the two-phase region at 1 bar, and one above it.
    {350.0, 1.0, 0.6},
    {360.0, 1.0, 0.6},
    {366.0, 1.0, 0.6},
    {372.0, 1.0, 0.6},
    {356.0, 0.5, 0.6},
    {356.0, 2.0, 0.6},
    {356.0, 5.0, 0.6},
    {-1, -1, -1},
  };

  private CpaSweep() {}

  /**
   * One key's declaration: the NeqSim method it comes from, the azoth field it corresponds
   * to, and any conversion applied. Printed once, in key order.
   */
  private static final String[][] KEYS = {
    // --- the state -------------------------------------------------------------------
    {"T_K", "system.getTemperature()", "ReducedParameters.t_kelvin", ""},
    {"P_bara", "system.getPressure()", "reduced.pressure / 1e5", ""},
    {"n_water", "the driver's argument", "x[0] (methanol is 1 - x[0])", ""},
    // --- the fluid's pure constants --------------------------------------------------
    // `aT` and `sqrtAT` are public fields on ComponentEos, set in `seta`/`setm` from
    // `aT(temperature)`. For a CPA component with a fitted set this is `aCPA * alpha(T)`
    // with the Soave alpha carrying the *fitted* `mCPA` (ComponentSrkCPA:142,
    // `getAttractiveTerm().setm(mCPA)`), which is the whole of the substitution.
    {"aT[0]", "ComponentEos.aT / 1e5", "reduced.a[0] * R^2 T^2 / P", "water; SI Pa*m^6/mol^2"},
    {"aT[1]", "ComponentEos.aT / 1e5", "reduced.a[1] * R^2 T^2 / P", "methanol"},
    {"sqrtAT[0]", "ComponentEos.sqrtAT / sqrt(1e5)", "sqrt(reduced.a[0] * R^2 T^2 / P)", "water"},
    {"sqrtAT[1]", "ComponentEos.sqrtAT / sqrt(1e5)", "sqrt(reduced.a[1] * R^2 T^2 / P)", "methanol"},
    {"b[0]", "ComponentEos.getBi() / 1e5", "reduced.b[0] * R T / P", "water; m^3/mol"},
    {"b[1]", "ComponentEos.getBi() / 1e5", "reduced.b[1] * R T / P", "methanol"},
    {"calca[0]", "ComponentSrkCPA.calca() / 1e5", "the fitted aCPA, SI", "water; internal * 1e5"},
    {"calca[1]", "ComponentSrkCPA.calca() / 1e5", "the fitted aCPA, SI", "methanol"},
    {"calcb[0]", "ComponentSrkCPA.calcb() / 1e5", "the fitted bCPA, SI", "water"},
    {"calcb[1]", "ComponentSrkCPA.calcb() / 1e5", "the fitted bCPA, SI", "methanol"},
    {"scheme[0]", "ComponentSrkCPA.getAssociationScheme()", "SiteScheme", "water"},
    {"scheme[1]", "ComponentSrkCPA.getAssociationScheme()", "SiteScheme", "methanol"},
    {"energy[0]", "ComponentSrkCPA.getAssociationEnergy()", "AssociationComponent.energy", "J/mol"},
    {"energy[1]", "ComponentSrkCPA.getAssociationEnergy()", "AssociationComponent.energy", "J/mol"},
    {"volume[0]", "ComponentSrkCPA.getAssociationVolume()", "AssociationComponent.volume", "kappa_AB"},
    {"volume[1]", "ComponentSrkCPA.getAssociationVolume()", "AssociationComponent.volume", "kappa_AB"},
    // --- the mixture -----------------------------------------------------------------
    // `getA()` and `getB()` are the *mixture's* reduced-scale parameters in NeqSim's
    // internal scale, so dividing by 1e5 gives the dimensional `a_mix`.
    // **`getAi()` is NOT a pure component's `a`**: `EosMixingRuleHandler.calcAi` returns
    // `sum_j n_j (1 - k_ij) sqrt(a_i a_j)`, the interaction row sum. It is printed here so
    // that it is on the record as what it is.
    {"A_mix", "PhaseSrkEos.getA() / 1e5", "a_mix * R^2 T^2 / P", "SI Pa*m^6/mol^2"},
    {"B_mix", "PhaseSrkEos.getB() / 1e5", "b_mix * R T / P", "m^3/mol"},
    {"Ai[0]", "PhaseSrkEos.getAi() / 1e5", "abar[0] * R^2 T^2 / P",
        "ROW SUM, not a pure a; measured 2x azoth's abar at every state, unexplained"},
    {"Ai[1]", "PhaseSrkEos.getAi() / 1e5", "abar[1] * R^2 T^2 / P",
        "ROW SUM, not a pure a; measured 2x azoth's abar at every state, unexplained"},
    // --- the root --------------------------------------------------------------------
    {"Z", "PhaseInterface.getZ()", "PhaseState.z", ""},
    // **`getMolarVolume()` is NOT m^3/mol.** It is 1e5 x the SI value, the same internal
    // scale `b` is carried in. `PV_over_RT` is what proves it: `ComponentEos.fugcoef` is
    // `ln phi = dFdN - ln(getPressure() * getMolarVolume() / (R T))`, and that argument
    // comes out exactly Z, because the pressure is bara (1e-5 x SI) and the volume is
    // 1e5 x SI and the two factors cancel. Without this key the scale is something a
    // reader has to infer, and the first reader inferred it wrong.
    {"v_molar", "PhaseInterface.getMolarVolume()", "z * R T / P * 1e5",
        "**1e5 x m^3/mol** - NeqSim's internal volume scale"},
    {"PV_over_RT", "getPressure() * getMolarVolume() / (R T)", "z",
        "the argument inside ComponentEos.fugcoef's logarithm; = Z confirms both scales"},
    // --- the distribution function ---------------------------------------------------
    {"gcpa", "PhaseCPAInterface.getGcpa()", "Rdf.g", "at contact, dimensionless"},
    {"gcpav", "PhaseCPAInterface.getGcpav()", "Rdf.d_ln_g_dv", "d ln g / dV"},
    {"sites", "PhaseCPAInterface.getTotalNumberOfAccociationSites()", "Association.site_count()", ""},
    // --- the site fractions ----------------------------------------------------------
    {"xsite[0][0]", "ComponentSrkCPA.getXsite()[0]", "SiteState.fractions[0]", "water site 0"},
    {"xsite[0][1]", "ComponentSrkCPA.getXsite()[1]", "SiteState.fractions[1]", "water site 1"},
    {"xsite[0][2]", "ComponentSrkCPA.getXsite()[2]", "SiteState.fractions[2]", "water site 2"},
    {"xsite[0][3]", "ComponentSrkCPA.getXsite()[3]", "SiteState.fractions[3]", "water site 3"},
    {"xsite[1][0]", "ComponentSrkCPA.getXsite()[0]", "SiteState.fractions[4]", "methanol site 0"},
    {"xsite[1][1]", "ComponentSrkCPA.getXsite()[1]", "SiteState.fractions[5]", "methanol site 1"},
    // --- the Helmholtz energy --------------------------------------------------------
    {"FCPA", "PhaseSrkCPA.FCPA()", "SiteState.helmholtz_rt", "A_assoc / (R T)"},
    {"dFCPAdV", "PhaseSrkCPA.dFCPAdV()", "SiteState.d_helmholtz_dv", "d(A_assoc/(R T)) / dV"},
    {"dFCPAdT", "PhaseSrkCPA.dFCPAdT()", "SiteDerivatives.d_helmholtz_dt", "d(A_assoc/(R T)) / dT"},
    {"dFCPAdN[0]", "ComponentSrkCPA.dFCPAdN(phase, n, T, P)", "SiteState.ln_phi[0]", "water"},
    {"dFCPAdN[1]", "ComponentSrkCPA.dFCPAdN(phase, n, T, P)", "SiteState.ln_phi[1]", "methanol"},
    // The kernel's second composition derivative, which CpaProbe's finding is about:
    // `calc_lngij` is not the derivative of `calc_lngi`, so this matrix is NOT azoth's
    // `d_ln_phi_dn`. Printed so that finding stays reproducible from this driver too.
    {"dFCPAdNdN[0][0]", "ComponentSrkCPA.dFCPAdNdN(0, phase, n, T, P)", "NOT d_ln_phi_dn[0][0]",
        "see association.rs's note on calc_lngij"},
    {"dFCPAdNdN[0][1]", "ComponentSrkCPA.dFCPAdNdN(1, phase, n, T, P)", "NOT d_ln_phi_dn[0][1]", ""},
    {"dFCPAdNdN[1][0]", "ComponentSrkCPA.dFCPAdNdN(0, phase, n, T, P)", "NOT d_ln_phi_dn[1][0]", ""},
    {"dFCPAdNdN[1][1]", "ComponentSrkCPA.dFCPAdNdN(1, phase, n, T, P)", "NOT d_ln_phi_dn[1][1]", ""},
    // --- the fugacity ----------------------------------------------------------------
    {"lnPhi[0]", "log(ComponentSrkCPA.getFugacityCoefficient())", "PhaseState.ln_phi[0]", "water"},
    {"lnPhi[1]", "log(ComponentSrkCPA.getFugacityCoefficient())", "PhaseState.ln_phi[1]", "methanol"},
    // `calc_lngi` is d ln g / dn_i and is the derivative the fugacity term is weighted by.
    {"lngi[0]", "ComponentSrkCPA.calc_lngi(phase)", "Rdf.d_ln_g_dn[0]", "water"},
    {"lngi[1]", "ComponentSrkCPA.calc_lngi(phase)", "Rdf.d_ln_g_dn[1]", "methanol"},
    // --- the flash -------------------------------------------------------------------
    {"flashPhases", "SystemInterface.getNumberOfPhases() after TPflash()", "-",
        "0 means the flash threw; see error"},
    {"flash_beta[0]", "PhaseInterface.getBeta()", "-", "phase 0 (NeqSim orders it first)"},
    {"flash_Z[0]", "PhaseInterface.getZ()", "-", ""},
    {"flash_v[0]", "PhaseInterface.getMolarVolume()", "-", "m^3/mol"},
    {"flash_x[0][0]", "ComponentInterface.getx()", "-", "water"},
    {"flash_x[0][1]", "ComponentInterface.getx()", "-", "methanol"},
    {"flash_beta[1]", "PhaseInterface.getBeta()", "-", "phase 1, absent when single phase"},
    {"flash_Z[1]", "PhaseInterface.getZ()", "-", ""},
    {"flash_v[1]", "PhaseInterface.getMolarVolume()", "-", "m^3/mol"},
    {"flash_x[1][0]", "ComponentInterface.getx()", "-", "water"},
    {"flash_x[1][1]", "ComponentInterface.getx()", "-", "methanol"},
  };

  private static void header() {
    System.out.println("# azoth CpaSweep - NeqSim 3.20.0's CPA internals over a grid of states.");
    System.out.println("#");
    System.out.println("# The fluid is water/methanol, built as SystemSrkCPA with setMixingRule(10).");
    System.out.println("# Rows are emitted in build order, so the first key that differs is the");
    System.out.println("# layer the divergence is in. Every value is in SI unless the note says");
    System.out.println("# otherwise; the conversion from NeqSim's internal scale is in the `neqsim`");
    System.out.println("# column, because NeqSim carries `a` in Pa*m^6/mol^2 * 1e5 and `b` in");
    System.out.println("# m^3/mol * 1e5 (Component.java:526-531).");
    System.out.println("#");
    System.out.println("# key\tneqsim\tazoth\tnote");
    for (String[] key : KEYS) {
      System.out.println("# " + key[0] + "\t" + key[1] + "\t" + key[2] + "\t" + key[3]);
    }
  }

  private static void put(Map<String, String> row, String key, double value) {
    row.put(key, String.format("%.15g", value));
  }

  /** One state, filled in build order. A failure anywhere records itself and stops. */
  private static Map<String, String> one(double temperature, double pressure, double nWater) {
    Map<String, String> row = new LinkedHashMap<>();
    put(row, "T_K", temperature);
    put(row, "P_bara", pressure);
    put(row, "n_water", nWater);

    SystemInterface system = new SystemSrkCPA(temperature, pressure);
    system.addComponent("water", nWater);
    system.addComponent("methanol", 1.0 - nWater);
    system.setMixingRule(10);
    system.init(0);
    system.init(1);

    int n = system.getNumberOfComponents();
    PhaseInterface phase = system.getPhase(0);

    for (int i = 0; i < n; i++) {
      ComponentEos c = (ComponentEos) phase.getComponent(i);
      // The internal a-scale is 1e5, but sqrtAT is the square root of that scale, so the
      // divisor is sqrt(1e5) - 316.22776..., not 1e5. Stated here because it is exactly
      // the kind of factor this file exists to keep visible.
      put(row, "aT[" + i + "]", c.aT / 1e5);
      put(row, "sqrtAT[" + i + "]", c.sqrtAT / Math.sqrt(1e5));
      put(row, "b[" + i + "]", c.getBi() / 1e5);
      put(row, "Ai[" + i + "]", c.getAi() / 1e5);
      ComponentSrkCPA cpa = (ComponentSrkCPA) c;
      put(row, "calca[" + i + "]", cpa.calca() / 1e5);
      put(row, "calcb[" + i + "]", cpa.calcb() / 1e5);
      row.put("scheme[" + i + "]", cpa.getAssociationScheme());
      put(row, "energy[" + i + "]", cpa.getAssociationEnergy());
      put(row, "volume[" + i + "]", cpa.getAssociationVolume());
      put(row, "lngi[" + i + "]", cpa.calc_lngi(phase));
    }

    put(row, "A_mix", phase.getA() / 1e5);
    put(row, "B_mix", phase.getB() / 1e5);
    put(row, "Z", phase.getZ());
    put(row, "v_molar", phase.getMolarVolume());
    // The two consistency keys. `ComponentEos.fugcoef` is
    // `ln phi = dFdN - ln(P V / (R T))`, so the log's argument is Z only if the volume is
    // in the same scale the root is. Printing both sides makes the scale a datum rather
    // than something a reader has to infer - which is the whole reason this file exists.
    double r = 8.3144621;
    // `phase.getPressure()` and not the driver's argument, because this is the value
    // `ComponentEos.fugcoef` itself evaluates with.
    double p_internal = phase.getPressure();
    put(row, "PV_over_RT", p_internal * phase.getMolarVolume() / (r * temperature));
    put(row, "gcpa", ((PhaseCPAInterface) phase).getGcpa());
    put(row, "gcpav", ((PhaseCPAInterface) phase).getGcpav());
    put(row, "sites", ((PhaseCPAInterface) phase).getTotalNumberOfAccociationSites());

    int site = 0;
    for (int i = 0; i < n; i++) {
      ComponentSrkCPA cpa = (ComponentSrkCPA) phase.getComponent(i);
      for (int j = 0; j < cpa.getNumberOfAssociationSites(); j++) {
        put(row, "xsite[" + i + "][" + j + "]", cpa.getXsite()[j]);
        site++;
      }
    }
    row.put("xsiteTotal", String.valueOf(site));

    put(row, "FCPA", ((PhaseSrkCPA) phase).FCPA());
    put(row, "dFCPAdV", ((PhaseSrkCPA) phase).dFCPAdV());
    put(row, "dFCPAdT", ((PhaseSrkCPA) phase).dFCPAdT());
    for (int i = 0; i < n; i++) {
      ComponentSrkCPA cpa = (ComponentSrkCPA) phase.getComponent(i);
      put(row, "dFCPAdN[" + i + "]", cpa.dFCPAdN(phase, n, temperature, pressure));
      put(row, "lnPhi[" + i + "]", Math.log(cpa.getFugacityCoefficient()));
      for (int j = 0; j < n; j++) {
        put(row, "dFCPAdNdN[" + i + "][" + j + "]",
            cpa.dFCPAdNdN(j, phase, n, temperature, pressure));
      }
    }

    // The flash last, because it mutates the system: everything above describes the
    // single-phase state at the feed, which is what the kernel comparison needs.
    ThermodynamicOperations operations = new ThermodynamicOperations(system);
    operations.TPflash();
    row.put("flashPhases", String.valueOf(system.getNumberOfPhases()));
    for (int i = 0; i < system.getNumberOfPhases() && i < 2; i++) {
      PhaseInterface flashed = system.getPhase(i);
      put(row, "flash_beta[" + i + "]", flashed.getBeta());
      put(row, "flash_Z[" + i + "]", flashed.getZ());
      put(row, "flash_v[" + i + "]", flashed.getMolarVolume());
      for (int j = 0; j < n; j++) {
        put(row, "flash_x[" + i + "][" + j + "]", flashed.getComponent(j).getx());
      }
    }
    return row;
  }

  public static void main(String[] args) {
    header();

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

    for (int index = 0; index < cases.length; index++) {
      double[] c = cases[index];
      if (c[0] < 0) {
        continue;
      }
      Map<String, String> row;
      try {
        row = one(c[0], c[1], c[2]);
      } catch (Throwable thrown) {
        // One state that will not solve must not take the sweep with it, and it must say
        // which one it was.
        row = new LinkedHashMap<>();
        put(row, "T_K", c[0]);
        put(row, "P_bara", c[1]);
        put(row, "n_water", c[2]);
        row.put("error", thrown.getClass().getSimpleName());
      }
      StringBuilder line = new StringBuilder("case ").append(index);
      for (Map.Entry<String, String> entry : row.entrySet()) {
        line.append('\t').append(entry.getKey()).append('=').append(entry.getValue());
      }
      System.out.println(line);
    }
  }
}
