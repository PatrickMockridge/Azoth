import java.util.LinkedHashMap;
import java.util.Map;

import neqsim.thermo.component.ComponentEos;
import neqsim.thermo.component.ComponentSrkCPA;
import neqsim.thermo.phase.PhaseCPAInterface;
import neqsim.thermo.phase.PhaseEos;
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
 * <b>What the finite-difference keys found.</b> At 356 K and 1 bara, water/methanol
 * 0.6/0.4: {@code F_scale_two} is exactly {@code 2 * F_res_over_R}, so {@code getF()} is
 * extensive and {@code dF = sum_i dFdN_i dn_i + FV dV} applies to it. The two isolated
 * differences then give {@code dFdN_fd = [-0.0419683, -0.0743505]} against NeqSim's own
 * {@code dFdN = [-0.0901941, -0.1211108]} - factors that differ per component, so not a
 * scale slip. The weighted sum {@code sum_i n_i dFdN_fd_i = -0.0549213} agrees to four
 * digits with the scale derivative {@code F - FV*V = -0.0549212}, which is derived without
 * any of those differences, so the finite differences and {@code FV()} are consistent with
 * each other and with extensivity.
 *
 * <p>
 * So <b>{@code ComponentEos.dFdN} is not the derivative of {@code ComponentEos.getF}</b>,
 * while {@code ComponentSrkCPA.dFCPAdN} is - the association's own derivative is a true
 * one. {@code getF()}'s form is verified, not read: {@code F_reconstructed} reproduces it
 * to every printed digit on every row.
 *
 * <p>
 * <b>What that does not yet settle</b> is whether the same is true for a plain cubic, and
 * it matters: if it is, NeqSim's ordinary cubic {@code ln phi} is not {@code integral dFdN}
 * either, and its cubic flashes still agree with azoth's - so the error would have to
 * cancel in {@code ln phi_L - ln phi_V}. That is a {@code SystemPrEos} sweep, and it is the
 * next measurement rather than a conclusion.
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
    {"FCPA", "PhaseSrkCPA.FCPA()", "SiteState.helmholtz_rt",
        "A_assoc / (R T), from the cached field"},
    // **`FCPA()` is a cached field, and the sum its own commented-out body gives is not what
    // fills it.** `PhaseSrkCPA` computes `FCPA = m^T (u - 0.5 ksi .* udot)` with
    // `u_i = ln x_i - x_i + 1` and `udot_i = S_i`, which equals
    // `sum_i n_i sum_A (ln x_A - x_A/2 + 1/2)` **only where the solve has converged**, since
    // that needs `x_i (1 + S_i) = 1`. `FCPA_sum` is the commented-out body from
    // `getXsite()`, so the two keys say which expression the cached field holds.
    {"FCPA_sum", "the commented-out body of PhaseSrkCPA.FCPA(), from getXsite()",
        "SiteState.helmholtz_rt", "what azoth ports"},
    // `hCPA = sum_i n_i sum_A (1 - X_A)`, the unbonded site count. `FCPA` depends on it through
    // the `-x/2 + 1/2` terms, and it is carried on the solve rather than derived from a
    // fugacity derivative, so it is the one quantity `FCPA` turns on that comparing
    // `dFCPAdN` does not pin.
    {"hcpa", "PhaseCPAInterface.getHcpatot()", "SiteState.unbonded_sites", ""},
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
    // **How NeqSim assembles the cubic `ln phi` at all**, which is not the textbook
    // expression azoth writes: `ComponentEos.fugcoef` is `dFdN_i - ln Z` and
    // `ComponentEos.dFdN` is `Fn() + FB()*getBi() + FD()*getAi()`, so the cubic part is
    // these three phase functions contracted with the component's own `Bi` and `Ai`. A
    // disagreement in `lnPhi` that the root, `A`, `B` and `kij` do not explain has to be
    // in one of them, and until this row existed there was no way to see which.
    // The raw building blocks of `getF()`, so that the identity
    // `getF() == -n*getg() - (getA()/T)*getf_loc()` can be **checked numerically** rather
    // than assumed. Every scale question about `Fn`/`FB`/`FD` reduces to this one line:
    // if it holds with the printed `getA()`, then `A` in the Helmholtz energy is in the
    // same internal scale `getA()` reports and the contraction in `dFdN` follows.
    // **`PhaseSrkCPA.getF() = -n*getg() - (getA()/T)*getf_loc() + FCPA`, verified
    // numerically to seven digits at 356 K / 1 bar.** The `F_reconstructed` key below is
    // that same expression evaluated independently on every row, so the decomposition is
    // checked rather than asserted, and a scale error in `getA()` would show as a
    // difference between the two columns rather than as a wrong answer somewhere later.
    {"F_res_over_R", "PhaseEos.getF(), overridden by PhaseSrkCPA", "no single azoth counterpart",
        "the residual Helmholtz over R; includes the association"},
    {"FV", "PhaseEos.FV()", "no single azoth counterpart", "dF/dV at constant n and T"},
    {"V_total", "PhaseInterface.getTotalVolume()", "z R T / P", "NeqSim's internal volume scale"},
    {"F_reconstructed", "-n*getg() - (getA()/T)*getf_loc() + FCPA", "-",
        "must equal F_res_over_R; the check on getA()'s scale"},
    // **The one key that settles whether NeqSim's `dFdN` is the derivative of its own `F`.**
    //
    // `F` is a function of `(T, V, n)`, so along *any* path
    // `dF = sum_i (dF/dn_i)|_{T,V} dn_i + (dF/dV)|_{T,n} dV`. Perturbing ONE mole number at
    // constant `T` and `P` therefore isolates `(dF/dn_i)|_{T,V}`:
    //
    //     dFdN_fd[i] = (F(n + d e_i) - F(n)) / d  -  FV() * (V(n + d e_i) - V(n)) / d
    //
    // with `V` the *total* volume, because `F` is extensive. No prescribed-volume state is
    // needed, which is what makes this measurable from a driver that only knows how to run
    // `init(1)` at a temperature and a pressure.
    //
    // `dFdN[i]` is NeqSim's own `ComponentEos.dFdN`. If the two agree, its assembly
    // `Fn + FB*getBi() + FD*getAi()` is a true derivative and the residual Helmholtz
    // *functions* differ between it and azoth; if they disagree, `dFdN` is not the
    // derivative of its own `F` and the finding is upstream, of the `calc_lngij` kind.
    {"dFdN_fd[0]", "the identity above, at a perturbation of 1e-7 in n_water", "no counterpart",
        "= dFdN[0] iff NeqSim's dFdN is the derivative of its own F"},
    {"dFdN_fd[1]", "the same, perturbing n_methanol", "no counterpart", "see above"},
    // The finite difference's raw ingredients, so the arithmetic can be checked rather than
    // trusted. `F_scale_two` is also the extensivity test: doubling every mole number at
    // the same T and P doubles an extensive `F` and leaves a molar one alone.
    {"F_nwater_up", "getF() at n_water + 1e-7", "-", ""},
    {"F_nwater_down", "getF() at n_water - 1e-7", "-", ""},
    {"V_nwater_up", "getTotalVolume() at n_water + 1e-7", "-", "internal scale"},
    {"V_nwater_down", "getTotalVolume() at n_water - 1e-7", "-", "internal scale"},
    {"F_methanol_up", "getF() at n_methanol + 1e-7", "-", ""},
    {"F_methanol_down", "getF() at n_methanol - 1e-7", "-", ""},
    {"F_scale_two", "getF() with every mole number doubled", "no counterpart",
        "= 2 * F_res_over_R iff F is extensive"},
    {"g_helmholtz", "PhaseEos.getg()", "no single azoth counterpart", "the (Z - B) term"},
    {"f_loc", "PhaseEos.getf_loc()", "no single azoth counterpart", "the (Z + d1 B)/(Z + d2 B) term"},
    {"n_moles", "PhaseInterface.getNumberOfMolesInPhase()", "1.0", ""},
    {"b_mixture", "PhaseEos.getb()", "b_mix * 1e5", "internal scale"},
    {"delta1", "PhaseEos.delta1", "Cubic::delta1()", ""},
    {"delta2", "PhaseEos.delta2", "Cubic::delta2()", ""},
    {"Fn", "PhaseEos.Fn()", "no single azoth counterpart", "d(F/RT)/dn, the total-moles term"},
    {"FB", "PhaseEos.FB()", "no single azoth counterpart", "d(F/RT)/dB"},
    {"FD", "PhaseEos.FD()", "no single azoth counterpart", "d(F/RT)/dA"},
    {"dFdN[0]", "ComponentEos.dFdN(phase, n, T, P)", "cubic ln_phi[0] + ln Z", "water"},
    {"dFdN[1]", "ComponentEos.dFdN(phase, n, T, P)", "cubic ln_phi[1] + ln Z", "methanol"},
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
    // The flashed phases' *association* values, which is the only route to NeqSim's liquid
    // branch: `init(1)` leaves whichever single phase the volume solve landed on - the
    // vapour, at low pressure - so the liquid root's `FCPA`, `dFCPAdN` and `gcpa` were not
    // reachable before. A `TPflash` that splits puts a real liquid composition in phase 1.
    {"flash_gcpa[0]", "PhaseCPAInterface.getGcpa() on the flashed phase", "Rdf.g", ""},
    {"flash_FCPA[0]", "PhaseSrkCPA.FCPA() on the flashed phase", "SiteState.helmholtz_rt", ""},
    {"flash_FCPA_sum[0]", "the commented-out body, from getXsite() on the flashed phase", "-",
        "= flash_FCPA[0] iff the cached field is current"},
    {"flash_hcpa[0]", "PhaseCPAInterface.getHcpatot() on the flashed phase", "-", ""},
    {"flash_hcpa_direct[0]", "PhaseCPAInterface.calc_hCPA() run on the flashed phase", "-",
        "= flash_hcpa[0] iff the cached field is current"},
    {"flash_hcpa_from_xsite[0]", "calc_hCPA's own body, re-run here from getXsite()", "-",
        "= flash_hcpa_direct[0] is an identity, not a claim"},
    {"flash_xsite[1][0]", "ComponentSrkCPA.getXsite() on the flashed liquid", "SiteState.fractions", ""},
    {"flash_dFCPAdN[0][0]", "ComponentSrkCPA.dFCPAdN on the flashed phase",
        "SiteState.ln_phi[0]", "the liquid's, when the flash splits"},
    {"flash_dFCPAdN[0][1]", "ComponentSrkCPA.dFCPAdN on the flashed phase",
        "SiteState.ln_phi[1]", ""},
    {"flash_gcpa[1]", "PhaseCPAInterface.getGcpa() on phase 1", "Rdf.g", ""},
    {"flash_FCPA[1]", "PhaseSrkCPA.FCPA() on phase 1", "SiteState.helmholtz_rt", ""},
    {"flash_dFCPAdN[1][0]", "ComponentSrkCPA.dFCPAdN on phase 1", "SiteState.ln_phi[0]", ""},
    {"flash_dFCPAdN[1][1]", "ComponentSrkCPA.dFCPAdN on phase 1", "SiteState.ln_phi[1]", ""},
    {"flash_beta[1]", "PhaseInterface.getBeta()", "-", "phase 1, absent when single phase"},
    {"flash_Z[1]", "PhaseInterface.getZ()", "-", ""},
    {"flash_v[1]", "PhaseInterface.getMolarVolume()", "-", "m^3/mol"},
    {"flash_x[1][0]", "ComponentInterface.getx()", "-", "water"},
    {"flash_x[1][1]", "ComponentInterface.getx()", "-", "methanol"},
  };

  /** `sum_i n_i sum_A (ln x_A - x_A/2 + 1/2)` from `getXsite()` - the commented-out body. */
  private static double fcpa_sum(PhaseInterface phase, int n) {
    double total = 0.0;
    for (int i = 0; i < n; i++) {
      ComponentSrkCPA c = (ComponentSrkCPA) phase.getComponent(i);
      double perSite = 0.0;
      for (int j = 0; j < c.getNumberOfAssociationSites(); j++) {
        double x = c.getXsite()[j];
        perSite += Math.log(x) - x / 2.0 + 0.5;
      }
      total += c.getNumberOfMolesInPhase() * perSite;
    }
    return total;
  }

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

  /**
   * One state, filled in build order, at an explicit pair of mole numbers.
   *
   * Split from {@link #one} because the finite-difference keys need to perturb one mole
   * number while holding the other, which a `(T, P, x_water)` signature cannot express.
   * A failure anywhere records itself and stops.
   */
  private static Map<String, String> oneRaw(
      double temperature, double pressure, double nWater, double nMethanol) {
    Map<String, String> row = new LinkedHashMap<>();
    put(row, "T_K", temperature);
    put(row, "P_bara", pressure);
    put(row, "n_water", nWater);
    put(row, "n_methanol", nMethanol);

    SystemInterface system = new SystemSrkCPA(temperature, pressure);
    system.addComponent("water", nWater);
    system.addComponent("methanol", nMethanol);
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

    // The raw building blocks, then the three phase functions `ComponentEos.dFdN`
    // contracts with `Bi` and `Ai`.
    put(row, "F_res_over_R", ((PhaseEos) phase).getF());
    put(row, "FCPA_sum", fcpa_sum(phase, n));
    put(row, "hcpa", ((PhaseCPAInterface) phase).getHcpatot());
    put(row, "FV", phase.FV());
    put(row, "V_total", phase.getTotalVolume());
    put(row, "g_helmholtz", ((PhaseEos) phase).getg());
    put(row, "f_loc", ((PhaseEos) phase).getf_loc());
    put(row, "n_moles", phase.getNumberOfMolesInPhase());
    put(row, "b_mixture", ((PhaseEos) phase).getb(phase, temperature, pressure, n));
    put(row, "delta1", ((PhaseEos) phase).delta1);
    put(row, "delta2", ((PhaseEos) phase).delta2);
    put(row, "Fn", phase.Fn());
    put(row, "F_reconstructed",
        -phase.getNumberOfMolesInPhase() * ((PhaseEos) phase).getg()
            - phase.getA() / temperature * ((PhaseEos) phase).getf_loc()
            + ((PhaseSrkCPA) phase).FCPA());
    put(row, "FB", phase.FB());
    put(row, "FD", phase.FD());

    put(row, "FCPA", ((PhaseSrkCPA) phase).FCPA());
    put(row, "dFCPAdV", ((PhaseSrkCPA) phase).dFCPAdV());
    put(row, "dFCPAdT", ((PhaseSrkCPA) phase).dFCPAdT());
    for (int i = 0; i < n; i++) {
      ComponentSrkCPA cpa = (ComponentSrkCPA) phase.getComponent(i);
      put(row, "dFCPAdN[" + i + "]", cpa.dFCPAdN(phase, n, temperature, pressure));
      put(row, "lnPhi[" + i + "]", Math.log(cpa.getFugacityCoefficient()));
      put(row, "dFdN[" + i + "]", cpa.dFdN(phase, n, temperature, pressure));
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
      put(row, "flash_gcpa[" + i + "]", ((PhaseCPAInterface) flashed).getGcpa());
      put(row, "flash_FCPA[" + i + "]", ((PhaseSrkCPA) flashed).FCPA());
      put(row, "flash_FCPA_sum[" + i + "]", fcpa_sum(flashed, n));
      put(row, "flash_hcpa[" + i + "]", ((PhaseCPAInterface) flashed).getHcpatot());
      put(row, "flash_hcpa_direct[" + i + "]", ((PhaseCPAInterface) flashed).calc_hCPA());
      double fromXsite = 0.0;
      for (int j = 0; j < n; j++) {
        ComponentSrkCPA c = (ComponentSrkCPA) flashed.getComponent(j);
        double perSite = 0.0;
        for (int k = 0; k < c.getNumberOfAssociationSites(); k++) {
          perSite += 1.0 - c.getXsite()[k];
        }
        fromXsite += c.getNumberOfMolesInPhase() * perSite;
      }
      put(row, "flash_hcpa_from_xsite[" + i + "]", fromXsite);
      int flashSite = 0;
      for (int j = 0; j < n; j++) {
        ComponentSrkCPA c = (ComponentSrkCPA) flashed.getComponent(j);
        for (int k = 0; k < c.getNumberOfAssociationSites(); k++) {
          put(row, "flash_xsite[" + i + "][" + flashSite + "]", c.getXsite()[k]);
          flashSite++;
        }
      }
      for (int j = 0; j < n; j++) {
        put(row, "flash_dFCPAdN[" + i + "][" + j + "]",
            ((ComponentSrkCPA) flashed.getComponent(j))
                .dFCPAdN(flashed, n, temperature, pressure));
      }
      for (int j = 0; j < n; j++) {
        put(row, "flash_x[" + i + "][" + j + "]", flashed.getComponent(j).getx());
      }
    }
    return row;
  }

  /**
   * The state, plus the two finite-difference keys that say whether NeqSim's `dFdN` is the
   * derivative of its own `F`.
   *
   * `F` is a function of `(T, V, n)`, so along any path `dF = sum_i dFdN_i dn_i + FV dV`.
   * Perturbing one mole number at constant `T` and `P` therefore isolates `dFdN_i` once
   * the volume's share is subtracted - and no prescribed-volume state is needed, which is
   * what makes this measurable from `init(1)` alone.
   */
  private static Map<String, String> one(double temperature, double pressure, double nWater) {
    Map<String, String> row = oneRaw(temperature, pressure, nWater, 1.0 - nWater);
    double nMethanol = 1.0 - nWater;
    double d = 1.0e-7;
    double[] base = {
      Double.parseDouble(row.get("F_res_over_R")),
      0.0,
    };
    base[1] = Double.parseDouble(row.get("FV"));
    // Component 0: `n_water` moves, methanol held.
    double fUp = Double.parseDouble(oneRaw(temperature, pressure, nWater + d, nMethanol).get("F_res_over_R"));
    double fDown = Double.parseDouble(oneRaw(temperature, pressure, nWater - d, nMethanol).get("F_res_over_R"));
    double vUp = Double.parseDouble(oneRaw(temperature, pressure, nWater + d, nMethanol).get("V_total"));
    double vDown = Double.parseDouble(oneRaw(temperature, pressure, nWater - d, nMethanol).get("V_total"));
    put(row, "dFdN_fd[0]", (fUp - fDown) / (2.0 * d) - base[1] * (vUp - vDown) / (2.0 * d));
    // Component 1: methanol moves, water held.
    fUp = Double.parseDouble(oneRaw(temperature, pressure, nWater, nMethanol + d).get("F_res_over_R"));
    fDown = Double.parseDouble(oneRaw(temperature, pressure, nWater, nMethanol - d).get("F_res_over_R"));
    vUp = Double.parseDouble(oneRaw(temperature, pressure, nWater, nMethanol + d).get("V_total"));
    vDown = Double.parseDouble(oneRaw(temperature, pressure, nWater, nMethanol - d).get("V_total"));
    put(row, "dFdN_fd[1]", (fUp - fDown) / (2.0 * d) - base[1] * (vUp - vDown) / (2.0 * d));

    // The ingredients, so the arithmetic above can be checked by hand, and the extensivity
    // test: doubling every mole number at the same T and P doubles an extensive `F` and
    // leaves a molar one alone. Which of the two `getF()` is decides whether the identity
    // the finite difference rests on even applies to it.
    put(row, "F_nwater_up", Double.parseDouble(oneRaw(temperature, pressure, nWater + d, nMethanol).get("F_res_over_R")));
    put(row, "F_nwater_down", Double.parseDouble(oneRaw(temperature, pressure, nWater - d, nMethanol).get("F_res_over_R")));
    put(row, "V_nwater_up", Double.parseDouble(oneRaw(temperature, pressure, nWater + d, nMethanol).get("V_total")));
    put(row, "V_nwater_down", Double.parseDouble(oneRaw(temperature, pressure, nWater - d, nMethanol).get("V_total")));
    put(row, "F_methanol_up", Double.parseDouble(oneRaw(temperature, pressure, nWater, nMethanol + d).get("F_res_over_R")));
    put(row, "F_methanol_down", Double.parseDouble(oneRaw(temperature, pressure, nWater, nMethanol - d).get("F_res_over_R")));
    put(row, "F_scale_two", Double.parseDouble(oneRaw(temperature, pressure, 2.0 * nWater, 2.0 * nMethanol).get("F_res_over_R")));
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
