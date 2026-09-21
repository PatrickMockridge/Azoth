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
 * <b>What the finite-difference keys found, and what they first said.</b> At 356 K and
 * 1 bara, water/methanol 0.6/0.4: {@code F_scale_two} is exactly
 * {@code 2 * F_res_over_R}, so {@code getF()} is extensive and
 * {@code dF = sum_i dFdN_i dn_i + dFdV dV} applies to it; {@code dFdN_fd} then reproduces
 * NeqSim's own {@code dFdN = [-0.0901941, -0.1211108]} to ten digits. It did not, at first.
 * The chain rule's volume term was written {@code FV() * dV}, and for a {@code PhaseSrkCPA}
 * that is the wrong method: {@code PhaseSrkCPA.dFdV} is {@code super.dFdV() + cpaon *
 * dFCPAdV()} while {@code PhaseEos.FV()} is the cubic part alone, and the association's
 * {@code dFCPAdV} is seventeen times the cubic's {@code FV()} at this state. The omitted
 * term is the whole of the difference that was read as a NeqSim defect - and it was read as
 * one for a whole session, because a probe that is wrong about which quantity to subtract
 * reports the other side as wrong.
 *
 * <p>
 * Both keys are now printed side by side, and the {@code euler_*} keys are the same check
 * without a finite difference in it: {@code F} is extensive in {@code (V, n)}, so
 * {@code V * dFdV() + sum_i n_i dF/dn_i|_V = F}, and NeqSim's {@code dFdN} satisfies it.
 *
 * <p>
 * <b>What the flash keys found.</b> At 356 K and 1 bara the flash splits, and the liquid's
 * mixture parameters are not the feed's. <b>NeqSim's {@code getA}, {@code getB} and
 * {@code getAi} are extensive over the phase's moles</b>: at that liquid, {@code beta} =
 * 0.7916164, {@code getB()} is {@code beta} times {@code sum_i x_i b_i} to every printed
 * digit, {@code getA()} is {@code beta^2} times the double sum {@code sum_i sum_j x_i x_j
 * (1 - k_ij) sqrt(a_i a_j)}, and {@code getAi_i} is {@code 2 beta} times the row sum
 * {@code sum_j x_j (1 - k_ij) sqrt(a_i a_j)} - with {@code k_ij = -0.153}, the
 * {@code cpakij_SRK} column, which is what {@code getA()} and the {@code aT} keys alone
 * return, to a relative {@code 5e-15}. The factor of two that was recorded as unexplained
 * is this and nothing else: {@code getAi} is {@code dA/dn_i} of a degree-two {@code A}, so
 * it carries a 2 that {@code getBi = dB/dn_i} does not. At the feed, one mole total, every
 * one of these equals the intensive azoth value - which is why nothing before the flash
 * could tell the two conventions apart.
 *
 * <p>
 * {@code CpaProbe} remains the single-state kernel probe the Rust tests cite; this is the
 * grid driver, and it reproduces every {@code CpaProbe} reading as a one-case run.
 *
 * <p>
 * Usage:
 *
 * <pre>
 * javac -proc:none -cp neqsim-f0c7436.jar CpaSweep.java
 * java -cp .:neqsim-f0c7436.jar CpaSweep [T_K P_bara n_water]...
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
    // Every value below is dimensional SI: the internal scale is undone by the division in
    // the method column, so the counterpart names azoth's quantity and nothing more.
    {"aT[0]", "ComponentEos.aT / 1e5", "reduced.a[0]", "water; SI Pa*m^6/mol^2"},
    {"aT[1]", "ComponentEos.aT / 1e5", "reduced.a[1]", "methanol"},
    {"sqrtAT[0]", "ComponentEos.sqrtAT / sqrt(1e5)", "sqrt(reduced.a[0])", "water"},
    {"sqrtAT[1]", "ComponentEos.sqrtAT / sqrt(1e5)", "sqrt(reduced.a[1])", "methanol"},
    {"b[0]", "ComponentEos.getBi() / 1e5", "reduced.b[0]", "water; m^3/mol"},
    {"b[1]", "ComponentEos.getBi() / 1e5", "reduced.b[1]", "methanol"},
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
    // **`getA()` and `getB()` are EXTENSIVE over the phase's moles**, `getAi()` and
    // `getBi()` their composition derivatives. `ClassicSRK.calcB` is `sum_i n_i b_i` (when
    // `bmixType == 0`) and `ClassicSRK.calcA` is `sum_i sum_j n_i n_j (1 - k_ij) sqrt(a_i
    // a_j)`; `ClassicSRK.calcAi` returns `2 sum_j n_j sqrt(a_i a_j)(1 - k_ij)` and `calcBi`
    // returns `b_i`, so `Ai_i` is `dA/dn_i` and `Bi_i` is `dB/dn_i`. `A` being degree two
    // and `B` degree one is the whole of why `Ai` carries a `2` and a `beta` that `Bi`
    // does not - the asymmetry that read for a session as an unexplained factor of two.
    //
    // azoth carries the intensive (per-mole) parameters. At the feed, one mole total, the
    // two coincide to every digit, which is why nothing before the flash exposed the
    // convention; at a phase carrying `beta` of the feed's moles they are `beta^2 * a_mix`,
    // `beta * b_mix` and `2 * beta * abar_i`.
    {"A_mix", "PhaseSrkEos.getA() / 1e5", "beta^2 * a_mix", "SI Pa*m^6/mol^2; beta = 1 here"},
    {"B_mix", "PhaseSrkEos.getB() / 1e5", "beta * b_mix", "m^3/mol; beta = 1 here"},
    {"Ai[0]", "PhaseSrkEos.getAi() / 1e5", "2 * beta * abar[0]",
        "the ROW SUM, not a pure a; dA/dn_0; beta = 1 here"},
    {"Ai[1]", "PhaseSrkEos.getAi() / 1e5", "2 * beta * abar[1]",
        "the ROW SUM, not a pure a; dA/dn_1; beta = 1 here"},
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
    {"FV", "PhaseEos.FV()", "no single azoth counterpart",
        "dF/dV of the CUBIC part only; see dFdV"},
    {"dFdV", "PhaseInterface.dFdV()", "no single azoth counterpart",
        "the phase's own dF/dV; = FV() + dFCPAdV() for a CPA phase"},
    {"V_total", "PhaseInterface.getTotalVolume()", "z R T / P", "NeqSim's internal volume scale"},
    {"F_reconstructed", "-n*getg() - (getA()/T)*getf_loc() + FCPA", "-",
        "must equal F_res_over_R; the check on getA()'s scale"},
    // **The key that settles whether NeqSim's `dFdN` is the derivative of its own `F`.**
    //
    // `F` is a function of `(T, V, n)`, so along *any* path
    // `dF = sum_i (dF/dn_i)|_{T,V} dn_i + (dF/dV)|_{T,n} dV`. Perturbing ONE mole number at
    // constant `T` and `P` therefore isolates `(dF/dn_i)|_{T,V}`:
    //
    //     dFdN_fd[i] = (F(n + d e_i) - F(n)) / d  -  dFdV() * (V(n + d e_i) - V(n)) / d
    //
    // with `V` the *total* volume, because `F` is extensive. No prescribed-volume state is
    // needed, which is what makes this measurable from a driver that only knows how to run
    // `init(1)` at a temperature and a pressure.
    //
    // **It is `dFdV()` here and not `FV()`, and the first draft of this driver used `FV()`.**
    // `PhaseSrkCPA.dFdV` overrides `PhaseEos.dFdV` - which is `return FV();` - as
    // `super.dFdV() + cpaon * dFCPAdV()`. At 356 K and 1 bara the association's `dFCPAdV`
    // is `1.694e-5` against a cubic `FV()` of `9.75e-7`: seventeen times larger. Subtracting
    // `FV()` where the chain rule needs `dFdV()` leaves the association's share of `dV` in
    // the residual, and it reported `dFdN` as 38-54% away from the derivative of `F` when the
    // two agree to ten digits. Both keys are printed so the ambiguity is a reading.
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
    // `dFdN`'s own three terms, each printed as the *product* `dFdN` forms, and their sum.
    // The chain rule is `dFdN_i = Fn + FB*Bi + FD*Ai` with `Bi = dB/dn_i` and `Ai = dA/dn_i`,
    // so if the form is right `dFdN_sum[i]` equals `dFdN[i]` and the difference from
    // `dFdN_fd[i]` is not in the assembly. A term that is orders out names its own factor.
    {"dFdN_term_n[0]", "PhaseEos.Fn()", "-", "water; the total-moles term, no Bi"},
    {"dFdN_term_B[0]", "PhaseEos.FB() * ComponentEos.getBi()", "-", "water"},
    {"dFdN_term_A[0]", "PhaseEos.FD() * ComponentEos.getAi()", "-", "water"},
    {"dFdN_sum[0]", "the three above, summed", "dFdN[0]", "= dFdN[0] iff the form is right"},
    {"dFdN_term_n[1]", "PhaseEos.Fn()", "-", "methanol; the same for every component"},
    {"dFdN_term_B[1]", "PhaseEos.FB() * ComponentEos.getBi()", "-", "methanol"},
    {"dFdN_term_A[1]", "PhaseEos.FD() * ComponentEos.getAi()", "-", "methanol"},
    {"dFdN_sum[1]", "the three above, summed", "dFdN[1]", "= dFdN[1] iff the form is right"},
    // **Euler's theorem, the check that needs no finite difference.** `F` is extensive in
    // `(V, n)`, so `V*dFdV() + sum_i n_i dF/dn_i|_V = F`. Every term is a NeqSim reading, so
    // the sides either agree or `dFdN` is not `dF/dn_i` - and unlike the finite-difference
    // pair above, this has no second system, no step size and no volume solver tolerance in
    // it. Run once whole and once with the association's own term removed, so a failure
    // would name the cubic rather than the model.
    {"euler_dFdN", "sum_i n_i * ComponentEos.dFdN_i", "euler_scale",
        "= euler_scale iff dFdN is dF/dn_i"},
    {"euler_dFdN_cubic", "the same, with dFCPAdN_i subtracted", "euler_scale_cubic",
        "the cubic's alone"},
    {"euler_scale", "PhaseEos.getF() - getTotalVolume() * PhaseInterface.dFdV()", "-", ""},
    {"euler_scale_cubic", "the same, from `getF() - FCPA` and `FV()`", "-",
        "PhaseSrkCPA overrides getF and dFdV, not FV"},
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
    // The two facts that decide whether a flash can be compared at all: `lnPhi` is what azoth
    // produces, and it is `dFdN - ln Z` here - so a disagreement in it is either `dFdN` (the
    // finding measured by `CpaFdMismatch`) or the root.
    {"flash_lnPhi[1][0]", "log(ComponentSrkCPA.getFugacityCoefficient()) on the flashed liquid",
        "PhaseState.ln_phi[0]", "the comparison this key exists for"},
    // The per-component and per-mixture inputs the cubic `ln phi` is built from, at a flashed
    // phase. Every one of them has only ever been checked at the *feed*, and the cubic part
    // is the whole of what still differs there. This is the pair of rows that exposed the
    // extensivity above: `aT` and `b` are per component and carry no `beta`, while `A`, `B`
    // and `Ai` are built from the phase's mole numbers and carry `beta^2`, `beta`, `beta`.
    {"flash_aT[1][0]", "ComponentEos.aT / 1e5 on the flashed liquid", "reduced.a[0]",
        "water; SI Pa*m^6/mol^2; the fitted aCPA times alpha, no beta"},
    {"flash_b[1][0]", "ComponentEos.getBi() / 1e5 on the flashed liquid", "reduced.b[0]",
        "water; m^3/mol; the fitted bCPA, no beta"},
    {"flash_Ai[1][0]", "PhaseSrkEos.getAi() / 1e5 on the flashed liquid", "2 * beta * abar[0]",
        "the ROW SUM, not a pure a; dA/dn_0"},
    {"flash_A_mix[1]", "PhaseSrkEos.getA() / 1e5 on the flashed liquid", "beta^2 * a_mix",
        "beta = this phase's moles over the feed's"},
    {"flash_B_mix[1]", "PhaseSrkEos.getB() / 1e5 on the flashed liquid", "beta * b_mix", ""},
    {"flash_lnPhi[1][1]", "the same, component 1", "PhaseState.ln_phi[1]", ""},
    {"flash_lnPhi[0][0]", "the same, on the flashed gas", "PhaseState.ln_phi[0]", ""},
    {"flash_lnPhi[0][1]", "the same, component 1", "PhaseState.ln_phi[1]", ""},
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
    System.out.println("# azoth CpaSweep - NeqSim master's CPA internals over a grid of states.");
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
    put(row, "dFdV", phase.dFdV());
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

    // `ComponentEos.dFdN` is `Fn() + FB()*getBi() + FD()*getAi()` (ComponentEos:265), and the
    // three products are printed separately. `Fn`, `FB` and `FD` on their own have been on
    // the record for a while and they settled nothing, because a scale error in a *product*
    // is invisible when only the factors are printed - the same shape as the `getAi`
    // misreading. Which factor is in which scale is a question these four rows answer by
    // arithmetic instead of argument.
    for (int i = 0; i < n; i++) {
      ComponentEos eos = (ComponentEos) phase.getComponent(i);
      put(row, "dFdN_term_n[" + i + "]", phase.Fn());
      put(row, "dFdN_term_B[" + i + "]", phase.FB() * eos.getBi());
      put(row, "dFdN_term_A[" + i + "]", phase.FD() * eos.getAi());
      put(row, "dFdN_sum[" + i + "]",
          phase.Fn() + phase.FB() * eos.getBi() + phase.FD() * eos.getAi());
    }

    // Euler's theorem, which needs no finite difference, no second system and no solver
    // tolerance: `F` is extensive in `(V, n)`, so `V*FV + sum_i n_i dF/dn_i|_V = F`. Every
    // term is a NeqSim reading, so the two sides either agree or `dFdN` is not `dF/dn_i`.
    double euler = 0.0;
    double eulerCubic = 0.0;
    for (int i = 0; i < n; i++) {
      ComponentSrkCPA cpa = (ComponentSrkCPA) phase.getComponent(i);
      double ni = cpa.getNumberOfMolesInPhase();
      double dFdN = cpa.dFdN(phase, n, temperature, pressure);
      euler += ni * dFdN;
      eulerCubic += ni * (dFdN - cpa.dFCPAdN(phase, n, temperature, pressure));
    }
    put(row, "euler_dFdN", euler);
    put(row, "euler_dFdN_cubic", eulerCubic);
    put(row, "euler_scale", ((PhaseEos) phase).getF() - phase.getTotalVolume() * phase.dFdV());
    put(row, "euler_scale_cubic",
        ((PhaseEos) phase).getF() - ((PhaseSrkCPA) phase).FCPA()
            - phase.getTotalVolume() * phase.FV());

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
      for (int j = 0; j < n; j++) {
        put(row, "flash_lnPhi[" + i + "][" + j + "]",
            Math.log(flashed.getComponent(j).getFugacityCoefficient()));
      }
      put(row, "flash_A_mix[" + i + "]", flashed.getA() / 1e5);
      put(row, "flash_B_mix[" + i + "]", flashed.getB() / 1e5);
      for (int j = 0; j < n; j++) {
        ComponentEos ce = (ComponentEos) flashed.getComponent(j);
        put(row, "flash_aT[" + i + "][" + j + "]", ce.aT / 1e5);
        put(row, "flash_b[" + i + "][" + j + "]", ce.getBi() / 1e5);
        put(row, "flash_Ai[" + i + "][" + j + "]", ce.getAi() / 1e5);
      }
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
   * `F` is a function of `(T, V, n)`, so along any path `dF = sum_i dFdN_i dn_i + dFdV dV`.
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
    base[1] = Double.parseDouble(row.get("dFdV"));
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
