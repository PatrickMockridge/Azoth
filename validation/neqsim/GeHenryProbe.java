// `ComponentGE.fugcoef`'s two reference states, measured component by component.
//
//     javac -proc:none -cp neqsim-f0c7436.jar GeHenryProbe.java
//     java -cp .:neqsim-f0c7436.jar GeHenryProbe > captures/ge_henry_probe.tsv
//
// The four solvent-branch drivers (`GeNrtlPhi`, `GeWilsonPhi`, `GeUnifacPhi`) chose
// solvent-tagged components on purpose: `ComponentGE.fugcoef` gives anything else a
// Henry's-law coefficient instead of `gamma_i P0_i / P`. This probe measures that other
// branch, so a port of it can be held to the same numbers.
//
// **What the branch is.** `referenceStateType.equals("solvent")`, and the IAPWS aqueous
// override not applying, takes `phi_i = gamma_i P0_i / P`. Everything else takes
// `phi_i = (gamma_i / gamma_i^inf) H_i / P`, where `gamma^inf` is the phase's own
// activity model evaluated on a reference: the binary (solute, water) at
// x = [1e-11, 1-1e-11] when the phase carries water, otherwise the phase's composition
// with the solute at x = 1e-10 and the rest renormalised. `H_i` prefers the IAPWS G7-04
// table for its gas if it is usable at T, and otherwise reads the databank's four
// HenryCoef columns through `getHenryCoef`.
//
// **The IAPWS override is what makes a solvent tag insufficient.** `hydrogen`, `argon`
// and `helium` are tagged `solvent` in COMP.csv but are supported species, so in a phase
// carrying water they take the Henry branch - and their reference phase is built as a
// *one-component* phase, because `initRefPhases` gives a solvent-tagged component a pure
// phase. `gamma_inf_water_error` records what that does.
//
// **The UNIFAC reference phase is built without `checkGroups` or `calcaij`,** so its
// components' group arrays are empty and its interaction matrix is unset. The raw value
// is printed as `gamma_inf_water`; `gamma_inf_water_synced` is the same reference phase
// after `setUnifacGroups(getUnifacGroups2())` and `calcaij`, so the two say whether the
// raw path is usable as an oracle or an upstream defect to record.
//
// `henry_branch_bar = phi P / gammaRefCor` is the Henry coefficient the branch actually
// used - and for a solvent-arm component the same expression is its `P0`. The four
// columns `gamma`, `gammaRefCor`, `henry_branch_bar` and `phi` therefore let a reader
// check the arithmetic separately from the correlations.

import java.lang.reflect.Field;

import neqsim.thermo.component.ComponentGE;
import neqsim.thermo.component.ComponentGEInterface;
import neqsim.thermo.component.ComponentGEUnifac;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.component.IapwsHenryLaw;
import neqsim.thermo.phase.Phase;
import neqsim.thermo.phase.PhaseGE;
import neqsim.thermo.phase.PhaseGEInterface;
import neqsim.thermo.phase.PhaseGEUnifac;
import neqsim.thermo.phase.PhaseGEUniquac;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.EosGeFlashModel;
import neqsim.thermo.system.SystemGEWilson;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemNRTL;
import neqsim.thermo.system.SystemUNIFAC;
import neqsim.thermodynamicoperations.flashops.RachfordRice;

public class GeHenryProbe {

  enum Model { WILSON, NRTL, UNIFAC }

  public static void main(String[] args) throws Exception {
    one(Model.WILSON, "water_co2_aqueous", names("water", "CO2"), z(0.98, 0.02), 298.15, 1.0);
    one(Model.WILSON, "water_methane_iapws", names("water", "methane"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.WILSON, "water_nitrogen_iapws", names("water", "nitrogen"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.WILSON, "water_nitrogen_hot", names("water", "nitrogen"), z(0.99, 0.01), 640.0, 100.0);
    one(Model.WILSON, "water_ethylene_db", names("water", "ethylene"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.WILSON, "water_nbutane_zeros", names("water", "n-butane"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.WILSON, "water_argon_onecomponent", names("water", "argon"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.WILSON, "methanol_co2_dry", names("methanol", "CO2"), z(0.98, 0.02), 298.15, 10.0);
    one(Model.WILSON, "methanol_hexane_dry", names("methanol", "n-hexane"), z(0.5, 0.5), 298.15, 1.0);
    one(Model.WILSON, "methanol_ethanol_control", names("methanol", "ethanol"), z(0.5, 0.5), 298.15, 1.0);
    one(Model.WILSON, "water_methanol_co2_ternary", names("water", "methanol", "CO2"),
        z(0.60, 0.35, 0.05), 298.15, 5.0);

    one(Model.NRTL, "water_co2_aqueous", names("water", "CO2"), z(0.98, 0.02), 298.15, 1.0);
    one(Model.NRTL, "water_methane_iapws", names("water", "methane"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.NRTL, "water_nitrogen_iapws", names("water", "nitrogen"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.NRTL, "water_nitrogen_hot", names("water", "nitrogen"), z(0.99, 0.01), 640.0, 100.0);
    one(Model.NRTL, "water_ethylene_db", names("water", "ethylene"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.NRTL, "water_nbutane_zeros", names("water", "n-butane"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.NRTL, "water_argon_onecomponent", names("water", "argon"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.NRTL, "methanol_co2_dry", names("methanol", "CO2"), z(0.98, 0.02), 298.15, 10.0);
    one(Model.NRTL, "methanol_hexane_dry", names("methanol", "n-hexane"), z(0.5, 0.5), 298.15, 1.0);
    one(Model.NRTL, "methanol_ethanol_control", names("methanol", "ethanol"), z(0.5, 0.5), 298.15, 1.0);
    one(Model.NRTL, "water_methanol_co2_ternary", names("water", "methanol", "CO2"),
        z(0.60, 0.35, 0.05), 298.15, 5.0);
    one(Model.NRTL, "water_nacl_co2_ions", names("water", "Na+", "Cl-", "CO2"),
        z(0.951, 0.024, 0.024, 0.001), 298.15, 1.0);

    one(Model.UNIFAC, "water_co2_aqueous", names("water", "CO2"), z(0.98, 0.02), 298.15, 1.0);
    one(Model.UNIFAC, "water_methane_iapws", names("water", "methane"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.UNIFAC, "water_nitrogen_iapws", names("water", "nitrogen"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.UNIFAC, "water_nitrogen_hot", names("water", "nitrogen"), z(0.99, 0.01), 640.0, 100.0);
    one(Model.UNIFAC, "water_ethylene_db", names("water", "ethylene"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.UNIFAC, "water_nbutane_zeros", names("water", "n-butane"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.UNIFAC, "water_argon_onecomponent", names("water", "argon"), z(0.99, 0.01), 298.15, 1.0);
    one(Model.UNIFAC, "methanol_co2_dry", names("methanol", "CO2"), z(0.98, 0.02), 298.15, 10.0);
    one(Model.UNIFAC, "methanol_hexane_dry", names("methanol", "n-hexane"), z(0.5, 0.5), 298.15, 1.0);
    one(Model.UNIFAC, "methanol_ethanol_control", names("methanol", "ethanol"), z(0.5, 0.5), 298.15, 1.0);
    one(Model.UNIFAC, "water_methanol_co2_ternary", names("water", "methanol", "CO2"),
        z(0.60, 0.35, 0.05), 298.15, 5.0);

    // `PhaseGEUniquac` is withdrawn upstream: its constructor throws, which the block
    // records rather than skips.
    uniquacAttempt();

    flash(names("water", "methanol", "CO2"), z(0.55, 0.40, 0.05), 350.0, 1.0);
  }

  static String[] names(String... names) {
    return names;
  }

  static double[] z(double... values) {
    return values;
  }

  static void one(Model model, String id, String[] names, double[] z, double tK, double pBar) {
    System.out.println();
    System.out.println("=== " + id + " " + model);
    SystemInterface system;
    try {
      system = build(model, names, z, tK, pBar);
    } catch (Throwable e) {
      System.out.println("build_error = " + describe(e));
      return;
    }
    System.out.println("T_K = " + fmt(tK));
    System.out.println("P_bar = " + fmt(pBar));
    PhaseInterface ge = findGe(system, model);
    if (ge == null) {
      System.out.println("phase_error = no GE phase was built");
      return;
    }
    if (ge instanceof PhaseGEUnifac) {
      ((PhaseGEUnifac) ge).checkGroups();
    }
    int n = ge.getNumberOfComponents();
    System.out.println("phase_class = " + ge.getClass().getSimpleName());
    System.out.println("n = " + n);
    int water = ge.hasComponent("water")
        ? ge.getComponent("water").getComponentNumber() : -1;
    System.out.println("water_index = " + water);

    try {
      ((PhaseGEInterface) ge).getExcessGibbsEnergy(ge, n, tK, pBar, PhaseType.LIQUID);
    } catch (Throwable e) {
      System.out.println("excess_error = " + describe(e));
    }

    for (int i = 0; i < n; i++) {
      final int index = i;
      ComponentInterface c = ge.getComponent(i);
      String name = c.getName();
      String tag = c.getReferenceStateType();
      boolean iapws = IapwsHenryLaw.isSupportedSpecies(name);
      boolean usable = false;
      try {
        usable = IapwsHenryLaw.isUsable(name, tK);
      } catch (Throwable e) {
        // ``isUsable`` is total; a throw here is itself a finding and printed below.
      }
      boolean henry = !tag.equals("solvent")
          || (!name.equalsIgnoreCase("water") && iapws && water >= 0);

      String p = "component[" + i + "].";
      System.out.println(p + "name = " + name);
      System.out.println(p + "tag = " + tag);
      System.out.println(p + "ion = " + c.isIsIon());
      System.out.println(p + "branch = " + (henry ? "henry" : "solvent"));
      System.out.println(p + "iapws_supported = " + iapws);
      System.out.println(p + "iapws_usable = " + usable);
      System.out.println(p + "x = " + fmt(c.getx()));

      double gamma = Double.NaN;
      double ref = Double.NaN;
      double phi = Double.NaN;
      try {
        // `fugcoef` reads the cached `gamma`; the excess-Gibbs call above is what fills it.
        c.fugcoef(ge);
        gamma = ((ComponentGEInterface) c).getGamma();
        ref = ((ComponentGE) c).getGammaRefCor();
        phi = c.getFugacityCoefficient();
      } catch (Throwable e) {
        System.out.println(p + "fugcoef_error = " + describe(e));
      }
      System.out.println(p + "gamma = " + fmt(gamma));
      System.out.println(p + "gamma_ref_cor = " + fmt(ref));
      System.out.println(p + "activinf = " + fmt(ref == 0.0 ? Double.NaN : gamma / ref));
      System.out.println(p + "p0_bar = " + fmt(tryDouble(() -> c.getAntoineVaporPressure(tK))));
      System.out.println(p + "phi = " + fmt(phi));
      System.out.println(p + "ln_phi = " + fmt(phi > 0.0 ? Math.log(phi) : Double.NaN));
      // The coefficient the branch actually used: `H` on the Henry arm, `P0` on the
      // solvent arm. Derived from the reported values rather than read from a protected
      // method, so the check is the identity the port must reproduce.
      System.out.println(p + "coefficient_branch_bar = "
          + fmt(ref > 0.0 && Double.isFinite(phi) ? phi * pBar / ref : Double.NaN));
      System.out.println(p + "henry_db_bar = " + fmt(tryDouble(() -> c.getHenryCoef(tK))));
      if (iapws && usable) {
        System.out.println(p + "henry_iapws_bar = "
            + fmt(tryDouble(() -> IapwsHenryLaw.getHenryCoefficientBar(name, tK))));
      }

      if (henry && water >= 0) {
        System.out.println(p + "gamma_inf_water = "
            + attempt(() -> ((PhaseGE) ge).getActivityCoefficientInfDilWater(index, water)));
        if (model == Model.UNIFAC) {
          System.out.println(p + "gamma_inf_water_synced = "
              + attempt(() -> syncedReferenceGamma(ge, index, water, tK, pBar)));
        }
      }
      System.out.println(p + "gamma_inf_clone = "
          + attempt(() -> ((PhaseGE) ge).getActivityCoefficientInfDil(index)));
    }
  }

  static SystemInterface build(Model model, String[] names, double[] z, double tK, double pBar) {
    SystemInterface system;
    switch (model) {
      case WILSON:
        system = new SystemGEWilson(tK, pBar);
        break;
      case NRTL:
        system = new SystemNRTL(tK, pBar);
        break;
      default:
        system = new SystemUNIFAC(tK, pBar);
        break;
    }
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    // `ComponentGEUnifac`'s constructor fills each component's group *list*; the array
    // `getUnifacGroup(i)` reads is only rebuilt by `setUnifacGroups`. `GeUnifacPhi.java`
    // documents the sequence and why each `init(0)` needs it again.
    system.createDatabase(true);
    if (model == Model.UNIFAC) {
      synchronise(system);
    }
    system.init(0);
    if (model == Model.UNIFAC) {
      synchronise(system);
    }
    system.setMixingRule("classic");
    if (model == Model.UNIFAC) {
      synchronise(system);
    }
    system.init(0);
    if (model == Model.UNIFAC) {
      synchronise(system);
    }
    return system;
  }

  static void synchronise(SystemInterface system) {
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      if (!(system.getPhase(p) instanceof PhaseGEUnifac)) {
        continue;
      }
      PhaseInterface ge = system.getPhase(p);
      for (int i = 0; i < ge.getNumberOfComponents(); i++) {
        ComponentGEUnifac component = (ComponentGEUnifac) ge.getComponent(i);
        component.setUnifacGroups(component.getUnifacGroups2());
      }
    }
  }

  static PhaseInterface findGe(SystemInterface system, Model model) {
    PhaseInterface ge = null;
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      if (model == Model.UNIFAC) {
        if (system.getPhase(p) instanceof PhaseGEUnifac) {
          ge = system.getPhase(p);
        }
      } else if (system.getPhase(p) instanceof PhaseGEInterface) {
        ge = system.getPhase(p);
      }
    }
    return ge;
  }

  /**
   * The UNIFAC reference phase as the model intends it: groups synchronised and the
   * interaction matrix built, which `initRefPhases` does neither of.
   *
   * The raw value is `getActivityCoefficientInfDilWater`'s own output; this is the
   * repair, so the capture says whether the raw path is an oracle or a defect.
   */
  static double syncedReferenceGamma(PhaseInterface ge, int k, int water, double tK, double pBar)
      throws Exception {
    ((Phase) ge).initRefPhases(false, ge.getComponent(water).getName());
    Field field = Phase.class.getDeclaredField("refPhase");
    field.setAccessible(true);
    PhaseInterface[] reference = (PhaseInterface[]) field.get(ge);
    if (reference == null || reference.length <= k || reference[k] == null) {
      throw new IllegalStateException("no reference phase for component " + k);
    }
    PhaseInterface pair = reference[k];
    if (pair.getNumberOfComponents() < 2) {
      throw new IllegalStateException("the reference phase has "
          + pair.getNumberOfComponents() + " component(s), not the solute-water pair");
    }
    for (int i = 0; i < pair.getNumberOfComponents(); i++) {
      ComponentGEUnifac component = (ComponentGEUnifac) pair.getComponent(i);
      component.setUnifacGroups(component.getUnifacGroups2());
    }
    if (pair instanceof PhaseGEUnifac) {
      ((PhaseGEUnifac) pair).calcaij();
    }
    pair.setTemperature(tK);
    pair.setPressure(pBar);
    pair.init(pair.getNumberOfMolesInPhase(), pair.getNumberOfComponents(), 1, ge.getType(), 1.0);
    ((PhaseGEInterface) pair).getExcessGibbsEnergy(pair, pair.getNumberOfComponents(),
        pair.getTemperature(), pair.getPressure(), pair.getType());
    return ((ComponentGEInterface) pair.getComponent(0)).getGamma();
  }

  /** A bare UNIQUAC phase: upstream's constructors throw, which the capture records. */
  static void uniquacAttempt() {
    System.out.println();
    System.out.println("=== uniquac_attempt UNIQUAC");
    try {
      PhaseGEUniquac phase = new PhaseGEUniquac();
      phase.setTemperature(298.15);
      phase.setPressure(1.0);
      phase.addComponent("methanol", 0.5, 0.5, 0);
      phase.addComponent("water", 0.5, 0.5, 1);
      System.out.println("constructed = true");
      ((PhaseGEInterface) phase).getExcessGibbsEnergy(phase, 2, 298.15, 1.0, PhaseType.LIQUID);
      for (int i = 0; i < phase.getNumberOfComponents(); i++) {
        ComponentInterface c = phase.getComponent(i);
        c.fugcoef(phase);
        System.out.println("component[" + i + "].name = " + c.getName());
        System.out.println("component[" + i + "].gamma = "
            + fmt(((ComponentGEInterface) c).getGamma()));
        System.out.println("component[" + i + "].phi = " + fmt(c.getFugacityCoefficient()));
      }
    } catch (Throwable e) {
      System.out.println("build_error = " + describe(e));
    }
  }

  /**
   * The gamma-phi flash of an SRK vapour over a `PhaseGENRTL` liquid, for a solute-bearing
   * feed.
   *
   * The loop is `NrtlGammaPhiFlash.java`'s, transcribed from `TPflash.sucsSubsGammaPhi`,
   * because `SystemNRTL` does not opt into NeqSim's direct gamma-phi route: its own
   * `TPflash` runs ordinary EOS successive substitution and never activates the GE liquid.
   * What is NeqSim's is the model and the update rule; what is printed is their fixed point.
   */
  static void flash(String[] names, double[] z, double tK, double pBar) throws Exception {
    System.out.println();
    System.out.println("=== water_methanol_co2_flash NRTL");
    SystemNRTL system = new SystemNRTL(tK, pBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.createDatabase(true);
    system.setMixingRule("classic");
    system.init(0);

    PhaseInterface vapour = system.getPhase(0);
    PhaseInterface liquid = system.getPhase(1);
    EosGeFlashModel model = (EosGeFlashModel) system;
    RachfordRice rachfordRice = new RachfordRice();

    int n = vapour.getNumberOfComponents();
    double[] k = new double[n];
    for (int i = 0; i < n; i++) {
      k[i] = Math.exp(Math.log(vapour.getComponent(i).getPC() / pBar)
          + 5.373 * (1.0 + vapour.getComponent(i).getAcentricFactor())
          * (1.0 - vapour.getComponent(i).getTC() / tK));
      vapour.getComponent(i).setK(k[i]);
      liquid.getComponent(i).setK(k[i]);
    }

    double beta = 0.5;
    double residual = Double.NaN;
    int iteration = 0;
    // The vapour's coefficients are a cached field; `init(0)` leaves them unset for this
    // three-component state, and the first K update would divide by `NaN`.
    system.setBeta(beta);
    system.calc_x_y();
    system.init(1);
    for (iteration = 1; iteration <= 500; iteration++) {
      ((PhaseGEInterface) liquid).getExcessGibbsEnergy(liquid, n, tK, pBar, PhaseType.LIQUID);
      residual = 0.0;
      for (int i = 0; i < n; i++) {
        ComponentInterface inVapour = vapour.getComponent(i);
        ComponentInterface inLiquid = liquid.getComponent(i);
        double kOld = inVapour.getK();
        inLiquid.fugcoef(liquid);
        double phiVapour = model.getGammaPhiVapourFugacityCoefficient(inVapour, vapour);
        double targetK = inLiquid.getFugacityCoefficient() / phiVapour;
        double kNew = model.relaxGammaPhiKValue(kOld, targetK);
        inVapour.setK(kNew);
        inLiquid.setK(kNew);
        k[i] = kNew;
        residual += Math.pow(Math.log(kNew / kOld), 2.0);
      }
      residual = Math.sqrt(residual / n);
      try {
        beta = rachfordRice.calcBeta(k, system.getzvector());
      } catch (Throwable e) {
        // A feed outside the two-phase window has no Rachford-Rice root, and NeqSim
        // throws rather than returning one; the capture records which state that was.
        System.out.println("iteration = " + iteration);
        System.out.println("k_at_error = " + vector(k));
        System.out.println("flash_error = " + describe(e));
        return;
      }
      system.setBeta(beta);
      system.calc_x_y();
      system.init(1);
      if (residual <= 1.0e-14) {
        break;
      }
    }

    for (int i = 0; i < n; i++) {
      vapour.getComponent(i).setK(k[i]);
      liquid.getComponent(i).setK(k[i]);
    }
    system.setBeta(beta);
    system.calc_x_y();
    system.init(1);
    ((PhaseGEInterface) liquid).getExcessGibbsEnergy(liquid, n, tK, pBar, PhaseType.LIQUID);

    System.out.println("T_K = " + fmt(tK));
    System.out.println("P_bar = " + fmt(pBar));
    System.out.println("z = " + vector(z));
    System.out.println("iterations = " + iteration);
    System.out.println("residual = " + fmt(residual));
    System.out.println("beta = " + fmt(beta));
    for (int i = 0; i < n; i++) {
      ComponentInterface inVapour = vapour.getComponent(i);
      ComponentInterface inLiquid = liquid.getComponent(i);
      inLiquid.fugcoef(liquid);
      double phiVapour = model.getGammaPhiVapourFugacityCoefficient(inVapour, vapour);
      String p = "component[" + i + "].";
      System.out.println(p + "name = " + inVapour.getName());
      System.out.println(p + "x = " + fmt(inLiquid.getx()));
      System.out.println(p + "y = " + fmt(inVapour.getx()));
      System.out.println(p + "K = " + fmt(inVapour.getK()));
      System.out.println(p + "gamma = " + fmt(((ComponentGEInterface) inLiquid).getGamma()));
      System.out.println(p + "gamma_ref_cor = " + fmt(((ComponentGE) inLiquid).getGammaRefCor()));
      System.out.println(p + "ln_phi_liquid = "
          + fmt(Math.log(inLiquid.getFugacityCoefficient())));
      System.out.println(p + "ln_phi_vapour = " + fmt(Math.log(phiVapour)));
    }
  }

  static String fmt(double value) {
    return Double.toString(value);
  }

  static String vector(double[] values) {
    StringBuilder out = new StringBuilder("[");
    for (int i = 0; i < values.length; i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(fmt(values[i]));
    }
    return out.append("]").toString();
  }

  static String describe(Throwable e) {
    return e.getClass().getName() + ": " + String.valueOf(e.getMessage());
  }

  interface DoubleCall {
    double call() throws Exception;
  }

  static double tryDouble(DoubleCall call) {
    try {
      return call.call();
    } catch (Throwable e) {
      return Double.NaN;
    }
  }

  static String attempt(DoubleCall call) {
    try {
      return fmt(call.call());
    } catch (Throwable e) {
      return "error: " + describe(e);
    }
  }
}
