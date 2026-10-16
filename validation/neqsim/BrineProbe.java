// The electrolyte-CPA brine: the phase's own layers, and the four operations that gate on it.
//
//     javac -proc:none -cp neqsim-f0c7436.jar BrineProbe.java
//     java -cp .:neqsim-f0c7436.jar BrineProbe
//
// `SystemElectrolyteCPAstatoil` is `SystemFurstElectrolyteEos` with `PhaseElectrolyteCPAstatoil`
// swapped into the phase array, and `PhaseElectrolyteCPA` is `PhaseModifiedFurstElectrolyteEos`
// with the Wertheim association added:
//
//     getF()   = super.getF()   + cpaon * FCPA()
//     dFdV()   = super.dFdV()   + cpaon * dFCPAdV()
//
// and the same for `dFdT`, `dFdTdT`, `dFdTdV`, `dFdVdV`, `dFdVdVdV`. This prints both halves
// separately, so the sum is measured rather than inferred from `getF()`.
//
// Two things the source does not make obvious, which the rows below are here to measure:
//
//   1. **The statoil RDF is not the interface default.** `PhaseCPAInterface.calc_g` is the
//      simplified Carnahan-Starling `(2-eta)/(2(1-eta)^3)`; `PhaseElectrolyteCPAstatoil`
//      overrides it with `1/(1 - 1.9*B/(4V))`. `gcpa`, `gcpav`, `gcpavv` and `gcpavvv` are
//      printed for both, because a port that reads only one would match at `eta -> 0`.
//   2. **`cpaon` is a property of the phase's site count**, not of the model: it is zero when
//      `getTotalNumberOfAccociationSites() == 0`. `n_sites` is printed beside it so a run
//      that cannot exercise the off case says so rather than looking like it did.
//   3. **`ader` is the Huron-Vidal quantity, and this fluid does not run Huron-Vidal.**
//      `SystemElectrolyteCPAstatoil` is given `setMixingRule(10)`, which the javadoc of
//      `SystemInterface.setMixingRule(int)` defines as "classic mixing rule with temperature
//      and composition dependent kij of CPA" - not rule 4, which is Huron-Vidal and which
//      `SystemFurstElectrolyteEos`'s own test uses. Every `ader` row below is therefore
//      `0.0`, and that is a measurement rather than a gap: a classic-rule phase has no
//      `alphaMix` to report, so `A` comes from the sum over pairs instead.
//
// The rows are `key = value` and the blocks are `#`-opened, the shape
// `tools/neqsim_layer_diff.py` reads.

import java.lang.reflect.Field;
import neqsim.thermo.phase.PhaseCPAInterface;
import neqsim.thermo.phase.PhaseElectrolyteCPA;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemElectrolyteCPAstatoil;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.flashops.CO2BrinePhaseEquilibrium;
import neqsim.thermodynamicoperations.flashops.CalcIonicComposition;
import neqsim.thermodynamicoperations.flashops.ReactiveCO2BrinePhaseEquilibrium;
import neqsim.thermodynamicoperations.flashops.SaturateWithWater;

public class BrineProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  private static void text(String key, String value) {
    System.out.printf("%s = %s%n", key, value);
  }

  /** A failure, with the NeqSim frames that produced it. */
  private static void fail(Throwable error) {
    text("error", error.getClass().getSimpleName() + ": " + error.getMessage());
    for (StackTraceElement frame : error.getStackTrace()) {
      if (frame.getClassName().startsWith("neqsim.")) {
        System.out.println("#   at " + frame);
      }
    }
  }

  /** One declared field, by name, walking the class chain, as an object. */
  private static Object reflect(Object target, String name) {
    for (Class<?> k = target.getClass(); k != null; k = k.getSuperclass()) {
      try {
        Field f = k.getDeclaredField(name);
        f.setAccessible(true);
        return f.get(target);
      } catch (Throwable ignored) {
        // keep walking up
      }
    }
    return null;
  }

  /** One declared field, by name, walking the class chain. */
  private static double field(Object target, String name) {
    for (Class<?> k = target.getClass(); k != null; k = k.getSuperclass()) {
      try {
        Field f = k.getDeclaredField(name);
        f.setAccessible(true);
        return f.getDouble(target);
      } catch (Throwable ignored) {
        // keep walking up
      }
    }
    return Double.NaN;
  }

  /**
   * The fluid NeqSim's own `CO2BrineHydratePhaseStateTest.brine()` builds.
   *
   * The salt is a mole number derived from a mass percent, not a mole fraction:
   * `salt/(100-salt)/0.05844277` with `0.05844277` the molar mass of NaCl in kg/mol.
   */
  private static SystemInterface brine(double pressure, double saltPercent, double co2Moles) {
    SystemInterface fluid = new SystemElectrolyteCPAstatoil(283.15, pressure);
    fluid.addComponent("CO2", co2Moles);
    fluid.addComponent("water", 1.0 / 0.01801528);
    if (saltPercent > 0.0) {
      double saltMoles = saltPercent / (100.0 - saltPercent) / 0.05844277;
      fluid.addComponent("Na+", saltMoles);
      fluid.addComponent("Cl-", saltMoles);
    }
    fluid.setMixingRule(10);
    fluid.setMultiPhaseCheck(true);
    return fluid;
  }

  private static void layers(String label, SystemInterface fluid) {
    System.out.printf("# %s%n", label);
    try {
      // A flash first, so every phase's `x` below is an equilibrium composition and the
      // `lnPhi` beside it is the pair a case file states.
      new ThermodynamicOperations(fluid).TPflash();
      fluid.init(1);
      row("phases", fluid.getNumberOfPhases());
      text("phase_class", fluid.getPhase(0).getClass().getSimpleName());
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        System.out.printf("# %s phase %d%n", label, p);
        text("type", String.valueOf(phase.getType()));
        row("T", phase.getTemperature());
        row("P", phase.getPressure());
        row("V", phase.getMolarVolume());
        row("Z", phase.getZ());
        row("n_total", phase.getNumberOfMolesInPhase());
        row("B_phase", phase.getB());
        row("A_phase", ((neqsim.thermo.phase.PhaseEos) phase).getA());

        int n = phase.getNumberOfComponents();
        for (int i = 0; i < n; i++) {
          neqsim.thermo.component.ComponentEosInterface species =
              (neqsim.thermo.component.ComponentEosInterface) phase.getComponent(i);
          row("x[" + i + "]", phase.getComponent(i).getx());
          row("charge[" + i + "]", phase.getComponent(i).getIonicCharge());
          row("b[" + i + "]", species.getb());
          row("a[" + i + "]", species.geta());
          // `a` is temperature-independent and `aT` is `a alpha(T)`, so the pair reads out
          // the alpha the phase actually runs. `ComponentElectrolyteCPA.setAttractiveTerm`
          // calls `getAttractiveTerm().setm(mCPA)` for a component with a fitted set, so
          // this is Schwartzentruber's form with the *fitted* Soave coefficient rather than
          // the one `0.48508 + 1.55191 w - 0.15613 w^2` gives.
          row("aT[" + i + "]", species.getaT());
          row("alpha[" + i + "]", species.getAttractiveTerm().alpha(phase.getTemperature()));
          // **Which attractive term the phase installs**, which the alpha alone only hints at.
          // `SystemElectrolyteCPAstatoil`'s constructor sets `attractiveTermNumber = 15`, and
          // that number is what each component's term is built from.
          text("term[" + i + "]", species.getAttractiveTerm().getClass().getSimpleName());
          // **The Huron-Vidal rule's own per-component quantity.** The phase's `A` is
          // `n B R T alphaMix` with `alphaMix = sum_i x_i ader_i`, so `ader` is where a
          // divergence in `A` localises to a component. An ion's first term is `1e-35`
          // over `b_i R T`, so everything an ion contributes is its `-ln gamma_i/lambda`.
          row("ader[" + i + "]", species.getAder());
          row("n_sites_c[" + i + "]", phase.getComponent(i).getNumberOfAssociationSites());
          // **What the phase actually reads as the association parameters**, which is not
          // always the table's row: `COMP.csv` gives `CO2` a `2A` scheme with both
          // `associationenergy` and `associationboundingvolume_SRK` at `0.0`, and the
          // component's two sites still come out at `X = 0.628` rather than `1.0`.
          row("assoc_energy[" + i + "]", phase.getComponent(i).getAssociationEnergy());
          row("assoc_volume[" + i + "]", phase.getComponent(i).getAssociationVolume());
        }

        if (!(phase instanceof PhaseElectrolyteCPA)) {
          System.out.println("note = not an electrolyte CPA phase: " + phase.getClass().getSimpleName());
          continue;
        }
        PhaseElectrolyteCPA cpa = (PhaseElectrolyteCPA) phase;

        // The RDF the statiil override installs, and the values `PhaseSrkCPA` would give.
        row("n_sites", ((PhaseCPAInterface) cpa).getTotalNumberOfAccociationSites());
        row("cpaon", field(cpa, "cpaon"));
        row("gcpa", field(cpa, "gcpa"));
        row("gcpav", field(cpa, "gcpav"));
        row("gcpavv", field(cpa, "gcpavv"));
        row("gcpavvv", field(cpa, "gcpavvv"));
        row("calc_g", cpa.calc_g());
        row("calc_lngV", cpa.calc_lngV());
        row("calc_lngVV", cpa.calc_lngVV());
        row("calc_lngVVV", cpa.calc_lngVVV());
        row("hcpatot", field(cpa, "hcpatot"));

        // `alphaMix`, summed here from the rows above rather than read from a getter, so the
        // two can be checked against each other.
        double alphaMix = 0.0;
        for (int i = 0; i < n; i++) {
          alphaMix += phase.getComponent(i).getx()
              * ((neqsim.thermo.component.ComponentEosInterface) phase.getComponent(i)).getAder();
        }
        row("alpha_mix", alphaMix);

        // **The association strength matrix the phase built.** `calcDeltaNog` for an
        // Elliott pair is `sqrt(t1 t2)` with `t = (exp(eps/RT) - 1) b kappa` of the *self*
        // pair, so a component whose energy and volume are zero should have an all-zero row
        // and a site fraction of exactly `1.0`. Printing the matrix is what says whether
        // that is what the phase computed.
        Object deltaNog = reflect(cpa, "deltaNog");
        if (deltaNog instanceof double[][] matrix) {
          for (int i = 0; i < matrix.length; i++) {
            System.out.printf("deltaNog[%d] = %s%n", i, java.util.Arrays.toString(matrix[i]));
          }
        }
        // The identity the port checks `A` against: `calcA() = n B R T alphaMix`.
        row("A_from_alpha_mix", phase.getNumberOfMolesInPhase() * phase.getB()
            * 8.3144621 * phase.getTemperature() * alphaMix);

        // The two halves of the sum, separately.
        row("FCPA", cpa.FCPA());
        row("dFCPAdT", cpa.dFCPAdT());
        row("dFCPAdTdT", cpa.dFCPAdTdT());
        row("dFCPAdTdV", cpa.dFCPAdTdV());
        row("dFCPAdV", cpa.dFCPAdV());
        row("dFCPAdVdV", cpa.dFCPAdVdV());
        row("dFCPAdVdVdV", cpa.dFCPAdVdVdV());
        row("F", cpa.getF());
        row("F_super", cpa.getF() - field(cpa, "cpaon") * cpa.FCPA());
        row("dFdT", cpa.dFdT());
        row("dFdT_super", cpa.dFdT() - field(cpa, "cpaon") * cpa.dFCPAdT());
        row("dFdV", cpa.dFdV());
        row("dFdV_super", cpa.dFdV() - field(cpa, "cpaon") * cpa.dFCPAdV());
        row("dFdVdV", cpa.dFdVdV());
        row("dFdVdVdV", cpa.dFdVdVdV());
        row("dFdTdT", cpa.dFdTdT());
        row("dFdTdV", cpa.dFdTdV());

        for (int i = 0; i < n; i++) {
          row("lnPhi[" + i + "]", Math.log(cpa.getComponent(i).getFugacityCoefficient()));
          // The association's own per-site fractions, which the port's site solve reproduces.
          try {
            double[] xsite = ((neqsim.thermo.component.ComponentCPAInterface) cpa.getComponent(i)).getXsite();
            System.out.printf("xsite[%d] = %s%n", i, java.util.Arrays.toString(xsite));
          } catch (Throwable ignored) {
            // a component with no sites
          }
        }
      }
    } catch (Throwable error) {
      System.out.printf("# %s failed%n", label);
      text("error", error.getClass().getSimpleName() + ": " + error.getMessage());
    }
  }

  private static void report(String label, SystemInterface system) {
    System.out.printf("# %s%n", label);
    row("phases", system.getNumberOfPhases());
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      PhaseInterface phase = system.getPhase(p);
      System.out.printf("# %s phase %d%n", label, p);
      text("type", String.valueOf(phase.getType()));
      row("beta", system.getBeta(p));
      row("T", phase.getTemperature());
      row("P", phase.getPressure());
      for (int i = 0; i < phase.getNumberOfComponents(); i++) {
        row("x[" + i + "]", phase.getComponent(i).getx());
        row("n[" + i + "]", phase.getComponent(i).getNumberOfMolesInPhase());
      }
    }
  }

  /** The non-reactive constrained equilibrium, and the trial distance it decided on. */
  private static void phaseEquilibrium(String label, double pressure, double salt, double co2) {
    SystemInterface fluid = brine(pressure, salt, co2);
    System.out.printf("# %s%n", label);
    try {
      CO2BrinePhaseEquilibrium equilibrium = new CO2BrinePhaseEquilibrium(fluid);
      text("applicable", String.valueOf(CO2BrinePhaseEquilibrium.isApplicable(fluid)));
      equilibrium.run();
      row("minimum_trial_distance", equilibrium.getMinimumTrialDistance());
      report(label, fluid);
    } catch (Throwable error) {
      fail(error);
    }
  }

  /**
   * `ReactiveCO2BrineHydrateOrderTest.brine(order = 0, ...)`: five components, and
   * `chemicalReactionInit()` **before** the mixing rule.
   *
   * The order is not cosmetic. `chemicalReactionInit` adds the species its tables name, and
   * `setMixingRule` builds the electrolyte rule's `wijCalcOrFitted` at the component count of
   * the moment - so a fluid that sets the rule first throws `Index 4 out of bounds for length
   * 4` from `EosMixingRuleHandler$ElectrolyteMixRule.calcWij` at its first `init`.
   */
  private static SystemInterface reactiveBrine(double pressure, double saltScale) {
    String[] names = { "CO2", "water", "Na+", "K+", "Cl-" };
    double sodium = saltScale * (3.0 / 95.0) / 0.05844277;
    double potassium = saltScale * (2.0 / 95.0) / 0.0745513;
    double[] amounts = { 10.0, 1.0 / 0.01801528, sodium, potassium, sodium + potassium };
    SystemInterface fluid = new SystemElectrolyteCPAstatoil(283.15, pressure);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], amounts[i]);
    }
    fluid.chemicalReactionInit();
    fluid.setMixingRule(10);
    fluid.setMultiPhaseCheck(true);
    fluid.setHydrateCheck(true);
    return fluid;
  }

  /** The same, with the aqueous chemistry coupled in. */
  private static void reactiveEquilibrium(String label, double pressure, double salt, double co2) {
    SystemInterface fluid = reactiveBrine(pressure, salt);
    System.out.printf("# %s%n", label);
    try {
      // `chemicalReactionInit` is what the five-component brine grows to nine with.
      row("components", fluid.getNumberOfComponents());
      text("chemical_system", String.valueOf(fluid.isChemicalSystem()));
      text("applicable", String.valueOf(ReactiveCO2BrinePhaseEquilibrium.isApplicable(fluid)));
      ReactiveCO2BrinePhaseEquilibrium equilibrium = new ReactiveCO2BrinePhaseEquilibrium(fluid);
      equilibrium.run();
      row("iterations", equilibrium.getIterations());
      row("minimum_trial_distance", equilibrium.getMinimumTrialDistance());
      report(label, fluid);
    } catch (Throwable error) {
      fail(error);
    }
  }

  /** The molality report, which solves nothing. */
  private static void ionicComposition(String label, double pressure, double salt, double co2, int phase) {
    SystemInterface fluid = brine(pressure, salt, co2);
    System.out.printf("# %s%n", label);
    try {
      new ThermodynamicOperations(fluid).TPflash();
      // **The phase the table is about**, which the table itself does not state. The model
      // takes a phase and not a feed - the class solves nothing - so its composition is the
      // input a case needs, and the table gives it only as six-figure molalities.
      PhaseInterface reported = fluid.getPhase(phase);
      text("phase_type", String.valueOf(reported.getType()));
      row("phase_T", reported.getTemperature());
      row("phase_P", reported.getPressure());
      for (int i = 0; i < reported.getNumberOfComponents(); i++) {
        text("phase_name[" + i + "]", reported.getComponent(i).getComponentName());
        row("phase_x[" + i + "]", reported.getComponent(i).getx());
        row("phase_n[" + i + "]", reported.getComponent(i).getNumberOfMolesInPhase());
      }
      CalcIonicComposition calc = new CalcIonicComposition(fluid, phase);
      calc.run();
      String[][] table = calc.getResultTable();
      row("rows", table.length);
      for (int i = 0; i < table.length; i++) {
        System.out.printf("row[%d] = %s | %s | %s | %s%n", i, table[i][0], table[i][1], table[i][2], table[i][3]);
      }
    } catch (Throwable error) {
      fail(error);
    }
  }

  /** The water flood, which does need a flash. */
  private static void saturateWithWater(String label, SystemInterface fluid) {
    System.out.printf("# %s%n", label);
    try {
      text("system", fluid.getClass().getSimpleName());
      ThermodynamicOperations ops = new ThermodynamicOperations(fluid);
      ops.TPflash();
      ops.saturateWithWater();
      row("phases", fluid.getNumberOfPhases());
      for (int i = 0; i < fluid.getNumberOfComponents(); i++) {
        text("name[" + i + "]", fluid.getPhase(0).getComponent(i).getComponentName());
        row("n[" + i + "]", fluid.getPhase(0).getComponent(i).getNumberOfmoles());
        row("z[" + i + "]", fluid.getPhase(0).getComponent(i).getz());
      }
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        System.out.printf("# %s phase %d%n", label, p);
        text("type", String.valueOf(fluid.getPhase(p).getType()));
        row("beta", fluid.getBeta(p));
        for (int i = 0; i < fluid.getNumberOfComponents(); i++) {
          row("x[" + i + "]", fluid.getPhase(p).getComponent(i).getx());
        }
      }
    } catch (Throwable error) {
      fail(error);
    }
  }

  /** `SaturateWithWaterTest.testRun4`'s fluid: a plain cubic, which this library can flash. */
  private static SystemInterface prWater() {
    SystemInterface fluid = new neqsim.thermo.system.SystemPrEos(273.15 + 20.0, 150.0);
    fluid.addComponent("methane", 2.0);
    fluid.addComponent("n-heptane", 75.0);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    return fluid;
  }

  /** `SaturateWithWaterTest.testRun2`'s fluid: the CPA one the flood is written for. */
  private static SystemInterface srkCpaWater() {
    SystemInterface fluid = new neqsim.thermo.system.SystemSrkCPAstatoil(273.15 + 70.0, 150.0);
    fluid.addComponent("methane", 75.0);
    fluid.addComponent("ethane", 7.5);
    fluid.addComponent("propane", 4.0);
    fluid.setMixingRule(10);
    return fluid;
  }

  public static void main(String[] args) {
    // The layer stack at the state NeqSim's own test fixes first: 40 bara, 10 % salt.
    layers("electrolyte_cpa_statoil at 40 bara 10 pct salt", brine(40.0, 10.0, 10.0));
    layers("electrolyte_cpa_statoil at 100 bara 5 pct salt", brine(100.0, 5.0, 2.0));
    layers("electrolyte_cpa_statoil at 40 bara no salt", brine(40.0, 0.0, 10.0));

    phaseEquilibrium("co2_brine_phase_equilibrium at 40 bara 10 pct salt", 40.0, 10.0, 10.0);
    phaseEquilibrium("co2_brine_phase_equilibrium at 100 bara 5 pct salt 0.5 mol", 100.0, 5.0, 0.5);
    phaseEquilibrium("co2_brine_phase_equilibrium at 200 bara 10 pct salt", 200.0, 10.0, 10.0);

    reactiveEquilibrium("reactive_co2_brine_phase_equilibrium at 50 bara salt 1.0", 50.0, 1.0, 10.0);
    reactiveEquilibrium("reactive_co2_brine_phase_equilibrium at 40 bara salt 0.9", 40.0, 0.9, 10.0);

    ionicComposition("calc_ionic_composition at 40 bara 10 pct salt aqueous", 40.0, 10.0, 10.0, 1);
    ionicComposition("calc_ionic_composition at 100 bara 5 pct salt aqueous", 100.0, 5.0, 2.0, 1);

    saturateWithWater("saturate_with_water PR methane n-heptane at 20 C 150 bara", prWater());
    saturateWithWater("saturate_with_water SRK CPAstatoil methane ethane propane 70 C 150 bara", srkCpaWater());
  }
}
