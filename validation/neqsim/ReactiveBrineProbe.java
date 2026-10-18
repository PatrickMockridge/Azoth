// The reactive CO2/brine equilibrium: the species `chemicalReactionInit` adds, and the
// residuals its coupling loop is accepted on.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ReactiveBrineProbe.java
//     java -cp .:neqsim-f0c7436.jar ReactiveBrineProbe > captures/reactive_brine_probe.tsv
//
// `BrineProbe` runs the class and reports the phase inventory, but it prints the *count* of
// components and never their names - and the nine the five-component brine grows to are the
// input a port has to state. Its `reactive_co2_brine_phase_equilibrium` block also does not
// carry the three residuals `ReactiveCO2BrinePhaseEquilibrium.run` accepts on:
//
//     reaction residual  <= 2e-6      `getMaximumAbsoluteReactionLogResidual`
//     electroneutrality  <= 1e-8      `getReactivePhaseChargeMoles`
//     element balance    <= its own   `getMaximumAbsoluteElementBalanceResidual`
//
// The third has no constant in the class: it is checked against the element vector `verifyElements`
// computes from the *original feed*, at `max(1e-9, |expected| * 1e-10)`. Printing it is what says
// how much room a port has.
//
// The rows are `key = value` and the blocks are `#`-opened, the shape `tools/oracle_sweep.py`
// and `tools/neqsim_layer_diff.py` read.

import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemElectrolyteCPAstatoil;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.flashops.ReactiveCO2BrinePhaseEquilibrium;

/** The reactive brine, and the state it is accepted at. */
public final class ReactiveBrineProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  private static void text(String key, String value) {
    System.out.printf("%s = %s%n", key, value);
  }

  /** `BrineProbe.reactiveBrine`: the same five components, grown by the reaction init. */
  private static SystemInterface reactiveBrine(double pressure, double saltScale) {
    String[] names = {"CO2", "water", "Na+", "K+", "Cl-"};
    double sodium = saltScale * (3.0 / 95.0) / 0.05844277;
    double potassium = saltScale * (2.0 / 95.0) / 0.0745513;
    double[] amounts = {10.0, 1.0 / 0.01801528, sodium, potassium, sodium + potassium};
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

  /** One state: the species, the loop's count, and the residuals it accepted on. */
  private static void run(String label, double pressure, double saltScale) {
    System.out.printf("# %s%n", label);
    try {
      SystemInterface fluid = reactiveBrine(pressure, saltScale);
      row("components", fluid.getNumberOfComponents());
      for (int i = 0; i < fluid.getNumberOfComponents(); i++) {
        text("name[" + i + "]", fluid.getComponent(i).getName());
      }

      // **Before the solve**, because the class reads this accessor after the phase flash and
      // before the chemistry: a port that starts from a fresh solve has to know what the
      // initialized object answers.
      ChemicalReactionOperations before = fluid.getChemicalReactionOperations();
      row("initial_reaction_log_residual", before.getMaximumAbsoluteReactionLogResidual());
      row("initial_element_balance_residual",
          before.getMaximumAbsoluteElementBalanceResidual());
      row("initial_charge_moles", before.getReactivePhaseChargeMoles());

      ReactiveCO2BrinePhaseEquilibrium equilibrium = new ReactiveCO2BrinePhaseEquilibrium(fluid);
      equilibrium.run();

      row("iterations", equilibrium.getIterations());
      row("minimum_trial_distance", equilibrium.getMinimumTrialDistance());

      ChemicalReactionOperations after = fluid.getChemicalReactionOperations();
      row("final_reaction_log_residual", after.getMaximumAbsoluteReactionLogResidual());
      row("final_element_balance_residual",
          after.getMaximumAbsoluteElementBalanceResidual());
      row("final_charge_moles", after.getReactivePhaseChargeMoles());

      row("phases", fluid.getNumberOfPhases());
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        System.out.printf("# %s phase %d%n", label, p);
        text("type", String.valueOf(phase.getType()));
        row("beta", fluid.getBeta(p));
        for (int i = 0; i < phase.getNumberOfComponents(); i++) {
          row("x[" + i + "]", phase.getComponent(i).getx());
          row("n[" + i + "]", phase.getComponent(i).getNumberOfMolesInPhase());
        }
      }
    } catch (Throwable error) {
      text("error", String.valueOf(error));
    }
  }

  /** Two states, the ones `BrineProbe` already reports. */
  public static void main(String[] args) {
    run("reactive brine at 50 bara salt 1.0", 50.0, 1.0);
    run("reactive brine at 40 bara salt 0.9", 40.0, 0.9);
  }

  private ReactiveBrineProbe() {}
}
