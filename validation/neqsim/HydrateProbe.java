// The hydrate formation layers, printed one at a time so each can be ported and checked.
//
// `HydrateFormationPressureFlash` is a fixed point on pressure and this is the same
// equilibrium read the other way round: the hydrate's water fugacity is scaled against the
// fluid's until they match, with each guest's fugacity fed in from the gas phase by
// `setFug`. The layers below are what a port has to reproduce between those two ends - the
// structure the phase selects, each cavity's Langmuir constant and occupancy, and the two
// water fugacities at the answer.
//
// **A plain SRK fluid on purpose.** The hydrate model reads its guests' fugacities from the
// fluid it is attached to, so the arithmetic can be checked without also checking a mixing
// rule: NeqSim's own published case for this is `SystemSrkCPAstatoil` with
// `CLASSIC_TX_CPA`, which is a second thing to be right about rather than a first.
//
// A capture, not a test: run it against the pinned jar and commit what it prints.
//   javac -proc:none -cp neqsim-3.20.0.jar HydrateProbe.java
//   java -cp .:neqsim-3.20.0.jar HydrateProbe > captures/hydrate_probe.tsv

import neqsim.thermo.component.ComponentHydrate;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class HydrateProbe {

  /** One `key = value` row, which is the shape `tools/neqsim_layer_diff.py` reads. */
  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    // The composition NeqSim's own `hydrateEquilibriumTemperature` test uses - methane,
    // ethane, propane and water at 100 bara - on a plain cubic system.
    double[] pressures = {100.0, 50.0, 200.0};
    for (double pressureBara : pressures) {
      report(pressureBara);
    }
  }

  private static void report(double pressureBara) {
    System.out.printf("# methane/ethane/propane/water at P = %.15g bara%n", pressureBara);
    try {
      SystemInterface fluid = new SystemSrkEos(273.15, pressureBara);
      fluid.addComponent("methane", 79.0);
      fluid.addComponent("ethane", 10.10);
      fluid.addComponent("propane", 2.050);
      fluid.addComponent("water", 10.0);
      fluid.setMixingRule(2);
      fluid.setMultiPhaseCheck(true);
      fluid.setHydrateCheck(true);

      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      operations.hydrateFormationTemperature();

      row("temperature_K", fluid.getTemperature());
      row("pressure_bara", fluid.getPressure());
      row("phases", fluid.getNumberOfPhases());

      // The hydrate phase, wherever the system put it, and the layers inside it.
      PhaseInterface hydrate = hydratePhase(fluid);
      if (hydrate == null) {
        System.out.println("note = no hydrate phase in the system");
        System.out.println();
        return;
      }
      row("hydrate_phase_index", hydratePhaseIndex(fluid));
      // **The stable structure lives on the *water* component**, and reading it from whichever
      // component comes first reads a stale zero: NeqSim's `ComponentHydratePVTsim.fugcoef`
      // runs its two-structure selection inside the water branch and writes the winner into
      // the water component's own field. A capture that read methane's reported structure I
      // for a phase whose cages are structure II, and every occupancy printed from it was the
      // wrong structure's.
      int structure = ((ComponentHydrate) hydrate.getComponent("water")).getHydrateStructure();
      row("hydrate_structure", structure);

      // **Which phase the equilibrium partners, and every phase's water fugacity.** NeqSim's
      // temperature flash prefers the *aqueous* phase for this comparison and its pressure
      // flash takes phase 0, so which one is the partner is not a detail a port can assume -
      // and the phases are printed so the answer is read rather than guessed.
      row("hydrate_water_fugacity", hydrate.getFugacity("water"));
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        System.out.printf("phase[%d] = %s%n", p, phase.getType());
        row("phase[" + p + "].beta", phase.getBeta());
        row("phase[" + p + "].water_fugacity", phase.getFugacity("water"));
      }
      row("phase0_water_fugacity", fluid.getPhase(0).getFugacity("water"));
      // **The kernel's own input**: the gas phase's per-component fugacities, which are what
      // `setFug` copies into the hydrate's `reffug` before every occupancy evaluation. Printed
      // so the hydrate arithmetic can be checked against the capture *without* the fluid - a
      // divergence in an occupancy is then the kernel's and not the flash's.
      for (int i = 0; i < fluid.getPhase(0).getNumberOfComponents(); i++) {
        System.out.printf(
            "f[%s] = %.15g%n",
            fluid.getPhase(0).getComponent(i).getName(), fluid.getPhase(0).getFugacity(i));
      }

      for (int i = 0; i < hydrate.getNumberOfComponents(); i++) {
        String name = hydrate.getComponent(i).getName();
        if (!(hydrate.getComponent(i) instanceof ComponentHydrate component)) {
          continue;
        }
        for (int cavity = 0; cavity < 2; cavity++) {
          row("cavprwat[" + name + "," + cavity + "]", component.getCavprwat(structure, cavity));
          row(
              "theta[" + name + "," + cavity + "]",
              component.calcYKI(structure, cavity, hydrate));
        }
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }

  private static PhaseInterface hydratePhase(SystemInterface system) {
    for (PhaseInterface phase : system.getPhases()) {
      if (phase != null && phase.getType() == neqsim.thermo.phase.PhaseType.HYDRATE) {
        return phase;
      }
    }
    return null;
  }

  private static double hydratePhaseIndex(SystemInterface system) {
    for (int i = 0; i < system.getPhases().length; i++) {
      if (system.getPhases()[i] != null
          && system.getPhases()[i].getType() == neqsim.thermo.phase.PhaseType.HYDRATE) {
        return i;
      }
    }
    return Double.NaN;
  }
}
