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
      int structure = 0;
      for (int i = 0; i < hydrate.getNumberOfComponents(); i++) {
        if (hydrate.getComponent(i) instanceof ComponentHydrate) {
          structure = ((ComponentHydrate) hydrate.getComponent(i)).getHydrateStructure();
          break;
        }
      }
      row("hydrate_structure", structure);

      // The gas phase the guests' fugacities come from, and the two water fugacities the
      // equilibrium is the equality of.
      row("gas_water_fugacity", fluid.getPhase(0).getFugacity("water"));
      row("hydrate_water_fugacity", hydrate.getFugacity("water"));

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
