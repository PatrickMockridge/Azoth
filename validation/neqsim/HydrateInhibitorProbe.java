// The hydrate inhibitor concentration, printed so the secant and its inner solve can be ported.
//
// `HydrateInhibitorConcentrationFlash` is a **secant on the inhibitor's moles** whose residual
// is `T_hydrate - T_target`: each trial adds MEG or methanol and asks what the hydrate
// temperature is now. Nothing in it reads a phase type, so unlike its `...wtFlash` sibling it
// needs no aqueous phase to exist.
//
// **Two fluids, on purpose.** NeqSim's own `main` for this class uses `SystemSrkCPAstatoil`
// with `setMixingRule(9)`, which is a CPA flash; the second block is the same composition on a
// plain SRK one, which is the fluid a port without a CPA flash can build. The two disagree,
// and printing both is what makes the disagreement a measurement rather than an assumption -
// the CPA association moves water's fugacity, and the hydrate equilibrium is on water.
//
// **The secant's own branches are visible in the iteration count.** It adds `error * 0.01`
// moles while `iter < 4` and `-error/derrordC * 0.5` after, and its loop condition is
// `(|error| > 1e-3 && iter < 100) || iter < 3`, so three iterations always run.
//
// A capture, not a test: run it against the pinned jar and commit what it prints.
//   javac -proc:none -cp neqsim-f0c7436.jar HydrateInhibitorProbe.java
//   java -cp .:neqsim-f0c7436.jar HydrateInhibitorProbe > captures/hydrate_inhibitor_probe.tsv

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkCPAstatoil;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class HydrateInhibitorProbe {

  /** NeqSim's own `main` for this class, which is a CPA fluid. */
  private static final String[] NAMES = {"methane", "ethane", "propane", "i-butane", "MEG", "water"};

  private static final double[] MOLES = {1.0, 0.10, 0.050, 0.0050, 0.1, 1.0};

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    for (double target : new double[] {270.9, 265.0, 275.0}) {
      report("cpa", target, 100.0);
      report("srk", target, 100.0);
    }
  }

  private static SystemInterface build(String fluid, double pressureBara) {
    SystemInterface system = fluid.equals("cpa")
        ? new SystemSrkCPAstatoil(273.15, pressureBara)
        : new SystemSrkEos(273.15, pressureBara);
    for (int i = 0; i < NAMES.length; i++) {
      system.addComponent(NAMES[i], MOLES[i]);
    }
    system.createDatabase(true);
    // NeqSim's own `main` uses the CPA mixing rule; the plain fluid has to use a rule it has.
    system.setMixingRule(fluid.equals("cpa") ? 9 : 2);
    system.init(0);
    system.setMultiPhaseCheck(true);
    system.setHydrateCheck(true);
    return system;
  }

  private static void report(String fluid, double target, double pressureBara) {
    System.out.printf("# %s fluid, MEG to a hydrate temperature of %.15g K at %.15g bara%n",
        fluid, target, pressureBara);
    try {
      SystemInterface system = build(fluid, pressureBara);
      ThermodynamicOperations operations = new ThermodynamicOperations(system);
      operations.hydrateInhibitorConcentration("MEG", target);

      // **The answer is the system's own inventory**, which the flash added to in place - there
      // is no result object, and NeqSim's `main` reads it the same way.
      double inhibitorMoles = system.getPhase(0).getComponent("MEG").getNumberOfmoles();
      double waterMoles = system.getPhase(0).getComponent("water").getNumberOfmoles();
      double inhibitorMass = inhibitorMoles * system.getPhase(0).getComponent("MEG").getMolarMass();
      double waterMass = waterMoles * system.getPhase(0).getComponent("water").getMolarMass();

      row("inhibitor_moles", inhibitorMoles);
      row("water_moles", waterMoles);
      row("total_moles", system.getTotalNumberOfMoles());
      row("inhibitor_weight_fraction", inhibitorMass / (inhibitorMass + waterMass));
      row("weight_percent", 100.0 * inhibitorMass / (inhibitorMass + waterMass));
      // The temperature the secant left the system at, which is the hydrate temperature its
      // last trial produced - the residual it stopped on is this against the target.
      row("hydrate_temperature_reached", system.getTemperature());
      row("residual", system.getTemperature() - target);

      row("phases", system.getNumberOfPhases());
      for (int p = 0; p < system.getNumberOfPhases(); p++) {
        PhaseInterface phase = system.getPhase(p);
        row("phase[" + p + "].beta", phase.getBeta());
        row("phase[" + p + "].temperature", phase.getTemperature());
        if (phase.hasComponent("MEG")) {
          row("phase[" + p + "].x_MEG", phase.getComponent("MEG").getx());
          row("phase[" + p + "].x_water", phase.getComponent("water").getx());
        }
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
