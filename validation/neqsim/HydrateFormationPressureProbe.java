// The hydrate formation pressure, printed a layer at a time so the answer can be ported.
//
// `HydrateFormationPressureFlash` is a fixed point on pressure: it copies the fluid's
// per-component fugacities into the hydrate's `reffug`, recomputes the hydrate water's
// fugacity coefficient at `x_water = 1`, reflashes, and then **scales the pressure by the
// ratio of the two water fugacities** until that ratio is one to `1e-8`. The layers below
// are the two ends of it - the pressure it lands on and the two water fugacities that are
// equal there - so a port can be checked against the answer *and* against the equality it
// is the answer to.
//
// The same feed as `HydrateProbe`, at temperatures below its formation temperature, so
// there is a pressure to find rather than none.
//
//   javac -proc:none -cp neqsim-3.20.0.jar HydrateFormationPressureProbe.java
//   java -cp .:neqsim-3.20.0.jar HydrateFormationPressureProbe > captures/hydrate_formation_pressure_probe.tsv

import neqsim.thermo.component.ComponentHydrate;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class HydrateFormationPressureProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    // The three states the fraction probe drives, all below the feed's 293.23 K at 100 bara.
    double[] temperatures = {288.15, 283.15, 278.15};
    for (double temperatureK : temperatures) {
      report(temperatureK, 100.0);
    }
  }

  private static void report(double temperatureK, double startBara) {
    System.out.printf(
        "# methane/ethane/propane/water at T = %.15g K, started at %.15g bara%n",
        temperatureK, startBara);
    try {
      SystemInterface fluid = new SystemSrkEos(temperatureK, startBara);
      fluid.addComponent("methane", 79.0);
      fluid.addComponent("ethane", 10.10);
      fluid.addComponent("propane", 2.050);
      fluid.addComponent("water", 10.0);
      fluid.setMixingRule(2);
      fluid.setMultiPhaseCheck(true);
      fluid.setHydrateCheck(true);

      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      operations.hydrateFormationPressure();

      row("formation_pressure_bara", fluid.getPressure());
      row("temperature_K", fluid.getTemperature());
      row("phases", fluid.getNumberOfPhases());

      PhaseInterface hydrate = hydratePhase(fluid);
      if (hydrate == null) {
        System.out.println("note = no hydrate phase in the system");
        System.out.println();
        return;
      }
      // The stable structure lives on the *water* component, as `HydrateProbe` records.
      row(
          "hydrate_structure",
          ((ComponentHydrate) hydrate.getComponent("water")).getHydrateStructure());
      row("hydrate_water_fugacity", hydrate.getFugacity("water"));

      // **The phase the pressure flash partners**, which NeqSim's `setFug` hardcodes to
      // phase 0 - the temperature flash prefers the aqueous one, so this is not the same
      // choice and is printed rather than assumed.
      row("phase0_water_fugacity", fluid.getPhase(0).getFugacity("water"));
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        System.out.printf("phase[%d] = %s%n", p, phase.getType());
        row("phase[" + p + "].beta", phase.getBeta());
        row("phase[" + p + "].water_fugacity", phase.getFugacity("water"));
      }
      // The objective at the answer, which is the quantity the fixed point drives to zero.
      row(
          "objective",
          Math.log(hydrate.getFugacity("water") / fluid.getPhase(0).getFugacity("water")));
      for (int i = 0; i < fluid.getPhase(0).getNumberOfComponents(); i++) {
        System.out.printf(
            "f[%s] = %.15g%n",
            fluid.getPhase(0).getComponent(i).getName(), fluid.getPhase(0).getFugacity(i));
      }
      row("phase0_z_water", fluid.getPhase(0).getComponent("water").getx());
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
}
