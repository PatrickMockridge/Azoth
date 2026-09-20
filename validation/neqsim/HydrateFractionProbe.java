// The hydrate fraction's layers, printed one at a time so each can be ported and checked.
//
// The formation temperature is where a hydrate *first* appears and its fraction is zero;
// this is the state inside the region, where the amount is the answer. NeqSim's
// `TPHydrateFlash` solves it by secant on the fraction with a material balance at every
// trial: the existing phases are scaled by `(1 - beta)/sum(beta)`, the hydrate takes `beta`
// with a composition from the cavity occupancies, and the objective is the same
// `ln(f_w^hydrate/f_w^fluid)` the formation temperature drives to zero.
//
// The two bounds bound one another: the most hydrate the feed can make is its water divided
// by the water fraction in the hydrate - `46/54` for structure I, `136/160` for structure II.
//
//   javac -proc:none -cp neqsim-3.20.0.jar HydrateFractionProbe.java
//   java -cp .:neqsim-3.20.0.jar HydrateFractionProbe > captures/hydrate_fraction_probe.tsv

import neqsim.thermo.component.ComponentHydrate;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class HydrateFractionProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    // The same feed the formation probe uses, at temperatures *below* its 293.23 K at
    // 100 bara - so the hydrate is present and the fraction is what is solved for.
    double[] temperatures = {288.15, 283.15, 278.15};
    for (double temperatureK : temperatures) {
      report(temperatureK, 100.0);
    }
  }

  private static void report(double temperatureK, double pressureBara) {
    System.out.printf(
        "# methane/ethane/propane/water at T = %.15g K, P = %.15g bara%n",
        temperatureK, pressureBara);
    try {
      SystemInterface fluid = new SystemSrkEos(temperatureK, pressureBara);
      fluid.addComponent("methane", 79.0);
      fluid.addComponent("ethane", 10.10);
      fluid.addComponent("propane", 2.050);
      fluid.addComponent("water", 10.0);
      fluid.setMixingRule(2);
      fluid.setMultiPhaseCheck(true);
      fluid.setHydrateCheck(true);

      ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
      operations.hydrateTPflash();

      row("hydrate_fraction", fluid.getHydrateFraction());
      row("phases", fluid.getNumberOfPhases());
      row("temperature_K", fluid.getTemperature());
      row("pressure_bara", fluid.getPressure());

      // Every phase's fraction and what it is made of, because the fraction is a *material
      // balance*: a port that solved the fugacity equality without moving the moles would
      // land on a number and mean a different state.
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        PhaseInterface phase = fluid.getPhase(p);
        System.out.printf("phase[%d] = %s%n", p, phase.getType());
        row("phase[" + p + "].beta", phase.getBeta());
        for (int i = 0; i < phase.getNumberOfComponents(); i++) {
          row(
              "phase[" + p + "].x[" + phase.getComponent(i).getName() + "]",
              phase.getComponent(i).getx());
        }
      }

      // The hydrate's own constants and the objective at the answer.
      PhaseInterface hydrate = null;
      for (PhaseInterface phase : fluid.getPhases()) {
        if (phase != null && phase.getType() == neqsim.thermo.phase.PhaseType.HYDRATE) {
          hydrate = phase;
        }
      }
      if (hydrate == null) {
        System.out.println("note = no hydrate phase in the system");
        System.out.println();
        return;
      }
      row(
          "hydrate_structure",
          ((ComponentHydrate) hydrate.getComponent("water")).getHydrateStructure());
      row("hydrate_water_fugacity", hydrate.getFugacity("water"));

      // The water-bearing phase the equilibrium is read against, whichever it is.
      int waterPhase = 0;
      double most = -1.0;
      for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
        if (fluid.getPhase(p).hasComponent("water")) {
          double fraction = fluid.getPhase(p).getComponent("water").getx();
          if (fraction > most) {
            most = fraction;
            waterPhase = p;
          }
        }
      }
      row("water_phase_index", waterPhase);
      row("fluid_water_fugacity", fluid.getPhase(waterPhase).getFugacity("water"));
      row(
          "objective",
          Math.log(
              hydrate.getFugacity("water")
                  / fluid.getPhase(waterPhase).getFugacity("water")));
      // And the guests' fugacities the occupancies are built from.
      for (int i = 0; i < fluid.getPhase(0).getNumberOfComponents(); i++) {
        row("f[" + fluid.getPhase(0).getComponent(i).getName() + "]",
            fluid.getPhase(0).getFugacity(i));
      }

      // **The material balance, which is what the fraction *is*.** The hydrate's fraction is
      // bounded by the feed's water over the water fraction in the hydrate, so a fraction at
      // that bound says every water molecule is in the hydrate - and then no other phase may
      // hold any. Counted per component across the phases, against the feed.
      int n = fluid.getPhase(0).getNumberOfComponents();
      for (int i = 0; i < n; i++) {
        double inPhases = 0.0;
        for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
          inPhases += fluid.getPhase(p).getBeta() * fluid.getPhase(p).getComponent(i).getx();
        }
        String name = fluid.getPhase(0).getComponent(i).getName();
        row("balance[" + name + "]", inPhases);
        row("feed[" + name + "]", fluid.getPhase(0).getComponent(i).getz());
        row("balance_error[" + name + "]", inPhases - fluid.getPhase(0).getComponent(i).getz());
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
