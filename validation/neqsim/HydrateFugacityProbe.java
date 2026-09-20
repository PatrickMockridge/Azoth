// The terms of the hydrate's water fugacity coefficient, printed one at a time.
//
// `HydrateProbe` prints the occupancies and the answer; this prints what sits between them.
// The formula's terms are local variables inside `fugcoef`, so they are recomputed here from
// the public accessors and compared with the coefficient NeqSim itself produced - which is
// the only way to tell a misread term from a misread formula.
//
//   javac -proc:none -cp neqsim-3.20.0.jar HydrateFugacityProbe.java
//   java -cp .:neqsim-3.20.0.jar HydrateFugacityProbe > captures/hydrate_fugacity_probe.tsv

import neqsim.thermo.component.ComponentHydrate;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class HydrateFugacityProbe {

  /** A declared field of any class in the object's hierarchy, for the ones with no getter. */
  private static Object field(Object target, String name) {
    for (Class<?> k = target.getClass(); k != null; k = k.getSuperclass()) {
      try {
        java.lang.reflect.Field found = k.getDeclaredField(name);
        found.setAccessible(true);
        return found.get(target);
      } catch (Throwable ignored) {
        // keep walking up
      }
    }
    throw new IllegalStateException("no field `" + name + "`");
  }

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    double pressureBara = 100.0;
    SystemInterface fluid = new SystemSrkEos(273.15, pressureBara);
    fluid.addComponent("methane", 79.0);
    fluid.addComponent("ethane", 10.10);
    fluid.addComponent("propane", 2.050);
    fluid.addComponent("water", 10.0);
    fluid.setMixingRule(2);
    fluid.setMultiPhaseCheck(true);
    fluid.setHydrateCheck(true);

    ThermodynamicOperations operations = new ThermodynamicOperations(fluid);
    try {
      operations.hydrateFormationTemperature();
    } catch (Exception error) {
      System.out.println("note = the flash failed: " + error);
      return;
    }

    PhaseInterface hydrate = null;
    for (PhaseInterface phase : fluid.getPhases()) {
      if (phase != null && phase.getType() == neqsim.thermo.phase.PhaseType.HYDRATE) {
        hydrate = phase;
      }
    }
    if (hydrate == null) {
      System.out.println("note = no hydrate phase");
      return;
    }

    ComponentHydrate water = (ComponentHydrate) hydrate.getComponent("water");
    int structure = water.getHydrateStructure();
    double temperature = fluid.getTemperature();
    double pressure = fluid.getPressure();

    row("temperature_K", temperature);
    row("pressure_bara", pressure);
    row("structure", structure);
    row("coefficient_reported", water.getFugacityCoefficient());
    row("x_water", hydrate.getComponent("water").getx());
    row("hydrate_water_fugacity", hydrate.getFugacity("water"));
    row("fluid_water_fugacity", fluid.getPhase(0).getFugacity("water"));

    // The occupancies' contribution, recomputed from the public accessors.
    double val = 0.0;
    for (int cavity = 0; cavity < 2; cavity++) {
      double occupied = 0.0;
      for (int j = 0; j < hydrate.getNumberOfComponents(); j++) {
        occupied += ((ComponentHydrate) hydrate.getComponent(j)).calcYKI(structure, cavity, hydrate);
      }
      row("occupied[" + cavity + "]", occupied);
      val += water.getCavprwat(structure, cavity) * Math.log(1.0 - occupied);
    }
    row("val", val);

    // The chemical-potential change, which is public on the component.
    // `calcDeltaChemPot` is on the fitted subclass, not the base.
    neqsim.thermo.component.ComponentHydratePVTsim fitted =
        (neqsim.thermo.component.ComponentHydratePVTsim) water;
    row("delta_chem_pot", fitted.calcDeltaChemPot(hydrate, hydrate.getNumberOfComponents(),
        temperature, pressure, structure));

    // **The reference water phase**, which is what the coefficient is built against and the
    // one term a port cannot guess: `setSolidRefFluidPhase` builds a water-only phase of the
    // *host* class, so its water component is the host's - an SRK water here, not the
    // `ComponentWater` whose `fugcoef` returns one - and its fugacity is a real number rather
    // than the pressure.
    Object refPhase = field(water, "refPhase");
    if (refPhase instanceof PhaseInterface reference) {
      row("ref_phase_temperature", reference.getTemperature());
      row("ref_phase_pressure_bara", reference.getPressure());
      row("ref_phase_number_of_components", reference.getNumberOfComponents());
      System.out.printf("ref_phase_class = %s%n", reference.getClass().getSimpleName());
      reference.init(reference.getNumberOfMolesInPhase(), 1, 3,
          neqsim.thermo.phase.PhaseType.LIQUID, 1.0);
      System.out.printf("ref_water_component = %s%n",
          reference.getComponent("water").getClass().getSimpleName());
      row("ref_water_fugacity_coefficient",
          reference.getComponent("water").fugcoef(reference));
      row("ref_water_fugacity", reference.getFugacity("water"));
    } else {
      System.out.println("note = no reference phase reachable");
    }
    // `reffug` is a field with no accessor, so it is read by reflection - the same reach
    // P8's `FurstProbe` made for `WTT`.
    double[] refFug = (double[]) field(water, "reffug");
    row("alpha_water", refFug[hydrate.getComponent("water").getComponentNumber()]);
    for (int i = 0; i < hydrate.getNumberOfComponents(); i++) {
      row("reffug[" + hydrate.getComponent(i).getName() + "]", refFug[i]);
    }
    row("R", neqsim.thermo.ThermodynamicConstantsInterface.R);
  }
}
