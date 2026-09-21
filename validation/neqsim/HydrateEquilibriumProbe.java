// The hydrate equilibrium line's points, and whether the *seed* between them changes them.
//
// `HydrateEquilibriumLine.run()` walks a pressure grid and, from the second point on, sets the
// system's temperature to the previous point's answer before each solve - "hydrate T increases
// with P", so the previous answer is the next point's starting guess. That makes the result a
// property of the grid *and the path through it*, which a port has to decide whether to keep.
//
// So this prints two vectors: the line's own, and the same ten pressures solved **independently**
// on a fresh system, where no answer can seed the next. If the two agree the seeding is a
// speed-up and not a model; if they differ the seeding is part of the answer.
//
// **Two more things the walk hides.** `numberOfPoints` is a field initialised to ten with no
// setter, so the grid is always ten points whatever a caller asks for. And a point whose solve
// *throws* is caught and then reported anyway - `hydratePoints[0][i]` reads the system's
// temperature unconditionally - so a failed point silently repeats the previous temperature
// rather than leaving a hole.
//
// A capture, not a test: run it against the pinned jar and commit what it prints.
//   javac -proc:none -cp neqsim-f0c7436.jar HydrateEquilibriumProbe.java
//   java -cp .:neqsim-f0c7436.jar HydrateEquilibriumProbe > captures/hydrate_equilibrium_probe.tsv

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermodynamicoperations.flashops.saturationops.HydrateEquilibriumLine;

public class HydrateEquilibriumProbe {

  private static final double MINIMUM_PRESSURE_BARA = 1.0;

  private static final double MAXIMUM_PRESSURE_BARA = 200.0;

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    System.out.printf("# methane/ethane/propane/water, %.15g to %.15g bara%n",
        MINIMUM_PRESSURE_BARA, MAXIMUM_PRESSURE_BARA);
    try {
      SystemInterface fluid = build();
      fluid.setHydrateCheck(true);
      // The operation itself rather than `ThermodynamicOperations.hydrateEquilibriumLine`: the
      // points are on the operation, and the line is the thing this probe is about.
      HydrateEquilibriumLine line =
          new HydrateEquilibriumLine(fluid, MINIMUM_PRESSURE_BARA, MAXIMUM_PRESSURE_BARA);
      line.run();

      double[][] points = line.getPoints(0);
      int count = points[0].length;
      row("points", count);
      for (int i = 0; i < count; i++) {
        row("line_temperature[" + i + "]", points[0][i]);
        row("line_pressure[" + i + "]", points[1][i]);
      }

      // The same pressures, each on a system that has just been built, so nothing can carry
      // over. This is the vector the line's seeding is measured against.
      double worst = 0.0;
      for (int i = 0; i < count; i++) {
        double pressure = points[1][i];
        SystemInterface fresh = build();
        fresh.setHydrateCheck(true);
        fresh.setPressure(pressure);
        try {
          new ThermodynamicOperations(fresh).hydrateFormationTemperature();
          row("independent_temperature[" + i + "]", fresh.getTemperature());
          double relative =
              Math.abs(fresh.getTemperature() / points[0][i] - 1.0);
          row("relative_gap[" + i + "]", relative);
          worst = Math.max(worst, relative);
        } catch (Exception ex) {
          System.out.printf("independent_temperature[%d] = NaN%n", i);
          System.out.printf("# independent solve failed: %s%n", ex.getMessage());
        }
      }
      row("worst_relative_gap", worst);
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }

  /** The fluid `HydrateProbe` uses, so the two captures describe one state. */
  private static SystemInterface build() {
    SystemInterface fluid = new SystemSrkEos(273.15, MINIMUM_PRESSURE_BARA);
    fluid.addComponent("methane", 79.0);
    fluid.addComponent("ethane", 10.10);
    fluid.addComponent("propane", 2.050);
    fluid.addComponent("water", 10.0);
    fluid.setMixingRule(2);
    fluid.setMultiPhaseCheck(true);
    return fluid;
  }
}
