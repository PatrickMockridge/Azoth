package neqsim.thermodynamicoperations.flashops;

// Which construct puts a phase into NeqSim's multiphase flash, attributed state by state.
//
//     javac -proc:none -cp neqsim-f0c7436.jar -d . TpMultiflashSeedProbe.java
//     java -cp .:neqsim-f0c7436.jar neqsim.thermodynamicoperations.flashops.TpMultiflashSeedProbe \
//       > captures/tp_multiflash_seed_probe.tsv
//
// # Why this probe declares NeqSim's package
//
// `TPmultiflash.multiPhaseTest` and the two ablation flags are package-private, and the
// gas/oil seeds are private methods. A probe outside the package would have to read the
// chain rather than measure it - and this class has already cost this tree one wrong
// reading. The two private seeds are called through reflection, which is stated here
// rather than hidden.
//
// # What it asks
//
// For each state, the parent two-phase flash is run first (the state `TPmultiflash.run()`
// starts from), then each construct is measured on its own **fresh** fluid:
//
//   basic      TPmultiflash.stabilityAnalysis()         phase count before -> after
//   enhanced   TPmultiflash.stabilityAnalysisEnhanced()  phase count before -> after
//   gas seed   seedAdditionalPhaseFromFeed()             return value
//   oil seed   seedHydrocarbonLiquidFromFeed()           return value
//   inline     the full run, with `aqueousPhaseSeedAttempted` pre-set (ablated) and not
//   post       the full run, with `postFlashStabilityChecked` pre-set (ablated) and not
//
// A construct whose ablation changes the final phase inventory fired. The full run is also
// printed, because that is NeqSim's answer and what a spec case records.
//
// **Every measurement builds its own fluid.** `setMultiPhaseCheck(true)` is a system flag,
// and reusing a fluid after a run makes the next `TPflash` a multiflash rather than the
// two-phase seed this probe is measuring from - the first version of this probe did that and
// reported every construct as a no-op on an answer the earlier call had already found.
//
// **The two-argument constructor, always.** `TPmultiflash(system)` never initialises
// `multTerm`/`multTerm2` (the parent always passes `doSolidPhaseCheck`), so a call to
// `solveBeta` on the one-argument form throws a NullPointerException before any answer.

import java.lang.reflect.Method;

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class TpMultiflashSeedProbe {

  /** How a state is built: the `pr` mixture and the water-dew-point gas differ in EOS and rule. */
  static SystemInterface build(String kind, String[] names, double[] moles, double temperatureK,
      double pressureBar, double ppmWater, boolean multi) {
    SystemInterface fluid;
    if (kind.equals("pr")) {
      fluid = new SystemPrEos(temperatureK, pressureBar);
      for (int i = 0; i < names.length; i++) {
        fluid.addComponent(names[i], moles[i]);
      }
      fluid.setMixingRule("classic");
      fluid.setAttractiveTerm(1);
    } else {
      double x = ppmWater * 1e-6;
      fluid = new SystemSrkEos(temperatureK, pressureBar);
      fluid.addComponent("CO2", 0.02);
      fluid.addComponent("nitrogen", 0.01);
      fluid.addComponent("methane", 0.9 - x);
      fluid.addComponent("ethane", 0.05);
      fluid.addComponent("propane", 0.01);
      fluid.addComponent("i-butane", 0.005);
      fluid.addComponent("n-butane", 0.005);
      fluid.addComponent("water", x);
      fluid.createDatabase(true);
      fluid.setMixingRule(2);
    }
    fluid.setMultiPhaseCheck(multi);
    return fluid;
  }

  static SystemInterface fresh(String kind, String[] names, double[] moles, double temperatureK,
      double pressureBar, double ppmWater) {
    return build(kind, names, moles, temperatureK, pressureBar, ppmWater, false);
  }

  /** The ordinary two-phase flash `TPmultiflash.run()` starts from, on a fresh fluid. */
  static SystemInterface twoPhase(String kind, String[] names, double[] moles, double temperatureK,
      double pressureBar, double ppmWater) {
    SystemInterface fluid = fresh(kind, names, moles, temperatureK, pressureBar, ppmWater);
    new ThermodynamicOperations(fluid).TPflash();
    fluid.setMultiPhaseCheck(true);
    return fluid;
  }

  static TPmultiflash flash(SystemInterface fluid) {
    return new TPmultiflash(fluid, false);
  }

  private static void row(String key, double value) {
    System.out.printf("%s = %.17g%n", key, value);
  }

  static double waterX(SystemInterface fluid, int phase) {
    try {
      return fluid.getPhase(phase).getComponent("water").getx();
    } catch (RuntimeException ex) {
      return Double.NaN;
    }
  }

  static void phases(String label, SystemInterface fluid) {
    try {
      row(label + ".phases", fluid.getNumberOfPhases());
      for (int phase = 0; phase < fluid.getNumberOfPhases(); phase++) {
        System.out.printf("    %s.phase_%d_type = %s%n", label, phase,
            fluid.getPhase(phase).getType());
        row(String.format("    %s.phase_%d_beta", label, phase), fluid.getBeta(phase));
        row(String.format("    %s.phase_%d_z", label, phase), fluid.getPhase(phase).getZ());
        row(String.format("    %s.phase_%d_water_x", label, phase), waterX(fluid, phase));
        StringBuilder x = new StringBuilder();
        for (int comp = 0; comp < fluid.getPhase(phase).getNumberOfComponents(); comp++) {
          x.append(String.format("  %s=%.17g", fluid.getPhase(phase).getComponent(comp).getName(),
              fluid.getPhase(phase).getComponent(comp).getx()));
        }
        System.out.printf("    %s.phase_%d_x =%s%n", label, phase, x);
      }
    } catch (RuntimeException ex) {
      System.out.printf("  %s = unreadable: %s%n", label, ex);
    }
  }

  static boolean privateSeed(TPmultiflash flash, String name) {
    try {
      Method method = TPmultiflash.class.getDeclaredMethod(name);
      method.setAccessible(true);
      return (Boolean) method.invoke(flash);
    } catch (ReflectiveOperationException ex) {
      System.out.printf("  %s = unavailable: %s%n", name, ex);
      return false;
    }
  }

  /** The final answer of the full chain, and the two ablations of it. */
  static void full(String kind, String[] names, double[] moles, double temperatureK,
      double pressureBar, double ppmWater) {
    SystemInterface normal = twoPhase(kind, names, moles, temperatureK, pressureBar, ppmWater);
    TPmultiflash full = flash(normal);
    full.run();
    phases("full", normal);
    row("full.multiPhaseTest", full.multiPhaseTest ? 1.0 : 0.0);

    SystemInterface noInline =
        twoPhase(kind, names, moles, temperatureK, pressureBar, ppmWater);
    TPmultiflash ablated = flash(noInline);
    ablated.aqueousPhaseSeedAttempted = true;
    ablated.run();
    phases("full_no_inline_aqueous", noInline);

    SystemInterface noPost = twoPhase(kind, names, moles, temperatureK, pressureBar, ppmWater);
    TPmultiflash noPostFlash = flash(noPost);
    noPostFlash.postFlashStabilityChecked = true;
    noPostFlash.run();
    phases("full_no_post_enhanced", noPost);
  }

  /** Every construct on its own fresh fluid, then the full chain. */
  static void attribute(String label, String kind, String[] names, double[] moles,
      double temperatureK, double pressureBar, double ppmWater) {
    System.out.printf("%n=== %s%n", label);

    SystemInterface initial = twoPhase(kind, names, moles, temperatureK, pressureBar, ppmWater);
    phases("initial", initial);

    SystemInterface basicFluid = twoPhase(kind, names, moles, temperatureK, pressureBar, ppmWater);
    TPmultiflash basic = flash(basicFluid);
    int before = basicFluid.getNumberOfPhases();
    basic.stabilityAnalysis();
    row("basic.added", basicFluid.getNumberOfPhases() - before);
    row("basic.multiPhaseTest", basic.multiPhaseTest ? 1.0 : 0.0);
    phases("basic", basicFluid);

    SystemInterface enhancedFluid =
        twoPhase(kind, names, moles, temperatureK, pressureBar, ppmWater);
    TPmultiflash enhanced = flash(enhancedFluid);
    int enhancedBefore = enhancedFluid.getNumberOfPhases();
    enhanced.stabilityAnalysisEnhanced();
    row("enhanced.added", enhancedFluid.getNumberOfPhases() - enhancedBefore);
    row("enhanced.multiPhaseTest", enhanced.multiPhaseTest ? 1.0 : 0.0);
    phases("enhanced", enhancedFluid);

    SystemInterface gasFluid = twoPhase(kind, names, moles, temperatureK, pressureBar, ppmWater);
    TPmultiflash gasFlash = flash(gasFluid);
    gasFlash.stabilityAnalysis();
    row("gas_seed.fired", privateSeed(gasFlash, "seedAdditionalPhaseFromFeed") ? 1.0 : 0.0);
    phases("gas_seed", gasFluid);

    SystemInterface oilFluid = twoPhase(kind, names, moles, temperatureK, pressureBar, ppmWater);
    TPmultiflash oilFlash = flash(oilFluid);
    oilFlash.stabilityAnalysis();
    row("oil_seed.fired", privateSeed(oilFlash, "seedHydrocarbonLiquidFromFeed") ? 1.0 : 0.0);
    phases("oil_seed", oilFluid);

    full(kind, names, moles, temperatureK, pressureBar, ppmWater);
  }

  /** One line per state where a construct changes the answer, for finding case states. */
  static void scan(String label, String kind, String[] names, double[] moles, double[] pressures,
      double temperatureStart, int temperatureSteps, double ppmWater) {
    for (int ti = 0; ti < temperatureSteps; ti++) {
      double temperature = temperatureStart + 10.0 * ti;
      for (double pressure : pressures) {
        StringBuilder fired = new StringBuilder();
        int phases = -1;
        int fullPhases = -1;
        String fullTypes = "";
        try {
          SystemInterface initial = fresh(kind, names, moles, temperature, pressure, ppmWater);
          new ThermodynamicOperations(initial).TPflash();
          phases = initial.getNumberOfPhases();

          SystemInterface basicFluid =
              twoPhase(kind, names, moles, temperature, pressure, ppmWater);
          TPmultiflash basic = flash(basicFluid);
          int before = basicFluid.getNumberOfPhases();
          basic.stabilityAnalysis();
          if (basicFluid.getNumberOfPhases() > before) {
            fired.append(" basic");
          }

          SystemInterface enhancedFluid =
              twoPhase(kind, names, moles, temperature, pressure, ppmWater);
          TPmultiflash enhanced = flash(enhancedFluid);
          int enhancedBefore = enhancedFluid.getNumberOfPhases();
          enhanced.stabilityAnalysisEnhanced();
          if (enhancedFluid.getNumberOfPhases() > enhancedBefore) {
            fired.append(" enhanced");
          }

          SystemInterface gasFluid =
              twoPhase(kind, names, moles, temperature, pressure, ppmWater);
          TPmultiflash gasFlash = flash(gasFluid);
          gasFlash.stabilityAnalysis();
          if (privateSeed(gasFlash, "seedAdditionalPhaseFromFeed")) {
            fired.append(" gas");
          }

          SystemInterface oilFluid =
              twoPhase(kind, names, moles, temperature, pressure, ppmWater);
          TPmultiflash oilFlash = flash(oilFluid);
          oilFlash.stabilityAnalysis();
          if (privateSeed(oilFlash, "seedHydrocarbonLiquidFromFeed")) {
            fired.append(" oil");
          }

          SystemInterface normal = twoPhase(kind, names, moles, temperature, pressure, ppmWater);
          TPmultiflash normalFlash = flash(normal);
          normalFlash.run();
          fullPhases = normal.getNumberOfPhases();
          StringBuilder types = new StringBuilder();
          for (int phase = 0; phase < fullPhases; phase++) {
            types.append(phase == 0 ? "" : ",").append(normal.getPhase(phase).getType());
          }
          fullTypes = types.toString();

          SystemInterface noInline = twoPhase(kind, names, moles, temperature, pressure, ppmWater);
          TPmultiflash noInlineFlash = flash(noInline);
          noInlineFlash.aqueousPhaseSeedAttempted = true;
          noInlineFlash.run();
          if (normal.getNumberOfPhases() != noInline.getNumberOfPhases()) {
            fired.append(" inline_aqueous");
          }

          SystemInterface noPost = twoPhase(kind, names, moles, temperature, pressure, ppmWater);
          TPmultiflash noPostFlash = flash(noPost);
          noPostFlash.postFlashStabilityChecked = true;
          noPostFlash.run();
          if (normal.getNumberOfPhases() != noPost.getNumberOfPhases()) {
            fired.append(" post_enhanced");
          }
        } catch (RuntimeException ex) {
          continue;
        }
        if (fired.length() > 0) {
          System.out.printf(
              "<state> %s T=%.0f P=%.0f plain_phases=%d full_phases=%d full_types=%s fired=%s%n",
              label, temperature, pressure, phases, fullPhases, fullTypes, fired.toString().trim());
        }
      }
    }
  }

  public static void main(String[] args) {
    neqsim.util.database.NeqSimDataBase.createTemporaryTables();

    String[] n2co2c8 = {"nitrogen", "CO2", "n-octane"};
    double[] n2co2c8moles = {0.1, 0.4, 0.5};
    String[] co2c1c10 = {"CO2", "methane", "nc10"};
    double[] co2c1c10moles = {0.4, 0.3, 0.3};
    String[] waterHexane = {"water", "n-hexane"};
    double[] waterHexaneMoles = {0.5, 0.5};

    attribute("N2/CO2/nC8 10/40/50 180 K 10 bar", "pr", n2co2c8, n2co2c8moles, 180.0, 10.0, 0.0);
    attribute("CO2/C1/nC10 40/30/30 180 K 2 bar", "pr", co2c1c10, co2c1c10moles, 180.0, 2.0, 0.0);
    attribute("water/nC6 50/50 350 K 20 bar", "pr", waterHexane, waterHexaneMoles, 350.0, 20.0,
        0.0);
    attribute("8-component gas 35 ppm water 260.545 K 70 bar", "water", new String[] {"water"},
        new double[] {35.0}, 260.545463562012, 70.0, 35.0);
    attribute("8-component gas 35 ppm water 230 K 80 bar", "water", new String[] {"water"},
        new double[] {35.0}, 230.0, 80.0, 35.0);
    attribute("8-component gas 35 ppm water 224.2 K 1 bar", "water", new String[] {"water"},
        new double[] {35.0}, 224.2, 1.0, 35.0);

    String[] c1c4 = {"methane", "n-butane"};
    double[] c1c4moles = {0.5, 0.5};
    String[] c1c4water = {"methane", "n-butane", "water"};
    double[] c1c4waterMoles = {0.5, 0.5, 1.0e-4};

    attribute("C1/nC4 50/50 370 K 40 bar", "pr", c1c4, c1c4moles, 370.0, 40.0, 0.0);
    attribute("C1/nC4 50/50 + 1e-4 water 370 K 40 bar", "pr", c1c4water, c1c4waterMoles, 370.0,
        40.0, 0.0);

    double[] pressures = {1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0};
    System.out.println();
    System.out.println("# scan");
    scan("C1/nC4", "pr", c1c4, c1c4moles, pressures, 300.0, 10, 0.0);
    scan("C1/nC4 + 1e-4 water", "pr", c1c4water, c1c4waterMoles, pressures, 300.0, 10, 0.0);
    scan("C1/nC7 + 1e-4 water", "pr", new String[] {"methane", "n-heptane", "water"},
        new double[] {0.5, 0.5, 1.0e-4}, pressures, 300.0, 10, 0.0);
    scan("CO2/nC6 + 1e-4 water", "pr", new String[] {"CO2", "n-hexane", "water"},
        new double[] {0.7, 0.3, 1.0e-4}, pressures, 260.0, 10, 0.0);
    scan("N2/CO2/nC8 + 1e-4 water", "pr",
        new String[] {"nitrogen", "CO2", "n-octane", "water"},
        new double[] {0.1, 0.4, 0.5, 1.0e-4}, pressures, 180.0, 10, 0.0);
    scan("N2/CO2/nC8", "pr", n2co2c8, n2co2c8moles, pressures, 180.0, 8, 0.0);
    scan("CO2/C1/nC10", "pr", co2c1c10, co2c1c10moles, pressures, 180.0, 8, 0.0);
    scan("water/nC6", "pr", waterHexane, waterHexaneMoles, pressures, 260.0, 8, 0.0);
    scan("water/methane 40/60", "pr", new String[] {"water", "methane"}, new double[] {0.4, 0.6},
        pressures, 280.0, 10, 0.0);
    scan("8-component gas 35 ppm", "water", new String[] {"water"}, new double[] {35.0}, pressures,
        230.0, 8, 35.0);
  }
}
