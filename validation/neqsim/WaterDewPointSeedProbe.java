package neqsim.thermodynamicoperations.flashops;

// Which construct puts an aqueous phase into NeqSim's multiphase flash, attributed.
//
//     javac -proc:none -cp neqsim-f0c7436.jar -d . WaterDewPointSeedProbe.java
//     java -cp .:neqsim-f0c7436.jar neqsim.thermodynamicoperations.flashops.WaterDewPointSeedProbe \
//       > captures/water_dew_point_seed_probe.tsv
//
// # Why this probe declares NeqSim's own package
//
// `TPmultiflash.multiPhaseTest` is package-private, and it is the hinge of the whole chain:
// the aqueous seed is inline in `run()` **inside `if (multiPhaseTest)`**, so it can only fire
// once a stability analysis has put a phase in the set. A probe outside the package can see
// neither the field nor the two named seeds, so the attribution would be a reading of the
// source rather than a measurement - and this tree has already recorded one wrong reading of
// exactly this class. `WaxReferenceProbe` and `PcsaftCompositionProbe` are the precedent for
// declaring the package; the recipe is in `captures/README.md`.
//
// # What it asks
//
// Both stability analyses are **public** methods on `TPmultiflash`, so each is called
// directly on a fresh fluid and what it leaves is printed: the phase count, every phase's
// type, beta and water mole fraction, and - the point - `multiPhaseTest`. The full
// `TPflash` under `setMultiPhaseCheck(true)` is run beside them as the control that ties the
// answer to the capture the water-dew-point tranche already holds.
//
// # The gas
//
// NeqSim's own `WaterDewPointTest` fixture, as `WaterDewPointProbe` builds it: CO2 0.02,
// nitrogen 0.01, methane 0.9 - x, ethane 0.05, propane 0.01, i-butane 0.005, n-butane 0.005
// and water x, on a `SystemSrkEos` at 70 bara.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class WaterDewPointSeedProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  /** The `WaterDewPointTest` gas at 70 bara. */
  static SystemInterface gas(double ppmWater) {
    double x = ppmWater * 1e-6;
    SystemInterface system = new SystemSrkEos(270.15, 70.0);
    system.addComponent("CO2", 0.02);
    system.addComponent("nitrogen", 0.01);
    system.addComponent("methane", 0.9 - x);
    system.addComponent("ethane", 0.05);
    system.addComponent("propane", 0.01);
    system.addComponent("i-butane", 0.005);
    system.addComponent("n-butane", 0.005);
    system.addComponent("water", x);
    system.createDatabase(true);
    system.setMixingRule(2);
    return system;
  }

  static SystemInterface at(double ppmWater, double temperatureK) {
    SystemInterface system = gas(ppmWater);
    system.setTemperature(temperatureK);
    system.setPressure(70.0);
    // **The ceiling the analyses step past.** `setMultiPhaseCheck(true)` is what raises the
    // system's phase limit to three, and `stabilityAnalysis` calls `addPhase()` before it
    // indexes the new slot - so a probe that skipped this would crash on NeqSim's own guard
    // rather than read the answer, which is exactly what the first run of this did.
    system.setMultiPhaseCheck(true);
    system.init(0);
    system.init(1);
    return system;
  }

  /**
   * Every phase the fluid holds, and the water each of them carries.
   *
   * **Guarded, because a phase `addPhase()` has just grown is not always readable.** NeqSim's
   * own `restartMultiphaseFromFreshEstimate` javadoc says `addPhase()` "only exposes a stale
   * phase slot", and here that is not a warning but the observed behaviour - so a probe that
   * threw on it would hide the one reading it exists to take.
   */
  static void phases(String label, SystemInterface system) {
    try {
      row(String.format("%s.phases", label), system.getNumberOfPhases());
      row(String.format("%s.max_phases", label), system.getMaxNumberOfPhases());
      for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
        System.out.printf("    phase_%d_type = %s%n", phase, system.getPhase(phase).getType());
        row(String.format("    phase_%d_beta", phase), system.getBeta(phase));
        row(String.format("    phase_%d_water_x", phase),
            system.getPhase(phase).getComponent("water").getx());
      }
    } catch (Exception failure) {
      System.out.printf("  %s = unreadable: %s%n", label, failure);
    }
  }

  /** The state, what each public stability analysis leaves, and the aqueous seed's guard. */
  static void attribute(double ppmWater, double temperatureK) {
    System.out.printf("# %g ppm water, 70 bara, %g K%n", ppmWater, temperatureK);

    SystemInterface fluid = at(ppmWater, temperatureK);
    row("water_z", fluid.getComponent("water").getz());
    row("phases_before", fluid.getNumberOfPhases());
    row("has_aqueous_before", fluid.hasPhaseType("aqueous") ? 1.0 : 0.0);

    SystemInterface basicFluid = at(ppmWater, temperatureK);
    TPmultiflash basic = new TPmultiflash(basicFluid);
    basic.stabilityAnalysis();
    row("multiPhaseTest_after_stabilityAnalysis", basic.multiPhaseTest ? 1.0 : 0.0);
    phases("basic", basicFluid);

    SystemInterface enhancedFluid = at(ppmWater, temperatureK);
    TPmultiflash enhanced = new TPmultiflash(enhancedFluid);
    enhanced.stabilityAnalysisEnhanced();
    row("multiPhaseTest_after_enhanced", enhanced.multiPhaseTest ? 1.0 : 0.0);
    phases("enhanced", enhancedFluid);

    SystemInterface flashed = at(ppmWater, temperatureK);
    flashed.setMultiPhaseCheck(true);
    try {
      new ThermodynamicOperations(flashed).TPflash();
      row("full_flash_has_aqueous", flashed.hasPhaseType("aqueous") ? 1.0 : 0.0);
      phases("full_flash", flashed);
    } catch (Exception failure) {
      System.out.printf("full_flash = failed: %s%n", failure);
    }
    System.out.println();
  }

  public static void main(String[] args) {
    neqsim.util.database.NeqSimDataBase.createTemporaryTables();

    // The temperatures azoth's own `eos.stability_test` was measured at, and the two
    // states NeqSim's capture calls the multiphase water dew point.
    attribute(35.0, 270.0);
    attribute(35.0, 260.545463562012);
    attribute(35.0, 250.0);
    attribute(35.0, 240.0);
    attribute(35.0, 217.0);
    attribute(200.0, 284.302406311035);
  }
}
