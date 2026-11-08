// `ThermalDesignCalculator` driven directly, on states of its own.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar ThermalRatingProbe.java
//     java -cp .:neqsim-f0c7436.jar ThermalRatingProbe > captures/thermal_rating_probe.tsv
//
// **The equipment path is measured separately.** `process_heat_exchanger_rating.tsv` runs
// `HeatExchanger.runThermalHydraulicRating`, which feeds this calculator from two flashed
// streams. This probe states the calculator's own inputs so the arithmetic is observable
// without a flash: tube side Gnielinski/laminar, shell side Kern or Bell-Delaware, the
// resistance-model overall U, and both pressure drops.
//
// **Both entry points are measured.** `calculate()` leaves zeros where its guards return
// early; `calculateStrict()` refuses the same state with an exception. The port implements
// the strict semantics, so each refusal row prints both answers and the class's message.
//
// Every row prints the inputs it read; the getters are the class's, and the two correction
// blocks are read from `toMap()` because the J factors have no getters of their own.

import java.util.Map;
import neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator;
import neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator.ShellSideMethod;

public class ThermalRatingProbe {

  public static void main(String[] args) {
    // 1. The constructor's own geometry and both fluids. This is the class's default state.
    ThermalDesignCalculator defaults = new ThermalDesignCalculator();
    row("defaults_kern", defaults, null);

    // 2. The fouling test's geometry (ID 0.020, OD 0.025, pitch 0.032, 100 tubes, 2 passes,
    //    L 4.0) on the calculator's default fluids: a geometry-only move from block 1.
    ThermalDesignCalculator geometry = geometry();
    row("kern_fouling_test_geometry", geometry, null);

    // 3. Bell-Delaware on the same state, so every correction factor is visible.
    ThermalDesignCalculator bell = geometry();
    bell.setShellSideMethod(ShellSideMethod.BELL_DELAWARE);
    row("bell_delaware_fouling_test_geometry", bell, null);

    // 4. A sealed Bell-Delaware bundle, which moves Jb through the sealing-pair count.
    ThermalDesignCalculator sealed = geometry();
    sealed.setShellSideMethod(ShellSideMethod.BELL_DELAWARE);
    sealed.setHasSealing(true);
    sealed.setSealingPairs(2);
    row("bell_delaware_sealed", sealed, null);

    // 5. Square pitch, which moves the equivalent diameter and the crossflow area.
    ThermalDesignCalculator square = geometry();
    square.setTriangularPitch(false);
    row("kern_square_pitch", square, null);

    // 6. Laminar tube side: Re below 2300, so Nu is the class's 3.66.
    ThermalDesignCalculator laminar = geometry();
    laminar.setTubeSideFluid(800.0, 5.0e-4, 2500.0, 0.13, 0.5, true);
    row("laminar_tube_side", laminar, null);

    // 7. Transition tube side: Re between 2300 and 3000, so the class interpolates between
    //    3.66 and the Gnielinski value at Re = 3000.
    ThermalDesignCalculator transition = geometry();
    transition.setTubeSideFluid(800.0, 5.0e-4, 2500.0, 0.13, 1.0, true);
    row("transition_tube_side", transition, null);

    // 8. The fouling pair: clean against 0.0003 m2 K/W on both sides.
    ThermalDesignCalculator clean = geometry();
    clean.setFoulingTube(0.0);
    clean.setFoulingShell(0.0);
    row("fouling_clean", clean, null);
    ThermalDesignCalculator fouled = geometry();
    fouled.setFoulingTube(3.0e-4);
    fouled.setFoulingShell(3.0e-4);
    row("fouling_pair", fouled, null);

    // 9. A deposit layer inside the tube, which moves the effective ID and adds a resistance.
    ThermalDesignCalculator deposit = geometry();
    deposit.setTubeFoulingLayer(0.001, 0.2);
    row("deposit_layer", deposit, null);

    // 10. A wall viscosity stated below the bulk, which the classes' correlations read.
    ThermalDesignCalculator wallViscosity = geometry();
    wallViscosity.setShellViscosityWall(1.0e-4);
    row("shell_wall_viscosity", wallViscosity, null);

    // 11. A bypass area stated rather than estimated.
    ThermalDesignCalculator statedBypass = geometry();
    statedBypass.setShellSideMethod(ShellSideMethod.BELL_DELAWARE);
    statedBypass.setBypassArea(0.01);
    row("bell_delaware_stated_bypass", statedBypass, null);

    // 12. `calculate()` leaves zeros where `calculateStrict()` refuses: zero tubes.
    refusalRows();

    // 13. The remaining permissive guards: zero inner diameter and zero tube length.
    zeroGeometryRows();

    System.out.flush();
  }

  /// The fouling test's geometry with the calculator's default fluids.
  static ThermalDesignCalculator geometry() {
    ThermalDesignCalculator calculator = new ThermalDesignCalculator();
    calculator.setTubeIDm(0.020);
    calculator.setTubeODm(0.025);
    calculator.setTubePitchm(0.032);
    calculator.setTubeCount(100);
    calculator.setTubePasses(2);
    calculator.setTubeLengthm(4.0);
    return calculator;
  }

  /// The permissive and strict entry points on the same invalid state.
  static void refusalRows() {
    // Zero tubes: `calculate()` returns early and leaves every result at zero.
    ThermalDesignCalculator zeroTubes = geometry();
    zeroTubes.setTubeCount(0);
    row("strict_refusal_zero_tubes_calculate", zeroTubes, null);
    ThermalDesignCalculator zeroTubesStrict = geometry();
    zeroTubesStrict.setTubeCount(0);
    strictRow("strict_refusal_zero_tubes_strict", zeroTubesStrict);

    // A pitch no greater than the outside diameter: the cross-field guard.
    ThermalDesignCalculator badPitch = geometry();
    badPitch.setTubePitchm(0.025);
    row("strict_refusal_pitch_not_above_od_calculate", badPitch, null);
    ThermalDesignCalculator badPitchStrict = geometry();
    badPitchStrict.setTubePitchm(0.025);
    strictRow("strict_refusal_pitch_not_above_od_strict", badPitchStrict);

    // A baffle cut at the half boundary, which strict refuses (`>= 0.5`).
    ThermalDesignCalculator badCut = geometry();
    badCut.setBaffleCut(0.5);
    strictRow("strict_refusal_baffle_cut_half", badCut);

    // A deposit that consumes the bore: the setter validates on assignment, so the refusal
    // happens there and not in calculateStrict().
    System.out.println("strict_refusal_deposit_blocks_bore");
    try {
      ThermalDesignCalculator badDeposit = geometry();
      badDeposit.setTubeFoulingLayer(0.011, 0.2);
      System.out.println("setter_result=accepted");
    } catch (RuntimeException ex) {
      System.out.println("setter_exception=" + ex.getClass().getSimpleName());
      System.out.println("setter_message=" + ex.getMessage());
    }
    System.out.println();
  }

  /// The permissive guards that leave zeros rather than throwing.
  static void zeroGeometryRows() {
    ThermalDesignCalculator zeroId = geometry();
    zeroId.setTubeIDm(0.0);
    row("zero_tube_id_calculate", zeroId, null);

    ThermalDesignCalculator zeroLength = geometry();
    zeroLength.setTubeLengthm(0.0);
    row("zero_tube_length_calculate", zeroLength, null);
  }

  /// One full row: call `calculate()` and print every getter and the map.
  static void row(String label, ThermalDesignCalculator calculator, String note) {
    try {
      calculator.calculate();
    } catch (RuntimeException ex) {
      System.out.println(label);
      System.out.println("calculate_exception=" + ex.getClass().getSimpleName());
      System.out.println("calculate_message=" + ex.getMessage());
      System.out.println();
      return;
    }
    printRow(label, calculator, note);
  }

  /// One strict row: call `calculateStrict()` and print either the answers or the refusal.
  static void strictRow(String label, ThermalDesignCalculator calculator) {
    System.out.println(label);
    try {
      calculator.calculateStrict();
      printGetters(calculator);
      System.out.println("strict_result=accepted");
    } catch (RuntimeException ex) {
      System.out.println("strict_exception=" + ex.getClass().getSimpleName());
      System.out.println("strict_message=" + ex.getMessage());
    }
    System.out.println();
  }

  static void printRow(String label, ThermalDesignCalculator calculator, String note) {
    System.out.println(label);
    if (note != null) {
      System.out.println("note=" + note);
    }
    printInputs(calculator);
    printGetters(calculator);
    System.out.println();
  }

  static void printInputs(ThermalDesignCalculator calculator) {
    // No getters for the geometry inputs, so the row states the block's own values through
    // the derived getters it can read and leaves the rest to the label's construction.
    System.out.println("effective_tube_inside_diameter_m=" + calculator.getEffectiveTubeIDm());
    System.out.println("outside_heat_transfer_area_m2=" + calculator.getOutsideHeatTransferArea());
    System.out.println("shell_side_method=" + calculator.getShellSideMethod());
  }

  static void printGetters(ThermalDesignCalculator calculator) {
    System.out.println("tube_side_htc_W_per_m2K=" + calculator.getTubeSideHTC());
    System.out.println("shell_side_htc_W_per_m2K=" + calculator.getShellSideHTC());
    System.out.println("overall_u_W_per_m2K=" + calculator.getOverallU());
    System.out.println("tube_side_pressure_drop_Pa=" + calculator.getTubeSidePressureDrop());
    System.out.println("tube_side_pressure_drop_bar=" + calculator.getTubeSidePressureDropBar());
    System.out.println("shell_side_pressure_drop_Pa=" + calculator.getShellSidePressureDrop());
    System.out.println("shell_side_pressure_drop_bar=" + calculator.getShellSidePressureDropBar());
    System.out.println("tube_side_velocity_m_per_s=" + calculator.getTubeSideVelocity());
    System.out.println("shell_side_velocity_m_per_s=" + calculator.getShellSideVelocity());
    System.out.println("tube_side_reynolds=" + calculator.getTubeSideRe());
    System.out.println("shell_side_reynolds=" + calculator.getShellSideRe());
    System.out.println("deposit_thickness_m=" + calculator.getTubeFoulingLayerThicknessM());
    System.out.println("deposit_conductivity_W_per_mK=" + calculator.getTubeFoulingLayerConductivity());
    System.out.println(
        "deposit_resistance_outside_m2K_per_W=" + calculator.getTubeFoulingLayerResistanceOutside());
    printMap("map", calculator.toMap());
  }

  /// A stable key-order walk of the class's own `toMap()`, which is a LinkedHashMap.
  static void printMap(String prefix, Map<String, Object> map) {
    for (Map.Entry<String, Object> entry : map.entrySet()) {
      String key = prefix + "." + entry.getKey();
      Object value = entry.getValue();
      if (value instanceof Map<?, ?>) {
        @SuppressWarnings("unchecked")
        Map<String, Object> nested = (Map<String, Object>) value;
        printMap(key, nested);
      } else {
        System.out.println(key + "=" + value);
      }
    }
  }
}
