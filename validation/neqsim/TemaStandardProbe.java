// `TEMAStandard` driven directly, on its configuration table, enum tables and static methods.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar TemaStandardProbe.java
//     java -cp .:neqsim-f0c7436.jar TemaStandardProbe > captures/tema_standard_probe.tsv
//
// **The class speaks millimetres in every static method**; the port converts from metres. The
// tube-count result is an int truncation, and the material factor is a substring match on the
// upper-cased name.

import neqsim.process.mechanicaldesign.heatexchanger.TEMAStandard;
import neqsim.process.mechanicaldesign.heatexchanger.TEMAStandard.BaffleType;
import neqsim.process.mechanicaldesign.heatexchanger.TEMAStandard.StandardTubeSize;
import neqsim.process.mechanicaldesign.heatexchanger.TEMAStandard.TEMAClass;
import neqsim.process.mechanicaldesign.heatexchanger.TEMAStandard.TEMAConfiguration;
import neqsim.process.mechanicaldesign.heatexchanger.TEMAStandard.TubePitchPattern;

public class TemaStandardProbe {

  public static void main(String[] args) {
    // 1. The most common configuration at Class R.
    configurationRow("aes_r", "AES", TEMAClass.R, StandardTubeSize.TUBE_3_4_INCH,
        TubePitchPattern.TRIANGULAR_30, BaffleType.SINGLE_SEGMENTAL, 500.0, 23.81, 2,
        "CARBON_STEEL", false, false, false, false);

    // 2. A fixed tubesheet at Class C, one pass, stainless.
    configurationRow("bem_c", "BEM", TEMAClass.C, StandardTubeSize.TUBE_1_2_INCH,
        TubePitchPattern.SQUARE_90, BaffleType.DOUBLE_SEGMENTAL, 400.0, 16.0, 1,
        "STAINLESS_STEEL", true, false, false, false);

    // 3. A sealed floating head at Class B, four passes, copper, a large temperature range.
    configurationRow("ajw_b", "AJW", TEMAClass.B, StandardTubeSize.TUBE_1_INCH,
        TubePitchPattern.TRIANGULAR_60, BaffleType.NO_TUBES_IN_WINDOW, 800.0, 31.75, 4,
        "COPPER", false, true, false, false);

    // 4. A U-tube at Class R, two passes, titanium.
    configurationRow("beu_r", "BEU", TEMAClass.R, StandardTubeSize.TUBE_3_8_INCH,
        TubePitchPattern.SQUARE_45, BaffleType.DISC_AND_DOUGHNUT, 300.0, 12.0, 2,
        "TITANIUM", false, false, false, false);

    // 5. A kettle reboiler with rod baffles and no material name.
    configurationRow("akt_c", "AKT", TEMAClass.C, StandardTubeSize.TUBE_1_INCH,
        TubePitchPattern.TRIANGULAR_30, BaffleType.ROD_BAFFLES, 1000.0, 31.75, 1, null,
        false, false, false, false);

    // 6-7. A custom configuration and the min baffle spacing at a small shell.
    customRow("custom_ceu", 'C', 'E', 'U', TEMAClass.B, StandardTubeSize.TUBE_5_8_INCH,
        TubePitchPattern.TRIANGULAR_30, BaffleType.TRIPLE_SEGMENTAL, 200.0, 19.84, 2,
        "STAINLESS", false, false, false, false);
    customRow("small_shell_r", 'B', 'E', 'M', TEMAClass.R, StandardTubeSize.TUBE_3_4_INCH,
        TubePitchPattern.SQUARE_90, BaffleType.SINGLE_SEGMENTAL, 200.0, 23.81, 2,
        "CARBON_STEEL", false, false, false, false);

    // 8-15. Every recommendation branch.
    recommendationRow("recommend_hazardous_dt", false, true, false, true);
    recommendationRow("recommend_hazardous", false, false, false, true);
    recommendationRow("recommend_high_pressure_dt", false, true, true, false);
    recommendationRow("recommend_high_pressure", false, false, true, false);
    recommendationRow("recommend_dt_cleanable", true, true, false, false);
    recommendationRow("recommend_dt", false, true, false, false);
    recommendationRow("recommend_cleanable", true, false, false, false);
    recommendationRow("recommend_default", false, false, false, false);

    // 16-18. The tube count at one, two and four passes on the same geometry.
    countRow("tube_count_one_pass", 500.0, 19.05, 23.81, TubePitchPattern.TRIANGULAR_30, 1);
    countRow("tube_count_two_passes", 500.0, 19.05, 23.81, TubePitchPattern.TRIANGULAR_30, 2);
    countRow("tube_count_four_passes_square", 500.0, 19.05, 23.81, TubePitchPattern.SQUARE_90, 4);

    // 19-21. The unsupported span across the diameter bands and the material factors.
    spanRow("span_small_carbon", 12.7, "CARBON_STEEL");
    spanRow("span_mid_stainless", 19.05, "STAINLESS_STEEL");
    spanRow("span_large_titanium", 28.0, "TITANIUM");

    System.out.flush();
  }

  static void configurationRow(String label, String designation, TEMAClass temaClass,
      StandardTubeSize tubeSize, TubePitchPattern pattern, BaffleType baffle, double shellIdMm,
      double pitchMm, int tubePasses, String material, boolean cleanable, boolean largeDt,
      boolean highPressure, boolean hazardous) {
    TEMAConfiguration config = TEMAStandard.getConfiguration(designation);
    row(label, config, temaClass, tubeSize, pattern, baffle, shellIdMm, pitchMm, tubePasses,
        material, cleanable, largeDt, highPressure, hazardous);
  }

  static void customRow(String label, char front, char shell, char rear, TEMAClass temaClass,
      StandardTubeSize tubeSize, TubePitchPattern pattern, BaffleType baffle, double shellIdMm,
      double pitchMm, int tubePasses, String material, boolean cleanable, boolean largeDt,
      boolean highPressure, boolean hazardous) {
    TEMAConfiguration config = TEMAStandard.createConfiguration(front, shell, rear);
    row(label, config, temaClass, tubeSize, pattern, baffle, shellIdMm, pitchMm, tubePasses,
        material, cleanable, largeDt, highPressure, hazardous);
  }

  static void row(String label, TEMAConfiguration config, TEMAClass temaClass,
      StandardTubeSize tubeSize, TubePitchPattern pattern, BaffleType baffle, double shellIdMm,
      double pitchMm, int tubePasses, String material, boolean cleanable, boolean largeDt,
      boolean highPressure, boolean hazardous) {
    double tubeOdMm = tubeSize.getOuterDiameterMm();
    System.out.println(label);
    System.out.println("designation=" + config.getDesignation());
    System.out.println("tema_class=" + temaClass.name());
    System.out.println("tube_size=" + tubeSize.name());
    System.out.println("pitch_pattern=" + pattern.name());
    System.out.println("baffle_type=" + baffle.name());
    System.out.println("shell_inside_diameter_mm=" + shellIdMm);
    System.out.println("tube_pitch_mm=" + pitchMm);
    System.out.println("tube_passes=" + tubePasses);
    System.out.println("tube_material=" + material);
    System.out.println("needs_cleanable=" + cleanable);
    System.out.println("large_temperature_difference=" + largeDt);
    System.out.println("high_pressure=" + highPressure);
    System.out.println("hazardous=" + hazardous);
    System.out.println("configuration_cost_factor=" + config.getCostFactor());
    System.out.println("shell_pressure_drop_factor="
        + config.getShell().getPressureDropFactor());
    System.out.println("bundle_removable=" + config.isBundleRemovable());
    System.out.println("good_thermal_expansion=" + config.hasGoodThermalExpansion());
    System.out.println("class_cost_factor=" + temaClass.getCostFactor());
    System.out.println("min_tube_wall_mm=" + temaClass.getMinTubeWallMm());
    System.out.println("min_corrosion_allowance_mm=" + temaClass.getMinCorrosionAllowanceMm());
    System.out.println("tube_outer_diameter_mm=" + tubeOdMm);
    System.out.println("tube_outer_diameter_inch=" + tubeSize.getOuterDiameterInch());
    System.out.println("layout_angle_deg=" + pattern.getLayoutAngle());
    System.out.println("min_pitch_ratio=" + pattern.getMinPitchRatio());
    System.out.println("pitch_heat_transfer_factor=" + pattern.getHeatTransferFactor());
    System.out.println("baffle_heat_transfer_factor=" + baffle.getHeatTransferFactor());
    System.out.println("baffle_pressure_drop_factor=" + baffle.getPressureDropFactor());
    System.out.println("min_tube_pitch_mm="
        + TEMAStandard.calculateMinTubePitch(tubeOdMm, pattern));
    System.out.println("estimated_tube_count="
        + TEMAStandard.estimateTubeCount(shellIdMm, tubeOdMm, pitchMm, pattern, tubePasses));
    System.out.println("max_unsupported_span_mm="
        + TEMAStandard.getMaxUnsupportedSpan(tubeOdMm, material));
    System.out.println("min_baffle_spacing_mm="
        + TEMAStandard.getMinBaffleSpacing(shellIdMm, temaClass));
    System.out.println("max_baffle_spacing_mm="
        + TEMAStandard.getMaxBaffleSpacing(shellIdMm));
    System.out.println("recommended_designation="
        + TEMAStandard.recommendConfiguration(cleanable, largeDt, highPressure, hazardous));
    System.out.println();
  }

  static void recommendationRow(String label, boolean cleanable, boolean largeDt,
      boolean highPressure, boolean hazardous) {
    System.out.println(label);
    System.out.println("needs_cleanable=" + cleanable);
    System.out.println("large_temperature_difference=" + largeDt);
    System.out.println("high_pressure=" + highPressure);
    System.out.println("hazardous=" + hazardous);
    System.out.println("recommended_designation="
        + TEMAStandard.recommendConfiguration(cleanable, largeDt, highPressure, hazardous));
    System.out.println();
  }

  static void countRow(String label, double shellIdMm, double tubeOdMm, double pitchMm,
      TubePitchPattern pattern, int tubePasses) {
    System.out.println(label);
    System.out.println("shell_inside_diameter_mm=" + shellIdMm);
    System.out.println("tube_outer_diameter_mm=" + tubeOdMm);
    System.out.println("tube_pitch_mm=" + pitchMm);
    System.out.println("layout_angle_deg=" + pattern.getLayoutAngle());
    System.out.println("tube_passes=" + tubePasses);
    System.out.println("estimated_tube_count="
        + TEMAStandard.estimateTubeCount(shellIdMm, tubeOdMm, pitchMm, pattern, tubePasses));
    System.out.println();
  }

  static void spanRow(String label, double tubeOdMm, String material) {
    System.out.println(label);
    System.out.println("tube_outer_diameter_mm=" + tubeOdMm);
    System.out.println("tube_material=" + material);
    System.out.println("max_unsupported_span_mm="
        + TEMAStandard.getMaxUnsupportedSpan(tubeOdMm, material));
    System.out.println();
  }
}
