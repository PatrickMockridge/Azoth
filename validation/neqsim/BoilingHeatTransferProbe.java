// `BoilingHeatTransfer` driven directly, on states of its own.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar BoilingHeatTransferProbe.java
//     java -cp .:neqsim-f0c7436.jar BoilingHeatTransferProbe > captures/boiling_heat_transfer_probe.tsv
//
// **The zone analysis's path is `calcAverageHTC`** (Gungor-Winterton, Simpson), and the
// Chen local, the two enhancement factors and the local Gungor-Winterton entries are the
// class's other public surface. The boiling quality clamp is `[0.001, 0.999]`, unlike the
// condensation correlation's `[0, 1]`.

import neqsim.process.mechanicaldesign.heatexchanger.BoilingHeatTransfer;

public class BoilingHeatTransferProbe {

  public static void main(String[] args) {
    // The class's test fluid, R134a-like at ~40 C and 10 bar.
    double liquidDensity = 1147.0;
    double vaporDensity = 50.1;
    double liquidViscosity = 1.65e-4;
    double vaporViscosity = 1.25e-5;
    double liquidCp = 1480.0;
    double liquidConductivity = 0.078;
    double surfaceTension = 0.006;
    double heatOfVaporization = 163000.0;
    double tubeId = 0.01905;
    double massFlux = 300.0;

    // 1-3. The Chen local coefficient: the test state, no wall superheat, and a zero flux.
    chenRow("chen_mid_quality", massFlux, 0.3, tubeId, liquidDensity, vaporDensity, liquidViscosity,
        vaporViscosity, liquidCp, liquidConductivity, surfaceTension, heatOfVaporization, 5.0, 50000.0);
    chenRow("chen_no_wall_superheat", massFlux, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 0.0, 0.0);
    chenRow("chen_zero_flux", 0.0, 0.3, tubeId, liquidDensity, vaporDensity, liquidViscosity,
        vaporViscosity, liquidCp, liquidConductivity, surfaceTension, heatOfVaporization, 5.0,
        50000.0);
    chenRow("chen_quality_zero", massFlux, 0.0, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 5.0, 50000.0);
    chenRow("chen_quality_one", massFlux, 1.0, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 5.0, 50000.0);

    // 4-7. The Gungor-Winterton local coefficient and its horizontal correction.
    gungorRow("gungor_mid_quality", massFlux, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization);
    gungorRow("gungor_low_quality", massFlux, 0.1, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization);
    gungorRow("gungor_high_quality", massFlux, 0.5, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization);
    gungorRow("gungor_quality_zero", massFlux, 0.0, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization);
    gungorRow("gungor_quality_one", massFlux, 1.0, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization);
    correctedRow("corrected_horizontal_low_flux", 5.0, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization, true);
    correctedRow("corrected_horizontal_high_flux", massFlux, 0.3, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization,
        true);
    correctedRow("corrected_vertical", massFlux, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization, false);

    // 8-10. The Martinelli parameter and the two Chen factors.
    martinelliRow("martinelli_mid_quality", 0.5, liquidDensity, vaporDensity, liquidViscosity,
        vaporViscosity);
    martinelliRow("martinelli_quality_zero", 0.0, liquidDensity, vaporDensity, liquidViscosity,
        vaporViscosity);
    martinelliRow("martinelli_quality_one", 1.0, liquidDensity, vaporDensity, liquidViscosity,
        vaporViscosity);
    enhancementRow("enhancement_xtt_half", 0.5);
    enhancementRow("enhancement_xtt_ten", 10.0);
    enhancementRow("enhancement_xtt_zero", 0.0);
    suppressionRow("suppression_re_50000", 50000.0);
    suppressionRow("suppression_re_zero", 0.0);

    // 11-12. The average over a quality range, including the zone analysis's 10 intervals.
    averageRow("average_high_to_low", massFlux, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization, 0.1, 0.9, 20);
    averageRow("average_low_to_high", massFlux, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization, 0.1, 0.9, 20);
    averageRow("average_equal_endpoints", massFlux, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization, 0.5, 0.5, 20);
    averageRow("average_odd_intervals", massFlux, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization, 0.1, 0.9, 7);
    averageRow("average_below_minimum", massFlux, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization, 0.1, 0.9, 2);
    averageRow("average_zone_intervals", massFlux, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, 20000.0, heatOfVaporization, 0.1, 0.9, 10);

    // 13. The composed state the calc publishes, at each branch the cases pin.
    composedRow("composed_test_fluid", massFlux, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 5.0, 50000.0, 20000.0, 0.1, 0.9, 20, false);
    composedRow("composed_no_wall_superheat", massFlux, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 0.0, 0.0, 20000.0, 0.1, 0.9, 20, false);
    composedRow("composed_quality_zero", massFlux, 0.0, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 5.0, 50000.0, 20000.0, 0.1, 0.9, 20, false);
    composedRow("composed_quality_one", massFlux, 1.0, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 5.0, 50000.0, 20000.0, 0.1, 0.9, 20, false);
    composedRow("composed_horizontal_low_flux", 5.0, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 5.0, 50000.0, 20000.0, 0.1, 0.9, 20, true);
    composedRow("composed_horizontal_high_flux", massFlux, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 5.0, 50000.0, 20000.0, 0.1, 0.9, 20, true);
    composedRow("composed_odd_intervals", massFlux, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 5.0, 50000.0, 20000.0, 0.1, 0.9, 7, false);
    composedRow("composed_zone_intervals", massFlux, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 5.0, 50000.0, 20000.0, 0.1, 0.9, 10, false);
    composedRow("composed_zero_heat_flux", massFlux, 0.3, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, liquidCp, liquidConductivity, surfaceTension,
        heatOfVaporization, 5.0, 50000.0, 0.0, 0.1, 0.9, 20, false);

    System.out.flush();
  }

  static void chenRow(String label, double massFlux, double quality, double tubeId,
      double liquidDensity, double vaporDensity, double liquidViscosity, double vaporViscosity,
      double liquidCp, double liquidConductivity, double surfaceTension, double heatOfVaporization,
      double wallSuperheat, double satPressureDiff) {
    double h = BoilingHeatTransfer.calcChenHTC(massFlux, quality, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, liquidCp, liquidConductivity,
        surfaceTension, heatOfVaporization, wallSuperheat, satPressureDiff);
    System.out.println(label);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("vapor_quality=" + quality);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("liquid_viscosity_pa_s=" + liquidViscosity);
    System.out.println("vapor_viscosity_pa_s=" + vaporViscosity);
    System.out.println("liquid_cp_j_per_kgk=" + liquidCp);
    System.out.println("liquid_conductivity_w_per_mk=" + liquidConductivity);
    System.out.println("surface_tension_n_per_m=" + surfaceTension);
    System.out.println("heat_of_vaporization_j_per_kg=" + heatOfVaporization);
    System.out.println("wall_superheat_k=" + wallSuperheat);
    System.out.println("sat_pressure_diff_pa=" + satPressureDiff);
    System.out.println("chen_htc_w_per_m2k=" + h);
    System.out.println();
  }

  static void gungorRow(String label, double massFlux, double quality, double tubeId,
      double liquidDensity, double vaporDensity, double liquidViscosity, double liquidCp,
      double liquidConductivity, double heatFlux, double heatOfVaporization) {
    double h = BoilingHeatTransfer.calcGungorWintertonHTC(massFlux, quality, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, liquidCp, liquidConductivity, heatFlux, heatOfVaporization);
    System.out.println(label);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("vapor_quality=" + quality);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("liquid_viscosity_pa_s=" + liquidViscosity);
    System.out.println("liquid_cp_j_per_kgk=" + liquidCp);
    System.out.println("liquid_conductivity_w_per_mk=" + liquidConductivity);
    System.out.println("heat_flux_w_per_m2=" + heatFlux);
    System.out.println("heat_of_vaporization_j_per_kg=" + heatOfVaporization);
    System.out.println("gungor_winterton_htc_w_per_m2k=" + h);
    System.out.println();
  }

  static void correctedRow(String label, double massFlux, double quality, double tubeId,
      double liquidDensity, double vaporDensity, double liquidViscosity, double liquidCp,
      double liquidConductivity, double heatFlux, double heatOfVaporization,
      boolean horizontal) {
    double h = BoilingHeatTransfer.calcGungorWintertonCorrectedHTC(massFlux, quality, tubeId,
        liquidDensity, vaporDensity, liquidViscosity, liquidCp, liquidConductivity, heatFlux,
        heatOfVaporization, horizontal);
    System.out.println(label);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("vapor_quality=" + quality);
    System.out.println("horizontal=" + horizontal);
    System.out.println("corrected_htc_w_per_m2k=" + h);
    System.out.println();
  }

  static void martinelliRow(String label, double quality, double liquidDensity,
      double vaporDensity, double liquidViscosity, double vaporViscosity) {
    double x = BoilingHeatTransfer.calcMartinelliParameter(quality, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity);
    System.out.println(label);
    System.out.println("vapor_quality=" + quality);
    System.out.println("martinelli_parameter=" + x);
    System.out.println();
  }

  static void enhancementRow(String label, double xtt) {
    double f = BoilingHeatTransfer.calcChenEnhancementFactor(xtt);
    System.out.println(label);
    System.out.println("martinelli_parameter=" + xtt);
    System.out.println("chen_enhancement_factor=" + f);
    System.out.println();
  }

  static void suppressionRow(String label, double reTp) {
    double s = BoilingHeatTransfer.calcChenSuppressionFactor(reTp);
    System.out.println(label);
    System.out.println("two_phase_reynolds=" + reTp);
    System.out.println("chen_suppression_factor=" + s);
    System.out.println();
  }

  static void composedRow(String label, double massFlux, double quality, double tubeId,
      double liquidDensity, double vaporDensity, double liquidViscosity, double vaporViscosity,
      double liquidCp, double liquidConductivity, double surfaceTension, double heatOfVaporization,
      double wallSuperheat, double satPressureDiff, double heatFlux, double qualityIn,
      double qualityOut, int intervals, boolean horizontal) {
    double chen = BoilingHeatTransfer.calcChenHTC(massFlux, quality, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, liquidCp, liquidConductivity,
        surfaceTension, heatOfVaporization, wallSuperheat, satPressureDiff);
    double gungor = BoilingHeatTransfer.calcGungorWintertonHTC(massFlux, quality, tubeId,
        liquidDensity, vaporDensity, liquidViscosity, liquidCp, liquidConductivity, heatFlux,
        heatOfVaporization);
    double corrected = BoilingHeatTransfer.calcGungorWintertonCorrectedHTC(massFlux, quality,
        tubeId, liquidDensity, vaporDensity, liquidViscosity, liquidCp, liquidConductivity,
        heatFlux, heatOfVaporization, horizontal);
    double average = BoilingHeatTransfer.calcAverageHTC(massFlux, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, liquidCp, liquidConductivity, heatFlux, heatOfVaporization,
        qualityIn, qualityOut, intervals);
    double xtt = BoilingHeatTransfer.calcMartinelliParameter(quality, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity);
    double enhancement = BoilingHeatTransfer.calcChenEnhancementFactor(xtt);
    double clamped = Math.max(0.001, Math.min(0.999, quality));
    double reLiquid = massFlux * (1.0 - clamped) * tubeId / liquidViscosity;
    double suppression = BoilingHeatTransfer.calcChenSuppressionFactor(
        reLiquid * Math.pow(enhancement, 1.25));
    System.out.println(label);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("vapor_quality=" + quality);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("liquid_viscosity_pa_s=" + liquidViscosity);
    System.out.println("vapor_viscosity_pa_s=" + vaporViscosity);
    System.out.println("liquid_cp_j_per_kgk=" + liquidCp);
    System.out.println("liquid_conductivity_w_per_mk=" + liquidConductivity);
    System.out.println("surface_tension_n_per_m=" + surfaceTension);
    System.out.println("heat_of_vaporization_j_per_kg=" + heatOfVaporization);
    System.out.println("wall_superheat_k=" + wallSuperheat);
    System.out.println("sat_pressure_diff_pa=" + satPressureDiff);
    System.out.println("heat_flux_w_per_m2=" + heatFlux);
    System.out.println("quality_in=" + qualityIn);
    System.out.println("quality_out=" + qualityOut);
    System.out.println("intervals=" + intervals);
    System.out.println("horizontal=" + horizontal);
    System.out.println("chen_htc_w_per_m2k=" + chen);
    System.out.println("gungor_winterton_htc_w_per_m2k=" + gungor);
    System.out.println("gungor_winterton_corrected_htc_w_per_m2k=" + corrected);
    System.out.println("average_htc_w_per_m2k=" + average);
    System.out.println("martinelli_parameter=" + xtt);
    System.out.println("chen_enhancement_factor=" + enhancement);
    System.out.println("chen_suppression_factor=" + suppression);
    System.out.println();
  }

  static void averageRow(String label, double massFlux, double tubeId, double liquidDensity,
      double vaporDensity, double liquidViscosity, double liquidCp, double liquidConductivity,
      double heatFlux, double heatOfVaporization, double qualityIn, double qualityOut,
      int intervals) {
    double h = BoilingHeatTransfer.calcAverageHTC(massFlux, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, liquidCp, liquidConductivity, heatFlux, heatOfVaporization, qualityIn,
        qualityOut, intervals);
    System.out.println(label);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("liquid_viscosity_pa_s=" + liquidViscosity);
    System.out.println("liquid_cp_j_per_kgk=" + liquidCp);
    System.out.println("liquid_conductivity_w_per_mk=" + liquidConductivity);
    System.out.println("heat_flux_w_per_m2=" + heatFlux);
    System.out.println("heat_of_vaporization_j_per_kg=" + heatOfVaporization);
    System.out.println("quality_in=" + qualityIn);
    System.out.println("quality_out=" + qualityOut);
    System.out.println("intervals=" + intervals);
    System.out.println("average_htc_w_per_m2k=" + h);
    System.out.println();
  }
}
