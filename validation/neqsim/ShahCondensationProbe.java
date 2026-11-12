// `ShahCondensation` driven directly, on states of its own.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-118d9af.jar ShahCondensationProbe.java
//     java -cp .:neqsim-118d9af.jar ShahCondensationProbe > captures/shah_condensation_probe.tsv
//
// **The composed screen is the zone analysis's pair**: `calcLiquidOnlyHTC` feeds
// `calcAverageHTC`, and the local and vertical coefficients are the class's other two
// public entries. The reduced-pressure clamps are the class's own.

import neqsim.process.mechanicaldesign.heatexchanger.ShahCondensation;

public class ShahCondensationProbe {

  public static void main(String[] args) {
    // The class's test fluid, R134a-like at ~40 C and 10 bar.
    double liquidDensity = 1147.0;
    double vaporDensity = 50.1;
    double liquidViscosity = 1.65e-4;
    double liquidCp = 1480.0;
    double liquidConductivity = 0.078;
    double tubeId = 0.01905;
    double massFlux = 300.0;

    // 1. The liquid-only coefficient, the class's own test state.
    liquidOnlyRow("liquid_only_htc", massFlux, tubeId, liquidDensity, liquidViscosity, liquidCp,
        liquidConductivity);

    // 2-4. The local coefficient at mid, zero and high quality, and its clamps.
    localRow("local_mid_quality", 3000.0, 0.5, 0.1);
    localRow("local_zero_quality", 3000.0, 0.0, 0.1);
    localRow("local_low_quality", 3000.0, 0.2, 0.1);
    localRow("local_high_quality", 3000.0, 0.8, 0.1);
    localRow("local_quality_above_one", 3000.0, 1.4, 0.1);
    localRow("local_quality_below_zero", 3000.0, -0.4, 0.1);
    localRow("local_reduced_pressure_below_floor", 3000.0, 0.5, 1e-4);
    localRow("local_reduced_pressure_above_ceiling", 3000.0, 0.5, 0.9);
    localRow("local_negative_liquid_only", -1.0, 0.5, 0.1);
    localRow("local_negative_pressure", 3000.0, 0.5, -0.1);

    // 5. The average coefficient over a quality range, and its degenerate forms.
    averageRow("average_high_to_low", 3000.0, 0.1, 0.9, 0.1, 20);
    averageRow("average_low_to_high", 3000.0, 0.1, 0.1, 0.9, 20);
    averageRow("average_equal_endpoints", 3000.0, 0.1, 0.5, 0.5, 20);
    averageRow("average_odd_intervals", 3000.0, 0.1, 0.9, 0.1, 7);
    averageRow("average_below_minimum", 3000.0, 0.1, 0.9, 0.1, 2);

    // 6. The vertical-tube coefficient: shear-dominated, gravity-dominated, and a liquid
    //    lighter than the vapour (the class's fallback).
    verticalRow("vertical_shear", 3000.0, 0.5, 0.1, 300.0, 0.01905, 50.1, 1147.0);
    verticalRow("vertical_gravity", 3000.0, 0.1, 0.1, 5.0, 0.01905, 50.1, 1147.0);
    verticalRow("vertical_liquid_lighter", 3000.0, 0.5, 0.1, 300.0, 0.01905, 1200.0, 1147.0);

    // 7. The composed state the zone analysis takes: the liquid-only coefficient feeds
    //    the local, average and vertical entries.
    double hLo = ShahCondensation.calcLiquidOnlyHTC(massFlux, tubeId, liquidDensity,
        liquidViscosity, liquidCp, liquidConductivity);
    double reLo = massFlux * tubeId / liquidViscosity;
    System.out.println("composed_test_fluid");
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("liquid_viscosity_pa_s=" + liquidViscosity);
    System.out.println("liquid_cp_j_per_kgk=" + liquidCp);
    System.out.println("liquid_conductivity_w_per_mk=" + liquidConductivity);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("reduced_pressure=0.1");
    System.out.println("vapor_quality=0.5");
    System.out.println("quality_in=0.9");
    System.out.println("quality_out=0.1");
    System.out.println("intervals=20");
    System.out.println("liquid_only_htc_w_per_m2k=" + hLo);
    System.out.println("liquid_only_reynolds=" + reLo);
    System.out.println("local_htc_w_per_m2k="
        + ShahCondensation.calcLocalHTC(hLo, 0.5, 0.1));
    System.out.println("average_htc_w_per_m2k="
        + ShahCondensation.calcAverageHTC(hLo, 0.1, 0.9, 0.1, 20));
    System.out.println("vertical_htc_w_per_m2k="
        + ShahCondensation.calcVerticalTubeHTC(hLo, 0.5, 0.1, massFlux, tubeId, vaporDensity,
            liquidDensity));
    System.out.println("in_valid_range="
        + ShahCondensation.isInValidRange(0.1, reLo, massFlux));
    System.out.println();

    // 8. The composed state at each branch the calc's cases pin.
    composedRow("composed_quality_one", massFlux, tubeId, liquidDensity, liquidViscosity, liquidCp,
        liquidConductivity, vaporDensity, 0.1, 1.0, 0.9, 0.1, 20);
    composedRow("composed_quality_below_zero", massFlux, tubeId, liquidDensity, liquidViscosity,
        liquidCp, liquidConductivity, vaporDensity, 0.1, -0.4, 0.9, 0.1, 20);
    composedRow("composed_pressure_floor", massFlux, tubeId, liquidDensity, liquidViscosity,
        liquidCp, liquidConductivity, vaporDensity, 1e-4, 0.5, 0.9, 0.1, 20);
    composedRow("composed_pressure_ceiling", massFlux, tubeId, liquidDensity, liquidViscosity,
        liquidCp, liquidConductivity, vaporDensity, 0.9, 0.5, 0.9, 0.1, 20);
    composedRow("composed_odd_intervals", massFlux, tubeId, liquidDensity, liquidViscosity,
        liquidCp, liquidConductivity, vaporDensity, 0.1, 0.5, 0.9, 0.1, 7);
    composedRow("composed_gravity", 5.0, tubeId, liquidDensity, liquidViscosity, liquidCp,
        liquidConductivity, vaporDensity, 0.1, 0.1, 0.9, 0.1, 20);
    composedRow("composed_liquid_lighter", massFlux, tubeId, liquidDensity, liquidViscosity,
        liquidCp, liquidConductivity, 1200.0, 0.1, 0.5, 0.9, 0.1, 20);
    composedRow("composed_low_flux", 3.0, tubeId, liquidDensity, liquidViscosity, liquidCp,
        liquidConductivity, vaporDensity, 0.1, 0.5, 0.9, 0.1, 20);
    composedRow("composed_high_flux", 900.0, tubeId, liquidDensity, liquidViscosity, liquidCp,
        liquidConductivity, vaporDensity, 0.1, 0.5, 0.9, 0.1, 20);
    composedRow("composed_low_reynolds", 10.0, tubeId, liquidDensity, 1e-3, liquidCp,
        liquidConductivity, vaporDensity, 0.1, 0.5, 0.9, 0.1, 20);

    // 9. The validity predicate, valid and one state per boundary.
    rangeRow("range_valid", 0.1, 1000.0, 300.0);
    rangeRow("range_pressure_low", 1e-4, 1000.0, 300.0);
    rangeRow("range_pressure_high", 0.5, 1000.0, 300.0);
    rangeRow("range_re_low", 0.1, 300.0, 300.0);
    rangeRow("range_flux_low", 0.1, 1000.0, 3.0);
    rangeRow("range_flux_high", 0.1, 1000.0, 900.0);

    System.out.flush();
  }

  static void liquidOnlyRow(String label, double massFlux, double tubeId, double liquidDensity,
      double liquidViscosity, double liquidCp, double liquidConductivity) {
    double h = ShahCondensation.calcLiquidOnlyHTC(massFlux, tubeId, liquidDensity, liquidViscosity,
        liquidCp, liquidConductivity);
    System.out.println(label);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("liquid_viscosity_pa_s=" + liquidViscosity);
    System.out.println("liquid_cp_j_per_kgk=" + liquidCp);
    System.out.println("liquid_conductivity_w_per_mk=" + liquidConductivity);
    System.out.println("liquid_only_htc_w_per_m2k=" + h);
    System.out.println();
  }

  static void localRow(String label, double hLiquidOnly, double quality, double reducedPressure) {
    double h = ShahCondensation.calcLocalHTC(hLiquidOnly, quality, reducedPressure);
    System.out.println(label);
    System.out.println("liquid_only_htc_w_per_m2k=" + hLiquidOnly);
    System.out.println("vapor_quality=" + quality);
    System.out.println("reduced_pressure=" + reducedPressure);
    System.out.println("local_htc_w_per_m2k=" + h);
    System.out.println();
  }

  static void averageRow(String label, double hLiquidOnly, double reducedPressure, double qualityIn,
      double qualityOut, int intervals) {
    double h = ShahCondensation.calcAverageHTC(hLiquidOnly, reducedPressure, qualityIn, qualityOut,
        intervals);
    System.out.println(label);
    System.out.println("liquid_only_htc_w_per_m2k=" + hLiquidOnly);
    System.out.println("reduced_pressure=" + reducedPressure);
    System.out.println("quality_in=" + qualityIn);
    System.out.println("quality_out=" + qualityOut);
    System.out.println("intervals=" + intervals);
    System.out.println("average_htc_w_per_m2k=" + h);
    System.out.println();
  }

  static void verticalRow(String label, double hLiquidOnly, double quality, double reducedPressure,
      double massFlux, double tubeId, double vaporDensity, double liquidDensity) {
    double h = ShahCondensation.calcVerticalTubeHTC(hLiquidOnly, quality, reducedPressure, massFlux,
        tubeId, vaporDensity, liquidDensity);
    System.out.println(label);
    System.out.println("liquid_only_htc_w_per_m2k=" + hLiquidOnly);
    System.out.println("vapor_quality=" + quality);
    System.out.println("reduced_pressure=" + reducedPressure);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("vertical_htc_w_per_m2k=" + h);
    System.out.println();
  }

  static void composedRow(String label, double massFlux, double tubeId, double liquidDensity,
      double liquidViscosity, double liquidCp, double liquidConductivity, double vaporDensity,
      double reducedPressure, double quality, double qualityIn, double qualityOut,
      int intervals) {
    double hLo = ShahCondensation.calcLiquidOnlyHTC(massFlux, tubeId, liquidDensity,
        liquidViscosity, liquidCp, liquidConductivity);
    double reLo = massFlux * tubeId / liquidViscosity;
    System.out.println(label);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("liquid_viscosity_pa_s=" + liquidViscosity);
    System.out.println("liquid_cp_j_per_kgk=" + liquidCp);
    System.out.println("liquid_conductivity_w_per_mk=" + liquidConductivity);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("reduced_pressure=" + reducedPressure);
    System.out.println("vapor_quality=" + quality);
    System.out.println("quality_in=" + qualityIn);
    System.out.println("quality_out=" + qualityOut);
    System.out.println("intervals=" + intervals);
    System.out.println("liquid_only_htc_w_per_m2k=" + hLo);
    System.out.println("liquid_only_reynolds=" + reLo);
    System.out.println("local_htc_w_per_m2k="
        + ShahCondensation.calcLocalHTC(hLo, quality, reducedPressure));
    System.out.println("average_htc_w_per_m2k="
        + ShahCondensation.calcAverageHTC(hLo, reducedPressure, qualityIn, qualityOut, intervals));
    System.out.println("vertical_htc_w_per_m2k="
        + ShahCondensation.calcVerticalTubeHTC(hLo, quality, reducedPressure, massFlux, tubeId,
            vaporDensity, liquidDensity));
    System.out.println("in_valid_range="
        + ShahCondensation.isInValidRange(reducedPressure, reLo, massFlux));
    System.out.println();
  }

  static void rangeRow(String label, double reducedPressure, double reLiquidOnly,
      double massFlux) {
    System.out.println(label);
    System.out.println("reduced_pressure=" + reducedPressure);
    System.out.println("liquid_only_reynolds=" + reLiquidOnly);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("in_valid_range=" + ShahCondensation.isInValidRange(reducedPressure,
        reLiquidOnly, massFlux));
    System.out.println();
  }
}
