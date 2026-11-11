// `TwoPhasePressureDrop` driven directly, on states of its own.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar TwoPhasePressureDropProbe.java
//     java -cp .:neqsim-f0c7436.jar TwoPhasePressureDropProbe > captures/two_phase_pressure_drop_probe.tsv
//
// **The zone analysis takes `calcFriedelPressureDrop` at an average quality**; the gradient,
// the average over a quality range, the Muller-Steinhagen-Heck gradient, the gravitational
// term and the acceleration term are the class's other public entries.

import neqsim.process.mechanicaldesign.heatexchanger.TwoPhasePressureDrop;

public class TwoPhasePressureDropProbe {

  public static void main(String[] args) {
    // The class's test fluid, R134a-like at ~40 C and 10 bar.
    double liquidDensity = 1147.0;
    double vaporDensity = 50.1;
    double liquidViscosity = 1.65e-4;
    double vaporViscosity = 1.25e-5;
    double surfaceTension = 0.006;
    double tubeId = 0.01905;
    double massFlux = 300.0;

    // 1-4. The Friedel gradient across quality, at zero surface tension and with no flux.
    gradientRow("friedel_gradient_mid_quality", massFlux, 0.5, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension);
    gradientRow("friedel_gradient_quality_zero", massFlux, 0.0, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension);
    gradientRow("friedel_gradient_low_quality", massFlux, 0.1, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension);
    gradientRow("friedel_gradient_high_quality", massFlux, 0.8, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension);
    gradientRow("friedel_gradient_quality_one", massFlux, 1.0, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension);
    gradientRow("friedel_gradient_no_surface_tension", massFlux, 0.5, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, 0.0);
    gradientRow("friedel_gradient_zero_flux", 0.0, 0.5, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, surfaceTension);

    // 5-6. The total drop and the average over a quality range.
    dropRow("friedel_drop_length_5", massFlux, 0.5, tubeId, 5.0, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, surfaceTension);
    dropRow("friedel_drop_length_10", massFlux, 0.5, tubeId, 10.0, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, surfaceTension);
    averageRow("friedel_average_high_to_low", massFlux, 0.1, 0.9, tubeId, 5.0, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension, 20);
    averageRow("friedel_average_low_to_high", massFlux, 0.9, 0.1, tubeId, 5.0, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension, 20);
    averageRow("friedel_average_equal_endpoints", massFlux, 0.5, 0.5, tubeId, 5.0, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension, 20);
    averageRow("friedel_average_odd_intervals", massFlux, 0.1, 0.9, tubeId, 5.0, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension, 7);
    averageRow("friedel_average_below_minimum", massFlux, 0.1, 0.9, tubeId, 5.0, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension, 2);

    // 7. The Muller-Steinhagen-Heck gradient.
    mshRow("msh_gradient_mid_quality", massFlux, 0.5, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity);
    mshRow("msh_gradient_quality_zero", massFlux, 0.0, tubeId, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity);

    // 8-9. The gravitational and acceleration terms.
    gravityRow("gravitational_mid_quality", 0.3, liquidDensity, vaporDensity);
    gravityRow("gravitational_quality_zero", 0.0, liquidDensity, vaporDensity);
    gravityRow("gravitational_quality_one", 1.0, liquidDensity, vaporDensity);
    accelerationRow("acceleration_up", massFlux, 0.1, 0.9, liquidDensity, vaporDensity);
    accelerationRow("acceleration_down", massFlux, 0.9, 0.1, liquidDensity, vaporDensity);

    // 10. The composed state the calc publishes, at each branch the cases pin.
    composedRow("composed_test_fluid", massFlux, 0.5, tubeId, 5.0, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, surfaceTension, 0.1, 0.9, 20);
    composedRow("composed_quality_zero", massFlux, 0.0, tubeId, 5.0, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, surfaceTension, 0.1, 0.9, 20);
    composedRow("composed_no_surface_tension", massFlux, 0.5, tubeId, 5.0, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, 0.0, 0.1, 0.9, 20);
    composedRow("composed_odd_intervals", massFlux, 0.5, tubeId, 5.0, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, surfaceTension, 0.1, 0.9, 7);
    composedRow("composed_zero_flux", 0.0, 0.5, tubeId, 5.0, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, surfaceTension, 0.1, 0.9, 20);
    composedRow("composed_quality_one", massFlux, 1.0, tubeId, 5.0, liquidDensity, vaporDensity,
        liquidViscosity, vaporViscosity, surfaceTension, 0.1, 0.9, 20);

    System.out.flush();
  }

  static void gradientRow(String label, double massFlux, double quality, double tubeId,
      double liquidDensity, double vaporDensity, double liquidViscosity, double vaporViscosity,
      double surfaceTension) {
    double g = TwoPhasePressureDrop.calcFriedelGradient(massFlux, quality, tubeId, liquidDensity,
        vaporDensity, liquidViscosity, vaporViscosity, surfaceTension);
    System.out.println(label);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("vapor_quality=" + quality);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("liquid_viscosity_pa_s=" + liquidViscosity);
    System.out.println("vapor_viscosity_pa_s=" + vaporViscosity);
    System.out.println("surface_tension_n_per_m=" + surfaceTension);
    System.out.println("friedel_gradient_pa_per_m=" + g);
    System.out.println();
  }

  static void dropRow(String label, double massFlux, double quality, double tubeId,
      double tubeLength, double liquidDensity, double vaporDensity, double liquidViscosity,
      double vaporViscosity, double surfaceTension) {
    double dp = TwoPhasePressureDrop.calcFriedelPressureDrop(massFlux, quality, tubeId, tubeLength,
        liquidDensity, vaporDensity, liquidViscosity, vaporViscosity, surfaceTension);
    System.out.println(label);
    System.out.println("tube_length_m=" + tubeLength);
    System.out.println("friedel_pressure_drop_pa=" + dp);
    System.out.println();
  }

  static void averageRow(String label, double massFlux, double qualityIn, double qualityOut,
      double tubeId, double tubeLength, double liquidDensity, double vaporDensity,
      double liquidViscosity, double vaporViscosity, double surfaceTension, int intervals) {
    double dp = TwoPhasePressureDrop.calcFriedelAveragePressureDrop(massFlux, qualityIn, qualityOut,
        tubeId, tubeLength, liquidDensity, vaporDensity, liquidViscosity, vaporViscosity,
        surfaceTension, intervals);
    System.out.println(label);
    System.out.println("quality_in=" + qualityIn);
    System.out.println("quality_out=" + qualityOut);
    System.out.println("intervals=" + intervals);
    System.out.println("friedel_average_pressure_drop_pa=" + dp);
    System.out.println();
  }

  static void mshRow(String label, double massFlux, double quality, double tubeId,
      double liquidDensity, double vaporDensity, double liquidViscosity, double vaporViscosity) {
    double g = TwoPhasePressureDrop.calcMullerSteinhagenHeckGradient(massFlux, quality, tubeId,
        liquidDensity, vaporDensity, liquidViscosity, vaporViscosity);
    System.out.println(label);
    System.out.println("vapor_quality=" + quality);
    System.out.println("muller_steinhagen_heck_gradient_pa_per_m=" + g);
    System.out.println();
  }

  static void gravityRow(String label, double quality, double liquidDensity,
      double vaporDensity) {
    double g = TwoPhasePressureDrop.calcGravitationalGradient(quality, liquidDensity,
        vaporDensity);
    System.out.println(label);
    System.out.println("vapor_quality=" + quality);
    System.out.println("gravitational_gradient_pa_per_m=" + g);
    System.out.println();
  }

  static void accelerationRow(String label, double massFlux, double qualityIn, double qualityOut,
      double liquidDensity, double vaporDensity) {
    double dp = TwoPhasePressureDrop.calcAccelerationPressureDrop(massFlux, qualityIn, qualityOut,
        liquidDensity, vaporDensity);
    System.out.println(label);
    System.out.println("quality_in=" + qualityIn);
    System.out.println("quality_out=" + qualityOut);
    System.out.println("acceleration_pressure_drop_pa=" + dp);
    System.out.println();
  }

  static void composedRow(String label, double massFlux, double quality, double tubeId,
      double tubeLength, double liquidDensity, double vaporDensity, double liquidViscosity,
      double vaporViscosity, double surfaceTension, double qualityIn, double qualityOut,
      int intervals) {
    double gradient = TwoPhasePressureDrop.calcFriedelGradient(massFlux, quality, tubeId,
        liquidDensity, vaporDensity, liquidViscosity, vaporViscosity, surfaceTension);
    double drop = TwoPhasePressureDrop.calcFriedelPressureDrop(massFlux, quality, tubeId,
        tubeLength, liquidDensity, vaporDensity, liquidViscosity, vaporViscosity, surfaceTension);
    double average = TwoPhasePressureDrop.calcFriedelAveragePressureDrop(massFlux, qualityIn,
        qualityOut, tubeId, tubeLength, liquidDensity, vaporDensity, liquidViscosity,
        vaporViscosity, surfaceTension, intervals);
    double msh = TwoPhasePressureDrop.calcMullerSteinhagenHeckGradient(massFlux, quality, tubeId,
        liquidDensity, vaporDensity, liquidViscosity, vaporViscosity);
    double gravity = TwoPhasePressureDrop.calcGravitationalGradient(quality, liquidDensity,
        vaporDensity);
    double acceleration = TwoPhasePressureDrop.calcAccelerationPressureDrop(massFlux, qualityIn,
        qualityOut, liquidDensity, vaporDensity);
    System.out.println(label);
    System.out.println("mass_flux_kg_per_m2s=" + massFlux);
    System.out.println("vapor_quality=" + quality);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("tube_length_m=" + tubeLength);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("liquid_viscosity_pa_s=" + liquidViscosity);
    System.out.println("vapor_viscosity_pa_s=" + vaporViscosity);
    System.out.println("surface_tension_n_per_m=" + surfaceTension);
    System.out.println("quality_in=" + qualityIn);
    System.out.println("quality_out=" + qualityOut);
    System.out.println("intervals=" + intervals);
    System.out.println("friedel_gradient_pa_per_m=" + gradient);
    System.out.println("friedel_pressure_drop_pa=" + drop);
    System.out.println("friedel_average_pressure_drop_pa=" + average);
    System.out.println("muller_steinhagen_heck_gradient_pa_per_m=" + msh);
    System.out.println("gravitational_gradient_pa_per_m=" + gravity);
    System.out.println("acceleration_pressure_drop_pa=" + acceleration);
    System.out.println();
  }
}
