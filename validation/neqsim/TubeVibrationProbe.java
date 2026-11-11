// `VibrationAnalysis` driven directly, on states of its own.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar TubeVibrationProbe.java
//     java -cp .:neqsim-f0c7436.jar TubeVibrationProbe > captures/tube_vibration_probe.tsv
//
// **The screening is the class's own `performScreening`, and the four helpers are its
// intermediates.** The natural frequency fixes the end condition at "pinned" in the
// screening; the other two conditions are measured here as evidence rather than as cases.
// The class that calls the screening is `ShellAndTubeDesignCalculator`, so this is the
// first slice of the detailed design tree.

import neqsim.process.mechanicaldesign.heatexchanger.VibrationAnalysis;

public class TubeVibrationProbe {

  public static void main(String[] args) {
    // 1-3. The natural frequency on the class's own test state, one per end condition.
    naturalRow("natural_pinned", "pinned");
    naturalRow("natural_fixed", "fixed");
    naturalRow("natural_clamped_pinned", "clamped-pinned");

    // 4. Vortex shedding at the class's test state: St ~ 0.22, f = 0.22 * 2 / 0.01905.
    System.out.println("vortex_shedding");
    double fvs = VibrationAnalysis.calcVortexSheddingFrequency(2.0, 0.01905, 0.02381);
    System.out.println("tube_od_m=0.01905");
    System.out.println("tube_pitch_m=0.02381");
    System.out.println("crossflow_velocity_m_per_s=2.0");
    System.out.println("vortex_shedding_frequency_hz=" + fvs);
    System.out.println("derived_strouhal=" + (2.0 > 0 ? fvs * 0.01905 / 2.0 : 0.0));
    System.out.println();

    // 5. Connors' critical velocity, both layouts.
    connorsRow("connors_triangular", true);
    connorsRow("connors_square", false);

    // 6. The acoustic frequency, the first three modes.
    acousticRow("acoustic_mode_1", 1);
    acousticRow("acoustic_mode_2", 2);
    acousticRow("acoustic_mode_3", 3);

    // 7. The effective acoustic velocity helper, whose only caller is itself.
    System.out.println("effective_acoustic_velocity");
    System.out.println("sonic_velocity_m_per_s=340.0");
    System.out.println("tube_volume_fraction=0.3");
    System.out.println(
        "effective_acoustic_velocity_m_per_s="
            + VibrationAnalysis.calcEffectiveAcousticVelocity(340.0, 0.3));
    System.out.println();

    // 8. The class's own pass state.
    screeningRow("screening_pass", 0.01905, 0.01483, 0.4, 0.02381, 200e9, 7800.0, 0.5, 50.0,
        800.0, 0.5, 340.0, 0.03, true);

    // 9. A crossflow above the reduced critical velocity: fluid-elastic instability.
    screeningRow("screening_fluid_elastic_critical", 0.01905, 0.01483, 0.4, 0.02381, 200e9,
        7800.0, 29.4, 50.0, 800.0, 0.5, 340.0, 0.03, true);

    // 10. The vortex-shedding frequency inside 0.8-1.2 of the natural frequency, and the
    //     acoustic ratio outside its band.
    screeningRow("screening_vortex_critical", 0.01905, 0.01483, 0.4, 0.02381, 200e9, 7800.0,
        15.0, 50.0, 800.0, 0.5, 340.0, 0.03, true);

    // 11. The vortex-shedding frequency inside 0.8-1.2 of the first acoustic mode, with the
    //     fluid-elastic ratio below one: St * V / d ~ c / (2 * shell_id).
    screeningRow("screening_acoustic_critical", 0.01905, 0.01483, 0.4, 0.02381, 200e9, 7800.0,
        20.0, 50.0, 800.0, 0.5, 340.0, 0.03, true);

    // 12. Degenerate geometry: every helper's guard, and the screening's MAX_VALUE sentinel.
    degenerateRow("degenerate_zero_diameter", 0.0, 0.01483, 0.4);
    degenerateRow("degenerate_id_above_od", 0.01905, 0.02, 0.4);
    degenerateRow("degenerate_zero_span", 0.01905, 0.01483, 0.0);

    System.out.flush();
  }

  static void naturalRow(String label, String endCondition) {
    double fn = VibrationAnalysis.calcNaturalFrequency(0.01905, 0.01483, 0.5, 200e9, 7800.0, 800.0,
        50.0, endCondition);
    System.out.println(label);
    System.out.println("tube_od_m=0.01905");
    System.out.println("tube_id_m=0.01483");
    System.out.println("unsupported_span_m=0.5");
    System.out.println("tube_material_e_pa=2.0E11");
    System.out.println("tube_density_kg_per_m3=7800.0");
    System.out.println("tube_fluid_density_kg_per_m3=800.0");
    System.out.println("shell_fluid_density_kg_per_m3=50.0");
    System.out.println("end_condition=" + endCondition);
    System.out.println("natural_frequency_hz=" + fn);
    System.out.println();
  }

  static void connorsRow(String label, boolean triangular) {
    double vCrit = VibrationAnalysis.calcCriticalVelocityConnors(100.0, 0.01905, 0.03, 1.5, 50.0,
        triangular);
    System.out.println(label);
    System.out.println("natural_frequency_hz=100.0");
    System.out.println("tube_od_m=0.01905");
    System.out.println("damping_ratio=0.03");
    System.out.println("effective_mass_per_length_kg_per_m=1.5");
    System.out.println("shell_fluid_density_kg_per_m3=50.0");
    System.out.println("triangular_pitch=" + triangular);
    System.out.println("critical_velocity_m_per_s=" + vCrit);
    System.out.println();
  }

  static void acousticRow(String label, int mode) {
    double fac = VibrationAnalysis.calcAcousticFrequency(0.5, 340.0, mode);
    System.out.println(label);
    System.out.println("shell_id_m=0.5");
    System.out.println("sonic_velocity_m_per_s=340.0");
    System.out.println("mode=" + mode);
    System.out.println("acoustic_frequency_hz=" + fac);
    System.out.println();
  }

  static void screeningRow(String label, double tubeOd, double tubeId, double span,
      double tubePitch, double e, double tubeDensity, double crossflowVelocity,
      double shellFluidDensity, double tubeFluidDensity, double shellId, double sonicVelocity,
      double damping, boolean triangular) {
    VibrationAnalysis.VibrationResult result = VibrationAnalysis.performScreening(tubeOd, tubeId, span,
        tubePitch, e, tubeDensity, crossflowVelocity, shellFluidDensity, tubeFluidDensity, shellId,
        sonicVelocity, damping, triangular);
    System.out.println(label);
    System.out.println("tube_od_m=" + tubeOd);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("unsupported_span_m=" + span);
    System.out.println("tube_pitch_m=" + tubePitch);
    System.out.println("tube_material_e_pa=" + e);
    System.out.println("tube_density_kg_per_m3=" + tubeDensity);
    System.out.println("crossflow_velocity_m_per_s=" + crossflowVelocity);
    System.out.println("shell_fluid_density_kg_per_m3=" + shellFluidDensity);
    System.out.println("tube_fluid_density_kg_per_m3=" + tubeFluidDensity);
    System.out.println("shell_id_m=" + shellId);
    System.out.println("sonic_velocity_m_per_s=" + sonicVelocity);
    System.out.println("damping_ratio=" + damping);
    System.out.println("triangular_pitch=" + triangular);
    printResult(result);
    System.out.println();
  }

  static void degenerateRow(String label, double tubeOd, double tubeId, double span) {
    VibrationAnalysis.VibrationResult result = VibrationAnalysis.performScreening(tubeOd, tubeId, span,
        0.02381, 200e9, 7800.0, 0.5, 50.0, 800.0, 0.5, 340.0, 0.03, true);
    System.out.println(label);
    System.out.println("tube_od_m=" + tubeOd);
    System.out.println("tube_id_m=" + tubeId);
    System.out.println("unsupported_span_m=" + span);
    printResult(result);
    System.out.println();
  }

  static void printResult(VibrationAnalysis.VibrationResult result) {
    System.out.println("natural_frequency_hz=" + result.naturalFrequencyHz);
    System.out.println("vortex_shedding_frequency_hz=" + result.vortexSheddingFrequencyHz);
    System.out.println("vortex_shedding_ratio=" + result.vortexSheddingRatio);
    System.out.println("vortex_shedding_critical=" + result.vortexSheddingCritical);
    System.out.println("critical_velocity_m_per_s=" + result.criticalVelocityMs);
    System.out.println("velocity_ratio=" + result.velocityRatio);
    System.out.println("fluid_elastic_critical=" + result.fluidElasticCritical);
    System.out.println("acoustic_frequency_hz=" + result.acousticFrequencyHz);
    System.out.println("acoustic_ratio=" + result.acousticRatio);
    System.out.println("acoustic_critical=" + result.acousticCritical);
    System.out.println("passed=" + result.passed);
  }
}
