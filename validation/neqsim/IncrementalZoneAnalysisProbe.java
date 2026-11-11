// `IncrementalZoneAnalysis` driven directly, on its default geometry and the phase regimes it
// dispatches over.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar IncrementalZoneAnalysisProbe.java
//     java -cp .:neqsim-f0c7436.jar IncrementalZoneAnalysisProbe \
//         > captures/incremental_zone_analysis_probe.tsv
//
// The class is a single pass over the zones in insertion order: no subdivision, no iteration,
// no convergence test. Condensing and evaporating zones use Shah and Gungor-Winterton on the
// tube side and Friedel for the drop; every other regime - including TWO_PHASE - falls through
// to the single-phase path. The shell side is always Kern and ignores its phase regime. The
// probe prints every geometry and fouling value, the constructor's target count, every zone's
// inputs and every computed field, so a row is reconstructible without the Java source.

import java.util.List;
import java.util.Map;

import neqsim.process.mechanicaldesign.heatexchanger.IncrementalZoneAnalysis;
import neqsim.process.mechanicaldesign.heatexchanger.IncrementalZoneAnalysis.IncrementalZone;
import neqsim.process.mechanicaldesign.heatexchanger.IncrementalZoneAnalysis.PhaseRegime;

public class IncrementalZoneAnalysisProbe {

  public static void main(String[] args) {
    // 1. One liquid zone on the class's defaults, with a stated temperature set.
    IncrementalZone liquid = zone("liquid", 150000.0, 400.0, 380.0, 300.0, 320.0);
    row("single_liquid_default", new IncrementalZoneAnalysis(), liquid);

    // 2. Two zones, a liquid and a vapour.
    IncrementalZoneAnalysis two = new IncrementalZoneAnalysis();
    two.addZone(liquid);
    IncrementalZone vapor = zone("vapor", 120000.0, 380.0, 360.0, 300.0, 315.0);
    vapor.setTubeSideProperties(5.0, 1e-5, 2100.0, 0.025, 10.0, PhaseRegime.VAPOR);
    two.addZone(vapor);
    row("two_zones_liquid_vapor", two, null);

    // 3. The condensing / evaporating / vapour cycle.
    IncrementalZoneAnalysis cycle = new IncrementalZoneAnalysis(20);
    IncrementalZone condensing = zone("condensing", 200000.0, 400.0, 370.0, 300.0, 340.0);
    condensing.setTubeSideProperties(20.0, 1.2e-5, 2100.0, 0.05, 10.0, PhaseRegime.CONDENSING);
    condensing.tubeQualityIn = 1.0;
    condensing.tubeQualityOut = 0.0;
    condensing.tubeReducedPressure = 0.3;
    condensing.tubeSurfaceTension = 0.01;
    cycle.addZone(condensing);
    IncrementalZone evaporating = zone("evaporating", 250000.0, 370.0, 330.0, 300.0, 320.0);
    evaporating.setTubeSideProperties(20.0, 1.2e-5, 2100.0, 0.05, 10.0, PhaseRegime.EVAPORATING);
    evaporating.tubeQualityIn = 0.0;
    evaporating.tubeQualityOut = 1.0;
    evaporating.tubeReducedPressure = 0.2;
    evaporating.tubeSurfaceTension = 0.01;
    evaporating.tubeHeatOfVaporization = 250000.0;
    cycle.addZone(evaporating);
    cycle.addZone(vapor);
    row("three_zones_cycle", cycle, null);

    // 4. TWO_PHASE falls through to single-phase; the shell regime is ignored.
    IncrementalZoneAnalysis mixed = new IncrementalZoneAnalysis();
    IncrementalZone twoPhase =
        zone("two_phase", 180000.0, 400.0, 370.0, 300.0, 330.0);
    twoPhase.setTubeSideProperties(800.0, 5e-4, 2500.0, 0.13, 10.0, PhaseRegime.TWO_PHASE);
    twoPhase.setShellSideProperties(1000.0, 8e-4, 4180.0, 0.6, 20.0, PhaseRegime.EVAPORATING);
    mixed.addZone(twoPhase);
    row("two_phase_falls_through", mixed, null);

    // 5. The three LMTD branches: a cross, a near-equal pair, and a true log mean.
    IncrementalZoneAnalysis lmtd = new IncrementalZoneAnalysis();
    IncrementalZone cross = zone("temperature_cross", 100000.0, 400.0, 340.0, 350.0, 390.0);
    IncrementalZone nearEqual = zone("near_equal", 100000.0, 500.0, 450.0, 400.0, 450.0);
    IncrementalZone logMean = zone("log_mean", 100000.0, 500.0, 450.0, 400.0, 420.0);
    lmtd.addZone(cross);
    lmtd.addZone(nearEqual);
    lmtd.addZone(logMean);
    row("lmtd_branches", lmtd, null);

    // 6. A condensing zone whose quality range is a single point.
    IncrementalZoneAnalysis point = new IncrementalZoneAnalysis();
    IncrementalZone single = zone("single_point", 120000.0, 400.0, 370.0, 300.0, 330.0);
    single.setTubeSideProperties(20.0, 1.2e-5, 2100.0, 0.05, 10.0, PhaseRegime.CONDENSING);
    single.tubeQualityIn = 0.5;
    single.tubeQualityOut = 0.5;
    point.addZone(single);
    row("zero_quality_range", point, null);

    // 7. A zero tube-side mass flow.
    IncrementalZoneAnalysis stopped = new IncrementalZoneAnalysis();
    IncrementalZone noFlow = zone("zero_flow", 100000.0, 400.0, 370.0, 300.0, 330.0);
    noFlow.setTubeSideProperties(800.0, 5e-4, 2500.0, 0.13, 0.0, PhaseRegime.LIQUID);
    stopped.addZone(noFlow);
    row("zero_flow", stopped, null);

    // 8. No zones at all: the sentinels.
    IncrementalZoneAnalysis empty = new IncrementalZoneAnalysis();
    empty.calculate();
    System.out.println("empty_zones");
    geometry(empty);
    printMap("result.", empty.toMap());
    System.out.println();

    // 9. The constructor's target count clamps and changes nothing.
    for (int target : new int[] {0, 1, 3, 20, 100}) {
      IncrementalZoneAnalysis targetCount = new IncrementalZoneAnalysis(target);
      targetCount.addZone(condensing);
      targetCount.calculate();
      System.out.println("target_zone_count." + target);
      System.out.println("constructor_argument=" + target);
      System.out.println("target_zone_count=" + targetCount.getTargetZoneCount());
      System.out.println("total_required_area_m2=" + targetCount.getTotalRequiredArea());
      System.out.println("weighted_average_u_Wpm2K=" + targetCount.getWeightedAverageU());
      System.out.println("total_shell_pressure_drop_Pa="
          + targetCount.getTotalShellSidePressureDrop());
      System.out.println();
    }

    // 10. The shell drop's baffle count scales with the zone count.
    for (int count : new int[] {2, 4, 10}) {
      IncrementalZoneAnalysis scaled = new IncrementalZoneAnalysis();
      for (int i = 0; i < count; i++) {
        scaled.addZone(condensing);
      }
      scaled.calculate();
      System.out.println("baffle_scale." + count);
      System.out.println("zone_count=" + count);
      geometry(scaled);
      System.out.println("total_shell_pressure_drop_Pa="
          + scaled.getTotalShellSidePressureDrop());
      System.out.println("total_tube_pressure_drop_Pa="
          + scaled.getTotalTubeSidePressureDrop());
      System.out.println();
    }

    // 11. Passes, including a count that does not divide.
    for (int passes : new int[] {1, 2, 4, 3}) {
      IncrementalZoneAnalysis analysis = new IncrementalZoneAnalysis();
      analysis.setGeometry(0.01483, 0.01905, 6.096, 100, passes);
      analysis.addZone(liquid);
      analysis.calculate();
      System.out.println("geometry_passes." + passes);
      System.out.println("tube_count=100");
      System.out.println("tube_passes=" + passes);
      System.out.println("total_tube_pressure_drop_Pa="
          + analysis.getTotalTubeSidePressureDrop());
      System.out.println("tube_htc_Wpm2K="
          + analysis.getZones().get(0).tubeSideHTC);
      System.out.println();
    }

    System.out.flush();
  }

  static IncrementalZone zone(String name, double duty, double hotIn, double hotOut,
      double coldIn, double coldOut) {
    IncrementalZone zone = new IncrementalZone();
    zone.zoneName = name;
    zone.duty = duty;
    zone.setTemperatures(hotIn, hotOut, coldIn, coldOut);
    zone.setShellSideProperties(998.0, 0.001, 4180.0, 0.60, 20.0, PhaseRegime.LIQUID);
    return zone;
  }

  static void row(String label, IncrementalZoneAnalysis analysis, IncrementalZone single) {
    if (single != null) {
      analysis.addZone(single);
    }
    analysis.calculate();
    System.out.println(label);
    geometry(analysis);
    printMap("result.", analysis.toMap());
    System.out.println();
  }

  static void geometry(IncrementalZoneAnalysis analysis) {
    System.out.println("target_zone_count=" + analysis.getTargetZoneCount());
    System.out.println("tube_inside_diameter_m=0.01483");
    System.out.println("tube_outside_diameter_m=0.01905");
    System.out.println("tube_length_m=6.096");
    System.out.println("tube_count=100");
    System.out.println("tube_passes=2");
    System.out.println("tube_pitch_m=0.02381");
    System.out.println("triangular_pitch=true");
    System.out.println("shell_inside_diameter_m=0.5");
    System.out.println("baffle_spacing_m=0.2");
    System.out.println("tube_wall_conductivity_WpmK=52.0");
    System.out.println("fouling_tube=1.76E-4");
    System.out.println("fouling_shell=1.76E-4");
    System.out.println("zone_count=" + analysis.getZones().size());
    List<IncrementalZone> zones = analysis.getZones();
    for (int i = 0; i < zones.size(); i++) {
      IncrementalZone zone = zones.get(i);
      String p = "zones." + i + ".";
      System.out.println(p + "zone_name=" + zone.zoneName);
      System.out.println(p + "duty_W=" + zone.duty);
      System.out.println(p + "hot_inlet_temp_K=" + zone.hotInletTemp);
      System.out.println(p + "hot_outlet_temp_K=" + zone.hotOutletTemp);
      System.out.println(p + "cold_inlet_temp_K=" + zone.coldInletTemp);
      System.out.println(p + "cold_outlet_temp_K=" + zone.coldOutletTemp);
      System.out.println(p + "tube_phase_regime=" + zone.tubePhaseRegime.name());
      System.out.println(p + "tube_density=" + zone.tubeDensity);
      System.out.println(p + "tube_viscosity=" + zone.tubeViscosity);
      System.out.println(p + "tube_cp=" + zone.tubeCp);
      System.out.println(p + "tube_conductivity=" + zone.tubeConductivity);
      System.out.println(p + "tube_mass_flow_kgps=" + zone.tubeMassFlowRate);
      System.out.println(p + "tube_liquid_density=" + zone.tubeLiquidDensity);
      System.out.println(p + "tube_vapor_density=" + zone.tubeVaporDensity);
      System.out.println(p + "tube_liquid_viscosity=" + zone.tubeLiquidViscosity);
      System.out.println(p + "tube_vapor_viscosity=" + zone.tubeVaporViscosity);
      System.out.println(p + "tube_liquid_cp=" + zone.tubeLiquidCp);
      System.out.println(p + "tube_liquid_conductivity=" + zone.tubeLiquidConductivity);
      System.out.println(p + "tube_quality_in=" + zone.tubeQualityIn);
      System.out.println(p + "tube_quality_out=" + zone.tubeQualityOut);
      System.out.println(p + "tube_reduced_pressure=" + zone.tubeReducedPressure);
      System.out.println(p + "tube_surface_tension=" + zone.tubeSurfaceTension);
      System.out.println(p + "tube_heat_of_vaporization=" + zone.tubeHeatOfVaporization);
      System.out.println(p + "shell_phase_regime=" + zone.shellPhaseRegime.name());
      System.out.println(p + "shell_density=" + zone.shellDensity);
      System.out.println(p + "shell_viscosity=" + zone.shellViscosity);
      System.out.println(p + "shell_cp=" + zone.shellCp);
      System.out.println(p + "shell_conductivity=" + zone.shellConductivity);
      System.out.println(p + "shell_mass_flow_kgps=" + zone.shellMassFlowRate);
    }
  }

  @SuppressWarnings("unchecked")
  static void printMap(String prefix, Map<String, Object> map) {
    for (Map.Entry<String, Object> entry : map.entrySet()) {
      if (entry.getValue() instanceof Map) {
        printMap(prefix + entry.getKey() + ".", (Map<String, Object>) entry.getValue());
      } else if (entry.getValue() instanceof List) {
        List<Object> rows = (List<Object>) entry.getValue();
        for (int i = 0; i < rows.size(); i++) {
          if (rows.get(i) instanceof Map) {
            printMap(prefix + entry.getKey() + "." + i + ".",
                (Map<String, Object>) rows.get(i));
          } else {
            System.out.println(prefix + entry.getKey() + "." + i + "=" + rows.get(i));
          }
        }
      } else {
        System.out.println(prefix + entry.getKey() + "=" + entry.getValue());
      }
    }
  }
}
