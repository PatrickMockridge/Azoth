// `ThermalDesignCalculator.calculateZones` driven directly, on the class's documented zone
// example, a three-zone case and the states the permissive path leaves unanswered.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-118d9af.jar ZoneAnalysisProbe.java
//     java -cp .:neqsim-118d9af.jar ZoneAnalysisProbe > captures/zone_analysis_probe.tsv
//
// `calculateZones` calls the permissive `calculate()`, not `calculateStrict()`. Each zone
// overrides only the eight transport properties, runs the full rating, and restores them; the
// required area is `dutyFraction * totalDuty / (U * lmtd)`, zero when either denominator is
// non-positive. The probe prints the base state, every zone's twelve inputs and all seven
// results, so a row can be reconstructed without the Java source.

import java.util.Map;

import neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator;
import neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator.ShellSideMethod;
import neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator.ZoneDefinition;
import neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator.ZoneResult;

public class ZoneAnalysisProbe {

  public static void main(String[] args) {
    // 1. The documented zone example, on the class's default base state.
    ThermalDesignCalculator doc = new ThermalDesignCalculator();
    ZoneDefinition desuperheating = zone("desuperheating", 0.15, 150000.0, 40.0,
        5.0, 1.2e-5, 2100.0, 0.025, 998.0, 0.001, 4180.0, 0.60);
    row("doc_single_desuperheating", doc, new ZoneDefinition[] {desuperheating}, false);

    // 2. Three zones over the test fixture's liquid/gas base.
    ThermalDesignCalculator simple = simpleCalculator();
    ZoneDefinition[] three = {
        zone("desuperheating", 0.2, 1500000.0, 60.0,
            20.0, 1.2e-5, 2200.0, 0.05, 998.0, 0.001, 4180.0, 0.60),
        zone("condensing", 0.5, 1500000.0, 40.0,
            5.0, 1.2e-5, 2100.0, 0.025, 998.0, 0.001, 4180.0, 0.60),
        zone("subcooling", 0.3, 1500000.0, 25.0,
            800.0, 5e-4, 2500.0, 0.13, 998.0, 0.001, 4180.0, 0.60)};
    row("simple_three_zone", simple, three, false);

    // 3. Bell-Delaware with two sealing pairs, on the same base.
    ThermalDesignCalculator bell = simpleCalculator();
    bell.setShellSideMethod(ShellSideMethod.BELL_DELAWARE);
    bell.setHasSealing(true);
    bell.setSealingPairs(2);
    row("bell_delaware_sealed_base", bell, new ZoneDefinition[] {three[1]}, false);

    // 4. A uniform inner deposit on the same base.
    ThermalDesignCalculator deposit = simpleCalculator();
    deposit.setTubeFoulingLayer(0.001, 0.2);
    row("deposit_base", deposit, new ZoneDefinition[] {three[1]}, true);

    // 5. Degenerate zone properties: the permissive path's unanswered states.
    ThermalDesignCalculator degenerate = simpleCalculator();
    ZoneDefinition[] bad = {
        zone("zero_viscosity", 1.0, 1000000.0, 40.0,
            800.0, 0.0, 2500.0, 0.13, 998.0, 0.001, 4180.0, 0.60),
        zone("zero_conductivity", 1.0, 1000000.0, 40.0,
            800.0, 5e-4, 2500.0, 0.0, 998.0, 0.001, 4180.0, 0.60),
        zone("zero_density", 1.0, 1000000.0, 40.0,
            0.0, 5e-4, 2500.0, 0.13, 998.0, 0.001, 4180.0, 0.60)};
    row("degenerate_properties", degenerate, bad, false);

    // 6. A zero and a negative LMTD: U positive, required area zero.
    ThermalDesignCalculator lmtd = simpleCalculator();
    ZoneDefinition[] lmtdZones = {
        zone("zero_lmtd", 1.0, 1000000.0, 0.0,
            800.0, 5e-4, 2500.0, 0.13, 998.0, 0.001, 4180.0, 0.60),
        zone("negative_lmtd", 1.0, 1000000.0, -25.0,
            800.0, 5e-4, 2500.0, 0.13, 998.0, 0.001, 4180.0, 0.60)};
    row("zero_and_negative_lmtd", lmtd, lmtdZones, false);

    // 7. Empty and null zone arrays.
    ThermalDesignCalculator empty = simpleCalculator();
    System.out.println("empty_and_null_zones");
    base(empty, "simple", ShellSideMethod.KERN, false, 0, 0.0, 1.0);
    System.out.println("empty_count=" + empty.calculateZones(new ZoneDefinition[0]).length);
    System.out.println("null_count=" + empty.calculateZones(null).length);
    System.out.println();

    // 8. The base state before, during and after the zone loop.
    ThermalDesignCalculator restore = simpleCalculator();
    restore.calculate();
    double before = restore.getOverallU();
    ZoneResult[] restored = restore.calculateZones(new ZoneDefinition[] {three[1]});
    double afterZones = restore.getOverallU();
    restore.calculate();
    double recalculated = restore.getOverallU();
    System.out.println("restore_after_zones");
    base(restore, "simple", ShellSideMethod.KERN, false, 0, 0.0, 1.0);
    System.out.println("base_u_before=" + before);
    System.out.println("base_u_during_last_zone=" + afterZones);
    System.out.println("zone_u=" + restored[0].overallU);
    System.out.println("base_u_recalculated=" + recalculated);
    System.out.println();

    // 9. The strict refusals the permissive path does not make.
    strictRefusal("tube_count_zero", calc -> calc.setTubeCount(0));
    strictRefusal("passes_do_not_divide", calc -> calc.setTubePasses(3));
    strictRefusal("id_not_below_od", calc -> calc.setTubeIDm(0.02));
    strictRefusal("pitch_not_above_od", calc -> calc.setTubePitchm(0.019));
    strictRefusal("baffle_cut_half", calc -> calc.setBaffleCut(0.5));
    strictRefusal("zero_tube_density",
        calc -> calc.setTubeSideFluid(0.0, 0.001, 4180.0, 0.60, 5.0, true));
    strictRefusal("deposit_blocks_bore", calc -> calc.setTubeFoulingLayer(0.008, 0.2));
    strictRefusal("sealing_without_pairs", calc -> {
      calc.setShellSideMethod(ShellSideMethod.BELL_DELAWARE);
      calc.setHasSealing(true);
      calc.setSealingPairs(0);
    });

    System.out.flush();
  }

  interface Mutation {
    void apply(ThermalDesignCalculator calc);
  }

  static ThermalDesignCalculator simpleCalculator() {
    ThermalDesignCalculator calc = new ThermalDesignCalculator();
    calc.setTubeODm(0.01905);
    calc.setTubeIDm(0.01483);
    calc.setTubeLengthm(6.096);
    calc.setTubeCount(100);
    calc.setTubePasses(2);
    calc.setTubePitchm(0.02381);
    calc.setTriangularPitch(true);
    calc.setShellIDm(0.5);
    calc.setBaffleSpacingm(0.2);
    calc.setBaffleCount(10);
    calc.setBaffleCut(0.25);
    calc.setTubeWallConductivity(52.0);
    calc.setTubeSideFluid(998.0, 0.001, 4180.0, 0.60, 5.0, true);
    calc.setShellSideFluid(50.0, 1.5e-5, 2200.0, 0.03, 3.0);
    calc.setFoulingTube(0.00018);
    calc.setFoulingShell(0.00035);
    return calc;
  }

  static ZoneDefinition zone(String name, double dutyFraction, double totalDuty, double lmtd,
      double tubeDensity, double tubeViscosity, double tubeCp, double tubeConductivity,
      double shellDensity, double shellViscosity, double shellCp, double shellConductivity) {
    ZoneDefinition zone = new ZoneDefinition();
    zone.zoneName = name;
    zone.dutyFraction = dutyFraction;
    zone.totalDuty = totalDuty;
    zone.lmtd = lmtd;
    zone.tubeDensity = tubeDensity;
    zone.tubeViscosity = tubeViscosity;
    zone.tubeCp = tubeCp;
    zone.tubeConductivity = tubeConductivity;
    zone.shellDensity = shellDensity;
    zone.shellViscosity = shellViscosity;
    zone.shellCp = shellCp;
    zone.shellConductivity = shellConductivity;
    return zone;
  }

  static void row(String label, ThermalDesignCalculator calc, ZoneDefinition[] zones,
      boolean deposited) {
    System.out.println(label);
    if (label.startsWith("doc_")) {
      base(calc, "defaults", ShellSideMethod.KERN, false, 0, 0.0, 1.0);
    } else {
      ShellSideMethod method = calc.getShellSideMethod();
      base(calc, "simple", method, method == ShellSideMethod.BELL_DELAWARE,
          method == ShellSideMethod.BELL_DELAWARE ? 2 : 0,
          deposited ? 0.001 : 0.0, deposited ? 0.2 : 1.0);
    }
    System.out.println("zone_count=" + zones.length);
    for (int i = 0; i < zones.length; i++) {
      ZoneDefinition z = zones[i];
      String p = "zones." + i + ".";
      System.out.println(p + "zone_name=" + z.zoneName);
      System.out.println(p + "duty_fraction=" + z.dutyFraction);
      System.out.println(p + "total_duty_W=" + z.totalDuty);
      System.out.println(p + "lmtd_K=" + z.lmtd);
      System.out.println(p + "tube_density=" + z.tubeDensity);
      System.out.println(p + "tube_viscosity=" + z.tubeViscosity);
      System.out.println(p + "tube_cp=" + z.tubeCp);
      System.out.println(p + "tube_conductivity=" + z.tubeConductivity);
      System.out.println(p + "shell_density=" + z.shellDensity);
      System.out.println(p + "shell_viscosity=" + z.shellViscosity);
      System.out.println(p + "shell_cp=" + z.shellCp);
      System.out.println(p + "shell_conductivity=" + z.shellConductivity);
    }
    ZoneResult[] results = calc.calculateZones(zones);
    for (int i = 0; i < results.length; i++) {
      ZoneResult r = results[i];
      String p = "zones." + i + ".";
      System.out.println(p + "zone_name=" + r.zoneName);
      System.out.println(p + "duty_fraction=" + r.dutyFraction);
      System.out.println(p + "tube_side_htc_Wpm2K=" + r.tubeSideHTC);
      System.out.println(p + "shell_side_htc_Wpm2K=" + r.shellSideHTC);
      System.out.println(p + "overall_u_Wpm2K=" + r.overallU);
      System.out.println(p + "lmtd_K=" + r.lmtd);
      System.out.println(p + "required_area_m2=" + r.requiredArea);
    }
    printMap("last_to_map.", calc.toMap());
    System.out.println();
  }

  static void base(ThermalDesignCalculator calc, String kind, ShellSideMethod method,
      boolean hasSealing, int sealingPairs, double depositThickness, double depositConductivity) {
    System.out.println("base.kind=" + kind);
    System.out.println("base.tube_outside_diameter_m=" + 0.01905);
    System.out.println("base.tube_inside_diameter_m=" + 0.01483);
    System.out.println("base.tube_length_m=" + 6.096);
    System.out.println("base.tube_count=" + 100);
    System.out.println("base.tube_passes=" + 2);
    System.out.println("base.tube_pitch_m=" + 0.02381);
    System.out.println("base.triangular_pitch=true");
    System.out.println("base.shell_inside_diameter_m=" + 0.5);
    System.out.println("base.baffle_spacing_m=" + 0.2);
    System.out.println("base.baffle_count=" + 10);
    System.out.println("base.baffle_cut=" + 0.25);
    System.out.println("base.tube_wall_conductivity_WpmK=" + 52.0);
    if ("defaults".equals(kind)) {
      System.out.println("base.tube_density=" + 800.0);
      System.out.println("base.tube_viscosity=" + 5e-4);
      System.out.println("base.tube_cp=" + 2500.0);
      System.out.println("base.tube_conductivity=" + 0.13);
      System.out.println("base.tube_mass_flow_kgps=" + 10.0);
      System.out.println("base.shell_density=" + 1.5);
      System.out.println("base.shell_viscosity=" + 1.5e-5);
      System.out.println("base.shell_cp=" + 1100.0);
      System.out.println("base.shell_conductivity=" + 0.025);
      System.out.println("base.shell_mass_flow_kgps=" + 5.0);
      System.out.println("base.fouling_tube=" + 0.000176);
      System.out.println("base.fouling_shell=" + 0.000176);
    } else {
      System.out.println("base.tube_density=" + 998.0);
      System.out.println("base.tube_viscosity=" + 0.001);
      System.out.println("base.tube_cp=" + 4180.0);
      System.out.println("base.tube_conductivity=" + 0.60);
      System.out.println("base.tube_mass_flow_kgps=" + 5.0);
      System.out.println("base.shell_density=" + 50.0);
      System.out.println("base.shell_viscosity=" + 1.5e-5);
      System.out.println("base.shell_cp=" + 2200.0);
      System.out.println("base.shell_conductivity=" + 0.03);
      System.out.println("base.shell_mass_flow_kgps=" + 3.0);
      System.out.println("base.fouling_tube=" + 0.00018);
      System.out.println("base.fouling_shell=" + 0.00035);
    }
    System.out.println("base.tube_fouling_layer_thickness_m=" + depositThickness);
    System.out.println("base.tube_fouling_layer_conductivity_WpmK=" + depositConductivity);
    System.out.println("base.shell_viscosity_wall=-1.0");
    System.out.println("base.has_sealing=" + hasSealing);
    System.out.println("base.sealing_pairs=" + sealingPairs);
    System.out.println("base.tube_to_baffle_clearance_m=4.0E-4");
    System.out.println("base.shell_to_baffle_clearance_m=0.003");
    System.out.println("base.bypass_area_m2=0.0");
    System.out.println("base.shell_side_method=" + method.name());
    System.out.println("base.effective_tube_inside_diameter_m=" + calc.getEffectiveTubeIDm());
    System.out.println("base.tube_fouling_layer_resistance_outside="
        + calc.getTubeFoulingLayerResistanceOutside());
    System.out.println("base.outside_heat_transfer_area_m2=" + calc.getOutsideHeatTransferArea());
  }

  static void strictRefusal(String label, Mutation mutation) {
    ThermalDesignCalculator calc = simpleCalculator();
    String thrown;
    try {
      mutation.apply(calc);
      calc.calculateStrict();
      thrown = "none";
    } catch (RuntimeException e) {
      thrown = e.getClass().getSimpleName() + ": " + e.getMessage();
    }
    System.out.println("strict_refusals." + label + "=" + thrown);
  }

  @SuppressWarnings("unchecked")
  static void printMap(String prefix, Map<String, Object> map) {
    for (Map.Entry<String, Object> entry : map.entrySet()) {
      if (entry.getValue() instanceof Map) {
        printMap(prefix + entry.getKey() + ".", (Map<String, Object>) entry.getValue());
      } else {
        System.out.println(prefix + entry.getKey() + "=" + entry.getValue());
      }
    }
  }
}
