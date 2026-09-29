// `TrayHydraulicsCalculator` driven directly, on states of its own.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar TrayHydraulicsProbe.java
//     java -cp .:neqsim-f0c7436.jar TrayHydraulicsProbe > captures/tray_hydraulics_probe.tsv
//
// **This is the tray hydraulics of `absorption_column`, `stripping_column` and
// `distillation_column`, and none of their captures holds its output.** All three declare
// `max_allowable_gas_load_factor` and refuse it rather than reading it, and `getFsFactor()` is the
// quantity that would answer it - so nothing in this directory measures the arithmetic itself.
// This probe states the geometry and the load and prints both sides, which is what makes the port's
// case set an oracle rather than a self-consistency check.
//
// **The class branches on the tray type**, and three of them are named in its own bytecode: `sieve`
// is the constructor's default and `valve` and `bubble-cap` are the other two.
//
// **Not every field has a getter.** `activeAreaFraction`, `holeDiameter`, `holeAreaFraction` and
// `relativeVolatility` are settable and unreadable, so the values printed for them are the ones
// this probe set - the class's own constructor initialisers - rather than values read back.

import neqsim.process.equipment.distillation.internals.TrayHydraulicsCalculator;

public class TrayHydraulicsProbe {

  public static void main(String[] args) {
    // The class's own geometry and its default state. Every other row moves one thing from this.
    row("sieve_default", "sieve", 1.0, 1.0, 3.0, 2.0, -1.0);
    // The other two tray types the class names, at the same state.
    row("valve_tray", "valve", 1.0, 1.0, 3.0, 2.0, -1.0);
    row("bubble_cap_tray", "bubble-cap", 1.0, 1.0, 3.0, 2.0, -1.0);
    // A heavier vapour load, which is what moves the flood and the pressure drop.
    row("sieve_high_vapor", "sieve", 1.0, 2.4, 3.0, 2.0, -1.0);
    // A light vapour load, which is where the weeping verdict turns: the actual velocity falls
    // below the minimum this tray needs, and `isWeepingOk` is the only verdict that goes false.
    row("sieve_low_vapor", "sieve", 1.0, 0.03, 3.0, 2.0, -1.0);
    // The same load on a valve tray, whose minimum velocity is higher - so the type is what
    // decides, at one state.
    row("valve_low_vapor", "valve", 1.0, 0.03, 3.0, 2.0, -1.0);
    // A stated weir length, where the class's own is `-1.0` and derived from the diameter.
    row("sieve_stated_weir", "sieve", 1.0, 1.0, 3.0, 2.0, 0.7);
    // A wider column at a proportionally larger load, so the areas are not all one diameter's.
    row("sieve_wide", "sieve", 1.8, 3.2, 9.7, 2.0, -1.0);
  }

  static void row(String label, String trayType, double diameterM, double vaporMassFlowKgPerS,
      double liquidMassFlowKgPerS, double vaporDensity, double weirLengthM) {
    // The class's own constructor initialisers, stated here so the capture records them.
    double traySpacingM = 0.6;
    double weirHeightM = 0.05;
    double downcommerAreaFraction = 0.1;
    double holeDiameterMm = 12.7;
    double holeAreaFraction = 0.1;
    double activeAreaFraction = 0.8;
    double designFloodFraction = 0.8;
    double liquidDensity = 800.0;
    double liquidViscosity = 1.0e-3;
    double surfaceTension = 0.02;
    double relativeVolatility = 2.0;

    TrayHydraulicsCalculator tray = new TrayHydraulicsCalculator();
    tray.setTrayType(trayType);
    tray.setColumnDiameter(diameterM);
    tray.setTraySpacing(traySpacingM);
    tray.setWeirHeight(weirHeightM);
    tray.setWeirLength(weirLengthM);
    tray.setDowncommerAreaFraction(downcommerAreaFraction);
    tray.setHoleDiameter(holeDiameterMm);
    tray.setHoleAreaFraction(holeAreaFraction);
    tray.setDesignFloodFraction(designFloodFraction);
    tray.setVaporMassFlow(vaporMassFlowKgPerS);
    tray.setLiquidMassFlow(liquidMassFlowKgPerS);
    tray.setVaporDensity(vaporDensity);
    tray.setLiquidDensity(liquidDensity);
    tray.setLiquidViscosity(liquidViscosity);
    tray.setSurfaceTension(surfaceTension);
    tray.setRelativeVolatility(relativeVolatility);
    tray.calculate();

    System.out.println(label);
    System.out.println("tray_type=" + trayType);
    System.out.println("column_diameter_m=" + diameterM);
    System.out.println("tray_spacing_m=" + traySpacingM);
    System.out.println("weir_height_m=" + weirHeightM);
    System.out.println("weir_length_m=" + weirLengthM);
    System.out.println("downcommer_area_fraction=" + downcommerAreaFraction);
    System.out.println("hole_diameter_mm=" + holeDiameterMm);
    System.out.println("hole_area_fraction=" + holeAreaFraction);
    System.out.println("active_area_fraction=" + activeAreaFraction);
    System.out.println("design_flood_fraction=" + designFloodFraction);
    System.out.println("vapor_mass_flow_kg_per_s=" + vaporMassFlowKgPerS);
    System.out.println("liquid_mass_flow_kg_per_s=" + liquidMassFlowKgPerS);
    System.out.println("vapor_density_kg_per_m3=" + vaporDensity);
    System.out.println("liquid_density_kg_per_m3=" + liquidDensity);
    System.out.println("liquid_viscosity_Pa_s=" + liquidViscosity);
    System.out.println("surface_tension_N_per_m=" + surfaceTension);
    System.out.println("relative_volatility=" + relativeVolatility);
    // ---- The class's own answers.
    System.out.println("flooding_velocity_m_per_s=" + tray.getFloodingVelocity());
    System.out.println("actual_vapor_velocity_m_per_s=" + tray.getActualVaporVelocity());
    System.out.println("percent_flood=" + tray.getPercentFlood());
    System.out.println("minimum_vapor_velocity_m_per_s=" + tray.getMinimumVaporVelocity());
    System.out.println("fs_factor=" + tray.getFsFactor());
    System.out.println("weeping_ok=" + tray.isWeepingOk());
    System.out.println("entrainment=" + tray.getEntrainment());
    System.out.println("entrainment_ok=" + tray.isEntrainmentOk());
    System.out.println("downcommer_backup_m=" + tray.getDowncommerBackup());
    System.out.println("downcommer_backup_fraction=" + tray.getDowncommerBackupFraction());
    System.out.println("downcommer_backup_ok=" + tray.isDowncommerBackupOk());
    System.out.println("total_pressure_drop_Pa=" + tray.getTotalTrayPressureDrop());
    System.out.println("total_pressure_drop_mbar=" + tray.getTotalTrayPressureDropMbar());
    System.out.println("dry_pressure_drop_Pa=" + tray.getDryTrayPressureDrop());
    System.out.println("liquid_head_pressure_drop_Pa=" + tray.getLiquidHeadPressureDrop());
    System.out.println("residual_head_pressure_drop_Pa=" + tray.getResidualHeadPressureDrop());
    System.out.println("tray_efficiency=" + tray.getTrayEfficiency());
    System.out.println("turndown_ratio=" + tray.getTurndownRatio());
    System.out.println("calculated_weir_length_m=" + tray.getCalculatedWeirLength());
    System.out.println("active_area_m2=" + tray.getActiveArea());
    System.out.println("total_area_m2=" + tray.getTotalArea());
    System.out.println("hole_area_m2=" + tray.getHoleArea());
    System.out.println("downcommer_area_m2=" + tray.getDowncommerArea());
    System.out.println("design_ok=" + tray.isDesignOk());
    // **Last, because it is the one call that changes the object it is read from.**
    // `sizeColumnDiameter` writes the trial `1.0` m into `columnDiameter`, re-derives the areas
    // and the flooding velocity there, and leaves the *sized* value behind - so every line above
    // is read at the diameter this row stated and this line is the only one that is not.
    System.out.println("sized_column_diameter_m=" + tray.sizeColumnDiameter());
    System.out.println();
  }
}
