// The scale family's other five classes: the barite/celestite solid solution, the mass
// calculator, the multi-mineral equilibrium in its three activity modes, the flowline
// profile (with its inert pH argument) and the water-compatibility screener (with its two
// inert pH fields).
//
//     javac -proc:none -cp neqsim-f0c7436.jar ScaleFamilyProbe.java
//     java -cp .:neqsim-f0c7436.jar ScaleFamilyProbe > captures/scale_family_probe.tsv

import java.util.Map;
import neqsim.pvtsimulation.flowassurance.BariteCelestiteSolidSolution;
import neqsim.pvtsimulation.flowassurance.FlowlineScaleProfile;
import neqsim.pvtsimulation.flowassurance.MultiMineralScaleEquilibrium;
import neqsim.pvtsimulation.flowassurance.ScaleMassCalculator;
import neqsim.pvtsimulation.flowassurance.ScalePredictionCalculator;
import neqsim.pvtsimulation.flowassurance.WaterCompatibilityScreener;

public class ScaleFamilyProbe {

  public static void main(String[] args) {
    solidSolution();
    masses();
    equilibrium();
    profile();
    compatibility();
  }

  // ---------------------------------------------------------------- solid solution
  static void solidSolution() {
    System.out.println("solid_solution");
    solidRow("pure_barite", 2.3, 298.15, 0.01, 0.0, 0.05, 1.08e-10, 3.44e-7);
    solidRow("pure_celestite", 2.3, 298.15, 0.0, 0.01, 0.05, 1.08e-10, 3.44e-7);
    solidRow("both_zero", 2.3, 298.15, 0.0, 0.0, 0.05, 1.08e-10, 3.44e-7);
    solidRow("strontium_trace", 2.3, 298.15, 0.01, 1.0e-31, 0.05, 1.08e-10, 3.44e-7);
    solidRow("mixed", 2.3, 298.15, 0.001, 0.005, 0.01, 1.08e-10, 3.44e-7);
    solidRow("mixed_at_350K", 2.3, 350.0, 0.001, 0.005, 0.01, 1.08e-10, 3.44e-7);
    solidRow("mixed_margules_zero", 0.0, 298.15, 0.001, 0.005, 0.01, 1.08e-10, 3.44e-7);
    solidRow("mixed_margules_five", 5.0, 298.15, 0.001, 0.005, 0.01, 1.08e-10, 3.44e-7);
    solidRow("mixed_custom_ksp", 2.3, 298.15, 0.001, 0.005, 0.01, 1.0e-9, 1.0e-6);
    System.out.println();
  }

  static void solidRow(String label, double margules, double temperature, double aBa, double aSr,
      double aSo4, double kspBa, double kspSr) {
    BariteCelestiteSolidSolution solution = new BariteCelestiteSolidSolution();
    solution.setMargules(margules);
    solution.setTemperature(temperature);
    solution.setAqueousActivities(aBa, aSr, aSo4);
    solution.setEndMemberKsp(kspBa, kspSr);
    solution.calculate();
    System.out.println("  " + label + " margules=" + margules + " temperature=" + temperature
        + " a_ba=" + aBa + " a_sr=" + aSr + " a_so4=" + aSo4 + " ksp_ba=" + kspBa + " ksp_sr="
        + kspSr + " x_baso4=" + solution.getBaSO4MoleFraction() + " x_srso4="
        + solution.getSrSO4MoleFraction() + " total_si=" + solution.getTotalSaturationIndex());
  }

  // ---------------------------------------------------------------- masses
  static void masses() {
    System.out.println("scale_mass");
    ScaleMassCalculator calculator = new ScaleMassCalculator(new ScalePredictionCalculator());
    calculator.setWaterVolume(1.0);
    double[][] states = {{2.5, 0.01, 0.02}, {-1.0, 0.01, 0.02}, {5.0, 0.01, 0.02},
        {0.0, 0.001, 0.001}};
    for (double[] state : states) {
      double si = state[0];
      double cation = state[1];
      double anion = state[2];
      System.out.println("  si=" + si + " cation=" + cation + " anion=" + anion
          + " caco3=" + calculator.calcCaCO3Mass(cation, anion, si)
          + " baso4=" + calculator.calcBaSO4Mass(cation, anion, si)
          + " srso4=" + calculator.calcSrSO4Mass(cation, anion, si)
          + " caso4=" + calculator.calcCaSO4Mass(cation, anion, si)
          + " feco3=" + calculator.calcFeCO3Mass(cation, anion, si));
    }
    calculator.setWaterVolume(2.0);
    System.out.println("  volume_2l baso4_mass_mg_per_l=" + calculator.calcBaSO4Mass(0.01, 0.02, 2.5));
    System.out.println();
  }

  // ---------------------------------------------------------------- equilibrium
  static void equilibrium() {
    System.out.println("multi_mineral_equilibrium");
    equilibriumRow("sulphate_limited_davies", MultiMineralScaleEquilibrium.ActivityModel.DAVIES,
        1000.0, 1000.0, 3000.0, 200.0, 1500.0, 60000.0, 70.0, 100.0, 1.0, 0.0, 0.0);
    equilibriumRow("sulphate_limited_bdot", MultiMineralScaleEquilibrium.ActivityModel.BDOT,
        1000.0, 1000.0, 3000.0, 200.0, 1500.0, 60000.0, 70.0, 100.0, 1.0, 0.0, 0.0);
    equilibriumRow("sulphate_limited_pitzer", MultiMineralScaleEquilibrium.ActivityModel.PITZER_BINARY,
        1000.0, 1000.0, 3000.0, 200.0, 1500.0, 60000.0, 70.0, 100.0, 1.0, 0.0, 0.0);
    equilibriumRow("high_salinity_davies", MultiMineralScaleEquilibrium.ActivityModel.DAVIES,
        15000.0, 400.0, 800.0, 200.0, 600.0, 200000.0, 90.0, 200.0, 2.0, 0.0, 0.0);
    equilibriumRow("high_salinity_bdot", MultiMineralScaleEquilibrium.ActivityModel.BDOT,
        15000.0, 400.0, 800.0, 200.0, 600.0, 200000.0, 90.0, 200.0, 2.0, 0.0, 0.0);
    equilibriumRow("high_salinity_pitzer", MultiMineralScaleEquilibrium.ActivityModel.PITZER_BINARY,
        15000.0, 400.0, 800.0, 200.0, 600.0, 200000.0, 90.0, 200.0, 2.0, 0.0, 0.0);
    equilibriumRow("co_precipitation", MultiMineralScaleEquilibrium.ActivityModel.DAVIES,
        500.0, 100.0, 2500.0, 100.0, 4000.0, 60000.0, 70.0, 100.0, 1.0, 0.0, 0.0);
    equilibriumRow("single_mineral_barite_pure_water",
        MultiMineralScaleEquilibrium.ActivityModel.DAVIES,
        0.0, 1000.0, 0.0, 0.0, 700.0, 20.0, 25.0, 1.013, 1.0, 7.0, 0.0);
    equilibriumRow("no_minerals", MultiMineralScaleEquilibrium.ActivityModel.DAVIES,
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 25.0, 1.013, 1.0, 7.0, 0.0);
    equilibriumControls("iteration_limit_one_step", 1, 0.5, 1.0e-30);
    equilibriumWaterVolume("water_volume_1000l", 1000.0);
    equilibriumPitzerStrength("pitzer_explicit_i_1p5", 1.5);
    System.out.println();
  }

  static ScalePredictionCalculator equilibriumPredictor(double calcium, double barium,
      double strontium, double bicarbonate, double sulphate, double tds, double temperatureC,
      double pressureBara, double co2pp, double ph, double automaticPh) {
    ScalePredictionCalculator predictor = new ScalePredictionCalculator();
    predictor.setCalciumConcentration(calcium);
    predictor.setBariumConcentration(barium);
    predictor.setStrontiumConcentration(strontium);
    predictor.setIronConcentration(0.0);
    predictor.setBicarbonateConcentration(bicarbonate);
    predictor.setSulphateConcentration(sulphate);
    predictor.setTotalDissolvedSolids(tds);
    predictor.setTemperatureCelsius(temperatureC);
    predictor.setPressureBara(pressureBara);
    predictor.setCO2PartialPressure(co2pp);
    if (automaticPh > 0.5) {
      predictor.enableAutoPH();
    } else {
      predictor.setPH(ph);
    }
    return predictor;
  }

  static void equilibriumRow(String label, MultiMineralScaleEquilibrium.ActivityModel model,
      double calcium, double barium, double strontium, double bicarbonate, double sulphate,
      double tds, double temperatureC, double pressureBara, double co2pp, double ph,
      double automaticPh) {
    ScalePredictionCalculator predictor = equilibriumPredictor(calcium, barium, strontium,
        bicarbonate, sulphate, tds, temperatureC, pressureBara, co2pp, ph, automaticPh);
    MultiMineralScaleEquilibrium equilibrium = new MultiMineralScaleEquilibrium(predictor);
    equilibrium.setActivityModel(model);
    printEquilibrium(label, equilibrium.solve());
  }

  static void equilibriumControls(String label, int maxIterations, double relaxation,
      double threshold) {
    ScalePredictionCalculator predictor = equilibriumPredictor(1000.0, 1000.0, 3000.0, 200.0,
        1500.0, 60000.0, 70.0, 100.0, 1.0, 0.0, 1.0);
    MultiMineralScaleEquilibrium equilibrium = new MultiMineralScaleEquilibrium(predictor);
    equilibrium.setSolverControls(maxIterations, relaxation, threshold);
    printEquilibrium(label, equilibrium.solve());
  }

  static void equilibriumWaterVolume(String label, double litres) {
    ScalePredictionCalculator predictor = equilibriumPredictor(1000.0, 1000.0, 3000.0, 200.0,
        1500.0, 60000.0, 70.0, 100.0, 1.0, 0.0, 1.0);
    MultiMineralScaleEquilibrium equilibrium = new MultiMineralScaleEquilibrium(predictor);
    equilibrium.setWaterVolume(litres);
    printEquilibrium(label, equilibrium.solve());
  }

  static void equilibriumPitzerStrength(String label, double strength) {
    ScalePredictionCalculator predictor = equilibriumPredictor(1000.0, 1000.0, 3000.0, 200.0,
        1500.0, 60000.0, 70.0, 100.0, 1.0, 0.0, 1.0);
    MultiMineralScaleEquilibrium equilibrium = new MultiMineralScaleEquilibrium(predictor);
    equilibrium.setActivityModel(MultiMineralScaleEquilibrium.ActivityModel.PITZER_BINARY);
    equilibrium.setPitzerIonicStrengthMolal(strength);
    printEquilibrium(label, equilibrium.solve());
  }

  static void printEquilibrium(String label, MultiMineralScaleEquilibrium equilibrium) {
    System.out.println("  " + label + " iterations=" + equilibrium.getIterationCount()
        + " limit_reached=" + equilibrium.hasReachedIterationLimit()
        + " complementarity=" + equilibrium.getMaximumComplementarityViolation()
        + " ion_balance=" + equilibrium.getMaximumIonBalanceResidualMolPerL()
        + " divalent_gamma=" + equilibrium.getDivalentActivityCoefficientUsed()
        + " total_mass_mg_per_l=" + equilibrium.getTotalScaleMassMgPerL());
    for (Map.Entry<String, MultiMineralScaleEquilibrium.MineralResult> entry : equilibrium
        .getResults().entrySet()) {
      MultiMineralScaleEquilibrium.MineralResult mineral = entry.getValue();
      System.out.println("    mineral=" + mineral.getName() + " common=" + mineral.getFormula()
          + " initial_si=" + mineral.getInitialSI() + " final_si="
          + mineral.getFinalSI() + " product="
          + mineral.getActivityCoefficientProduct() + " mol_per_l="
          + mineral.getPrecipitatedMolPerL() + " mg_per_l="
          + mineral.getPrecipitatedMassMgPerL() + " grams=" + mineral.getPrecipitatedMassGrams()
          + " limiting=" + mineral.isLimitingIonExhausted());
    }
    for (String ion : new String[] {"Ca++", "Ba++", "Sr++", "Fe++", "SO4--", "CO3--"}) {
      System.out.println("    residual_" + ion + "=" + equilibrium.getResidualFreeIonMolPerL(ion));
    }
    for (Map.Entry<String, Double> entry : equilibrium.getActivityCoefficientsUsed().entrySet()) {
      System.out.println("    activity_" + entry.getKey() + "=" + entry.getValue());
    }
  }

  // ---------------------------------------------------------------- flowline profile
  static void profile() {
    System.out.println("flowline_scale_profile");
    profileRow("upstream_state", 6.5, 10);
    profileRow("one_segment_clamped", 6.5, 1);
    profileRow("ph_argument_four", 4.0, 10);
    profileRow("high_sulphate_120c", 6.5, 10);
    System.out.println();
  }

  static void profileRow(String label, double ph, int segments) {
    FlowlineScaleProfile profile = new FlowlineScaleProfile();
    if (label.equals("high_sulphate_120c")) {
      profile.setWaterChemistry(2000.0, 0.0, 0.0, 0.0, 100.0, 3000.0, 60000.0, 2.0, ph);
      profile.setInletConditions(120.0, 200.0);
      profile.setOutletConditions(40.0, 50.0);
    } else {
      profile.setWaterChemistry(400.0, 10.0, 5.0, 0.0, 150.0, 10.0, 40000.0, 2.0, ph);
      profile.setInletConditions(80.0, 200.0);
      profile.setOutletConditions(20.0, 50.0);
    }
    profile.setMgNaConcentrations(150.0, 10.0);
    profile.setNumberOfSegments(segments);
    profile.setAutoPH(false);
    profile.calculate();
    System.out.println("  " + label + " ph_argument=" + ph + " segments=" + segments);
    for (FlowlineScaleProfile.SegmentResult segment : profile.getResults()) {
      System.out.println("    fraction=" + segment.distanceFraction + " t_c="
          + segment.temperatureC + " p_bar=" + segment.pressureBar + " si_caco3="
          + segment.siCaCO3 + " si_baso4=" + segment.siBaSO4 + " si_srso4=" + segment.siSrSO4
          + " si_caso4=" + segment.siCaSO4 + " si_feco3=" + segment.siFeCO3);
    }
    for (String mineral : new String[] {"CaCO3", "BaSO4", "SrSO4", "CaSO4", "FeCO3"}) {
      System.out.println("    max_" + mineral + "=" + profile.getMaxSI(mineral));
    }
  }

  // ---------------------------------------------------------------- compatibility
  static void compatibility() {
    System.out.println("water_compatibility");
    compatibilityRow("upstream_seawater", new double[] {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0,
        70.0, 80.0, 90.0, 100.0}, 6.5);
    compatibilityRow("single_ratio", new double[] {0.0}, 6.5);
    compatibilityRow("ph_argument_four", new double[] {0.0, 50.0, 100.0}, 4.0);
    compatibilityRow("ph_argument_eleven", new double[] {0.0, 50.0, 100.0}, 11.0);
    compatibilityRow("five_ratios", new double[] {0.0, 25.0, 50.0, 75.0, 100.0}, 6.5);
    System.out.println();
  }

  static void compatibilityRow(String label, double[] ratios, double ph) {
    WaterCompatibilityScreener screener = new WaterCompatibilityScreener();
    screener.setFormationWater(400.0, 200.0, 5.0, 0.0, 150.0, 10.0, 50000.0, 80.0, 100.0, 2.0, ph);
    screener.setInjectionWater(20.0, 0.0, 0.0, 0.0, 100.0, 2700.0, 35000.0, 20.0, 100.0, 0.5, ph);
    screener.setFormationWaterMgNa(150.0, 10000.0);
    screener.setInjectionWaterMgNa(100.0, 11000.0);
    screener.setMixingRatios(ratios);
    screener.calculate();
    System.out.println("  " + label + " ph_argument=" + ph);
    for (WaterCompatibilityScreener.MixingResult result : screener.getResults()) {
      System.out.println("    fraction_pct=" + result.injectionWaterPct + " si_caco3="
          + result.siCaCO3 + " si_baso4=" + result.siBaSO4 + " si_srso4=" + result.siSrSO4
          + " si_caso4=" + result.siCaSO4 + " si_feco3=" + result.siFeCO3);
    }
    System.out.println("    worst_ratio=" + screener.getWorstCaseRatio() + " worst_scale="
        + screener.getWorstCaseScale() + " worst_si=" + screener.getWorstCaseSI());
  }
}
