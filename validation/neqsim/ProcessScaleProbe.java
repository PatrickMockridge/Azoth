// ProcessScaleProbe - the process-chemistry scale calculators at the pinned NeqSim commit.
//
// Five classes that share nothing but the saturation index they consume: the activity-corrected
// `ElectrolyteScaleCalculator`, the two-brine `BrineMixingScaleEvaluator`, the empirical
// `ScaleInhibitorPerformance` model, the `ScaleKinetics` layer and the `ScaleControlAssessor`
// that wraps a `ScalePredictionCalculator`. The states are the upstream tests' own (the brucite
// seawater, the advanced-models brine, the two brine-mixing pairs, the kinetics rows and the RCA
// assessor state) plus the corners each class handles specially: the sentinel for an absent ion,
// the Debye-Hueckel floor, the steps floor, the temperature knees and warnings, the clamped
// setters, the efficiency cap and the NaN that passes through the assessor.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar ProcessScaleProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ProcessScaleProbe > captures/process_scale_probe.tsv

import java.util.Map;
import neqsim.process.chemistry.scale.BrineMixingScaleEvaluator;
import neqsim.process.chemistry.scale.BrineMixingScaleEvaluator.BrineComposition;
import neqsim.process.chemistry.scale.BrineMixingScaleEvaluator.MixPoint;
import neqsim.process.chemistry.scale.ElectrolyteScaleCalculator;
import neqsim.process.chemistry.scale.ScaleControlAssessor;
import neqsim.process.chemistry.scale.ScaleInhibitorPerformance;
import neqsim.process.chemistry.scale.ScaleInhibitorPerformance.InhibitorChemistry;
import neqsim.process.chemistry.scale.ScaleInhibitorPerformance.ScaleType;
import neqsim.process.chemistry.scale.ScaleKinetics;
import neqsim.pvtsimulation.flowassurance.ScalePredictionCalculator;

public class ProcessScaleProbe {

  static void electrolyteRow(String label, ElectrolyteScaleCalculator calc) {
    Map<String, Double> gamma = calc.getActivityCoefficients();
    System.out.println("  " + label + " ionic_strength=" + calc.getIonicStrength() + " debye_hueckel_a="
        + calc.getDebyeHueckelA() + " gamma_ca=" + gamma.get("Ca2+") + " gamma_na=" + gamma.get("Na+")
        + " gamma_so4=" + gamma.get("SO42-") + " gamma_co3=" + gamma.get("CO32-") + " caco3_si="
        + calc.getCaCO3SaturationIndex() + " baso4_si=" + calc.getBaSO4SaturationIndex() + " caso4_si="
        + calc.getCaSO4SaturationIndex() + " srso4_si=" + calc.getSrSO4SaturationIndex() + " brucite_si="
        + calc.getBruciteSaturationIndex() + " hydroxide_activity=" + calc.getHydroxideActivity() + " evaluated="
        + calc.isEvaluated());
  }

  static ElectrolyteScaleCalculator seawater(double tC, double pH) {
    return new ElectrolyteScaleCalculator().setTemperatureCelsius(tC).setPressureBara(1.013).setPH(pH)
        .setCations(412.0, 0.02, 7.9, 1290.0, 10780.0, 399.0, 0.002)
        .setAnions(19350.0, 2710.0, 142.0, 0.0);
  }

  static void brineMixingRow(String label, BrineComposition a, BrineComposition b, double tC, double p, double pH,
      int steps) {
    BrineMixingScaleEvaluator evaluator = new BrineMixingScaleEvaluator(a, b);
    evaluator.setConditions(tC, p).setPH(pH).setSteps(steps).evaluate();
    int index = 0;
    for (MixPoint point : evaluator.getPoints()) {
      System.out.println("  " + label + " point=" + index + " fraction_a=" + point.getFractionA() + " caco3_si="
          + point.getCaCO3SI() + " baso4_si=" + point.getBaSO4SI() + " caso4_si=" + point.getCaSO4SI()
          + " srso4_si=" + point.getSrSO4SI());
      index++;
    }
    System.out.println("  " + label + " worst fraction_a=" + evaluator.getWorstFractionA() + " mineral="
        + evaluator.getWorstMineral() + " si=" + evaluator.getWorstSI() + " points="
        + evaluator.getPoints().size());
  }

  static void inhibitorRow(String label, ScaleInhibitorPerformance perf) {
    perf.evaluate();
    Map<?, ?> inputs = (Map<?, ?>) perf.toMap().get("inputs");
    String warnings = String.join(",", perf.getWarnings().keySet());
    System.out.println("  " + label + " scale=" + inputs.get("scaleType") + " chemistry="
        + inputs.get("inhibitorChemistry") + " temperature_c=" + inputs.get("temperatureC") + " sr="
        + inputs.get("saturationRatio") + " tds_mg_l=" + inputs.get("tdsMgL") + " calcium_mg_l="
        + inputs.get("calciumMgL") + " dose_mg_l=" + inputs.get("availableDoseMgL") + " mic="
        + perf.getMinimumInhibitorConcentrationMgL() + " recommended=" + perf.getRecommendedDoseMgL()
        + " efficiency=" + perf.getEfficiency() + " adequate=" + perf.isAdequate() + " warnings="
        + (warnings.isEmpty() ? "-" : warnings));
  }

  static void kineticsRow(String label, ScaleKinetics kinetics) {
    kinetics.evaluate();
    System.out.println("  " + label + " s=" + kinetics.getSupersaturationRatio() + " induction_h="
        + kinetics.getInductionTimeHours() + " surface_mm_yr=" + kinetics.getSurfaceLimitedGrowthMmYr()
        + " transport_mm_yr=" + kinetics.getTransportLimitedGrowthMmYr() + " effective_mm_yr="
        + kinetics.getEffectiveGrowthRateMmYr() + " regime=" + kinetics.getLimitingRegime());
  }

  static void assessorRow(String label, ScaleControlAssessor assessor) {
    System.out.println("  " + label + " thermodynamic_caco3=" + assessor.getThermodynamicSaturationIndex(ScaleType.CACO3)
        + " risk_caco3=" + assessor.getKineticRiskIndex(ScaleType.CACO3) + " risk_baso4="
        + assessor.getKineticRiskIndex(ScaleType.BASO4) + " risk_srso4="
        + assessor.getKineticRiskIndex(ScaleType.SRSO4) + " risk_caso4="
        + assessor.getKineticRiskIndex(ScaleType.CASO4) + " risk_feco3="
        + assessor.getKineticRiskIndex(ScaleType.FECO3) + " worst_risk=" + assessor.getWorstKineticRiskIndex()
        + " controlled_0.5=" + assessor.isKineticallyControlled(0.5) + " evaluated=" + assessor.isEvaluated());
  }

  static ScalePredictionCalculator rcaPredictor() {
    ScalePredictionCalculator predictor = new ScalePredictionCalculator();
    predictor.setTemperatureCelsius(80.0);
    predictor.setPH(7.0);
    predictor.setCalciumConcentration(2000.0);
    predictor.setBicarbonateConcentration(800.0);
    predictor.setSulphateConcentration(500.0);
    predictor.setBariumConcentration(50.0);
    predictor.setSodiumConcentration(20000.0);
    predictor.setTotalDissolvedSolids(80000.0);
    predictor.calculate();
    return predictor;
  }

  static ScaleInhibitorPerformance rcaInhibitor(ScalePredictionCalculator predictor, double dose) {
    ScaleInhibitorPerformance inhibitor = new ScaleInhibitorPerformance();
    inhibitor.setScaleType(ScaleType.CACO3);
    inhibitor.setInhibitorChemistry(InhibitorChemistry.PHOSPHONATE);
    inhibitor.setTemperatureCelsius(80.0);
    inhibitor.setCalciumMgL(2000.0);
    inhibitor.setSaturationRatio(Math.pow(10.0, Math.max(0.0, predictor.getCaCO3SaturationIndex())));
    inhibitor.setTdsMgL(80000.0);
    inhibitor.setAvailableDoseMgL(dose);
    inhibitor.evaluate();
    return inhibitor;
  }

  static ScalePredictionCalculator skillPredictor() {
    ScalePredictionCalculator predictor = new ScalePredictionCalculator();
    predictor.setTemperatureCelsius(95.0);
    predictor.setPressureBara(150.0);
    predictor.setCalciumConcentration(1200.0);
    predictor.setBariumConcentration(250.0);
    predictor.setSulphateConcentration(1800.0);
    predictor.setBicarbonateConcentration(300.0);
    predictor.setTotalDissolvedSolids(90000.0);
    predictor.setPH(6.2);
    return predictor;
  }

  static ScaleInhibitorPerformance skillInhibitor() {
    ScaleInhibitorPerformance inhibitor = new ScaleInhibitorPerformance();
    inhibitor.setScaleType(ScaleType.BASO4);
    inhibitor.setInhibitorChemistry(InhibitorChemistry.PHOSPHONATE);
    inhibitor.setTemperatureCelsius(95.0);
    inhibitor.setSaturationRatio(12.0);
    inhibitor.setTdsMgL(90000.0);
    inhibitor.setCalciumMgL(1200.0);
    inhibitor.setAvailableDoseMgL(15.0);
    inhibitor.evaluate();
    return inhibitor;
  }

  public static void main(String[] args) {
    System.out.println("process_scale");

    System.out.println(" electrolyte_calculator");
    electrolyteRow("seawater_10c_ph8.1", seawater(10.0, 8.1).calculate());
    electrolyteRow("seawater_10c_ph9.0", seawater(10.0, 9.0).calculate());
    electrolyteRow("seawater_10c_ph10.0", seawater(10.0, 10.0).calculate());
    electrolyteRow("seawater_10c_ph10.5", seawater(10.0, 10.5).calculate());
    electrolyteRow("seawater_10c_ph14.0", seawater(10.0, 14.0).calculate());
    electrolyteRow("seawater_25c_ph14.0", seawater(25.0, 14.0).calculate());
    electrolyteRow("seawater_5c_ph10.0", seawater(5.0, 10.0).calculate());
    electrolyteRow("seawater_40c_ph10.0", seawater(40.0, 10.0).calculate());
    electrolyteRow("seawater_0c_ph10.0", seawater(0.0, 10.0).calculate());
    electrolyteRow("seawater_120c_ph10.0", seawater(120.0, 10.0).calculate());
    electrolyteRow("advanced_60c_50bara_ph6.5", new ElectrolyteScaleCalculator().setTemperatureCelsius(60.0)
        .setPressureBara(50.0).setPH(6.5).setCO2PartialPressureBar(2.0)
        .setCations(420.0, 0.05, 8.0, 1280.0, 10800.0, 400.0, 0.0).setAnions(19400.0, 2700.0, 142.0, 0.0)
        .calculate());
    electrolyteRow("pure_water_25c_ph7", new ElectrolyteScaleCalculator().setTemperatureCelsius(25.0)
        .setPressureBara(1.013).setPH(7.0).setCations(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        .setAnions(0.0, 0.0, 0.0, 0.0).calculate());
    electrolyteRow("absent_ions_25c_ph7", new ElectrolyteScaleCalculator().setTemperatureCelsius(25.0)
        .setPressureBara(1.013).setPH(7.0).setCations(400.0, 0.0, 0.0, 0.0, 10000.0, 0.0, 0.0)
        .setAnions(19000.0, 0.0, 0.0, 0.0).calculate());

    System.out.println(" brine_mixing");
    BrineComposition formation = new BrineComposition(5000.0, 300.0, 800.0, 500.0, 30000.0, 400.0, 55000.0, 0.0,
        400.0);
    BrineComposition injection = new BrineComposition(400.0, 0.0, 8.0, 1300.0, 11000.0, 400.0, 19000.0, 2700.0, 140.0);
    brineMixingRow("mixing_11", formation, injection, 60.0, 100.0, 6.0, 11);
    BrineComposition formation2 = new BrineComposition(2000.0, 500.0, 200.0, 300.0, 25000.0, 300.0, 45000.0, 0.0,
        300.0);
    BrineComposition injection2 = new BrineComposition(400.0, 0.0, 8.0, 1300.0, 11000.0, 400.0, 19000.0, 2800.0, 140.0);
    brineMixingRow("mixing_21", formation2, injection2, 80.0, 100.0, 6.2, 21);
    brineMixingRow("steps_floor", formation, injection, 60.0, 100.0, 6.0, 1);

    System.out.println(" inhibitor");
    ScaleInhibitorPerformance baso4Low = new ScaleInhibitorPerformance();
    baso4Low.setScaleType(ScaleType.BASO4);
    baso4Low.setInhibitorChemistry(InhibitorChemistry.PHOSPHONATE);
    baso4Low.setSaturationRatio(5.0);
    inhibitorRow("baso4_phosphonate_sr5", baso4Low);
    ScaleInhibitorPerformance baso4High = new ScaleInhibitorPerformance();
    baso4High.setScaleType(ScaleType.BASO4);
    baso4High.setInhibitorChemistry(InhibitorChemistry.PHOSPHONATE);
    baso4High.setSaturationRatio(50.0);
    inhibitorRow("baso4_phosphonate_sr50", baso4High);
    inhibitorRow("caco3_phosphonate_default", new ScaleInhibitorPerformance());
    ScaleInhibitorPerformance highCalcium = new ScaleInhibitorPerformance();
    highCalcium.setInhibitorChemistry(InhibitorChemistry.PHOSPHONATE);
    highCalcium.setCalciumMgL(5000.0);
    inhibitorRow("caco3_phosphonate_ca5000", highCalcium);
    ScaleInhibitorPerformance recommended = new ScaleInhibitorPerformance();
    recommended.setSaturationRatio(10.0);
    recommended.evaluate();
    recommended.setAvailableDoseMgL(recommended.getRecommendedDoseMgL());
    inhibitorRow("caco3_at_recommended_dose", recommended);
    inhibitorRow("barite_95c", skillInhibitor());
    ScaleInhibitorPerformance hot = new ScaleInhibitorPerformance();
    hot.setInhibitorChemistry(InhibitorChemistry.PHOSPHONATE);
    hot.setTemperatureCelsius(180.0);
    inhibitorRow("phosphonate_180c", hot);
    ScaleInhibitorPerformance coldEster = new ScaleInhibitorPerformance();
    coldEster.setInhibitorChemistry(InhibitorChemistry.PHOSPHATE_ESTER);
    coldEster.setTemperatureCelsius(30.0);
    inhibitorRow("phosphate_ester_30c", coldEster);
    ScaleInhibitorPerformance extreme = new ScaleInhibitorPerformance();
    extreme.setSaturationRatio(200.0);
    inhibitorRow("caco3_sr200", extreme);
    ScaleInhibitorPerformance polymaleate = new ScaleInhibitorPerformance();
    polymaleate.setInhibitorChemistry(InhibitorChemistry.POLYMALEATE);
    inhibitorRow("caco3_polymaleate", polymaleate);
    ScaleInhibitorPerformance polyacrylate = new ScaleInhibitorPerformance();
    polyacrylate.setScaleType(ScaleType.BASO4);
    polyacrylate.setInhibitorChemistry(InhibitorChemistry.POLYACRYLATE);
    polyacrylate.setCalciumMgL(5000.0);
    inhibitorRow("baso4_polyacrylate_ca5000", polyacrylate);
    ScaleInhibitorPerformance vinyl = new ScaleInhibitorPerformance();
    vinyl.setInhibitorChemistry(InhibitorChemistry.VINYL_SULPHONATE);
    vinyl.setTemperatureCelsius(130.0);
    inhibitorRow("caco3_vinyl_130c", vinyl);
    ScaleInhibitorPerformance saline = new ScaleInhibitorPerformance();
    saline.setTdsMgL(200000.0);
    inhibitorRow("caco3_tds200000", saline);
    ScaleInhibitorPerformance fresh = new ScaleInhibitorPerformance();
    fresh.setTdsMgL(1.0);
    inhibitorRow("caco3_tds1", fresh);
    ScaleInhibitorPerformance lowSr = new ScaleInhibitorPerformance();
    lowSr.setSaturationRatio(0.5);
    inhibitorRow("caco3_sr0.5", lowSr);
    ScaleInhibitorPerformance halfDose = new ScaleInhibitorPerformance();
    halfDose.setSaturationRatio(10.0);
    halfDose.evaluate();
    halfDose.setAvailableDoseMgL(0.5 * halfDose.getMinimumInhibitorConcentrationMgL());
    inhibitorRow("caco3_half_mic", halfDose);
    ScaleInhibitorPerformance overdosed = new ScaleInhibitorPerformance();
    overdosed.setSaturationRatio(10.0);
    overdosed.setAvailableDoseMgL(1.0e6);
    inhibitorRow("caco3_huge_dose", overdosed);

    System.out.println(" kinetics");
    kineticsRow("si_minus_0.5", new ScaleKinetics().setSaturationIndex(-0.5));
    kineticsRow("si_0.0", new ScaleKinetics().setSaturationIndex(0.0));
    kineticsRow("si_0.3", new ScaleKinetics().setSaturationIndex(0.3));
    kineticsRow("si_1.0", new ScaleKinetics().setSaturationIndex(1.0));
    kineticsRow("si_1.5", new ScaleKinetics().setSaturationIndex(1.5));
    kineticsRow("si_2.0", new ScaleKinetics().setSaturationIndex(2.0));
    kineticsRow("transport_limited", new ScaleKinetics().setSaturationIndex(2.0)
        .setSurfaceReaction(50.0, 2.0).setMassTransfer(1.0e-6, 1.0e-3, 0.0));
    kineticsRow("reaction_limited", new ScaleKinetics().setSaturationIndex(1.0)
        .setSurfaceReaction(5.0, 2.0).setMassTransfer(1.0e-3, 1.0e-3, 0.0));
    kineticsRow("custom_nucleation", new ScaleKinetics().setSaturationIndex(0.5).setNucleation(0.01, 5.0));
    kineticsRow("custom_order", new ScaleKinetics().setSaturationIndex(1.2).setSurfaceReaction(2.0, 3.0));
    kineticsRow("clamped_setters", new ScaleKinetics().setSaturationIndex(1.0).setSurfaceReaction(-5.0, -1.0)
        .setMassTransfer(-1.0, -1.0, -1.0).setNucleation(-1.0, -1.0));

    System.out.println(" control_assessor");
    ScalePredictionCalculator rca = rcaPredictor();
    ScaleControlAssessor rcaAssessor = new ScaleControlAssessor(rca);
    rcaAssessor.addInhibitor(ScaleType.CACO3, rcaInhibitor(rca, 20.0));
    rcaAssessor.evaluate();
    assessorRow("rca_state", rcaAssessor);

    ScalePredictionCalculator skill = skillPredictor();
    ScaleControlAssessor skillAssessor = new ScaleControlAssessor(skill);
    skillAssessor.addInhibitor(ScaleType.BASO4, skillInhibitor());
    skillAssessor.evaluate();
    assessorRow("skill_state", skillAssessor);

    ScalePredictionCalculator bare = rcaPredictor();
    ScaleControlAssessor bareAssessor = new ScaleControlAssessor(bare);
    bareAssessor.evaluate();
    assessorRow("no_inhibitor", bareAssessor);

    ScalePredictionCalculator zeroDose = rcaPredictor();
    ScaleControlAssessor zeroDoseAssessor = new ScaleControlAssessor(zeroDose);
    zeroDoseAssessor.addInhibitor(ScaleType.CACO3, rcaInhibitor(zeroDose, 0.0));
    zeroDoseAssessor.evaluate();
    assessorRow("zero_dose", zeroDoseAssessor);

    ScalePredictionCalculator hugeDose = rcaPredictor();
    ScaleControlAssessor hugeDoseAssessor = new ScaleControlAssessor(hugeDose);
    hugeDoseAssessor.addInhibitor(ScaleType.CACO3, rcaInhibitor(hugeDose, 1.0e6));
    hugeDoseAssessor.evaluate();
    assessorRow("huge_dose", hugeDoseAssessor);

    ScalePredictionCalculator notEvaluated = rcaPredictor();
    ScaleControlAssessor notEvaluatedAssessor = new ScaleControlAssessor(notEvaluated);
    System.out.println("  not_evaluated evaluated=" + notEvaluatedAssessor.isEvaluated() + " worst_risk="
        + notEvaluatedAssessor.getWorstKineticRiskIndex());
  }
}
