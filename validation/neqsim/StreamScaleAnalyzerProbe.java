// StreamScaleAnalyzerProbe - the stream scale analyser at the pinned NeqSim commit.
//
// The class is `StreamScaleAnalyzer` over a configured `ScalePredictionCalculator`: the
// coupled `MultiMineralScaleEquilibrium` at one litre of water, and the per-litre
// precipitated masses converted to a kg/day rate by the aqueous throughput. The states are
// the equilibrium capture's brines so the two layers can be compared directly, plus the
// zero-flow, no-scale, auto-pH and absent-ion corners.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar StreamScaleAnalyzerProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar StreamScaleAnalyzerProbe > captures/stream_scale_analyzer_probe.tsv

import java.util.Map;
import neqsim.process.chemistry.scale.StreamScaleAnalyzer;
import neqsim.pvtsimulation.flowassurance.MultiMineralScaleEquilibrium;
import neqsim.pvtsimulation.flowassurance.MultiMineralScaleEquilibrium.MineralResult;
import neqsim.pvtsimulation.flowassurance.ScalePredictionCalculator;

public class StreamScaleAnalyzerProbe {

  static ScalePredictionCalculator predictor(double calcium, double barium, double strontium,
      double bicarbonate, double sulphate, double tds, double temperatureC, double pressureBara,
      double co2pp, double ph, boolean autoPh) {
    ScalePredictionCalculator p = new ScalePredictionCalculator();
    p.setCalciumConcentration(calcium);
    p.setBariumConcentration(barium);
    p.setStrontiumConcentration(strontium);
    p.setBicarbonateConcentration(bicarbonate);
    p.setSulphateConcentration(sulphate);
    p.setTotalDissolvedSolids(tds);
    p.setTemperatureCelsius(temperatureC);
    p.setPressureBara(pressureBara);
    p.setCO2PartialPressure(co2pp);
    p.setPH(ph);
    if (autoPh) {
      p.enableAutoPH();
    }
    return p;
  }

  static void row(String label, ScalePredictionCalculator predictor,
      MultiMineralScaleEquilibrium.ActivityModel model, double flowLitrePerDay) {
    StreamScaleAnalyzer analyzer = new StreamScaleAnalyzer(predictor);
    analyzer.setActivityModel(model);
    analyzer.setWaterFlowLitrePerDay(flowLitrePerDay);
    analyzer.analyze();
    System.out.println("  " + label + " model=" + model + " flow_l_per_day=" + flowLitrePerDay
        + " total_mg_per_l=" + analyzer.getTotalScaleMassMgPerL() + " total_kg_per_day="
        + analyzer.getTotalScaleRateKgPerDay() + " dominant=" + analyzer.getDominantScale()
        + " has_scaling=" + analyzer.hasScaling());
    for (Map.Entry<String, MineralResult> entry : analyzer.getEquilibrium().getResults().entrySet()) {
      MineralResult result = entry.getValue();
      System.out.println("    " + entry.getKey() + " initial_si=" + result.getInitialSI()
          + " final_si=" + result.getFinalSI() + " mg_per_l=" + result.getPrecipitatedMassMgPerL()
          + " kg_per_day=" + analyzer.getScaleRateKgPerDay(entry.getKey()));
    }
  }

  public static void main(String[] args) {
    System.out.println("stream_scale_analyzer");

    System.out.println(" sulphate_limited");
    row("davies", predictor(1000.0, 1000.0, 3000.0, 200.0, 1500.0, 60000.0, 70.0, 100.0, 1.0, 0.0,
        false), MultiMineralScaleEquilibrium.ActivityModel.DAVIES, 100000.0);
    row("bdot", predictor(1000.0, 1000.0, 3000.0, 200.0, 1500.0, 60000.0, 70.0, 100.0, 1.0, 0.0,
        false), MultiMineralScaleEquilibrium.ActivityModel.BDOT, 100000.0);

    System.out.println(" co_precipitation");
    row("davies", predictor(500.0, 100.0, 2500.0, 100.0, 4000.0, 60000.0, 70.0, 100.0, 1.0, 0.0,
        false), MultiMineralScaleEquilibrium.ActivityModel.DAVIES, 250000.0);

    System.out.println(" corners");
    row("zero_flow", predictor(1000.0, 1000.0, 3000.0, 200.0, 1500.0, 60000.0, 70.0, 100.0, 1.0,
        0.0, false), MultiMineralScaleEquilibrium.ActivityModel.DAVIES, 0.0);
    row("no_minerals", predictor(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 25.0, 1.013, 1.0, 7.0, false),
        MultiMineralScaleEquilibrium.ActivityModel.DAVIES, 1000.0);
    row("absent_barium", predictor(1000.0, 0.0, 0.0, 200.0, 1500.0, 60000.0, 70.0, 100.0, 1.0,
        0.0, false), MultiMineralScaleEquilibrium.ActivityModel.DAVIES, 100000.0);
    row("auto_ph", predictor(1000.0, 1000.0, 3000.0, 200.0, 1500.0, 60000.0, 70.0, 100.0, 1.0,
        6.5, true), MultiMineralScaleEquilibrium.ActivityModel.DAVIES, 100000.0);
  }
}
