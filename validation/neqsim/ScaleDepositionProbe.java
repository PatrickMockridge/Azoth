// ScaleDepositionProbe - the segmented deposition accumulator at the pinned NeqSim commit.
//
// `ScaleDepositionAccumulator` walks a `PipeBeggsAndBrills` segment by segment, evaluates the
// scale predictor's CaCO3 index at each segment's temperature and pressure, and integrates a
// deposit mass and wall thickness from `max(0, S - 1) * base rate`, then extrapolates the
// time to half-block the worst segment. The probe first runs a real Beggs-and-Brills pipe and
// prints its profile - those segment temperatures and pressures are the state the port's
// model takes directly - then prints the accumulator's per-segment indices, masses and
// thicknesses and its summaries, across the inhibitor, service-time, base-rate and
// no-scaling corners.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar ScaleDepositionProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ScaleDepositionProbe > captures/scale_deposition_probe.tsv

import java.util.List;
import neqsim.process.chemistry.scale.ScaleDepositionAccumulator;
import neqsim.process.equipment.pipeline.PipeBeggsAndBrills;
import neqsim.process.equipment.stream.Stream;
import neqsim.thermo.system.SystemSrkEos;

public class ScaleDepositionProbe {

  static PipeBeggsAndBrills cooledPipe() {
    SystemSrkEos fluid = new SystemSrkEos(333.15, 100.0);
    fluid.addComponent("methane", 0.85);
    fluid.addComponent("CO2", 0.05);
    fluid.addComponent("water", 0.10);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    Stream inlet = new Stream("inlet", fluid);
    inlet.setFlowRate(5000.0, "kg/hr");
    inlet.setTemperature(60.0, "C");
    inlet.setPressure(100.0, "bara");
    inlet.run();
    PipeBeggsAndBrills pipe = new PipeBeggsAndBrills("cooled", inlet);
    pipe.setLength(1000.0);
    pipe.setDiameter(0.1);
    pipe.setNumberOfIncrements(10);
    pipe.setElevation(0.0);
    pipe.setConstantSurfaceTemperature(20.0, "C");
    pipe.run();
    return pipe;
  }

  static PipeBeggsAndBrills pipe() {
    SystemSrkEos fluid = new SystemSrkEos(333.15, 100.0);
    fluid.addComponent("methane", 0.85);
    fluid.addComponent("CO2", 0.05);
    fluid.addComponent("water", 0.10);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    Stream inlet = new Stream("inlet", fluid);
    inlet.setFlowRate(5000.0, "kg/hr");
    inlet.setTemperature(60.0, "C");
    inlet.setPressure(100.0, "bara");
    inlet.run();
    PipeBeggsAndBrills pipe = new PipeBeggsAndBrills("pipe", inlet);
    pipe.setLength(1000.0);
    pipe.setDiameter(0.1);
    pipe.setNumberOfIncrements(10);
    pipe.setElevation(0.0);
    pipe.run();
    return pipe;
  }

  static ScaleDepositionAccumulator accumulator(PipeBeggsAndBrills pipe) {
    ScaleDepositionAccumulator accumulator = new ScaleDepositionAccumulator(pipe);
    accumulator.setBrineChemistry(1000.0, 500.0, 100.0, 10.0, 20000.0, 50000.0);
    accumulator.setpHAndCo2(6.5, 2.0);
    accumulator.setInhibitorEfficiency(0.0);
    accumulator.setServiceYears(1.0);
    accumulator.setScaleDensity(2700.0);
    accumulator.setBaseRate(0.5);
    return accumulator;
  }

  static void row(String label, ScaleDepositionAccumulator accumulator) {
    accumulator.evaluate();
    System.out.println("  " + label + " total_mass_kg=" + accumulator.getTotalDepositionMassKg()
        + " max_thickness_mm=" + accumulator.getMaxThicknessMm() + " time_to_blockage_years="
        + accumulator.getTimeToBlockageYears() + " segments="
        + accumulator.getSegmentSaturationIndex().size());
    List<Double> si = accumulator.getSegmentSaturationIndex();
    List<Double> mass = accumulator.getSegmentDepositionMassKg();
    List<Double> thickness = accumulator.getSegmentThicknessMm();
    for (int i = 0; i < si.size(); i++) {
      System.out.println("    segment=" + i + " si=" + si.get(i) + " mass_kg=" + mass.get(i)
          + " thickness_mm=" + thickness.get(i));
    }
  }

  public static void main(String[] args) {
    System.out.println("scale_deposition");
    PipeBeggsAndBrills pipe = pipe();

    System.out.println(" profile");
    List<Double> temperatures = pipe.getTemperatureProfileList();
    List<Double> pressures = pipe.getPressureProfileList();
    System.out.println("  pipe diameter_m=" + pipe.getDiameter() + " length_m=" + pipe.getLength()
        + " increments=" + pipe.getNumberOfIncrements() + " profile_size=" + temperatures.size()
        + " pressure_size=" + pressures.size());
    for (int i = 0; i < temperatures.size(); i++) {
      System.out.println("    segment=" + i + " temperature=" + temperatures.get(i)
          + " pressure=" + pressures.get(i));
    }

    System.out.println(" accumulation");
    row("base", accumulator(pipe));

    ScaleDepositionAccumulator inhibited = accumulator(pipe);
    inhibited.setInhibitorEfficiency(0.5);
    row("inhibitor_0.5", inhibited);

    ScaleDepositionAccumulator longService = accumulator(pipe);
    longService.setServiceYears(5.0);
    row("service_5yr", longService);

    ScaleDepositionAccumulator faster = accumulator(pipe);
    faster.setBaseRate(1.0);
    row("base_rate_1.0", faster);

    System.out.println(" cooled_profile");
    PipeBeggsAndBrills cooled = cooledPipe();
    List<Double> cooledTemperatures = cooled.getTemperatureProfileList();
    List<Double> cooledPressures = cooled.getPressureProfileList();
    System.out.println("  profile_size=" + cooledTemperatures.size() + " pressure_size="
        + cooledPressures.size());
    for (int i = 0; i < cooledTemperatures.size(); i++) {
      System.out.println("    segment=" + i + " temperature=" + cooledTemperatures.get(i)
          + " pressure=" + (i < cooledPressures.size() ? cooledPressures.get(i) : "absent"));
    }
    row("cooled", accumulator(cooled));

    ScaleDepositionAccumulator noScale = new ScaleDepositionAccumulator(pipe);
    noScale.setBrineChemistry(10.0, 10.0, 100.0, 0.0, 20000.0, 50000.0);
    noScale.setpHAndCo2(4.0, 0.0);
    row("no_calcite", noScale);
  }
}
