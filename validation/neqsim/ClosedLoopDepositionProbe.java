// ClosedLoopDepositionProbe - the coupled deposition-flow loop at the pinned NeqSim commit.
//
// `ClosedLoopDepositionSolver` re-runs a `PipeBeggsAndBrills` at a shrinking effective diameter,
// evaluates a bound `ScaleDepositionAccumulator` after each run, and sets the next diameter from
// the worst segment's thickness, stopping on a diameter tolerance, the iteration cap or a
// blocked line. The probe prints each solve's iteration count, convergence, final diameter and
// the three histories across the corners: a one-year service where the first step already meets
// the tolerance, a long service that iterates, the cap reached before convergence, a blocked
// line, a brine that deposits nothing, an inhibited line, a tolerance wider than the step, and
// a Peng-Robinson fluid.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar ClosedLoopDepositionProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ClosedLoopDepositionProbe > captures/closed_loop_deposition_probe.tsv

import java.util.List;
import neqsim.process.chemistry.scale.ClosedLoopDepositionSolver;
import neqsim.process.chemistry.scale.ScaleDepositionAccumulator;
import neqsim.process.equipment.pipeline.PipeBeggsAndBrills;
import neqsim.process.equipment.stream.Stream;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;

public class ClosedLoopDepositionProbe {

  static PipeBeggsAndBrills pipe(String eos) {
    SystemInterface fluid = eos.equals("pr") ? new SystemPrEos(333.15, 100.0) : new SystemSrkEos(333.15, 100.0);
    fluid.addComponent("methane", 0.85);
    fluid.addComponent("CO2", 0.05);
    fluid.addComponent("water", 0.10);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    Stream inlet = new Stream("inlet", fluid);
    inlet.setFlowRate(80.0, "mol/sec");
    inlet.setTemperature(60.0, "C");
    inlet.setPressure(100.0, "bara");
    inlet.run();
    PipeBeggsAndBrills pipe = new PipeBeggsAndBrills("pipe", inlet);
    pipe.setLength(1000.0);
    pipe.setDiameter(0.1);
    pipe.setNumberOfIncrements(10);
    pipe.setElevation(0.0);
    pipe.setPipeWallRoughness(1.0e-5);
    return pipe;
  }

  static void row(String label, String eos, double serviceYears, double inhibitorEfficiency,
      double calciumMgL, double ph, double toleranceM, int maxIterations) {
    PipeBeggsAndBrills pipe = pipe(eos);
    ScaleDepositionAccumulator accumulator = new ScaleDepositionAccumulator(pipe);
    accumulator.setBrineChemistry(calciumMgL, 500.0, 100.0, 10.0, 20000.0, 50000.0);
    accumulator.setpHAndCo2(ph, 2.0);
    accumulator.setInhibitorEfficiency(inhibitorEfficiency);
    accumulator.setServiceYears(serviceYears);
    accumulator.setScaleDensity(2700.0);
    accumulator.setBaseRate(0.5);
    ClosedLoopDepositionSolver solver = new ClosedLoopDepositionSolver(pipe, accumulator)
        .setToleranceM(toleranceM)
        .setMaxIterations(maxIterations);
    solver.solve();
    System.out.println("state=" + label
        + " iterations=" + solver.getIterationsTaken()
        + " converged=" + solver.isConverged()
        + " final_diameter_m=" + solver.getFinalEffectiveDiameterM()
        + " tolerance_m=" + toleranceM
        + " max_iterations=" + maxIterations
        + " last_max_thickness_mm=" + accumulator.getMaxThicknessMm()
        + " last_time_to_blockage_years=" + accumulator.getTimeToBlockageYears());
    System.out.println("  diameter_history_m=" + join(solver.getDiameterHistoryM())
        + " thickness_history_mm=" + join(solver.getMaxThicknessHistoryMm())
        + " velocity_history_m_s=" + join(solver.getVelocityHistoryMs()));
    List<Double> pressures = pipe.getPressureProfileList();
    List<Double> mixtureVelocity = pipe.getMixtureSuperficialVelocityProfile();
    System.out.println("  restored_diameter_m=" + pipe.getDiameter()
        + " profile_size=" + Math.min(pipe.getTemperatureProfileList().size(), pressures.size())
        + " velocity_profile_entries=" + mixtureVelocity.size()
        + " segment0_pressure_bara=" + pipe.getSegmentPressure(0)
        + " segment0_temperature_k=" + pipe.getSegmentTemperature(0));
  }

  static String join(List<Double> values) {
    StringBuilder text = new StringBuilder("[");
    for (int i = 0; i < values.size(); i++) {
      if (i > 0) {
        text.append(",");
      }
      text.append(values.get(i));
    }
    return text.append("]").toString();
  }

  public static void main(String[] args) {
    System.out.println("ClosedLoopDepositionSolver at f0c7436c69 - 80 mol/sec, 100 bara, 60 C, 1000 m, 0.1 m, 10 increments");
    row("base_1yr", "srk", 1.0, 0.0, 1000.0, 6.5, 1.0e-4, 10);
    row("service_1000yr", "srk", 1000.0, 0.0, 1000.0, 6.5, 1.0e-4, 10);
    row("max_iterations_2", "srk", 1000.0, 0.0, 1000.0, 6.5, 1.0e-4, 2);
    row("blocked", "srk", 1.0e7, 0.0, 1000.0, 6.5, 1.0e-4, 10);
    row("no_scale", "srk", 1000.0, 0.0, 10.0, 4.0, 1.0e-4, 10);
    row("inhibited_0_5", "srk", 1000.0, 0.5, 1000.0, 6.5, 1.0e-4, 10);
    row("large_tolerance", "srk", 1000.0, 0.0, 1000.0, 6.5, 0.02, 10);
    row("pr_eos", "pr", 1000.0, 0.0, 1000.0, 6.5, 1.0e-4, 10);
  }
}
