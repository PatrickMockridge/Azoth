// PipeBeggsAndBrillsProbe - the Beggs-and-Brill pipe's per-segment state at the pinned commit.
//
// The class marches a pipeline in increments, flashing the fluid at each segment's pressure
// and temperature, and at every segment classifies the flow regime, computes the liquid
// holdup with its inclination correction, the no-slip density and viscosity, the Reynolds
// number, the friction factor and the pressure drop. The probe prints that per-segment
// state, plus the outlet and the total drop, over the states the correlation's branches
// need: a horizontal two-phase line, uphill and downhill, a wider pipe, a liquid-rich line,
// a three-phase oil/water/gas line, a single-phase gas and a single-phase liquid, and a
// rough wall.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar PipeBeggsAndBrillsProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PipeBeggsAndBrillsProbe > captures/pipe_beggs_and_brills_probe.tsv

import java.util.List;
import neqsim.process.equipment.pipeline.PipeBeggsAndBrills;
import neqsim.process.equipment.stream.Stream;
import neqsim.thermo.system.SystemSrkEos;

public class PipeBeggsAndBrillsProbe {

  static PipeBeggsAndBrills pipe(String[] names, double[] fractions, double flowMolSec,
      double diameter, double length, double elevation, double roughness) {
    SystemSrkEos fluid = new SystemSrkEos(333.15, 100.0);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], fractions[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    Stream inlet = new Stream("inlet", fluid);
    inlet.setFlowRate(flowMolSec, "mol/sec");
    inlet.setTemperature(60.0, "C");
    inlet.setPressure(100.0, "bara");
    inlet.run();
    PipeBeggsAndBrills pipe = new PipeBeggsAndBrills("pipe", inlet);
    pipe.setLength(length);
    pipe.setDiameter(diameter);
    pipe.setNumberOfIncrements(10);
    pipe.setElevation(elevation);
    pipe.setPipeWallRoughness(roughness);
    return pipe;
  }

  static void row(String label, PipeBeggsAndBrills pipe) {
    try {
      pipe.run();
    } catch (RuntimeException failure) {
      System.out.println("  " + label + " status=refused reason="
          + failure.getClass().getSimpleName() + ":" + failure.getMessage());
      return;
    }
    System.out.println("  " + label + " outlet_pressure_bara=" + pipe.getOutletPressure()
        + " outlet_temperature=" + pipe.getOutletStream().getTemperature("K")
        + " inlet_minus_outlet_bara=" + (pipe.getInletStream().getPressure("bara")
            - pipe.getOutletStream().getPressure("bara"))
        + " outlet_flow=" + pipe.getOutletStream().getFlowRate("mol/sec"));
    List<Double> pressures = pipe.getPressureProfileList();
    List<Double> temperatures = pipe.getTemperatureProfileList();
    List<Double> drops = pipe.getPressureDropProfile();
    List<Double> holdup = pipe.getLiquidHoldupProfileList();
    List<Double> liquidVelocity = pipe.getLiquidSuperficialVelocityProfile();
    List<Double> gasVelocity = pipe.getGasSuperficialVelocityProfile();
    List<Double> mixtureVelocity = pipe.getMixtureSuperficialVelocityProfile();
    List<Double> density = pipe.getMixtureDensityProfile();
    List<Double> viscosity = pipe.getMixtureViscosityProfile();
    List<Double> reynolds = pipe.getMixtureReynoldsNumber();
    for (int i = 0; i < pressures.size(); i++) {
      System.out.println("    segment=" + i + " pressure_bara=" + pressures.get(i)
          + " temperature_k=" + (i < temperatures.size() ? temperatures.get(i) : "absent")
          + " regime=" + (i < pressures.size() ? pipe.getSegmentFlowRegime(i) : "absent")
          + " holdup=" + (i < holdup.size() ? holdup.get(i) : "absent")
          + " sup_liquid_m_s=" + (i < liquidVelocity.size() ? liquidVelocity.get(i) : "absent")
          + " sup_gas_m_s=" + (i < gasVelocity.size() ? gasVelocity.get(i) : "absent")
          + " sup_mix_m_s=" + (i < mixtureVelocity.size() ? mixtureVelocity.get(i) : "absent")
          + " mixture_density_kg_m3=" + (i < density.size() ? density.get(i) : "absent")
          + " mixture_viscosity_cp=" + (i < viscosity.size() ? viscosity.get(i) : "absent")
          + " reynolds=" + (i < reynolds.size() ? reynolds.get(i) : "absent")
          + " pressure_drop=" + (i < drops.size() ? drops.get(i) : "absent"));
    }
  }

  static void correlationRow(String label, PipeBeggsAndBrills pipe) {
    try {
      pipe.run();
    } catch (RuntimeException failure) {
      System.out.println("  " + label + " status=refused reason="
          + failure.getClass().getSimpleName());
      return;
    }
    var system = pipe.getThermoSystem();
    System.out.println("  " + label + " phases=" + system.getNumberOfPhases() + " temperature="
        + system.getTemperature("K") + " pressure=" + system.getPressure("Pa") + " phase0_density="
        + system.getPhase(0).getDensity("kg/m3") + " phase0_viscosity_cp="
        + system.getPhase(0).getViscosity("cP") + " phase0_flow="
        + system.getPhase(0).getFlowRate("m3/sec")
        + (system.getNumberOfPhases() > 1
            ? " phase1_density=" + system.getPhase(1).getDensity("kg/m3") + " phase1_viscosity_cp="
                + system.getPhase(1).getViscosity("cP") + " phase1_flow="
                + system.getPhase(1).getFlowRate("m3/sec")
            : ""));
    List<Double> pressures = pipe.getPressureProfileList();
    int last = pressures.size() - 1;
    List<Double> holdup = pipe.getLiquidHoldupProfileList();
    List<Double> density = pipe.getMixtureDensityProfile();
    List<Double> viscosity = pipe.getMixtureViscosityProfile();
    List<Double> reynolds = pipe.getMixtureReynoldsNumber();
    List<Double> drops = pipe.getPressureDropProfile();
    List<Double> liquidVelocity = pipe.getLiquidSuperficialVelocityProfile();
    List<Double> gasVelocity = pipe.getGasSuperficialVelocityProfile();
    List<Double> mixtureVelocity = pipe.getMixtureSuperficialVelocityProfile();
    System.out.println("    outlet regime=" + pipe.getSegmentFlowRegime(last) + " holdup="
        + holdup.get(last) + " mixture_density=" + density.get(last) + " mixture_viscosity_cp="
        + viscosity.get(last) + " reynolds=" + reynolds.get(last) + " pressure_drop="
        + drops.get(last) + " sup_liquid_m_s=" + liquidVelocity.get(last) + " sup_gas_m_s="
        + gasVelocity.get(last) + " sup_mix_m_s=" + mixtureVelocity.get(last)
        + " diameter_m=" + pipe.getDiameter() + " length_m=" + pipe.getLength());
  }

  static PipeBeggsAndBrills segment(String[] names, double[] fractions, double flowMolSec,
      double diameter, double elevation, double roughness) {
    PipeBeggsAndBrills pipe = pipe(names, fractions, flowMolSec, diameter, 1.0, elevation, roughness);
    return pipe;
  }

  public static void main(String[] args) {
    System.out.println("pipe_beggs_and_brills");

    String[] wetGas = {"methane", "CO2", "water"};
    double[] wetFractions = {0.85, 0.05, 0.10};

    System.out.println(" horizontal");
    row("two_phase_0.1m", pipe(wetGas, wetFractions, 80.0, 0.1, 1000.0, 0.0, 1e-5));
    row("two_phase_0.2m", pipe(wetGas, wetFractions, 80.0, 0.2, 1000.0, 0.0, 1e-5));

    System.out.println(" inclined");
    row("uphill_100m", pipe(wetGas, wetFractions, 80.0, 0.1, 1000.0, 100.0, 1e-5));
    row("downhill_100m", pipe(wetGas, wetFractions, 80.0, 0.1, 1000.0, -100.0, 1e-5));

    System.out.println(" liquid_rich");
    row("half_water", pipe(wetGas, new double[] {0.5, 0.05, 0.45}, 80.0, 0.1, 1000.0, 0.0, 1e-5));

    System.out.println(" three_phase");
    row("oil_water_gas", pipe(new String[] {"methane", "n-octane", "water"},
        new double[] {0.5, 0.25, 0.25}, 80.0, 0.1, 1000.0, 0.0, 1e-5));

    System.out.println(" single_phase");
    row("gas_only", pipe(new String[] {"methane"}, new double[] {1.0}, 80.0, 0.1, 1000.0, 0.0, 1e-5));
    row("water_only", pipe(new String[] {"water"}, new double[] {1.0}, 80.0, 0.1, 1000.0, 0.0, 1e-5));

    System.out.println(" rough_wall");
    row("rough_1e-3", pipe(wetGas, wetFractions, 80.0, 0.1, 1000.0, 0.0, 1e-3));

    System.out.println(" correlation");
    correlationRow("segment_two_phase", segment(wetGas, wetFractions, 80.0, 0.1, 0.0, 1e-5));
    correlationRow("segment_uphill", segment(wetGas, wetFractions, 80.0, 0.1, 0.1, 1e-5));
    correlationRow("segment_downhill", segment(wetGas, wetFractions, 80.0, 0.1, -0.1, 1e-5));
    correlationRow("segment_gas_only", segment(new String[] {"methane"}, new double[] {1.0},
        80.0, 0.1, 0.0, 1e-5));
    correlationRow("segment_water_only", segment(new String[] {"water"}, new double[] {1.0},
        80.0, 0.1, 0.0, 1e-5));
    correlationRow("segment_three_phase", segment(new String[] {"methane", "n-octane", "water"},
        new double[] {0.5, 0.25, 0.25}, 80.0, 0.1, 0.0, 1e-5));
    correlationRow("segment_liquid_rich", segment(wetGas, new double[] {0.2, 0.05, 0.75},
        320.0, 0.1, 0.0, 1e-5));

    System.out.println(" regime_branches");
    row("distributed_high_flow",
        pipe(wetGas, wetFractions, 1600.0, 0.1, 1000.0, 0.0, 1e-5));
    row("intermittent_liquid_rich",
        pipe(wetGas, new double[] {0.5, 0.05, 0.45}, 800.0, 0.1, 1000.0, 0.0, 1e-5));
    row("slug_liquid_dominant",
        pipe(wetGas, new double[] {0.2, 0.05, 0.75}, 320.0, 0.1, 1000.0, 0.0, 1e-5));
  }
}
