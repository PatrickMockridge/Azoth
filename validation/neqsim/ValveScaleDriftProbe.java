// ValveScaleDriftProbe - the valve trim-fouling drift at the pinned NeqSim commit.
//
// `ValveScaleDrift` binds a `ThrottlingValve` and grows a uniform radial deposit on its
// trim: the open-area fraction is `((d0 - 2t)/d0)^2`, the fouling fraction its complement,
// and the effective Kv the clean Kv times the open area. Every getter is a function of the
// accumulated thickness, so the probe prints a state per row: the thickness, open area,
// fouling, effective Kv, time to plug, and the opening and time-to-pin for a stated clean
// opening. The valve is constructed without a stream and given a Kv directly; nothing in
// the drift arithmetic reads a fluid.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar ValveScaleDriftProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar ValveScaleDriftProbe > captures/valve_scale_drift_probe.tsv

import neqsim.process.chemistry.scale.ScaleKinetics;
import neqsim.process.chemistry.scale.ValveScaleDrift;
import neqsim.process.equipment.valve.ThrottlingValve;

public class ValveScaleDriftProbe {

  static ThrottlingValve valve(double kv) {
    ThrottlingValve valve = new ThrottlingValve("trim");
    valve.setKv(kv);
    return valve;
  }

  static ValveScaleDrift drift(double portDiameterM, double growthMmPerYear, double kv) {
    ValveScaleDrift drift = new ValveScaleDrift(valve(kv));
    drift.setPortDiameter(portDiameterM);
    drift.setGrowthRateMmPerYear(growthMmPerYear);
    return drift;
  }

  static void row(String label, ValveScaleDrift drift, double portDiameterM,
      double cleanOpeningPercent) {
    System.out.println("  " + label + " port_diameter_m=" + portDiameterM
        + " deposit_thickness_mm=" + drift.getDepositThicknessMm() + " elapsed_days="
        + drift.getElapsedDays() + " open_area=" + drift.getOpenAreaFraction() + " fouling="
        + drift.getFoulingFraction() + " effective_kv=" + drift.getEffectiveKv() + " time_to_plug_days="
        + drift.getTimeToPlugDays() + " clean_opening_percent=" + cleanOpeningPercent
        + " opening_percent=" + drift.predictOpeningPercent(cleanOpeningPercent) + " time_to_pin_days="
        + drift.predictTimeToPinDays(cleanOpeningPercent));
  }

  public static void main(String[] args) {
    System.out.println("valve_scale_drift");

    System.out.println(" growth");
    ValveScaleDrift thirty = drift(0.05, 1.0, 100.0);
    thirty.advance(30.0);
    row("thirty_days", thirty, 0.05, 50.0);

    ValveScaleDrift year = drift(0.05, 1.0, 100.0);
    year.advance(365.25);
    row("one_year", year, 0.05, 50.0);

    ValveScaleDrift pinned = drift(0.05, 1.0, 100.0);
    pinned.advance(10000.0);
    row("pinned", pinned, 0.05, 50.0);

    System.out.println(" corners");
    ValveScaleDrift clean = drift(0.05, 1.0, 100.0);
    row("not_advanced", clean, 0.05, 50.0);

    ValveScaleDrift noGrowth = drift(0.05, 0.0, 100.0);
    noGrowth.advance(365.25);
    row("zero_growth", noGrowth, 0.05, 50.0);

    ValveScaleDrift noPort = drift(0.0, 1.0, 100.0);
    noPort.advance(365.25);
    row("no_port", noPort, 0.0, 50.0);

    ValveScaleDrift nearlyPinned = drift(0.05, 1.0, 100.0);
    nearlyPinned.advance(10000.0);
    row("clean_opening_100", nearlyPinned, 0.05, 100.0);

    ValveScaleDrift negative = drift(0.05, 1.0, 100.0);
    negative.advance(30.0);
    negative.advance(-10.0);
    row("negative_advance_ignored", negative, 0.05, 50.0);

    ValveScaleDrift reset = drift(0.05, 1.0, 100.0);
    reset.advance(30.0);
    reset.reset();
    row("reset", reset, 0.05, 50.0);

    System.out.println(" from_kinetics");
    ScaleKinetics kinetics = new ScaleKinetics().setSaturationIndex(1.0).evaluate();
    ValveScaleDrift fromKinetics = new ValveScaleDrift(valve(100.0));
    fromKinetics.setPortDiameter(0.05);
    fromKinetics.setKinetics(kinetics);
    fromKinetics.advance(30.0);
    row("kinetics_si_1.0", fromKinetics, 0.05, 50.0);
  }
}
