// What a `Heater` reports about itself when two specifications are set, in either order.
//
// `ProcessProbe heat_exchanger` covers `Heater` one specification at a time, which is the
// shape the port's cases are pinned to. It cannot show the defect equinor/neqsim#4087
// reports, because that defect is about what the object *says* after the second setter has
// silently dropped the first - and a single-specification row never asks.
//
// Each row here sets two specifications and then reads back three things: the state the run
// actually reached, the duty it actually applied, and `getSpecifiedOutletTemperature()`. The
// row is a defect when the readback names a specification the run did not honour.
//
// The two temperature overloads are both exercised, and the zero-Celsius case is the one
// that would catch a `0.0` sentinel for "no temperature specification" - 0 degC is a
// legitimate setpoint and 273.15 K is not a value a guard should refuse.
//
//     javac -proc:none -cp neqsim-f0c7436.jar HeaterSpecProbe.java
//     java -Xmx2g -cp .:neqsim-f0c7436.jar HeaterSpecProbe > captures/heater_spec_probe.tsv

import neqsim.process.equipment.heatexchanger.Heater;
import neqsim.process.equipment.stream.Stream;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

public class HeaterSpecProbe {

  /** The feed the `process_heater.tsv` capture uses, so the states are comparable. */
  static Stream feed() {
    SystemInterface fluid = new SystemPrEos(320.0, 30.0);
    fluid.addComponent("methane", 0.9);
    fluid.addComponent("n-butane", 0.1);
    fluid.setMixingRule(2);
    Stream stream = new Stream("feed", fluid);
    stream.setFlowRate(1.0, "mol/sec");
    stream.run();
    return stream;
  }

  interface Spec {
    void apply(Heater heater);
  }

  static void row(String label, Spec spec) {
    Heater heater = new Heater("h1", feed());
    spec.apply(heater);
    heater.run();
    System.out.println(label);
    System.out.printf("outlet_T_K=%.15g%n", heater.getOutletStream().getTemperature());
    System.out.printf("duty_W=%.15g%n", heater.getDuty());
    System.out.printf("reported_specified_T_K=%.15g%n", heater.getSpecifiedOutletTemperature());
    System.out.printf("reported_specified_T_unit=%s%n", heater.getSpecifiedOutletTemperatureUnit());
    System.out.printf("isSetEnergyInput=%b%n", heater.isSetEnergyInput());
    System.out.println();
  }

  /** The same row run twice, because a dropped specification must not come back. */
  static void twice(String label, Spec spec) {
    Heater heater = new Heater("h1", feed());
    spec.apply(heater);
    heater.run();
    heater.run();
    System.out.println(label);
    System.out.printf("outlet_T_K=%.15g%n", heater.getOutletStream().getTemperature());
    System.out.printf("reported_specified_T_K=%.15g%n", heater.getSpecifiedOutletTemperature());
    System.out.println();
  }

  public static void main(String[] args) {
    System.out.println("=== one specification at a time (the capture's rows) ===");
    row("outlet_temperature_380",
        h -> h.setOutletTemperature(380.0));
    row("duty_5000_W",
        h -> h.setDuty(5000.0));

    System.out.println("=== two specifications, both orders ===");
    row("temperature_380_then_duty_5000",
        h -> {
          h.setOutletTemperature(380.0);
          h.setDuty(5000.0);
        });
    row("duty_5000_then_temperature_380",
        h -> {
          h.setDuty(5000.0);
          h.setOutletTemperature(380.0);
        });

    System.out.println("=== the energy-input spelling, both orders ===");
    row("temperature_380_then_energy_input_5000",
        h -> {
          h.setOutletTemperature(380.0);
          h.setEnergyInput(5000.0);
        });
    row("energy_input_5000_then_temperature_380",
        h -> {
          h.setEnergyInput(5000.0);
          h.setOutletTemperature(380.0);
        });

    System.out.println("=== the unit-carrying overload and the zero-Celsius setpoint ===");
    row("temperature_380K_by_unit_then_duty_5000",
        h -> {
          h.setOutletTemperature(380.0, "K");
          h.setDuty(5000.0);
        });
    row("temperature_0C_by_unit_then_duty_5000",
        h -> {
          h.setOutletTemperature(0.0, "C");
          h.setDuty(5000.0);
        });
    row("temperature_0C_by_unit_alone",
        h -> h.setOutletTemperature(0.0, "C"));
    row("temperature_273_15K_alone",
        h -> h.setOutletTemperature(273.15));

    System.out.println("=== the dropped specification after a repeated run ===");
    twice("temperature_380_then_duty_5000, run twice",
        h -> {
          h.setOutletTemperature(380.0);
          h.setDuty(5000.0);
        });
    twice("duty_5000_then_temperature_380, run twice",
        h -> {
          h.setDuty(5000.0);
          h.setOutletTemperature(380.0);
        });
  }
}
