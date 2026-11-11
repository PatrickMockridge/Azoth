// `FoulingModel` driven directly, with `updateConditions` called on every row.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-f0c7436.jar FoulingModelProbe.java
//     java -cp .:neqsim-f0c7436.jar FoulingModelProbe > captures/fouling_model_probe.tsv
//
// **The Ebert-Panchal resistance is the class's analytical path** (`getFoulingResistance(time)`),
// not the stateful `advanceTime` Euler step, which is a separate row. The threshold temperature
// answers `Double.MAX_VALUE` when the target rate exceeds the deposition rate at zero fouling.

import java.util.Map;
import neqsim.process.mechanicaldesign.heatexchanger.FoulingModel;

public class FoulingModelProbe {

  // The stated condition set every row starts from.
  static final double VELOCITY = 1.5;
  static final double DENSITY = 800.0;
  static final double VISCOSITY = 5e-4;
  static final double WALL_TEMPERATURE = 473.15;
  static final double TUBE_ID = 0.01905;

  public static void main(String[] args) {
    // 1-3. The fixed model, the default.
    row("fixed_default", new FoulingModel(), 1000.0, 0.0001, 1e-10);
    FoulingModel fixedCustom = new FoulingModel();
    fixedCustom.setFixedFoulingResistance(0.0003);
    row("fixed_custom", fixedCustom, 1000.0, 0.0001, 1e-10);
    row("fixed_predict_beyond", fixedCustom, 1000.0, 0.001, 1e-10);

    // 4-7. Kern-Seaton cooling water: the time constant, zero time, a zero maximum and late.
    row("kern_seaton_cooling_water", FoulingModel.createCoolingWaterModel(0.0002, 720.0), 720.0,
        0.0001, 1e-10);
    row("kern_seaton_time_zero", FoulingModel.createCoolingWaterModel(0.0002, 720.0), 0.0, 0.0001,
        1e-10);
    FoulingModel kernZero = FoulingModel.createCoolingWaterModel(0.0002, 720.0);
    kernZero.setRfMax(0.0);
    row("kern_seaton_rfmax_zero", kernZero, 720.0, 0.0001, 1e-10);
    row("kern_seaton_late", FoulingModel.createCoolingWaterModel(0.0002, 720.0), 5000.0, 0.0001,
        1e-10);

    // 8-11. The crude standard across the stated flow and wall.
    row("ebert_crude", FoulingModel.createCrudeOilModel(), 1000.0, 0.0001, 1e-10);
    rowWithConditions("ebert_crude_slow_flow", FoulingModel.createCrudeOilModel(), 0.3, WALL_TEMPERATURE,
        1000.0, 0.0001, 1e-10);
    rowWithConditions("ebert_crude_hot_wall", FoulingModel.createCrudeOilModel(), VELOCITY, 523.15,
        1000.0, 0.0001, 1e-10);
    // A removal constant high enough that the asymptotic crosses the fixed resistance inside
    // the binary search's bracket, so the threshold velocity is an interior answer.
    FoulingModel suppressed = FoulingModel.createCrudeOilModel();
    suppressed.setGamma(5e-8);
    rowWithConditions("ebert_suppressed_threshold_velocity", suppressed, 0.3, WALL_TEMPERATURE,
        1000.0, 0.0001, 1e-10);

    // 12-15. No removal, the threshold sentinel and an interior threshold, a target past the
    // asymptote, and the heavy crude factory.
    FoulingModel noRemoval = FoulingModel.createCrudeOilModel();
    noRemoval.setGamma(0.0);
    row("ebert_no_removal", noRemoval, 1000.0, 0.0001, 1e-10);
    row("ebert_threshold_sentinel", FoulingModel.createCrudeOilModel(), 1000.0, 0.0001, 0.01);
    row("ebert_threshold_mid", FoulingModel.createCrudeOilModel(), 1000.0, 0.0001, 1e-6);
    FoulingModel crude = FoulingModel.createCrudeOilModel();
    crude.updateConditions(VELOCITY, DENSITY, VISCOSITY, WALL_TEMPERATURE, TUBE_ID);
    row("ebert_predict_beyond_asymptote", crude, 1000.0,
        crude.getAsymptoticFoulingResistance() + 1e-6, 1e-10);
    row("heavy_crude", FoulingModel.createHeavyCrudeModel(), 1000.0, 0.0001, 1e-10);

    // 16. The stateful Euler path: three one-hour steps, then a reset.
    FoulingModel stepped = FoulingModel.createCrudeOilModel();
    stepped.updateConditions(VELOCITY, DENSITY, VISCOSITY, WALL_TEMPERATURE, TUBE_ID);
    for (int i = 0; i < 3; i++) {
      stepped.advanceTime(1.0);
    }
    row("ebert_advance_three_steps", stepped, 1000.0, 0.0001, 1e-10);

    System.out.flush();
  }

  static void row(String label, FoulingModel model, double timeHours, double targetRf,
      double thresholdRate) {
    rowWithConditions(label, model, VELOCITY, WALL_TEMPERATURE, timeHours, targetRf, thresholdRate);
  }

  static void rowWithConditions(String label, FoulingModel model, double velocity,
      double wallTemperature, double timeHours, double targetRf, double thresholdRate) {
    model.updateConditions(velocity, DENSITY, VISCOSITY, wallTemperature, TUBE_ID);
    System.out.println(label);
    System.out.println("model=" + model.getModelType().name());
    System.out.println("velocity_m_per_s=" + velocity);
    System.out.println("fluid_density_kg_per_m3=" + DENSITY);
    System.out.println("fluid_viscosity_pa_s=" + VISCOSITY);
    System.out.println("wall_temperature_k=" + wallTemperature);
    System.out.println("tube_inner_diameter_m=" + TUBE_ID);
    System.out.println("time_hours=" + timeHours);
    System.out.println("target_fouling_resistance=" + targetRf);
    System.out.println("threshold_rate=" + thresholdRate);
    System.out.println("fouling_resistance_at_time=" + model.getFoulingResistance(timeHours));
    System.out.println("current_fouling_resistance=" + model.getFoulingResistance());
    System.out.println("asymptotic_fouling_resistance=" + model.getAsymptoticFoulingResistance());
    System.out.println(
        "threshold_wall_temperature_k=" + model.calcThresholdTemperature(thresholdRate));
    System.out.println("threshold_velocity_m_per_s=" + model.calcThresholdVelocity());
    System.out.println("time_to_target_hours=" + model.predictTimeToFouling(targetRf));
    System.out.println("operating_time_hours=" + model.getOperatingTimeHours());
    System.out.println(
        "rate_at_current_rf=" + model.calcEbertPanchalFoulingRate(model.getFoulingResistance()));
    Map<String, Object> map = model.toMap();
    if (map.containsKey("ebertPanchalParameters")) {
      @SuppressWarnings("unchecked")
      Map<String, Object> params = (Map<String, Object>) map.get("ebertPanchalParameters");
      System.out.println("reynolds_number=" + params.get("reynoldsNumber"));
      System.out.println("wall_shear_stress_pa=" + params.get("wallShearStress_Pa"));
    }
    System.out.println();
  }
}
