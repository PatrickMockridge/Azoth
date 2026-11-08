// `FoulingModel` driven directly, on the class's own factories and conditions.
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
import neqsim.process.mechanicaldesign.heatexchanger.FoulingModel.FoulingModelType;

public class FoulingModelProbe {

  public static void main(String[] args) {
    // 1. The fixed model, the default.
    FoulingModel fixed = new FoulingModel();
    row("fixed_default", fixed, 1000.0, 0.0001, 1e-10);

    // 2. A custom fixed resistance.
    FoulingModel fixedCustom = new FoulingModel();
    fixedCustom.setFixedFoulingResistance(0.0003);
    row("fixed_custom", fixedCustom, 1000.0, 0.0001, 1e-10);

    // 3-5. Kern-Seaton cooling water: at the time constant, at zero time, and at a zero maximum.
    FoulingModel kern = FoulingModel.createCoolingWaterModel(0.0002, 720.0);
    row("kern_seaton_cooling_water", kern, 720.0, 0.0001, 1e-10);
    row("kern_seaton_time_zero", kern, 0.0, 0.0001, 1e-10);
    FoulingModel kernZero = FoulingModel.createCoolingWaterModel(0.0002, 720.0);
    kernZero.setRfMax(0.0);
    row("kern_seaton_rfmax_zero", kernZero, 720.0, 0.0001, 1e-10);
    row("kern_seaton_late", FoulingModel.createCoolingWaterModel(0.0002, 720.0), 5000.0, 0.0001, 1e-10);

    // 6-8. The crude standard: the factory's own conditions, then a stated flow and wall.
    row("ebert_crude_default_conditions", FoulingModel.createCrudeOilModel(), 1000.0, 0.0001, 1e-10);
    FoulingModel crude = FoulingModel.createCrudeOilModel();
    crude.updateConditions(1.5, 800.0, 5e-4, 473.15, 0.01905);
    row("ebert_crude_updated_conditions", crude, 1000.0, 0.0001, 1e-10);
    FoulingModel slow = FoulingModel.createCrudeOilModel();
    slow.updateConditions(0.3, 800.0, 5e-4, 473.15, 0.01905);
    row("ebert_crude_slow_flow", slow, 1000.0, 0.0001, 1e-10);
    FoulingModel hot = FoulingModel.createCrudeOilModel();
    hot.updateConditions(1.5, 800.0, 5e-4, 523.15, 0.01905);
    row("ebert_crude_hot_wall", hot, 1000.0, 0.0001, 1e-10);
    // A removal constant high enough that the asymptotic crosses the fixed resistance inside
    // the binary search's bracket, so the threshold velocity is an interior answer.
    FoulingModel suppressed = FoulingModel.createCrudeOilModel();
    suppressed.setGamma(5e-8);
    suppressed.updateConditions(0.3, 800.0, 5e-4, 473.15, 0.01905);
    row("ebert_suppressed_threshold_velocity", suppressed, 1000.0, 0.0001, 1e-10);

    // 9-10. No removal, and a threshold rate above the deposition rate.
    FoulingModel noRemoval = FoulingModel.createCrudeOilModel();
    noRemoval.setGamma(0.0);
    row("ebert_no_removal", noRemoval, 1000.0, 0.0001, 1e-10);
    FoulingModel sentinel = FoulingModel.createCrudeOilModel();
    row("ebert_threshold_sentinel", sentinel, 1000.0, 0.0001, 0.01);
    row("ebert_threshold_mid", sentinel, 1000.0, 0.0001, 1e-6);
    row("ebert_predict_beyond_asymptote", crude, 1000.0,
        crude.getAsymptoticFoulingResistance() + 1e-6, 1e-10);

    // 11. The stateful Euler path: three one-hour steps.
    FoulingModel stepped = FoulingModel.createCrudeOilModel();
    stepped.updateConditions(1.5, 800.0, 5e-4, 473.15, 0.01905);
    for (int i = 0; i < 3; i++) {
      stepped.advanceTime(1.0);
    }
    row("ebert_advance_three_steps", stepped, 1000.0, 0.0001, 1e-10);
    stepped.reset();
    row("ebert_after_reset", stepped, 1000.0, 0.0001, 1e-10);

    // 12. The heavy crude factory.
    FoulingModel heavy = FoulingModel.createHeavyCrudeModel();
    heavy.updateConditions(1.5, 800.0, 5e-4, 473.15, 0.01905);
    row("heavy_crude_updated_conditions", heavy, 1000.0, 0.0001, 1e-10);

    System.out.flush();
  }

  static void row(String label, FoulingModel model, double timeHours, double targetRf,
      double thresholdRate) {
    System.out.println(label);
    System.out.println("model=" + model.getModelType().name());
    System.out.println("time_hours=" + timeHours);
    System.out.println("target_fouling_resistance=" + targetRf);
    System.out.println("threshold_rate=" + thresholdRate);
    System.out.println("fouling_resistance_at_time=" + model.getFoulingResistance(timeHours));
    System.out.println("current_fouling_resistance=" + model.getFoulingResistance());
    System.out.println("asymptotic_fouling_resistance=" + model.getAsymptoticFoulingResistance());
    System.out.println("threshold_wall_temperature_k=" + model.calcThresholdTemperature(thresholdRate));
    System.out.println("threshold_velocity_m_per_s=" + model.calcThresholdVelocity());
    System.out.println("time_to_target_hours=" + model.predictTimeToFouling(targetRf));
    System.out.println("operating_time_hours=" + model.getOperatingTimeHours());
    System.out.println("fouling_rate=" + model.calcEbertPanchalFoulingRate(0.0));
    Map<String, Object> map = model.toMap();
    if (map.containsKey("ebertPanchalParameters")) {
      @SuppressWarnings("unchecked")
      Map<String, Object> params = (Map<String, Object>) map.get("ebertPanchalParameters");
      System.out.println("reynolds_number=" + params.get("reynoldsNumber"));
      System.out.println("wall_shear_stress_pa=" + params.get("wallShearStress_Pa"));
      System.out.println("rate_at_current_rf="
          + model.calcEbertPanchalFoulingRate(model.getFoulingResistance()));
    }
    System.out.println();
  }
}
