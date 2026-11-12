// `ShellAndTubeDesignCalculator` driven directly, on its defaults and every settable surface.
//
// Compile and run from this directory:
//
//     javac -proc:none -cp neqsim-118d9af.jar ShellAndTubeDesignProbe.java
//     java -cp .:neqsim-118d9af.jar ShellAndTubeDesignProbe \
//         > captures/shell_and_tube_design_probe.tsv
//
// The calculator is a single-pass pipeline: TEMA configuration, material lookup from the H2
// design database, tube bundle, shell dimensions, baffles, UHX-13, UG-37, UG-27, UG-99, the
// optional NACE assessment, the optional thermal and vibration block, weights and costs. The
// probe prints every setter it used and every `toMap()` leaf, so a row is reconstructible.
//
// **A failed material lookup is visible**: the class leaves `shell.material` at its
// "Carbon Steel ..." default. A row whose material string is a default rather than a
// database specification is marked invalid in the captures README.

import java.util.List;
import java.util.Map;

import neqsim.process.mechanicaldesign.heatexchanger.ShellAndTubeDesignCalculator;
import neqsim.process.mechanicaldesign.heatexchanger.TEMAStandard.StandardTubeSize;
import neqsim.process.mechanicaldesign.heatexchanger.ThermalDesignCalculator;

public class ShellAndTubeDesignProbe {

  public static void main(String[] args) {
    // 1. The class's own defaults.
    row("defaults", calculator(null, null, null, null, null, null));

    // 2. The rear-head clearance branches and the classes.
    for (String designation : new String[] {"AES", "BEM", "AEU", "AET", "AKT", "CEW"}) {
      ShellAndTubeDesignCalculator calc = calculator(designation, null, null, null, null, null);
      row("tema." + designation, calc);
    }

    // 3. The standard-size table's ends and the pass rounding.
    for (double area : new double[] {1.0, 50.0, 100.0, 150.0, 500.0, 5000.0}) {
      row("area." + area, calculator(null, area, null, null, null, null));
    }
    for (int passes : new int[] {1, 2, 4}) {
      row("passes." + passes, calculator(null, 150.0, null, null, passes, null));
    }

    // 4. Every standard tube size.
    for (StandardTubeSize size : StandardTubeSize.values()) {
      row("tube_size." + size.name(), calculator(null, null, null, size, null, null));
    }

    // 5. Pressures and the allowable-stress derating branches.
    for (double pressure : new double[] {1.0, 30.0, 100.0}) {
      row("pressure." + pressure,
          calculator(null, null, pressure, null, null, null));
    }
    for (double temperature : new double[] {199.0, 200.0, 201.0, 299.0, 300.0, 301.0, 399.0,
        400.0, 401.0}) {
      row("temperature." + temperature,
          calculator(null, null, null, null, null, temperature));
    }

    // 6. The material table, including a grade that is not in it.
    materialRow("carbon_carbon", "SA-516-70", "SA-179");
    materialRow("carbon_ss", "SA-516-70", "SA-213-TP316");
    materialRow("ss_ss", "SA-240-Tp316", "SA-213-TP316");
    materialRow("duplex_duplex", "SA-240-S31803", "SA-789-S31803");
    materialRow("unknown_grade", "SA-999", "SA-999");

    // 7. The NACE assessment across its threshold and both material flags.
    naceRow("nace_below_threshold_carbon", 0.001, "SA-516-70", "SA-179");
    naceRow("nace_at_threshold_carbon", 0.003, "SA-516-70", "SA-179");
    naceRow("nace_above_threshold_carbon", 0.0031, "SA-516-70", "SA-179");
    naceRow("nace_above_threshold_ss", 0.05, "SA-240-Tp316", "SA-213-TP316");
    naceRow("nace_off_above_threshold", -1.0, "SA-516-70", "SA-179");

    // 8. Joint efficiency and nozzle reinforcement.
    for (double efficiency : new double[] {0.65, 0.85, 1.0}) {
      ShellAndTubeDesignCalculator calc = calculator(null, null, null, null, null, null);
      calc.setShellJointEfficiency(efficiency);
      row("joint_efficiency." + efficiency, calc);
    }
    for (double diameter : new double[] {50.0, 150.0, 300.0}) {
      ShellAndTubeDesignCalculator calc = calculator(null, null, null, null, null, null);
      calc.setNozzleDiameter(diameter);
      calc.setNozzleWallThickness(diameter <= 50.0 ? 5.54 : (diameter <= 150.0 ? 7.11 : 12.7));
      row("nozzle." + diameter, calc);
    }

    // 9. The optional thermal and vibration block, on both shell-side methods and one side only.
    thermalRow("thermal_kern", ThermalDesignCalculator.ShellSideMethod.KERN, true, true);
    thermalRow("thermal_bell_delaware", ThermalDesignCalculator.ShellSideMethod.BELL_DELAWARE,
        true, true);
    thermalRow("thermal_shell_only", ThermalDesignCalculator.ShellSideMethod.KERN, false, true);
    thermalRow("thermal_tube_only", ThermalDesignCalculator.ShellSideMethod.KERN, true, false);

    // 10. Two calls on one instance: the standards list.
    ShellAndTubeDesignCalculator twice =
        calculator(null, null, null, null, null, null);
    twice.calculate();
    twice.calculate();
    System.out.println("calculate_twice");
    printMap("result.", twice.toMap());
    System.out.println();

    // 11. The refusals.
    refusalRow("designation_short", calc -> calc.setTemaDesignation("AE"));
    refusalRow("designation_unknown_shell", calc -> calc.setTemaDesignation("AZE"));
    refusalRow("zero_passes", calc -> calc.setTubePasses(0));
    refusalRow("zero_area", calc -> calc.setRequiredArea(0.0));

    System.out.flush();
  }

  interface Mutation {
    void apply(ShellAndTubeDesignCalculator calc);
  }

  static ShellAndTubeDesignCalculator calculator(String designation, Double area, Double pressure,
      StandardTubeSize size, Integer passes, Double temperature) {
    ShellAndTubeDesignCalculator calc = new ShellAndTubeDesignCalculator();
    if (designation != null) {
      calc.setTemaDesignation(designation);
    }
    if (area != null) {
      calc.setRequiredArea(area);
    }
    if (pressure != null) {
      calc.setShellSidePressure(pressure);
      calc.setTubeSidePressure(pressure);
    }
    if (size != null) {
      calc.setTubeSize(size);
    }
    if (passes != null) {
      calc.setTubePasses(passes);
    }
    if (temperature != null) {
      calc.setDesignTemperature(temperature);
    }
    return calc;
  }

  static void row(String label, ShellAndTubeDesignCalculator calc) {
    System.out.println(label);
    inputs(calc);
    calc.calculate();
    printMap("result.", calc.toMap());
    System.out.println();
  }

  static void materialRow(String label, String shellGrade, String tubeGrade) {
    ShellAndTubeDesignCalculator calc = calculator(null, null, null, null, null, null);
    calc.setShellMaterialGrade(shellGrade);
    calc.setTubeMaterialGrade(tubeGrade);
    System.out.println(label);
    inputs(calc);
    try {
      calc.calculate();
      printMap("result.", calc.toMap());
    } catch (RuntimeException refused) {
      System.out.println("refused=" + refused.getMessage());
    }
    System.out.println();
  }

  static void naceRow(String label, double h2s, String shellGrade, String tubeGrade) {
    ShellAndTubeDesignCalculator calc = calculator(null, null, null, null, null, null);
    calc.setShellMaterialGrade(shellGrade);
    calc.setTubeMaterialGrade(tubeGrade);
    if (h2s >= 0.0) {
      calc.setSourServiceAssessment(true);
      calc.setH2sPartialPressure(h2s);
    }
    System.out.println(label);
    inputs(calc);
    calc.calculate();
    printMap("result.", calc.toMap());
    System.out.println();
  }

  static void thermalRow(String label, ThermalDesignCalculator.ShellSideMethod method,
      boolean tube, boolean shell) {
    ShellAndTubeDesignCalculator calc = calculator(null, null, null, null, null, null);
    if (tube) {
      calc.setTubeSideFluidProperties(998.0, 0.001, 4180.0, 0.60, 5.0, true);
    }
    if (shell) {
      calc.setShellSideFluidProperties(50.0, 1.5e-5, 2200.0, 0.03, 3.0);
    }
    calc.setFoulingTube(0.00018);
    calc.setFoulingShell(0.00035);
    calc.setShellSideMethod(method);
    System.out.println(label);
    inputs(calc);
    calc.calculate();
    printMap("result.", calc.toMap());
    System.out.println();
  }

  static void refusalRow(String label, Mutation mutation) {
    ShellAndTubeDesignCalculator calc = new ShellAndTubeDesignCalculator();
    String thrown;
    try {
      mutation.apply(calc);
      calc.calculate();
      thrown = "none";
    } catch (RuntimeException e) {
      thrown = e.getClass().getSimpleName() + ": " + e.getMessage();
    }
    System.out.println("refusals." + label + "=" + thrown);
  }

  static void inputs(ShellAndTubeDesignCalculator calc) {
    System.out.println("input.tema_designation=" + calc.getTemaDesignation());
    System.out.println("input.required_area_m2=" + calc.getRequiredArea());
    System.out.println("input.shell_material_grade=" + calc.getShellMaterialGrade());
    System.out.println("input.tube_material_grade=" + calc.getTubeMaterialGrade());
    System.out.println("input.tube_thermal_conductivity_WpmK="
        + calc.getTubeThermalConductivity());
    System.out.println("input.tube_length_mm=6096.0");
    System.out.println("input.tube_passes=2");
    System.out.println("input.tube_size=TUBE_3_4_INCH");
    System.out.println("input.tube_wall_thickness_mm=2.108");
    System.out.println("input.tube_pitch_ratio=1.25");
    System.out.println("input.pitch_pattern=TRIANGULAR_30");
    System.out.println("input.baffle_type=SINGLE_SEGMENTAL");
    System.out.println("input.baffle_cut=0.25");
    System.out.println("input.corrosion_allowance_shell_mm=3.175");
    System.out.println("input.corrosion_allowance_tube_mm=0.0");
    System.out.println("input.shell_side_pressure_bara=10.0");
    System.out.println("input.tube_side_pressure_bara=10.0");
    System.out.println("input.design_temperature_C=150.0");
    System.out.println("input.shell_joint_efficiency=" + calc.getShellJointEfficiency());
    System.out.println("input.nozzle_diameter_mm=150.0");
    System.out.println("input.nozzle_wall_thickness_mm=7.11");
    System.out.println("input.sour_service_assessment=" + calc.isSourServiceAssessment());
    System.out.println("input.h2s_partial_pressure_bar=" + calc.getH2sPartialPressure());
    System.out.println("input.shell_max_hardness_HV=248");
    System.out.println("input.tube_max_hardness_HV=248");
  }

  @SuppressWarnings("unchecked")
  static void printMap(String prefix, Map<String, Object> map) {
    for (Map.Entry<String, Object> entry : map.entrySet()) {
      if (entry.getValue() instanceof Map) {
        printMap(prefix + entry.getKey() + ".", (Map<String, Object>) entry.getValue());
      } else if (entry.getValue() instanceof List) {
        List<Object> rows = (List<Object>) entry.getValue();
        for (int i = 0; i < rows.size(); i++) {
          if (rows.get(i) instanceof Map) {
            printMap(prefix + entry.getKey() + "." + i + ".",
                (Map<String, Object>) rows.get(i));
          } else {
            System.out.println(prefix + entry.getKey() + "." + i + "=" + rows.get(i));
          }
        }
      } else {
        System.out.println(prefix + entry.getKey() + "=" + entry.getValue());
      }
    }
  }
}
