// CharacterizationFluidRebuildProbe - the fluid-rebuild half of PseudoComponentCombiner at the pin.
//
// `combineReservoirFluids` pools N fluids' pseudo components, places equal-mass cut points on the
// pooled mass, distributes each fluid's own list onto them, and merges the profiles; the Tb, Tc,
// Pc and omega of each cut are weighted by a fluid's mass fraction times its mole fraction within
// that fluid, while the rest are mass-weighted. `characterizeToCommonSlate` does the same onto a
// slate inferred from the fluids or stated, with optional per-fluid weights, and shares the pseudo
// properties across the returned clones while keeping each fluid's own pseudo moles.
// `transferBinaryInteractionParameters` copies a reference's nonzero pairs onto a target, bases by
// name and pseudo components by position.
//
// Every source below is a real TBP fraction (`addTBPfraction` with a specific gravity), so every
// property profile.applyTo sets is finite - the TBP model's own defaults cannot leak into a row.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar CharacterizationFluidRebuildProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar CharacterizationFluidRebuildProbe > captures/characterization_fluid_rebuild_probe.tsv

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import neqsim.thermo.characterization.PseudoComponentCombiner;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseEos;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class CharacterizationFluidRebuildProbe {

  private static final String COLUMNS = "fluid|name|pseudo|moles|molar_mass|density|"
      + "normal_boiling_point|critical_temperature|critical_pressure|acentric_factor|"
      + "critical_volume|racket_z|racket_z_cpa|parachor|critical_viscosity|"
      + "triple_point_temperature|heat_of_fusion|ideal_gas_enthalpy_of_formation|"
      + "cp_a|cp_b|cp_c|cp_d|attractive_m";

  private static final String UNITS = "index|string|boolean|mol|kg/mol|kg/m3|K|K|Pa|"
      + "1|m3/mol|1|1|1|Pa*s|K|J/mol|J/mol|J/(mol*K)|J/(mol*K2)|J/(mol*K3)|J/(mol*K4)|1";

  /** One fluid: base components by name, then TBP fractions as {name, moles, MW, SG}. */
  private static SystemInterface fluid(String[] baseNames, double[] baseMoles,
      String[] pcNames, double[] pcMoles, double[] pcMass, double[] pcSg) {
    SystemInterface system = new SystemSrkEos(298.15, 1.0);
    for (int i = 0; i < baseNames.length; i++) {
      system.addComponent(baseNames[i], baseMoles[i]);
    }
    for (int i = 0; i < pcNames.length; i++) {
      system.addTBPfraction(pcNames[i], pcMoles[i], pcMass[i], pcSg[i]);
    }
    system.setMixingRule(2);
    system.init(0);
    return system;
  }

  private static SystemInterface fluidA() {
    return fluid(new String[] {"methane"}, new double[] {10.0},
        new String[] {"C7", "C8", "C9"}, new double[] {3.0, 2.0, 1.5},
        new double[] {0.095, 0.107, 0.121}, new double[] {0.73, 0.76, 0.79});
  }

  private static SystemInterface fluidB() {
    return fluid(new String[] {"methane", "CO2"}, new double[] {4.0, 1.5},
        new String[] {"C10", "C12", "C15"}, new double[] {4.0, 3.0, 2.0},
        new double[] {0.134, 0.161, 0.206}, new double[] {0.81, 0.84, 0.86});
  }

  private static SystemInterface fluidC() {
    return fluid(new String[] {"ethane"}, new double[] {5.0},
        new String[] {"C11", "C14"}, new double[] {2.0, 1.5},
        new double[] {0.147, 0.190}, new double[] {0.82, 0.85});
  }

  private static SystemInterface baseOnlyD() {
    return fluid(new String[] {"methane", "ethane"}, new double[] {5.0, 3.0},
        new String[] {}, new double[] {}, new double[] {}, new double[] {});
  }

  private static SystemInterface baseOnlyE() {
    return fluid(new String[] {"methane", "CO2", "propane"}, new double[] {2.0, 1.0, 0.5},
        new String[] {}, new double[] {}, new double[] {}, new double[] {});
  }

  private static String number(double value) {
    return Double.isFinite(value) ? String.format(Locale.ROOT, "%.17g", value) : "-";
  }

  private static boolean pseudo(ComponentInterface component) {
    return component.isIsTBPfraction() || component.isIsPlusFraction();
  }

  /** One component as a row; pseudo rows carry every property, base rows carry none. */
  private static String row(int fluid, ComponentInterface component, boolean withFluid) {
    StringBuilder text = new StringBuilder();
    if (withFluid) {
      text.append(fluid).append('|');
    }
    boolean isPseudo = pseudo(component);
    text.append(component.getComponentName()).append('|').append(isPseudo).append('|')
        .append(number(component.getNumberOfmoles()));
    if (!isPseudo) {
      for (int i = 0; i < 19; i++) {
        text.append("|-");
      }
      return text.toString();
    }
    text.append('|').append(number(component.getMolarMass()));
    text.append('|').append(number(component.getNormalLiquidDensity() * 1000.0));
    text.append('|').append(number(component.getNormalBoilingPoint()));
    text.append('|').append(number(component.getTC()));
    text.append('|').append(number(component.getPC() * 1.0e5));
    text.append('|').append(number(component.getAcentricFactor()));
    text.append('|').append(number(component.getCriticalVolume() / 1.0e6));
    text.append('|').append(number(component.getRacketZ()));
    text.append('|').append(number(component.getRacketZCPA()));
    text.append('|').append(number(component.getParachorParameter()));
    text.append('|').append(number(component.getCriticalViscosity()));
    text.append('|').append(number(component.getTriplePointTemperature()));
    text.append('|').append(number(component.getHeatOfFusion()));
    text.append('|').append(number(component.getIdealGasEnthalpyOfFormation()));
    text.append('|').append(number(component.getCpA()));
    text.append('|').append(number(component.getCpB()));
    text.append('|').append(number(component.getCpC()));
    text.append('|').append(number(component.getCpD()));
    text.append('|').append(
        component.getAttractiveTerm() == null ? "-" : number(component.getAttractiveTerm().getm()));
    return text.toString();
  }

  private static void rows(String prefix, SystemInterface system, int fluidIndex,
      boolean withFluid) {
    for (int i = 0; i < system.getNumberOfComponents(); i++) {
      ComponentInterface component = system.getComponent(i);
      if (!(component.getNumberOfmoles() > 0.0)) {
        continue;
      }
      System.out.println(prefix + i + " = " + row(fluidIndex, component, withFluid));
    }
  }

  private static double totalMoles(SystemInterface system) {
    double total = 0.0;
    for (int i = 0; i < system.getNumberOfComponents(); i++) {
      total += system.getComponent(i).getNumberOfmoles();
    }
    return total;
  }

  private static void combine(String label, int target, SystemInterface... fluids) {
    System.out.println("case=" + label + " kind=combine target=" + target);
    List<SystemInterface> list = new ArrayList<>();
    for (int j = 0; j < fluids.length; j++) {
      list.add(fluids[j]);
      System.out.println("fluid." + j + " total_moles = " + number(totalMoles(fluids[j])));
      rows("input." + j + ".", fluids[j], j, true);
    }
    SystemInterface combined = PseudoComponentCombiner.combineReservoirFluids(target, list);
    rows("output.", combined, 0, false);
    System.out.println("total_moles = " + number(totalMoles(combined)));
    System.out.println();
  }

  private static void slate(String label, double[] weights, int target,
      SystemInterface... fluids) {
    System.out.println("case=" + label + " kind=common_slate weights="
        + (weights == null ? "equal" : java.util.Arrays.toString(weights))
        + " target=" + (target <= 0 ? "inferred" : Integer.toString(target)));
    List<SystemInterface> list = new ArrayList<>();
    for (int j = 0; j < fluids.length; j++) {
      list.add(fluids[j]);
      System.out.println("fluid." + j + " total_moles = " + number(totalMoles(fluids[j])));
      rows("input." + j + ".", fluids[j], j, true);
    }
    List<SystemInterface> out = target <= 0
        ? PseudoComponentCombiner.characterizeToCommonSlate(list, weights)
        : PseudoComponentCombiner.characterizeToCommonSlate(list, weights, target);
    for (int j = 0; j < out.size(); j++) {
      rows("output." + j + ".", out.get(j), j, true);
      System.out.println("fluid_moles." + j + " = " + number(totalMoles(out.get(j))));
    }
    System.out.println();
  }

  private static double kij(SystemInterface system, int i, int j) {
    if (system.getPhase(0) instanceof PhaseEos phase) {
      return phase.getMixingRule().getBinaryInteractionParameter(i, j);
    }
    return 0.0;
  }

  private static void matrix(String prefix, SystemInterface system) {
    int n = system.getNumberOfComponents();
    StringBuilder text = new StringBuilder(prefix).append(" = [");
    for (int i = 0; i < n; i++) {
      if (i > 0) {
        text.append(',');
      }
      text.append('[');
      for (int j = 0; j < n; j++) {
        if (j > 0) {
          text.append(',');
        }
        text.append(number(kij(system, i, j)));
      }
      text.append(']');
    }
    text.append(']');
    System.out.println(text);
  }

  private static void bip(String label, SystemInterface reference, SystemInterface target,
      boolean refuse) {
    System.out.println("case=" + label + " kind=bip");
    for (int i = 0; i < reference.getNumberOfComponents(); i++) {
      ComponentInterface component = reference.getComponent(i);
      System.out.println("reference." + i + " = " + component.getComponentName() + "|"
          + pseudo(component));
    }
    for (int i = 0; i < target.getNumberOfComponents(); i++) {
      ComponentInterface component = target.getComponent(i);
      System.out.println("target." + i + " = " + component.getComponentName() + "|"
          + pseudo(component));
    }
    matrix("reference_kij", reference);
    matrix("target_initial_kij", target);
    if (refuse) {
      try {
        PseudoComponentCombiner.transferBinaryInteractionParameters(reference, target);
        System.out.println("status=answered");
      } catch (RuntimeException failure) {
        System.out.println("status=refused reason=" + failure.getClass().getSimpleName() + ":"
            + failure.getMessage());
      }
      System.out.println();
      return;
    }
    PseudoComponentCombiner.transferBinaryInteractionParameters(reference, target);
    matrix("target_final_kij", target);
    System.out.println();
  }

  public static void main(String[] args) {
    System.out.println("characterization fluid rebuild at f0c7436c69");
    System.out.println("# columns = " + COLUMNS);
    System.out.println("# units = " + UNITS);
    System.out.println();

    combine("combine_two_fluids_target_4", 4, fluidA(), fluidB());
    combine("combine_three_fluids_target_2", 2, fluidA(), fluidB(), fluidC());
    combine("combine_target_beyond_pool", 9, fluidA(), fluidB());
    combine("combine_no_pseudo_components", 3, baseOnlyD(), baseOnlyE());
    combine("combine_one_fluid_has_no_pseudo", 3, fluidA(), baseOnlyD());

    slate("common_slate_two_fluids_equal_inferred", null, 0, fluidA(), fluidB());
    slate("common_slate_two_fluids_weighted_target_3", new double[] {2.0, 1.0}, 3,
        fluidA(), fluidB());
    slate("common_slate_three_fluids_weighted_inferred", new double[] {1.0, 0.5, 2.0}, 0,
        fluidA(), fluidB(), fluidC());
    slate("common_slate_zero_weight_fluid", new double[] {1.0, 0.0}, 3, fluidA(), fluidB());
    slate("common_slate_no_pseudo_components", null, 0, baseOnlyD(), baseOnlyE());

    SystemInterface reference = fluid(new String[] {"methane", "CO2"}, new double[] {10.0, 2.0},
        new String[] {"C7", "C8", "C9"}, new double[] {3.0, 2.0, 1.5},
        new double[] {0.095, 0.107, 0.121}, new double[] {0.73, 0.76, 0.79});
    reference.setBinaryInteractionParameter("methane", "CO2", 0.09);
    reference.setBinaryInteractionParameter("methane", "C8_PC", 0.04);
    reference.setBinaryInteractionParameter("CO2", "C9_PC", 0.02);
    reference.setBinaryInteractionParameter("C7_PC", "C8_PC", 0.01);
    SystemInterface target = fluid(new String[] {"methane", "CO2"}, new double[] {5.0, 1.0},
        new String[] {"C10", "C11"}, new double[] {2.0, 1.0},
        new double[] {0.134, 0.147}, new double[] {0.81, 0.82});
    bip("bip_base_and_pseudo_pairs", reference, target, false);

    SystemInterface fewTarget = fluid(new String[] {"methane", "CO2"}, new double[] {5.0, 1.0},
        new String[] {"C10"}, new double[] {2.0},
        new double[] {0.134}, new double[] {0.81});
    bip("bip_target_has_fewer_pcs", reference, fewTarget, false);

    SystemInterface missing = fluid(new String[] {"methane"}, new double[] {5.0},
        new String[] {"C7"}, new double[] {2.0},
        new double[] {0.095}, new double[] {0.73});
    bip("bip_missing_base_name", reference, missing, true);
  }
}
