// `PlusCharacterize.setHeavyTBPtoPlus` - the TBP row whose NAME carries the largest carbon
// number, flagged as the fluid's plus end.
//
// **The carbon number is read out of the component name**, not out of a stored property:
//
//   name.substring(3, 4).equals("_") ? name.substring(1, 3) : name.substring(1, 2)
//
// so `C7_PC` reads 7 and `C10_PC` reads 10. `addTBPfraction("C7", ...)` is what produces those
// names - it appends `_PC` - so the walk's answer depends on the spelling the caller used, and a
// row whose name does not parse is caught, logged and skipped rather than stopping the walk.
//
// `densLastTBP` is the density of the row **below the winner in index order**, not the row with
// the next carbon number: the class indexes `compNumber - 1` after picking `compNumber`, so a
// fluid whose names are not in ascending carbon-number order takes the neighbour by position.
//
// Each case prints the fluid as the walk sees it (names, molar masses, mole fractions, normal
// liquid densities) and then what the walk set on the characteriser.
//
//   javac -proc:none -cp /path/to/neqsim-f0c7436.jar SetHeavyTbpProbe.java
//   java -cp .:/path/to/neqsim-f0c7436.jar SetHeavyTbpProbe > captures/set_heavy_tbp_probe.tsv

import java.util.function.ToDoubleFunction;
import neqsim.thermo.characterization.PlusCharacterize;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class SetHeavyTbpProbe {

  /** One fluid: its TBP rows in the order they are added, and what the walk did with them. */
  private static void sweep(String label, String[] names, double[] moles, double[] molarMass, double[] specificGravity) {
    SystemInterface system = new SystemSrkEos(298.15, 1.0);
    for (int i = 0; i < names.length; i++) {
      system.addTBPfraction(names[i], moles[i], molarMass[i], specificGravity[i]);
    }
    system.setMixingRule(2);
    system.init(0);
    report(label, system);
  }

  /** The walk, and everything it set, over a fluid someone else built. */
  private static void report(String label, SystemInterface system) {
    System.out.printf("# fluid = %s%n", label);
    System.out.printf("component_count = %d%n", system.getPhase(0).getNumberOfComponents());
    System.out.printf("names = [%s]%n", names(system));
    System.out.printf("molar_mass = [%s]%n", molarMasses(system));
    System.out.printf("mole_fraction = [%s]%n", moleFractions(system));
    System.out.printf("normal_liquid_density = [%s]%n", normalLiquidDensities(system));

    PlusCharacterize characterise = new PlusCharacterize(system);
    try {
      characterise.setHeavyTBPtoPlus();
    } catch (Exception ex) {
      System.out.printf("# refused: %s%n", ex);
      return;
    }

    System.out.printf("plus_index = %d%n", plusIndex(system));
    System.out.printf("plus_name = %s%n", plusName(system));
    System.out.printf("m_plus = %.15g%n", characterise.getMPlus());
    System.out.printf("z_plus = %.15g%n", characterise.getZPlus());
    System.out.printf("dens_plus = %.15g%n", characterise.getDensPlus());
    System.out.printf("dens_last_tbp = %.15g%n", characterise.getDensLastTBP());
    System.out.printf("coefs_2 = %.15g%n", characterise.getCoef(2));
    System.out.printf("is_plus = [%s]%n", plusFlags(system));
  }

  /** One fluid with a light end that is not a TBP fraction ahead of its TBP rows. */
  private static void sweepWithLightEnd(
      String label, String[] names, double[] moles, double[] molarMass, double[] specificGravity) {
    SystemInterface system = new SystemSrkEos(298.15, 1.0);
    system.addComponent("methane", 1.0);
    for (int i = 0; i < names.length; i++) {
      system.addTBPfraction(names[i], moles[i], molarMass[i], specificGravity[i]);
    }
    system.setMixingRule(2);
    system.init(0);
    report(label, system);
  }

  private static int plusIndex(SystemInterface system) {
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      if (system.getPhase(0).getComponent(i).isIsPlusFraction()) {
        return i;
      }
    }
    return -1;
  }

  private static String plusName(SystemInterface system) {
    int index = plusIndex(system);
    return index < 0 ? "" : system.getPhase(0).getComponent(index).getComponentName();
  }

  private static String plusFlags(SystemInterface system) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(system.getPhase(0).getComponent(i).isIsPlusFraction() ? "true" : "false");
    }
    return out.toString();
  }

  private static String names(SystemInterface system) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(system.getPhase(0).getComponent(i).getComponentName());
    }
    return out.toString();
  }

  private static String molarMasses(SystemInterface system) {
    return vector(system, component -> component.getMolarMass());
  }

  private static String moleFractions(SystemInterface system) {
    return vector(system, component -> component.getz());
  }

  private static String normalLiquidDensities(SystemInterface system) {
    return vector(system, component -> component.getNormalLiquidDensity());
  }

  /** One component property, over every component in the phase, in index order. */
  private static String vector(SystemInterface system, ToDoubleFunction<ComponentInterface> property) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(String.format("%.15g", property.applyAsDouble(system.getPhase(0).getComponent(i))));
    }
    return out.toString();
  }

  public static void main(String[] args) {
    // The winner is the eleventh row by index but the largest carbon number by name, so a walk
    // that took the last row, or the row with the most moles, answers differently.
    sweep(
        "the_plus_end_is_not_the_last_row",
        new String[] {"C7", "C20", "C12", "C9"},
        new double[] {1.0, 1.0, 1.0, 1.0},
        new double[] {0.10, 0.35, 0.20, 0.14},
        new double[] {0.72, 0.90, 0.80, 0.75});

    // Two digits against one, which is the branch the substring test selects.
    sweep(
        "one_and_two_digit_names",
        new String[] {"C6", "C9", "C10", "C11"},
        new double[] {1.0, 1.0, 1.0, 1.0},
        new double[] {0.086, 0.13, 0.15, 0.17},
        new double[] {0.68, 0.75, 0.78, 0.80});

    // A name that does not parse: `AB_PC` reads `B` where the carbon number should be, the
    // `NumberFormatException` is caught, and the walk carries on with the rest.
    sweep(
        "a_name_that_does_not_parse",
        new String[] {"C7", "AB", "C13"},
        new double[] {1.0, 1.0, 1.0},
        new double[] {0.10, 0.20, 0.22},
        new double[] {0.72, 0.80, 0.84});

    // One row: the winner is index 0, and `component(compNumber - 1)` is `component(-1)`.
    sweep(
        "one_row",
        new String[] {"C9"},
        new double[] {1.0},
        new double[] {0.15},
        new double[] {0.78});

    // **The first fluid's four rows in another order.** The plus end is the same row by name and
    // a different position, and `densLastTBP` is a different row's density, so a walk that sorted
    // the rows or took `compNumber + 1` answers differently.
    sweep(
        "the_same_rows_in_another_order",
        new String[] {"C7", "C9", "C20", "C12"},
        new double[] {1.0, 1.0, 1.0, 1.0},
        new double[] {0.10, 0.14, 0.35, 0.20},
        new double[] {0.72, 0.75, 0.90, 0.80});

    // **A component that is not a TBP fraction is not a candidate**, even when its name carries a
    // carbon number: the walk's first test is `isIsTBPfraction`. `addTBPfraction` appends `_PC`
    // itself, so the names here are the fluid's own.
    sweepWithLightEnd(
        "a_regular_component_in_the_mix",
        new String[] {"C7", "C10"},
        new double[] {1.0, 1.0},
        new double[] {0.10, 0.15},
        new double[] {0.72, 0.78});
  }
}
