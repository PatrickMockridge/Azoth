// `PseudoComponentCombiner.delumpContributions` and `determineReferenceEqualMassBoundaries` -
// the two private statics, read out of the class by reflection.
//
// Both are pure: they take a list of the class's own `PseudoComponentContribution` row and an
// int, and touch no `SystemInterface`. The class exposes neither, and its only caller selects
// them from `CharacterizationOptions` inside `characterizeToReferenceCore` - so a run through
// the public API answers with the redistributed fluid and never with the delumped rows. This
// is the probe that reaches them, and it builds the row the same way `extractComponents` does:
// every component that is a TBP fraction or a plus fraction and carries a positive mole amount.
//
//   javac -proc:none -cp /path/to/neqsim-f0c7436.jar PseudoComponentDelumpingProbe.java
//   java -cp .:/path/to/neqsim-f0c7436.jar PseudoComponentDelumpingProbe > captures/pseudo_component_delumping_probe.tsv

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class PseudoComponentDelumpingProbe {

  private static final String COMBINER = "neqsim.thermo.characterization.PseudoComponentCombiner";
  private static final String CONTRIBUTION = COMBINER + "$PseudoComponentContribution";

  /** A fluid of TBP fractions, which is what the combiner reads as pseudo components. */
  private static SystemInterface fluid(double[][] cuts) {
    SystemInterface system = new SystemSrkEos(298.15, 1.0);
    for (int i = 0; i < cuts.length; i++) {
      system.addTBPfraction("C" + (7 + i), cuts[i][0], cuts[i][1], cuts[i][2]);
    }
    system.setMixingRule(2);
    system.init(0);
    return system;
  }

  /** The combiner's own row, built through the private constructor `extractComponents` uses. */
  private static List<Object> contributions(SystemInterface system) throws Exception {
    Class<?> row = Class.forName(CONTRIBUTION);
    Constructor<?> ctor =
        row.getDeclaredConstructor(String.class, Class.forName("neqsim.thermo.component.ComponentInterface"));
    ctor.setAccessible(true);
    List<Object> rows = new ArrayList<>();
    for (String name : system.getComponentNames()) {
      ComponentInterface component = system.getComponent(name);
      if (component == null || !(component.getNumberOfmoles() > 0.0)) {
        continue;
      }
      if (component.isIsTBPfraction() || component.isIsPlusFraction()) {
        rows.add(ctor.newInstance(name, component));
      }
    }
    return rows;
  }

  /** One private field of a row, read reflectively. */
  private static Object field(Object row, String name) throws Exception {
    Field found = Class.forName(CONTRIBUTION).getDeclaredField(name);
    found.setAccessible(true);
    return found.get(row);
  }

  private static double number(Object row, String name) throws Exception {
    return ((Number) field(row, name)).doubleValue();
  }

  private static List<Object> delump(List<Object> rows, int resolution) throws Exception {
    Method method = Class.forName(COMBINER).getDeclaredMethod("delumpContributions", List.class, int.class);
    method.setAccessible(true);
    @SuppressWarnings("unchecked")
    List<Object> out = (List<Object>) method.invoke(null, rows, resolution);
    return out;
  }

  private static List<Double> boundaries(List<Object> rows, int resolution) throws Exception {
    Method method =
        Class.forName(COMBINER)
            .getDeclaredMethod("determineReferenceEqualMassBoundaries", List.class, int.class);
    method.setAccessible(true);
    @SuppressWarnings("unchecked")
    List<Double> out = (List<Double>) method.invoke(null, rows, resolution);
    return out;
  }

  private static void rows(String prefix, List<Object> table) throws Exception {
    StringBuilder name = new StringBuilder();
    StringBuilder moles = new StringBuilder();
    StringBuilder mass = new StringBuilder();
    StringBuilder density = new StringBuilder();
    StringBuilder boiling = new StringBuilder();
    for (int i = 0; i < table.size(); i++) {
      if (i > 0) {
        name.append(", ");
        moles.append(", ");
        mass.append(", ");
        density.append(", ");
        boiling.append(", ");
      }
      name.append(field(table.get(i), "name"));
      moles.append(String.format("%.15g", number(table.get(i), "moles")));
      mass.append(String.format("%.15g", number(table.get(i), "molarMass")));
      density.append(String.format("%.15g", number(table.get(i), "density")));
      boiling.append(String.format("%.15g", number(table.get(i), "normalBoilingPoint")));
    }
    System.out.printf("%s_count = %d%n", prefix, table.size());
    System.out.printf("%s_name = [%s]%n", prefix, name);
    System.out.printf("%s_moles = [%s]%n", prefix, moles);
    System.out.printf("%s_molar_mass = [%s]%n", prefix, mass);
    System.out.printf("%s_density = [%s]%n", prefix, density);
    System.out.printf("%s_boiling_point = [%s]%n", prefix, boiling);
  }

  private static void boundaries(String prefix, List<Double> values) {
    StringBuilder out = new StringBuilder();
    for (int i = 0; i < values.size(); i++) {
      if (i > 0) {
        out.append(", ");
      }
      out.append(String.format("%.15g", values.get(i)));
    }
    System.out.printf("%s_count = %d%n", prefix, values.size());
    System.out.printf("%s = [%s]%n", prefix, out);
  }

  private static void sweep(String label, SystemInterface source, SystemInterface reference, int resolution)
      throws Exception {
    System.out.printf("# case = %s%n", label);
    List<Object> lumps = contributions(source);
    System.out.printf("resolution = %d%n", resolution);
    rows("lump", lumps);
    rows("fraction", delump(lumps, resolution));
    boundaries("boundary", boundaries(contributions(reference), resolution));
    System.out.println();
  }

  public static void main(String[] args) throws Exception {
    // Three lumps of unequal mass and molar mass, delumped four ways - the case that exercises
    // the rescale and the light-end bias at once.
    sweep("three_lumps_resolution_4",
        fluid(new double[][] {{1.0, 0.100, 0.75}, {2.0, 0.180, 0.82}, {0.5, 0.320, 0.88}}),
        fluid(new double[][] {{1.0, 0.100, 0.75}, {2.0, 0.180, 0.82}, {0.5, 0.320, 0.88}}),
        4);

    // A resolution of one is the class's own identity: the input list comes back unchanged, and
    // the boundaries are the boiling-point midpoints rather than the equal-mass cuts.
    sweep("a_resolution_of_one_is_the_identity",
        fluid(new double[][] {{1.0, 0.100, 0.75}, {2.0, 0.180, 0.82}}),
        fluid(new double[][] {{1.0, 0.100, 0.75}, {2.0, 0.180, 0.82}}),
        1);

    // A lump whose normal boiling point was never set: the sub-fractions carry the non-finite
    // value through and the sort falls back to the molar mass.
    SystemInterface unset = fluid(new double[][] {{1.0, 0.100, 0.75}, {1.0, 0.200, 0.84}});
    unset.getComponent(unset.getComponentNames()[0]).setNormalBoilingPoint(Double.NaN);
    sweep("a_lump_with_no_boiling_point", unset, unset, 3);

    // A single lump, which leaves no adjacent pair and therefore no boundary at all.
    sweep("one_lump_leaves_no_boundary",
        fluid(new double[][] {{1.0, 0.150, 0.80}}),
        fluid(new double[][] {{1.0, 0.150, 0.80}}),
        2);
  }
}
