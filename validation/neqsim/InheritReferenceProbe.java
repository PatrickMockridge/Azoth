// `PseudoComponentCombiner.characterizeToReference` with `inheritReferenceProperties` on.
//
// `characterizeToReference(source, reference)` passes the flag `false` and takes each group's own
// molar mass and density - its accumulated mass over its accumulated moles, and that mass over the
// volume its rows occupy. The public overload that takes a `CharacterizationOptions` passes it
// through, and with it **true** every group takes the *reference* row's own molar mass and
// density instead. The walk is identical either way: the same boundaries from the same reference
// keys, the same binning, the same groups. Only the two numbers each group reports change.
//
//   javac -proc:none -cp /path/to/neqsim-f0c7436.jar InheritReferenceProbe.java
//   java -cp .:/path/to/neqsim-f0c7436.jar InheritReferenceProbe > captures/inherit_reference_probe.tsv

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import neqsim.thermo.characterization.CharacterizationOptions;
import neqsim.thermo.characterization.PseudoComponentCombiner;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class InheritReferenceProbe {

  private record Row(String name, double moles, double molarMass, double density, double tb) {}

  private static SystemInterface fluid(double[][] cuts) {
    SystemInterface system = new SystemSrkEos(298.15, 1.0);
    for (int i = 0; i < cuts.length; i++) {
      system.addTBPfraction("C" + (7 + i), cuts[i][0], cuts[i][1], cuts[i][2]);
    }
    system.setMixingRule(2);
    system.init(0);
    return system;
  }

  private static List<Row> table(SystemInterface system) {
    List<Row> rows = new ArrayList<>();
    for (String name : system.getComponentNames()) {
      var component = system.getComponent(name);
      if (component == null || !(component.getNumberOfmoles() > 0.0)) {
        continue;
      }
      if (component.isIsTBPfraction() || component.isIsPlusFraction()) {
        rows.add(new Row(name, component.getNumberOfmoles(), component.getMolarMass(),
            component.getNormalLiquidDensity(), component.getNormalBoilingPoint()));
      }
    }
    rows.sort(Comparator.comparingDouble(InheritReferenceProbe::sortingKey));
    return rows;
  }

  private static double sortingKey(Row row) {
    double key = row.tb();
    if (!Double.isFinite(key) || key <= 0.0) {
      key = row.molarMass();
    }
    return key;
  }

  private static void rows(String prefix, List<Row> table) {
    System.out.printf("%s_count = %d%n", prefix, table.size());
    StringBuilder moles = new StringBuilder();
    StringBuilder mass = new StringBuilder();
    StringBuilder density = new StringBuilder();
    for (int i = 0; i < table.size(); i++) {
      if (i > 0) {
        moles.append(", ");
        mass.append(", ");
        density.append(", ");
      }
      moles.append(String.format("%.15g", table.get(i).moles()));
      mass.append(String.format("%.15g", table.get(i).molarMass()));
      density.append(String.format("%.15g", table.get(i).density()));
    }
    System.out.printf("%s_moles = [%s]%n", prefix, moles);
    System.out.printf("%s_molar_mass = [%s]%n", prefix, mass);
    System.out.printf("%s_density = [%s]%n", prefix, density);
  }

  private static void run(String label, SystemInterface source, SystemInterface reference, boolean inherit) {
    System.out.printf("# case = %s%n", label);
    System.out.printf("inherit_reference_properties = %s%n", inherit);
    rows("source", table(source));
    rows("reference", table(reference));
    CharacterizationOptions options =
        CharacterizationOptions.builder().inheritReferenceProperties(inherit).build();
    try {
      SystemInterface out = PseudoComponentCombiner.characterizeToReference(source, reference, options);
      rows("result", table(out));
    } catch (RuntimeException e) {
      System.out.printf("# refused: %s%n",
          e.getMessage() == null ? e.toString() : e.getMessage().split("\n")[0]);
    }
    System.out.println();
  }

  /** The same five-cut source and three-cut reference the combiner probe uses. */
  private static SystemInterface fiveCutSource() {
    return fluid(new double[][] {{1.0, 0.10, 0.70}, {1.5, 0.14, 0.76}, {2.0, 0.18, 0.80},
        {1.0, 0.24, 0.84}, {0.5, 0.32, 0.88}});
  }

  private static SystemInterface threeCutReference() {
    return fluid(new double[][] {{2.0, 0.12, 0.74}, {2.0, 0.20, 0.82}, {1.0, 0.30, 0.87}});
  }

  public static void main(String[] args) {
    // The same fluid twice, so the two tables differ in exactly the two numbers the flag changes.
    run("five_cuts_onto_three", fiveCutSource(), threeCutReference(), false);
    run("five_cuts_onto_three", fiveCutSource(), threeCutReference(), true);

    // A source entirely below the reference's first cut: one group, and the flag still takes the
    // *first* reference row's properties for it.
    run("a_source_below_the_first_cut",
        fluid(new double[][] {{3.0, 0.09, 0.68}, {1.0, 0.10, 0.70}}),
        fluid(new double[][] {{1.0, 0.20, 0.82}, {1.0, 0.30, 0.87}}),
        true);
  }
}
