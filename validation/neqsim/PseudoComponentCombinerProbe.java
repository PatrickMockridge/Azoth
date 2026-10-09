// `PseudoComponentCombiner.characterizeToReference` - one fluid re-cut onto another's slate.
//
// The combiner is 1841 lines and every public method takes or returns a `SystemInterface`, so
// the arithmetic has to be read out of the middle: `extractComponents` lifts each fluid's pseudo
// components into a table keyed by a boiling point that falls back to the molar mass,
// `determineReferenceBoundaries` puts a cut between each adjacent pair of the *reference's* keys,
// and `distributeToProfiles` walks the source table binning each row into the group its key falls
// in. A group's molar mass is its accumulated mass over its accumulated moles, and its density is
// the mass over the volume the rows would occupy.
//
//   javac -proc:none -cp /path/to/neqsim-f0c7436.jar PseudoComponentCombinerProbe.java
//   java -cp .:/path/to/neqsim-f0c7436.jar PseudoComponentCombinerProbe > captures/pseudo_component_combiner_probe.tsv

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import neqsim.thermo.characterization.PseudoComponentCombiner;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class PseudoComponentCombinerProbe {

  private record Row(String name, double moles, double molarMass, double density, double tb) {}

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

  /** The combiner's own row, read the way `extractComponents` reads it. */
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
    rows.sort(Comparator.comparingDouble(PseudoComponentCombinerProbe::sortingKey));
    return rows;
  }

  /** `PseudoComponentContribution.sortingKey`: the boiling point, or the molar mass without one. */
  private static double sortingKey(Row row) {
    double key = row.tb();
    if (!Double.isFinite(key) || key <= 0.0) {
      key = row.molarMass();
    }
    return key;
  }

  private static void rows(String prefix, List<Row> table) {
    System.out.printf("%s_count = %d%n", prefix, table.size());
    StringBuilder key = new StringBuilder();
    StringBuilder moles = new StringBuilder();
    StringBuilder mass = new StringBuilder();
    StringBuilder density = new StringBuilder();
    for (int i = 0; i < table.size(); i++) {
      if (i > 0) {
        key.append(", ");
        moles.append(", ");
        mass.append(", ");
        density.append(", ");
      }
      key.append(String.format("%.15g", sortingKey(table.get(i))));
      moles.append(String.format("%.15g", table.get(i).moles()));
      mass.append(String.format("%.15g", table.get(i).molarMass()));
      density.append(String.format("%.15g", table.get(i).density()));
    }
    System.out.printf("%s_sorting_key = [%s]%n", prefix, key);
    System.out.printf("%s_moles = [%s]%n", prefix, moles);
    System.out.printf("%s_molar_mass = [%s]%n", prefix, mass);
    System.out.printf("%s_density = [%s]%n", prefix, density);
  }

  private static void sweep(String label, SystemInterface source, SystemInterface reference) {
    System.out.printf("# case = %s%n", label);
    rows("source", table(source));
    rows("reference", table(reference));
    try {
      SystemInterface out = PseudoComponentCombiner.characterizeToReference(source, reference);
      System.out.printf("result_count = %d%n", table(out).size());
      List<Row> result = table(out);
      StringBuilder moles = new StringBuilder();
      StringBuilder mass = new StringBuilder();
      StringBuilder density = new StringBuilder();
      for (int i = 0; i < result.size(); i++) {
        if (i > 0) {
          moles.append(", ");
          mass.append(", ");
          density.append(", ");
        }
        moles.append(String.format("%.15g", result.get(i).moles()));
        mass.append(String.format("%.15g", result.get(i).molarMass()));
        density.append(String.format("%.15g", result.get(i).density()));
      }
      System.out.printf("result_moles = [%s]%n", moles);
      System.out.printf("result_molar_mass = [%s]%n", mass);
      System.out.printf("result_density = [%s]%n", density);
    } catch (RuntimeException e) {
      System.out.printf("# refused: %s%n",
          e.getMessage() == null ? e.toString() : e.getMessage().split("\n")[0]);
    }
    System.out.println();
  }

  public static void main(String[] args) {
    // A five-cut source onto a three-cut reference, which is the case the method exists for.
    sweep("five_cuts_onto_three",
        fluid(new double[][] {{1.0, 0.10, 0.70}, {1.5, 0.14, 0.76}, {2.0, 0.18, 0.80},
            {1.0, 0.24, 0.84}, {0.5, 0.32, 0.88}}),
        fluid(new double[][] {{2.0, 0.12, 0.74}, {2.0, 0.20, 0.82}, {1.0, 0.30, 0.87}}));

    // A source lighter than every reference cut, which lands entirely in the first group.
    sweep("a_source_below_the_first_cut",
        fluid(new double[][] {{3.0, 0.09, 0.68}, {1.0, 0.10, 0.70}}),
        fluid(new double[][] {{1.0, 0.20, 0.82}, {1.0, 0.30, 0.87}}));

    // A single reference cut, which leaves nothing to cut on and one group holding everything.
    sweep("one_reference_cut",
        fluid(new double[][] {{1.0, 0.10, 0.70}, {1.0, 0.20, 0.82}}),
        fluid(new double[][] {{1.0, 0.25, 0.85}}));
  }
}
