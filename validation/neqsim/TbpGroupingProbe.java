// `Phase.groupTBPfractions` - a phase's components binned by normal boiling point.
//
// Twenty bins, indexed 6 to 19, each opened by a boiling-point threshold in `C` and never by a
// carbon number. A component below 69.2 `C` reaches no bin at all, so the light end of a fluid is
// not merely in bin zero - it is absent from the answer.
//
// **Nothing in NeqSim calls this, and it is not on the interface its own accessor returns.**
// `SystemInterface.getPhase` hands back a `PhaseInterface`, which declares no such method, so
// reaching it takes a cast to the concrete `Phase`. The other two implementations of the name
// are worse: `PlusCharacterize`'s returns `true` and computes nothing, and `TBPCharacterize`'s
// is reachable only from inside its own package, where nothing calls it.
//
//   javac -proc:none -cp /path/to/neqsim-f0c7436.jar TbpGroupingProbe.java
//   java -cp .:/path/to/neqsim-f0c7436.jar TbpGroupingProbe > captures/tbp_grouping_probe.tsv

import neqsim.thermo.phase.Phase;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class TbpGroupingProbe {

  private static final double[] ZERO_CELSIUS = {273.15};

  /** One fluid, its components' boiling points and mole fractions, and its twenty bins. */
  private static void sweep(String label, String[] names, double[] amounts) {
    SystemInterface system = new SystemSrkEos(298.15, 1.0);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], amounts[i]);
    }
    system.setMixingRule(2);
    system.init(0);

    System.out.printf("# grouping = %s%n", label);
    System.out.printf("component_count = %d%n", names.length);
    System.out.printf("boiling_point = [%s]%n", boilings(system));
    System.out.printf("mole_fraction = [%s]%n", moleFractions(system));
    double[] bins = ((Phase) system.getPhase(0)).groupTBPfractions();
    System.out.printf("group_fraction = [%s]%n", vector(bins));
    for (int i = 0; i < bins.length; i++) {
      if (bins[i] != 0.0) {
        System.out.printf("bin[%d] = %.15g%n", i, bins[i]);
      }
    }
    System.out.println();
  }

  /**
   * The stored boiling points, **in kelvin**, because that is what the class subtracts 273.15
   * from. Printing the `C` value instead would make the case state a number the kernel then has to
   * round-trip back, which is a last-bit difference at every threshold.
   */
  private static String boilings(SystemInterface system) {
    StringBuilder line = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      if (i > 0) {
        line.append(", ");
      }
      line.append(String.format("%.15g",
          system.getPhase(0).getComponent(i).getNormalBoilingPoint()));
    }
    return line.toString();
  }

  private static String moleFractions(SystemInterface system) {
    StringBuilder line = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      if (i > 0) {
        line.append(", ");
      }
      line.append(String.format("%.15g", system.getPhase(0).getComponent(i).getx()));
    }
    return line.toString();
  }

  private static String vector(double[] values) {
    StringBuilder line = new StringBuilder();
    for (int i = 0; i < values.length; i++) {
      if (i > 0) {
        line.append(", ");
      }
      line.append(String.format("%.15g", values[i]));
    }
    return line.toString();
  }

  public static void main(String[] args) {
    if (ZERO_CELSIUS[0] != 273.15) {
      throw new IllegalStateException("the probe's constants moved");
    }
    // A spread that lands in eleven of the fourteen reachable bins, with the light end left
    // where the class leaves it: methane to n-hexane are below the lowest threshold.
    sweep("a_spread_of_n_alkanes",
        new String[] {"methane", "propane", "n-butane", "n-pentane", "n-hexane", "n-heptane",
            "n-octane", "n-nonane", "n-decane", "n-undecane", "n-dodecane", "n-tridecane",
            "n-tetradecane", "n-pentadecane", "n-hexadecane"},
        new double[] {10.0, 8.0, 7.0, 6.0, 5.0, 5.0, 4.0, 4.0, 4.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0});
    // Two components inside one bin, which is what makes the bin a sum rather than a bucket: the
    // class has no carbon number to separate them by.
    sweep("two_inside_one_bin",
        new String[] {"n-heptane", "n-octane", "n-nonane"},
        new double[] {2.0, 3.0, 5.0});
    // One component only, and a light one: every bin stays zero, and the answer is zeros rather
    // than a refusal.
    sweep("a_fluid_below_the_lowest_threshold",
        new String[] {"methane", "propane"},
        new double[] {1.0, 1.0});
  }
}
