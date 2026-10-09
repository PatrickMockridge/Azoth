// `TbpClosure`'s four members, over a grid of boiling points and specific gravities.
//
// The closure is the inverse of a cut's property correlation: given a normal boiling point and a
// specific gravity, what molar mass does the cut have? Three of the four members answer by
// bisection over `[0.010, 0.800]` kg/mol and one answers in closed form, and all four refuse the
// *other* direction except the 1980 Riazi-Daubert pair - so the refusals are captured here too,
// as the measurement the port's `[[unported]]` rows are held to.
//
//   javac -proc:none -cp neqsim-f0c7436.jar TbpClosureProbe.java
//   java -cp .:neqsim-f0c7436.jar TbpClosureProbe > captures/tbp_closure_probe.tsv

import neqsim.thermo.characterization.TBPfractionModel;
import neqsim.thermo.characterization.TBPModelInterface;
import neqsim.thermo.characterization.TbpClosure;

public class TbpClosureProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    // The PedersenSRK model, which is what `Characterise` selects, for the `TBP_MODEL` closure.
    TBPModelInterface model = new TBPfractionModel().getModel("PedersenSRK");

    // K, and specific gravity in g/cm3. The boiling points straddle 540 g/mol's equivalent, and
    // the gravities straddle the 0.80 kg/mol upper bracket in the other direction.
    double[] boilingPoints = {300.0, 400.0, 500.0, 600.0, 700.0, 800.0};
    double[] gravities = {0.65, 0.75, 0.85, 0.95};

    for (TbpClosure closure : TbpClosure.values()) {
      System.out.printf("# closure = %s%n", closure.name());
      System.out.printf("# supports_density_from_molar_mass = %b%n", closure.supportsDensityFromMolarMass());
      for (int i = 0; i < boilingPoints.length; i++) {
        for (int j = 0; j < gravities.length; j++) {
          double tb = boilingPoints[i];
          double sg = gravities[j];
          System.out.printf("# tb = %.15g K, sg = %.15g%n", tb, sg);
          try {
            row("molar_mass", closure.calcMolarMass(tb, sg, model));
          } catch (RuntimeException e) {
            System.out.printf("molar_mass = # refused: %s%n", firstLine(e));
          }
          try {
            row("density", closure.calcDensity(tb, 0.2, model));
          } catch (RuntimeException e) {
            System.out.printf("density = # refused: %s%n", firstLine(e));
          }
          System.out.println();
        }
      }
      // The forward direction, which every member shares and which the bisection inverts.
      for (double molarMass : new double[] {0.05, 0.1, 0.2, 0.4, 0.7}) {
        for (double sg : gravities) {
          row("soreide_tb[" + molarMass + "," + sg + "]", TbpClosure.calcBoilingPointSoreide(molarMass, sg));
        }
      }
      System.out.println();
    }
  }

  /** The first line of an exception's message: the rest is a suggestion, not the measurement. */
  private static String firstLine(RuntimeException e) {
    String message = e.getMessage();
    if (message == null) {
      message = e.toString();
    }
    return message.split("\n")[0];
  }
}
