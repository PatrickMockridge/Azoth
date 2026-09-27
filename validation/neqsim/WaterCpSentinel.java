import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;
import neqsim.thermo.system.SystemSrkEos;

/**
 * How many of NeqSim's own components come back with <em>water's</em> ideal-gas heat capacity?
 *
 * <p>`COMP.csv` gives 131 of its 389 rows the whole of water's ideal-gas Cp polynomial
 * (`36.54003, -0.034802404, 0.000116811, -1.3e-07, 5.254448e-11` - the same five numbers
 * `devtools/generate_water_caloric_alpha_reference.py` fits for water), and 130 of those rows are
 * substances other than water. The critical constants are not shared - only water's row also
 * carries water's `Tc` and `Pc` - so it is the caloric column alone.
 *
 * <p>This probe asks whether the table is what the class actually returns, or whether something
 * downstream overrides it. It builds one single-component `SystemSrkEos` per name at 298.15 K and
 * 1 bara and compares `getCp0` against water's, to the last bit. **Two different molecules cannot
 * share an ideal-gas heat capacity**, so an exact match is not a coincidence and is not a
 * tolerance question.
 *
 * <p>Reads the component table's `NAME` column, which is quoted and carries commas inside names
 * (`"1,2,4-trimethylbenzene"`), so the split is quote-aware rather than a `split(",")`.
 *
 * <pre>
 * cd validation/neqsim
 * javac -proc:none -cp neqsim-f0c7436.jar WaterCpSentinel.java
 * java -cp .:neqsim-f0c7436.jar WaterCpSentinel &gt; captures/water_cp_sentinel.tsv
 * </pre>
 *
 * <p>The table path defaults to the vendored copy and may be given as the first argument.
 */
public class WaterCpSentinel {

  public static void main(String[] args) throws Exception {
    Path table =
        Path.of(args.length > 0 ? args[0] : "../../databank/sources/neqsim/COMP.csv");
    double temperature = 298.15;

    List<String> names = new ArrayList<>();
    for (String line : Files.readAllLines(table)) {
      List<String> fields = splitCsv(line);
      if (fields.size() < 2 || "NAME".equals(fields.get(1).trim())) {
        continue;
      }
      String name = fields.get(1).trim();
      if (!name.isEmpty()) {
        names.add(name);
      }
    }

    double water = cp0("water", temperature);
    System.out.printf("table = %s%n", table);
    System.out.printf("temperature_K = %.2f%n", temperature);
    System.out.printf("water_cp0 = %.17g%n", water);
    System.out.printf("names_in_table = %d%n", names.size());

    TreeSet<String> matching = new TreeSet<>();
    int unusable = 0;
    for (String name : names) {
      try {
        if (cp0(name, temperature) == water) {
          matching.add(name);
        }
      } catch (Throwable failure) {
        unusable++;
      }
    }

    System.out.printf("sharing_water_cp0 = %d%n", matching.size());
    System.out.printf("sharing_and_not_water = %d%n", matching.size() - (matching.contains("water") ? 1 : 0));
    System.out.printf("unusable = %d%n", unusable);
    for (String name : matching) {
      System.out.printf("  %s%n", name);
    }
  }

  /** Water's ideal-gas heat capacity at the probe's own temperature, for one component. */
  private static double cp0(String name, double temperature) {
    SystemSrkEos system = new SystemSrkEos(temperature, 1.0);
    system.addComponent(name, 1.0);
    system.createDatabase(true);
    return system.getPhase(0).getComponent(0).getCp0(temperature);
  }

  /**
   * Split one CSV line on commas that are not inside quotes, dropping the quotes.
   *
   * <p>Not a general CSV reader: the table quotes the columns whose values carry a comma and
   * leaves the rest bare, which is all this has to survive.
   */
  private static List<String> splitCsv(String line) {
    List<String> fields = new ArrayList<>();
    StringBuilder field = new StringBuilder();
    boolean quoted = false;
    for (int i = 0; i < line.length(); i++) {
      char c = line.charAt(i);
      if (c == '"') {
        quoted = !quoted;
      } else if (c == ',' && !quoted) {
        fields.add(field.toString());
        field.setLength(0);
      } else {
        field.append(c);
      }
    }
    fields.add(field.toString());
    return fields;
  }
}
