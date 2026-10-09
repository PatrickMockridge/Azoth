// The plus-fraction split, as `Characterise.characterisePlusFraction` runs it.
//
// A cut's C7+ end is one pseudo-component described by three numbers - its mole fraction, its
// molar mass and its density - and characterisation divides it into a row of carbon-number cuts.
// `PedersenPlusModel` does that by solving two 2x2 Newton systems: an `AB` pair for the two
// coefficients of `z = exp(a + b*CN)`, and a `CD` pair for the two of `density = c + d*ln(CN)`.
//
// NeqSim's solver damps its steps by `iter/(iter+50)` for AB and `iter/(iter+5)` for CD, and its
// CD Jacobian leaves the `(1,1)` entry at zero. All three are reproduced rather than improved,
// so they are printed where the port can be held to them.
//
// Two sections per case:
//
//   * the **raw split** - the four coefficients and the model's own `z`, `M` and `dens` arrays,
//     over the half-open carbon-number range `[first, last)` the class fills. This is what a
//     characterisation kernel returns, and it is what `characterization.pedersen_plus_split` is
//     held to.
//   * the **lumped fluid** - what `Characterise` then puts into the system, which is the split
//     after `LumpingModel.generateLumpedComposition` has grouped it. That is a different id's
//     number, and it is here because the two are easy to confuse.
//
//   javac -proc:none -cp /path/to/neqsim-f0c7436.jar PlusFractionProbe.java
//   java -cp .:/path/to/neqsim-f0c7436.jar PlusFractionProbe > captures/plus_fraction_probe.tsv

import neqsim.thermo.characterization.Characterise;
import neqsim.thermo.characterization.LumpingModelInterface;
import neqsim.thermo.characterization.PlusFractionModelInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class PlusFractionProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  /** A bracketed vector, one `%.15g` per entry, over the half-open index range. */
  private static String vector(double[] values, int from, int to) {
    StringBuilder line = new StringBuilder();
    for (int i = from; i < to; i++) {
      if (i > from) {
        line.append(", ");
      }
      line.append(String.format("%.15g", values[i]));
    }
    return line.toString();
  }

  /** The carbon numbers the split was filled at, which is the half-open range itself. */
  private static String carbonNumbers(int from, int to) {
    StringBuilder line = new StringBuilder();
    for (int i = from; i < to; i++) {
      if (i > from) {
        line.append(", ");
      }
      line.append(i);
    }
    return line.toString();
  }

  /**
   * The split before it is lumped: the four solved coefficients, and the model's own per-cut
   * `z`, molar mass and specific gravity. A model that refused prints nothing but the refusal.
   */
  private static void rawSplit(Characterise characterise) {
    PlusFractionModelInterface model = characterise.getPlusFractionModel();
    System.out.printf("# raw split, model = %s%n", model.getName());
    System.out.printf("# first = %d last = %d%n",
        model.getFirstPlusFractionNumber(), model.getLastPlusFractionNumber());
    double[] coefs = model.getCoefs();
    System.out.printf("coefficient = [%s]%n", vector(coefs, 0, coefs.length));
    double[] z = model.getZ();
    double[] m = model.getM();
    double[] dens = model.getDens();
    if (z == null || m == null || dens == null) {
      System.out.println("# refused: the model returned no split");
      return;
    }
    int first = model.getFirstPlusFractionNumber();
    int last = model.getLastPlusFractionNumber();
    System.out.printf("number_of_plus_pseudocomponents = %.15g%n",
        model.getNumberOfPlusPseudocomponents());
    System.out.printf("carbon_number = [%s]%n", carbonNumbers(first, last));
    System.out.printf("cut_z = [%s]%n", vector(z, first, last));
    System.out.printf("cut_molar_mass = [%s]%n", vector(m, first, last));
    System.out.printf("cut_density = [%s]%n", vector(dens, first, last));
  }

  /**
   * One of `WhitsonGammaModel`'s switches, run over the same feed.
   *
   * The gamma model is the only plus model with settable shape parameters and a choice of density
   * correlation, so its defaults are the only configuration the `sweep` rows carry. Each variant
   * below moves one switch, so the port's handling of each is an oracle row rather than a default
   * that happens to be right.
   */
  private static void whitsonVariant(
      String label, double mPlus, double densityPlus, java.util.function.Consumer<Characterise> configure) {
    SystemInterface system = feed(mPlus, densityPlus);
    Characterise characterise = new Characterise(system);
    characterise.setPlusFractionModel("Whitson Gamma Model");
    configure.accept(characterise);
    System.out.printf("# gamma variant = %s, m_plus = %.15g kg/mol, density_plus = %.15g%n",
        label, mPlus, densityPlus);
    try {
      characterise.characterisePlusFraction();
    } catch (RuntimeException e) {
      System.out.printf("# refused: %s%n",
          e.getMessage() == null ? e.toString() : e.getMessage().split("\n")[0]);
      System.out.println();
      return;
    }
    System.out.printf("# selected_model = %s%n", characterise.getPlusFractionModel().getName());
    rawSplit(characterise);
    System.out.println();
  }

  /** A fluid with a C7+ end, described by the three numbers a plus fraction carries. */
  private static SystemInterface feed(double mPlus, double densityPlus) {
    SystemInterface system = new SystemSrkEos(298.15, 50.0);
    system.addComponent("methane", 60.0);
    system.addComponent("propane", 20.0);
    system.addComponent("n-butane", 10.0);
    // One pseudo-component stands for everything above it: 10 mol of a C20-scale cut.
    system.addPlusFraction("C20", 10.0, mPlus, densityPlus);
    system.setMixingRule(2);
    return system;
  }

  private static void sweep(String plusModel, double mPlus, double densityPlus) {
    SystemInterface system = feed(mPlus, densityPlus);
    Characterise characterise = new Characterise(system);
    characterise.setPlusFractionModel(plusModel);
    System.out.printf("# plus_model = %s, m_plus = %.15g kg/mol, density_plus = %.15g%n",
        plusModel, mPlus, densityPlus);
    try {
      characterise.characterisePlusFraction();
    } catch (RuntimeException e) {
      System.out.printf("# refused: %s%n", e.getMessage() == null ? e.toString() : e.getMessage().split("\n")[0]);
      System.out.println();
      return;
    }
    System.out.printf("# selected_model = %s%n", characterise.getPlusFractionModel().getName());
    rawSplit(characterise);
    // `LumpingModel.generateLumpedComposition`'s own state, which is what it decided rather
    // than what the system then reports.
    LumpingModelInterface lumping = characterise.getLumpingModel();
    System.out.printf("# lumping_model = %s%n", lumping.getName());
    // A plus model that returned without a split leaves the lumping model untouched, so its
    // arrays are null and `Characterise` never called `generateLumpedComposition` at all.
    if (lumping.getLumpedComponentNames() == null) {
      System.out.println("# no lumping: the plus model returned no split");
    } else {
      System.out.printf("# lump_names = %s%n", String.join(",", lumping.getLumpedComponentNames()));
      row("number_of_lumped_components", lumping.getNumberOfLumpedComponents());
      row("number_of_pseudo_components", lumping.getNumberOfPseudoComponents());
      for (int k = 0; k < lumping.getNumberOfLumpedComponents(); k++) {
        row("fraction_of_heavy_end[" + k + "]", lumping.getFractionOfHeavyEnd(k));
      }
    }
    int n = system.getPhase(0).getNumberOfComponents();
    row("component_count", n);
    row("total_moles", system.getNumberOfMoles());
    for (int i = 0; i < n; i++) {
      var component = system.getPhase(0).getComponent(i);
      String name = component.getComponentName();
      System.out.printf("# cut %s%n", name);
      row(name + ".z", component.getz());
      // **`getz()` is not the whole story for a component added after the last `init`.**
      // `addTBPfraction` adds *moles* - `totalNumberOfMoles * zPlus[k]` - and NeqSim derives `z`
      // from them, so a reading taken before the next initialisation says nothing about what was
      // added. Both are printed so the row is not read as either.
      row(name + ".moles", component.getNumberOfMolesInPhase());
      row(name + ".molar_mass", component.getMolarMass());
      row(name + ".density", component.getNormalLiquidDensity());
      row(name + ".tc", component.getTC());
      row(name + ".pc", component.getPC());
      row(name + ".acentric", component.getAcentricFactor());
      row(name + ".tb", component.getNormalBoilingPoint());
    }
    System.out.println();
  }

  public static void main(String[] args) {
    // kg/mol and kg/m3, which is the boundary `addTBPfraction` takes.
    double[][] feeds = {
      {0.4, 850.0},
      {0.2, 780.0},
      {0.65, 900.0},
    };
    for (double[] f : feeds) {
      for (String model : new String[] {"Pedersen", "Pedersen Heavy Oil", "Whitson Gamma Model"}) {
        sweep(model, f[0], f[1]);
      }
    }
    // The gamma model's switches, each over the feed the rows above use. `0.4` kg/mol is below
    // every model's `maxPlusMolarMass`, so none of these is swept away by `Characterise`'s swap.
    double[][] gammaFeeds = {{0.4, 850.0}, {0.2, 780.0}};
    for (double[] f : gammaFeeds) {
      whitsonVariant("default", f[0], f[1], c -> {});
      whitsonVariant("shape=0.7", f[0], f[1], c -> c.setGammaShapeParameter(0.7));
      whitsonVariant("minMW=84", f[0], f[1], c -> c.setGammaMinMW(84.0));
      whitsonVariant("density=soreide", f[0], f[1], c -> c.setGammaDensityModel("Soreide"));
      whitsonVariant("autoAlpha", f[0], f[1], c -> c.setAutoEstimateGammaAlpha(true));
    }
  }
}
