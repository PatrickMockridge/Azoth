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
//   javac -proc:none -cp neqsim-f0c7436.jar PlusFractionProbe.java
//   java -cp .:neqsim-f0c7436.jar PlusFractionProbe > captures/plus_fraction_probe.tsv

import neqsim.thermo.characterization.Characterise;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class PlusFractionProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
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
    int n = system.getPhase(0).getNumberOfComponents();
    row("component_count", n);
    for (int i = 0; i < n; i++) {
      var component = system.getPhase(0).getComponent(i);
      String name = component.getComponentName();
      System.out.printf("# cut %s%n", name);
      row(name + ".z", component.getz());
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
  }
}
