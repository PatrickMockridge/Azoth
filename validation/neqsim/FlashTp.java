// Reproduces the numbers in validation/eos/*_against_neqsim.json.
//
//     javac -proc:none -cp neqsim-3.20.0.jar FlashTp.java
//     java -cp .:neqsim-3.20.0.jar FlashTp
//
// This is the only file in the repository that is not part of either
// implementation, and it is here on purpose: it is the ground truth the port is
// measured against, and a measurement nobody can repeat is a claim rather than a
// measurement. It is not built by CI - it needs a JVM and NeqSim's jar, neither of
// which belongs in a build - so the numbers it prints are recorded in the case
// files beside it.
//
// The jar is NeqSim 3.20.0, the version `databank/manifest.toml` names as the source
// of every vendored table. `databank/sources/neqsim/` holds its data files; the jar
// is not vendored, because it is 20 MB of compiled Java that nothing here links
// against and its licence is carried by the attribution in NOTICE.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermo.system.SystemRKEos;
import neqsim.thermo.system.SystemPrEosvolcor;
import neqsim.thermo.system.SystemSrkPenelouxEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class FlashTp {
  /** The cubic `cubic` selects: "pr", "srk" or "rk". */
  static SystemInterface system(String cubic, double temperatureK, double pressureBar) {
    switch (cubic) {
      case "srk":
        return new SystemSrkEos(temperatureK, pressureBar);
      case "rk":
        return new SystemRKEos(temperatureK, pressureBar);
      default:
        return new SystemPrEos(temperatureK, pressureBar);
    }
  }

  /** One PT flash, printed in the fields azoth's `eos.pt_flash` reports. */
  static void flash(String label, double temperatureK, double pressureBar, String[] names,
      double[] moles, String cubic, int alphaTerm) {
    SystemInterface fluid = system(cubic, temperatureK, pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    // "classic" is the van der Waals one-fluid mixing rule with the database's
    // interaction parameters. `setAttractiveTerm` overrides the cubic's default alpha
    // correlation - the `attractiveTermNumber` NeqSim decouples from the cubic shape.
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(alphaTerm);

    new ThermodynamicOperations(fluid).TPflash();
    fluid.initProperties();

    StringBuilder k = new StringBuilder();
    StringBuilder x = new StringBuilder();
    StringBuilder y = new StringBuilder();
    for (int i = 0; i < names.length; i++) {
      if (i > 0) {
        k.append(", ");
        x.append(", ");
        y.append(", ");
      }
      k.append(fluid.getPhase(0).getComponent(i).getK());
      x.append(fluid.getPhase(1).getComponent(i).getx());
      y.append(fluid.getPhase(0).getComponent(i).getx());
    }
    System.out.println(label);
    System.out.println("  beta      " + fluid.getBeta());
    System.out.println("  k         [" + k + "]");
    System.out.println("  x         [" + x + "]");
    System.out.println("  y         [" + y + "]");
    System.out.println("  z_liquid  " + fluid.getPhase(1).getZ());
    System.out.println("  z_vapour  " + fluid.getPhase(0).getZ());
  }

  /** One PT flash with its molar enthalpy and entropy, for `eos.molar_enthalpy_entropy`. */
  static void enthalpy(String label, double temperatureK, double pressureBar, String[] names,
      double[] moles, String cubic, int alphaTerm) {
    SystemInterface fluid = system(cubic, temperatureK, pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(alphaTerm);

    new ThermodynamicOperations(fluid).TPflash();
    fluid.initProperties();

    System.out.println(label);
    System.out.println("  H        " + fluid.getEnthalpy("J/mol"));
    System.out.println("  S        " + fluid.getEntropy("J/molK"));
    System.out.println("  z_vapour " + fluid.getPhase(0).getZ());
  }

  /** The per-component Peneloux volume-translation parameter, for `eos.*_peneloux_shift`. */
  static void volcorr() {
    SystemInterface pr = new SystemPrEosvolcor(300.0, 10.0);
    SystemInterface srk = new SystemSrkPenelouxEos(300.0, 10.0);
    String[] names = {"methane", "n-butane", "propane"};
    for (String name : names) {
      pr.addComponent(name, 1.0);
      srk.addComponent(name, 1.0);
    }
    System.out.println("Peneloux volume correction (m**3/mol; NeqSim prints 1e5 * this):");
    for (int i = 0; i < names.length; i++) {
      System.out.println("  PR  " + names[i] + "  c = " + pr.getPhase(0).getComponent(i).getVolumeCorrection());
      System.out.println("  SRK " + names[i] + "  c = " + srk.getPhase(0).getComponent(i).getVolumeCorrection());
    }
  }

  /** The pure-component correlations `eos.heat_of_vaporization` and `eos.liquid_heat_capacity` port. */
  static void corr() {
    SystemInterface fluid = new SystemPrEos(300.0, 10.0);
    fluid.addComponent("n-butane", 1.0);
    fluid.addComponent("water", 1.0);
    System.out.println("n-butane  hov(300 K) = "
        + fluid.getPhase(0).getComponent(0).getPureComponentHeatOfVaporization(300.0) + " J/mol");
    System.out.println("water     cpL(300 K) = "
        + fluid.getPhase(0).getComponent(1).getPureComponentCpLiquid(300.0) + " J/(mol*K)");
  }

  public static void main(String[] args) {
    volcorr();
    corr();
    flash("methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 1);
    flash("propane, 1.0, 300 K, 9 bar",
        300.0, 9.0, new String[] {"propane"}, new double[] {1.0}, "pr", 1);
    enthalpy("propane, 1.0, 300 K, 9 bar",
        300.0, 9.0, new String[] {"propane"}, new double[] {1.0}, "pr", 1);
    flash("SRK methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "srk", 0);
    flash("SRK propane, 1.0, 300 K, 9 bar",
        300.0, 9.0, new String[] {"propane"}, new double[] {1.0}, "srk", 0);
    flash("RK methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "rk", 5);
    flash("RK propane, 1.0, 300 K, 9 bar",
        300.0, 9.0, new String[] {"propane"}, new double[] {1.0}, "rk", 5);
    flash("SRK+Twu methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "srk", 14);
    flash("SRK+TwuCoon methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "srk", 11);
    flash("PR+Gassem2001 methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 8);
    flash("PR+LeeKesler methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 21);
    flash("PR+Danesh methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 9);
    flash("PR+Danesh nc12, 1.0, 500 K, 15 bar",
        500.0, 15.0, new String[] {"nc12"}, new double[] {1.0}, "pr", 9);
    flash("PR+Schwartzentruber methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 2);
    flash("PR+Mollerup methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 3);
    flash("PR+MatCop methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 4);
    flash("PR+MatCopPR methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 13);
    flash("PR+MatCopPRUMR methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 17);
    flash("PR+MatCop5PRUMR methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 22);
    flash("PR+Delft1998 methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4}, "pr", 7);
    flash("PR+PR78 nc12, 1.0, 500 K, 15 bar",
        500.0, 15.0, new String[] {"nc12"}, new double[] {1.0}, "pr", 6);
  }
}
