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
import neqsim.thermo.system.SystemGEWilson;
import neqsim.thermo.system.SystemNRTL;
import neqsim.thermo.system.SystemUNIFAC;
import neqsim.thermo.system.SystemUNIFACpsrk;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermo.system.SystemRKEos;
import neqsim.thermo.system.SystemPrEosvolcor;
import neqsim.thermo.system.SystemSrkPenelouxEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermo.util.empiric.NitricSulfuricAcidVaporPressure;
import neqsim.physicalproperties.methods.gasphysicalproperties.viscosity.ChungViscosityMethod;
import neqsim.physicalproperties.methods.gasphysicalproperties.conductivity.ChungConductivityMethod;
import neqsim.physicalproperties.system.PhysicalProperties;

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

  /** The four Antoine forms `eos.antoine_vapor_pressure` ports, one component each. */
  static void antoine() {
    SystemInterface fluid = new SystemPrEos(300.0, 10.0);
    fluid.addComponent("methane", 1.0);
    fluid.addComponent("ethylene", 1.0);
    fluid.addComponent("nc12", 1.0);
    fluid.addComponent("223-TM-C4", 1.0);
    String[] forms = {"pow10", "pow10KPa", "log", "loglog"};
    for (int i = 0; i < 4; i++) {
      System.out.println(fluid.getPhase(0).getComponent(i).getName() + "  (" + forms[i]
          + ")  psat(300 K) = "
          + fluid.getPhase(0).getComponent(i).getAntoineVaporPressure(300.0) + " bar");
    }
  }

  /** The pure-component Chung viscosity, for `eos.chung_viscosity`. */
  static void chung() {
    String[] names = {"methane", "propane", "MDEA"};
    double[] temps = {300.0, 300.0, 400.0};
    double[] press = {10.0, 9.0, 10.0};
    System.out.println("Chung pure-component viscosity (V in m**3/mol; NeqSim's molarVolume field is 1e5 * this):");
    for (int i = 0; i < names.length; i++) {
      SystemInterface fluid = new SystemPrEos(temps[i], press[i]);
      fluid.addComponent(names[i], 1.0);
      fluid.setMixingRule("classic");
      new ThermodynamicOperations(fluid).TPflash();
      fluid.initProperties();
      PhysicalProperties pp = fluid.getPhase(0).getPhysicalProperties();
      ChungViscosityMethod chung = new ChungViscosityMethod(pp);
      chung.initChungPureComponentViscosity();
      System.out.println("  " + names[i] + "  V = " + fluid.getPhase(0).getMolarVolume() * 1e-5
          + "  mu = " + chung.pureComponentViscosity[0] * 1e-7 + " Pa*s");
    }
  }

  /** The Wilke mixture viscosity, for `eos.wilke_viscosity`. */
  static void wilke() {
    String[][] names = {
      {"methane", "propane"}, {"methane", "n-butane", "propane"}, {"methane", "propane"}};
    double[][] moles = {{0.5, 0.5}, {0.5, 0.3, 0.2}, {0.8, 0.2}};
    double[] temps = {300.0, 350.0, 300.0};
    double[] press = {10.0, 15.0, 10.0};
    System.out.println("Wilke mixture viscosity (V in m**3/mol; NeqSim's molarVolume field is 1e5 * this):");
    for (int k = 0; k < names.length; k++) {
      SystemInterface fluid = new SystemPrEos(temps[k], press[k]);
      for (int i = 0; i < names[k].length; i++) {
        fluid.addComponent(names[k][i], moles[k][i]);
      }
      fluid.setMixingRule("classic");
      new ThermodynamicOperations(fluid).TPflash();
      fluid.initProperties();
      PhysicalProperties pp = fluid.getPhase(0).getPhysicalProperties();
      ChungViscosityMethod chung = new ChungViscosityMethod(pp);
      System.out.println("  " + String.join("+", names[k]) + "  V = "
          + fluid.getPhase(0).getMolarVolume() * 1e-5 + "  mu = " + chung.calcViscosity() + " Pa*s");
    }
  }

  /** The pure-component Chung conductivity, for `eos.chung_conductivity`. */
  static void chungCond() {
    String[] names = {"methane", "propane", "methanol"};
    double[] temps = {300.0, 300.0, 350.0};
    double[] press = {10.0, 9.0, 10.0};
    System.out.println("Chung pure-component conductivity (W/(m*K)):");
    for (int i = 0; i < names.length; i++) {
      SystemInterface fluid = new SystemPrEos(temps[i], press[i]);
      fluid.addComponent(names[i], 1.0);
      fluid.setMixingRule("classic");
      new ThermodynamicOperations(fluid).TPflash();
      fluid.initProperties();
      PhysicalProperties pp = fluid.getPhase(0).getPhysicalProperties();
      ChungConductivityMethod cond = new ChungConductivityMethod(pp);
      cond.calcPureComponentConductivity();
      System.out.println("  " + names[i] + "  Cv0 = "
          + fluid.getPhase(0).getComponent(0).getCv0(temps[i])
          + "  k = " + cond.pureComponentConductivity[0] + " W/(m*K)");
    }
  }

  /** The Mason-Saxena mixture conductivity, for `eos.mason_saxena_conductivity`. */
  static void masonSaxena() {
    String[][] names = {
      {"methane", "propane"}, {"methane", "n-butane", "propane"}, {"methane", "propane"}};
    double[][] moles = {{0.5, 0.5}, {0.5, 0.3, 0.2}, {0.8, 0.2}};
    double[] temps = {300.0, 350.0, 300.0};
    double[] press = {10.0, 15.0, 10.0};
    System.out.println("Mason-Saxena mixture conductivity (W/(m*K)):");
    for (int k = 0; k < names.length; k++) {
      SystemInterface fluid = new SystemPrEos(temps[k], press[k]);
      for (int i = 0; i < names[k].length; i++) {
        fluid.addComponent(names[k][i], moles[k][i]);
      }
      fluid.setMixingRule("classic");
      new ThermodynamicOperations(fluid).TPflash();
      fluid.initProperties();
      PhysicalProperties pp = fluid.getPhase(0).getPhysicalProperties();
      ChungConductivityMethod cond = new ChungConductivityMethod(pp);
      System.out.println("  " + String.join("+", names[k])
          + "  k = " + cond.calcConductivity() + " W/(m*K)");
    }
  }

  /** The NRTL activity coefficients, for `eos.nrtl_activity_coefficients`. */
  static void nrtl() {
    String[][] names = {{"methanol", "water"}, {"ethanol", "water"}};
    double[][] moles = {{0.5, 0.5}, {0.5, 0.5}};
    System.out.println("NRTL activity coefficients (liquid phase; gamma = exp(ln gamma)):");
    for (int k = 0; k < names.length; k++) {
      SystemNRTL system = new SystemNRTL(298.15, 1.0);
      for (int i = 0; i < names[k].length; i++) {
        system.addComponent(names[k][i], moles[k][i]);
      }
      system.createDatabase(true);
      system.setMixingRule("classic");
      new ThermodynamicOperations(system).TPflash();
      StringBuilder x = new StringBuilder();
      StringBuilder gamma = new StringBuilder();
      for (int i = 0; i < names[k].length; i++) {
        if (i > 0) {
          x.append(", ");
          gamma.append(", ");
        }
        x.append(system.getPhase(1).getComponent(i).getx());
        gamma.append(system.getPhase(1).getActivityCoefficient(i));
      }
      System.out.println("  " + String.join("+", names[k]));
      System.out.println("    x     [" + x + "]");
      System.out.println("    gamma [" + gamma + "]");
    }
  }

  /**
   * The UNIFAC activity coefficients `eos.unifac_activity_coefficients` would be checked
   * against, if NeqSim 3.20.0 could produce them.
   *
   * It cannot, and the two attempts below record why rather than asserting it. The
   * number in the bracket after each is the point of it: an oracle nobody can reach is
   * worth a measurement, not a claim.
   *
   * 1. Classic UNIFAC. `ComponentGEUnifac`'s constructor fills `unifacGroups` (an
   *    ArrayList) but never `unifacGroupsArray`, which only `addUNIFACgroup` and
   *    `setUnifacGroups` write. `getNumberOfUNIFACgroups()` reads the list and
   *    `getUnifacGroup(i)` reads the array, so `PhaseGEUnifac.checkGroups` - reached from
   *    `setMixingRule` - indexes a length-zero array and throws.
   * 2. UNIFAC-PSRK. `ComponentGEUnifac`'s constructor returns early for a subclass
   *    (`if (!this.getClass().equals(ComponentGEUnifac.class)) return;`), so a
   *    `ComponentGEUnifacPSRK` never reads its groups at all. Nothing throws - with zero
   *    groups the two loops in `checkGroups` never run - and the damage surfaces later
   *    and silently: `getR`/`getQ` sum over no groups, and every gamma comes back NaN
   *    while the flash reports a single phase. `ComponentGEUnifacUMRPRU` shares that
   *    early return, and `PhaseGEUnifacUMRPRU` is reached only through
   *    `SystemUMRPRUMCEosNew`.
   */
  static void unifac(String label, SystemInterface fluid, String[] names, double[] moles) {
    System.out.println("  " + label);
    try {
      for (int i = 0; i < names.length; i++) {
        fluid.addComponent(names[i], moles[i]);
      }
      fluid.createDatabase(true);
      fluid.setMixingRule("classic");
      new ThermodynamicOperations(fluid).TPflash();
      fluid.initProperties();

      System.out.println("    phases " + fluid.getNumberOfPhases() + "  phase 1 holds "
          + fluid.getPhase(1).getComponent(0).getClass().getSimpleName());

      StringBuilder x = new StringBuilder();
      StringBuilder gamma = new StringBuilder();
      for (int i = 0; i < names.length; i++) {
        if (i > 0) {
          x.append(", ");
          gamma.append(", ");
        }
        x.append(fluid.getPhase(1).getComponent(i).getx());
        gamma.append(fluid.getPhase(1).getActivityCoefficient(i));
      }
      System.out.println("    x     [" + x + "]");
      System.out.println("    gamma [" + gamma + "]");
    } catch (Throwable failure) {
      System.out.println("    THREW " + failure);
      System.out.println("      at " + failure.getStackTrace()[0]);
    }
  }

  /** The UNIFAC attempts, for `eos.unifac_activity_coefficients`. */
  static void unifac() {
    System.out.println("UNIFAC activity coefficients (the oracle that is not there):");
    String[] names = {"methanol", "water"};
    double[] moles = {0.5, 0.5};
    unifac("classic UNIFAC, methanol/water 0.5/0.5, 298.15 K", new SystemUNIFAC(298.15, 1.0), names, moles);
    unifac("UNIFAC-PSRK, methanol/water 0.5/0.5, 298.15 K", new SystemUNIFACpsrk(298.15, 1.0), names, moles);
  }

  /** The paraffin-wax Wilson activity coefficients, for `eos.wilson_activity_coefficients`. */
  static void wilson() {
    String[][] names = {{"n-butane", "nc12"}};
    double[][] moles = {{0.5, 0.5}};
    System.out.println("Wilson activity coefficients (gamma = exp(ln gamma)):");
    for (int k = 0; k < names.length; k++) {
      SystemGEWilson system = new SystemGEWilson(298.15, 1.0);
      for (int i = 0; i < names[k].length; i++) {
        system.addComponent(names[k][i], moles[k][i]);
      }
      system.createDatabase(true);
      system.setMixingRule("classic");
      new ThermodynamicOperations(system).TPflash();
      StringBuilder x = new StringBuilder();
      StringBuilder gamma = new StringBuilder();
      for (int i = 0; i < names[k].length; i++) {
        if (i > 0) {
          x.append(", ");
          gamma.append(", ");
        }
        x.append(system.getPhase(1).getComponent(i).getx());
        gamma.append(((neqsim.thermo.component.ComponentGEWilson) system.getPhase(1).getComponent(i))
            .getWilsonActivityCoefficient(system.getPhase(1)));
      }
      System.out.println("  " + String.join("+", names[k]));
      System.out.println("    x     [" + x + "]");
      System.out.println("    gamma [" + gamma + "]");
    }
  }

  /**
   * The Taleb-Ponche-Mirabel Van Laar activity coefficients, for
   * `eos.van_laar_acid_activity_coefficients`.
   *
   * Called on the static utility rather than through a phase: the three expressions
   * take the acid-basis mole fractions and the temperature directly, which is the same
   * arithmetic `ComponentGEVanLaarAcid.computeGamma` performs on a phase's composition.
   * Reaching them this way also reaches the one part of the model that a flash could
   * not - the ternary's own basis.
   */
  static void vanLaarAcid() {
    String[][] sets = {{"ternary", "0.5", "0.3", "0.2", "250.0"}};
    System.out.println("Van Laar acid activity coefficients (H2O/HNO3/H2SO4 basis):");
    for (String[] set : sets) {
      double x1 = Double.parseDouble(set[1]);
      double x2 = Double.parseDouble(set[2]);
      double x3 = Double.parseDouble(set[3]);
      double t = Double.parseDouble(set[4]);
      double water = NitricSulfuricAcidVaporPressure.activityCoefficientWater(x1, x2, x3, t);
      double nitric = NitricSulfuricAcidVaporPressure.activityCoefficientNitricAcid(x1, x2, x3, t);
      double sulfuric =
          NitricSulfuricAcidVaporPressure.activityCoefficientSulfuricAcid(x1, x2, x3, t);
      System.out.println("  " + set[0] + "  x [" + x1 + ", " + x2 + ", " + x3 + "]  T " + t);
      System.out.println("    lngamma [" + Math.log(water) + ", " + Math.log(nitric) + ", "
          + Math.log(sulfuric) + "]");
      System.out.println("    gamma   [" + water + ", " + nitric + ", " + sulfuric + "]");
    }
  }

  public static void main(String[] args) {
    volcorr();
    chung();
    wilke();
    chungCond();
    masonSaxena();
    corr();
    antoine();
    nrtl();
    unifac();
    wilson();
    vanLaarAcid();
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
