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
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class FlashTp {
  /** One PT flash, printed in the fields azoth's `eos.pt_flash` reports. */
  static void flash(String label, double temperatureK, double pressureBar, String[] names,
      double[] moles) {
    SystemInterface fluid = new SystemPrEos(temperatureK, pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    // "classic" is the van der Waals one-fluid mixing rule with the database's
    // Peng-Robinson `kij` - NeqSim's default for this system, and what `INTER.csv`
    // carries a value for.
    fluid.setMixingRule("classic");

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

  public static void main(String[] args) {
    flash("methane/n-butane, 0.6/0.4, 330 K, 25 bar",
        330.0, 25.0, new String[] {"methane", "n-butane"}, new double[] {0.6, 0.4});
    flash("propane, 1.0, 300 K, 9 bar",
        300.0, 9.0, new String[] {"propane"}, new double[] {1.0});
  }
}
