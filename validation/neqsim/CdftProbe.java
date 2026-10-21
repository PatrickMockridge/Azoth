// `CDFTSurfaceTension` on the states the port holds. The class consumes a flashed system
// and reads only each phase's molar volume and composition, so this records the component
// densities it derives (mol/m3) beside the answer. Pure states come from the bubble-point
// flash the class's own calibration test uses - a pure fluid is two-phase only there - and
// mixtures from TPflash. The class catches every exception and falls back to Parachor
// silently, so the parachor answer is printed beside it: equal values mean the fallback fired.
//
//     javac -proc:none -cp neqsim-f0c7436.jar CdftProbe.java
//     java -cp .:neqsim-f0c7436.jar CdftProbe > captures/cdft_probe.tsv

import neqsim.physicalproperties.interfaceproperties.surfacetension.CDFTSurfaceTension;
import neqsim.physicalproperties.interfaceproperties.surfacetension.ParachorSurfaceTension;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class CdftProbe {

  public static void main(String[] args) {
    pure("methane_150K_pr", "methane", 150.0, "pr");
    pure("methane_150K_srk", "methane", 150.0, "srk");
    pure("nbutane_320K_pr", "n-butane", 320.0, "pr");
    pure("nbutane_320K_srk", "n-butane", 320.0, "srk");
    mixture("methane_nbutane_300K_20bara_srk", new String[] { "methane", "n-butane" },
        new double[] { 0.5, 0.5 }, 300.0, 20.0, "srk");
    mixture("methane_nheptane_313K_30bara_srk", new String[] { "methane", "n-heptane" },
        new double[] { 0.5, 0.5 }, 313.15, 30.0, "srk");
    mixture("co2_water_313K_50bara_srk", new String[] { "CO2", "water" },
        new double[] { 0.02, 0.98 }, 313.15, 50.0, "srk");
    mixture("nheptane_water_300K_1bara_srk", new String[] { "n-heptane", "water" },
        new double[] { 0.5, 0.5 }, 300.0, 1.0, "srk");
  }

  static void pure(String label, String name, double temperatureK, String eos) {
    SystemInterface system = system(eos, temperatureK, 1.0);
    system.addComponent(name, 1.0);
    system.setMixingRule("classic");
    system.setMultiPhaseCheck(true);
    ThermodynamicOperations operations = new ThermodynamicOperations(system);
    try {
      operations.bubblePointPressureFlash(false);
    } catch (Exception exception) {
      throw new AssertionError("bubble-point flash failed", exception);
    }
    system.initProperties();
    report(label, system);
  }

  static void mixture(String label, String[] names, double[] z, double temperatureK,
      double pressureBara, String eos) {
    SystemInterface system = system(eos, temperatureK, pressureBara);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.createDatabase(true);
    system.setMixingRule(2);
    new ThermodynamicOperations(system).TPflash();
    system.initProperties();
    report(label, system);
  }

  static SystemInterface system(String eos, double temperatureK, double pressureBara) {
    return "pr".equals(eos) ? new SystemPrEos(temperatureK, pressureBara)
        : new SystemSrkEos(temperatureK, pressureBara);
  }

  static void report(String label, SystemInterface system) {
    System.out.println(label);
    System.out.println("eos_class=" + system.getClass().getSimpleName());
    System.out.println("temperature=" + system.getTemperature());
    System.out.println("pressure=" + system.getPressure());
    System.out.println("phase_count=" + system.getNumberOfPhases());
    int nc = system.getPhase(0).getNumberOfComponents();
    if (system.getNumberOfPhases() >= 2) {
      for (int p = 0; p < system.getNumberOfPhases(); p++) {
        double molarVolume = system.getPhase(p).getMolarVolume() * 1.0e-5;
        double total = 1.0 / molarVolume;
        double[] x = new double[nc];
        double[] rho = new double[nc];
        for (int c = 0; c < nc; c++) {
          x[c] = system.getPhase(p).getComponent(c).getx();
          rho[c] = x[c] * total;
        }
        System.out.println("phase" + p + "_type=" + system.getPhase(p).getType());
        System.out.println("phase" + p + "_molar_volume=" + molarVolume);
        System.out.println("phase" + p + "_x=" + arr(x));
        System.out.println("phase" + p + "_rho=" + arr(rho));
      }
      for (int c = 0; c < nc; c++) {
        System.out.println("component" + c + "_name=" + system.getPhase(0).getComponent(c).getName());
        System.out.println("component" + c + "_tc=" + system.getPhase(0).getComponent(c).getTC());
        System.out.println("component" + c + "_pc=" + system.getPhase(0).getComponent(c).getPC());
        System.out.println("component" + c + "_omega="
            + system.getPhase(0).getComponent(c).getAcentricFactor());
      }
      CDFTSurfaceTension cdft = new CDFTSurfaceTension(system);
      System.out.println("cdft_sigma=" + cdft.calcSurfaceTension(0, 1));
      ParachorSurfaceTension parachor = new ParachorSurfaceTension(system);
      System.out.println("parachor_sigma=" + parachor.calcSurfaceTension(0, 1));
    } else {
      System.out.println("no_interface=true");
    }
    System.out.println();
  }

  static String arr(double[] values) {
    StringBuilder sb = new StringBuilder("[");
    for (int i = 0; i < values.length; i++) {
      if (i > 0) {
        sb.append(", ");
      }
      sb.append(values[i]);
    }
    return sb.append("]").toString();
  }
}
