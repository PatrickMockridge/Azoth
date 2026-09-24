import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.process.equipment.stream.Stream;

public class VolCorrProbe {
  public static void main(String[] args) {
    row("water_liquid", new String[] { "water" }, new double[] { 1.0 }, 300.0, 1.0);
    row("butane_liquid", new String[] { "n-butane" }, new double[] { 1.0 }, 300.0, 10.0);
    row("methane_butane_gas", new String[] { "methane", "n-butane" }, new double[] { 0.9, 0.1 }, 320.0, 30.0);
  }

  static void row(String label, String[] names, double[] z, double t, double p) {
    SystemInterface f = new SystemPrEos(t, p);
    for (int i = 0; i < names.length; i++) f.addComponent(names[i], z[i]);
    f.setMixingRule(2);
    f.setTotalFlowRate(1.0, "mol/sec");
    f.init(0);
    new neqsim.thermodynamicoperations.ThermodynamicOperations(f).TPflash();
    f.initProperties();
    double vm = f.getPhase(0).getMolarVolume();
    double corr = 0.0;
    StringBuilder per = new StringBuilder();
    for (int i = 0; i < f.getPhase(0).getNumberOfComponents(); i++) {
      double c = f.getPhase(0).getComponent(i).getVolumeCorrection();
      double cT = f.getPhase(0).getComponent(i).getVolumeCorrectionT();
      double x = f.getPhase(0).getComponent(i).getx();
      per.append(String.format(" %s: c=%.6e cT=%.6e",
          f.getPhase(0).getComponent(i).getName(), c, cT));
      corr += x * (c + cT * (t - 288.15));
    }
    System.out.println(label);
    System.out.println("  phase_type=" + f.getPhase(0).getType());
    System.out.println("  molar_volume=" + vm + "  correction_term=" + corr);
    System.out.println("  density_translated=" + 1.0 / (vm - corr) * f.getPhase(0).getMolarMass() * 1.0e5);
    System.out.println("  density_cubic=" + f.getPhase(0).getDensity());
    System.out.println("  density_pp=" + f.getPhase(0).getPhysicalProperties().getDensity());
    System.out.println("  per-component:" + per);
    System.out.println();
  }
}
