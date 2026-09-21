import neqsim.thermo.phase.PhasePCSAFT;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPCSAFT;

/**
 * PC-SAFT's volume solve at one state, and the same state reached from another.
 *
 * <p>
 * `PhasePCSAFTRahmat.calcVolume` seeds its Newton loop from `cachedMolarVolume`, which is
 * set from the last converged volume and is never cleared when the temperature or the
 * pressure changes. So the same state can be answered two ways, and this prints both:
 * once on an object that has only ever seen this state, and once on an object that walked
 * to it from a neighbouring one.
 *
 * <pre>
 * javac -proc:none -cp neqsim-f0c7436.jar PcsaftVolumeProbe.java
 * java -cp .:neqsim-f0c7436.jar PcsaftVolumeProbe [T_K P_bara name:z ...]
 * </pre>
 */
public final class PcsaftVolumeProbe {

  private PcsaftVolumeProbe() {}

  private static void show(String label, SystemInterface system) {
    System.out.printf("%s: numberOfPhases=%d%n", label, system.getNumberOfPhases());
    for (int i = 0; i < system.getNumberOfPhases(); i++) {
      PhasePCSAFT phase = (PhasePCSAFT) system.getPhase(i);
      System.out.printf(
          "   phase[%d] type=%s beta=%.15g n=%.15g volumeSAFT=%.15g nSAFT=%.15g Z=%.15g phi=%.15g%n",
          i, phase.getType(), phase.getBeta(), phase.getNumberOfMolesInPhase(),
          phase.getVolumeSAFT(), phase.getNSAFT(), phase.getZ(),
          phase.getComponent(0).getFugacityCoefficient());
    }
  }

  private static SystemInterface build(double t, double pBar, String[] names, double[] z) {
    SystemInterface system = new SystemPCSAFT(t, pBar);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], z[i]);
    }
    system.setMixingRule("classic");
    return system;
  }

  private static void probe(double t, double pBar, String[] names, double[] z) {
    System.out.printf("%n# %s at T=%g K, P=%g bara%n", String.join("/", names), t, pBar);

    SystemInterface fresh = build(t, pBar, names, z);
    fresh.init(0);
    fresh.init(1);
    show("fresh", fresh);

    // The same state, on an object that reached it from 200% of the pressure.
    SystemInterface walked = build(t, 2.0 * pBar, names, z);
    walked.init(0);
    walked.init(1);
    show("first at twice the pressure", walked);
    walked.setPressure(pBar);
    walked.init(2);
    show("then set back to this pressure", walked);
  }

  public static void main(String[] args) {
    System.out.println("# azoth PcsaftVolumeProbe - NeqSim master's SystemPCSAFT volume solve.");
    if (args.length >= 5) {
      double t = Double.parseDouble(args[0]);
      double p = Double.parseDouble(args[1]);
      String[] names = new String[(args.length - 2) / 2];
      double[] z = new double[names.length];
      for (int i = 0; i < names.length; i++) {
        names[i] = args[2 + 2 * i];
        z[i] = Double.parseDouble(args[3 + 2 * i]);
      }
      probe(t, p, names, z);
      return;
    }
    probe(300.0, 50.0, new String[] {"n-butane"}, new double[] {1.0});
  }
}
