import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermo.system.SystemPrEos;

/** What the cubic's kij column does: SRK reads KIJSRK, only PhasePrEos reads KIJPR. */
public final class SrkKijProbe {
  private SrkKijProbe() {}

  private static void show(String label, SystemInterface s) {
    s.setMixingRule("classic");
    s.init(0);
    s.init(1);
    System.out.printf("   (%s leaves %d phase(s))%n", label, s.getNumberOfPhases());
    System.out.printf("%-10s phase=%s Z=%.15g lnPhi=[", label,
        s.getPhase(0).getClass().getSimpleName(), s.getPhase(0).getZ());
    for (int i = 0; i < s.getNumberOfComponents(); i++) {
      System.out.printf("%.15g%s", Math.log(s.getPhase(0).getComponent(i).getFugacityCoefficient()),
          i + 1 < s.getNumberOfComponents() ? ", " : "");
    }
    System.out.println("]");
  }

  private static void pair(String[] names, double[] z, double t, double pBar, String note) {
    System.out.printf("%n# %s at %g K, %g bara: %s%n", String.join("/", names), t, pBar, note);
    SystemInterface srk = new SystemSrkEos(t, pBar);
    SystemInterface pr = new SystemPrEos(t, pBar);
    for (int i = 0; i < names.length; i++) {
      srk.addComponent(names[i], z[i]);
      pr.addComponent(names[i], z[i]);
    }
    show("SRK", srk);
    show("PR", pr);
  }

  public static void main(String[] args) {
    // Pairs whose `KIJSRK` and `KIJPR` differ, both light enough to stay one phase.
    pair(new String[] {"methane", "h2s"}, new double[] {0.7, 0.3}, 300.0, 20.0,
        "kijsrk = 0.085, kijpr = 0.08");
    pair(new String[] {"propane", "co2"}, new double[] {0.5, 0.5}, 350.0, 30.0,
        "kijsrk = 0.1018, kijpr = 0.135");
  }
}
