import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/** NeqSim's own `VHflashQfunc` on the states `eos.vh_flash` carries cases for. */
public class VhFlashProbe {
  private static void one(double t, double pbar) {
    SystemInterface base = new SystemPrEos(t, pbar);
    base.addComponent("methane", 0.6);
    base.addComponent("n-butane", 0.4);
    base.setMixingRule("classic");
    base.setAttractiveTerm(1);
    new ThermodynamicOperations(base).TPflash();
    base.initProperties();
    double v = base.getVolume();
    double h = base.getEnthalpy("J/mol");
    System.out.printf("T=%.1f P=%.1f: V=%.12e H=%.9f%n", t, pbar, v, h);

    SystemInterface f = new SystemPrEos(300.0, pbar * 0.5);
    f.addComponent("methane", 0.6);
    f.addComponent("n-butane", 0.4);
    f.setMixingRule("classic");
    f.setAttractiveTerm(1);
    long t0 = System.nanoTime();
    try {
      new ThermodynamicOperations(f).VHflash(v, h);
      System.out.printf("   -> P=%.6f bar T=%.9f beta=%.9f  (%.1f ms)%n", f.getPressure(),
          f.getTemperature(), f.getBeta(), (System.nanoTime() - t0) / 1e6);
    } catch (Exception e) {
      System.out.printf("   -> THREW %s: %s%n", e.getClass().getSimpleName(), e.getMessage());
    }
  }

  public static void main(String[] args) {
    one(400.0, 10.0);
    one(330.0, 25.0);
    one(300.0, 50.0);
  }
}
