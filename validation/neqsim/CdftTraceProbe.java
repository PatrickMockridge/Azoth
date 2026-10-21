// The mixture branch's intermediates, read from the class's own private methods: the bulk
// component densities it builds, the equilibrium potentials its central difference returns,
// the bulk pressure, and the sigma at two fixed widths. The port has to reproduce these
// numbers, not a re-derivation of them.
//
//     javac -proc:none -cp neqsim-f0c7436.jar CdftTraceProbe.java
//     java -cp .:neqsim-f0c7436.jar CdftTraceProbe

import java.lang.reflect.Method;
import neqsim.physicalproperties.interfaceproperties.surfacetension.CDFTSurfaceTension;
import neqsim.thermo.component.ComponentEos;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class CdftTraceProbe {

  public static void main(String[] args) throws Exception {
    SystemInterface system = new SystemSrkEos(300.0, 20.0);
    system.addComponent("methane", 0.5);
    system.addComponent("n-butane", 0.5);
    system.createDatabase(true);
    system.setMixingRule(2);
    new ThermodynamicOperations(system).TPflash();
    system.initProperties();

    PhaseInterface liqPhase = system.getPhase(1);
    PhaseInterface vapPhase = system.getPhase(0);
    int nc = liqPhase.getNumberOfComponents();
    double temperature = system.getPhase(0).getTemperature();

    double[] aComp = new double[nc];
    double[] bComp = new double[nc];
    double[] rhoLiq = new double[nc];
    double[] rhoVap = new double[nc];
    double totalRhoLiq = 1.0 / (liqPhase.getMolarVolume() * 1.0e-5);
    double totalRhoVap = 1.0 / (vapPhase.getMolarVolume() * 1.0e-5);
    for (int i = 0; i < nc; i++) {
      aComp[i] = ((ComponentEos) liqPhase.getComponent(i)).getaT() * 1.0e-5;
      bComp[i] = ((ComponentEos) liqPhase.getComponent(i)).getb() * 1.0e-5;
      rhoLiq[i] = liqPhase.getComponent(i).getx() * totalRhoLiq;
      rhoVap[i] = vapPhase.getComponent(i).getx() * totalRhoVap;
    }
    double[][] aij = new double[nc][nc];
    for (int i = 0; i < nc; i++) {
      for (int j = 0; j < nc; j++) {
        aij[i][j] = Math.sqrt(aComp[i] * aComp[j]);
      }
    }
    double bAvg = 0.0;
    for (int i = 0; i < nc; i++) {
      bAvg += (rhoLiq[i] / totalRhoLiq) * bComp[i];
    }
    double dMol = Math.pow(6.0 * bAvg / (Math.PI * 6.023e23), 1.0 / 3.0);

    CDFTSurfaceTension cdft = new CDFTSurfaceTension(system);
    Method free = CDFTSurfaceTension.class.getDeclaredMethod("mixtureFreeEnergyDensity",
        double[].class, double[][].class, double[].class, double.class, double.class,
        double.class, int.class);
    free.setAccessible(true);
    Method press = CDFTSurfaceTension.class.getDeclaredMethod("mixturePressure",
        double[].class, double[][].class, double[].class, double.class, double.class,
        double.class, int.class);
    press.setAccessible(true);
    Method sigma = CDFTSurfaceTension.class.getDeclaredMethod("computeMixtureSigmaForDelta",
        double[].class, double[].class, double.class, double[][].class, double[].class,
        double.class, double.class, double.class, double.class, double[].class, double.class,
        int.class);
    sigma.setAccessible(true);
    double d1 = 1.0;
    double d2 = 0.0;

    System.out.println("totalRhoLiq=" + totalRhoLiq);
    System.out.println("totalRhoVap=" + totalRhoVap);
    System.out.println("rhoLiq=" + arr(rhoLiq));
    System.out.println("rhoVap=" + arr(rhoVap));
    System.out.println("aComp=" + arr(aComp));
    System.out.println("bComp=" + arr(bComp));
    System.out.println("dMol=" + dMol);

    double eps = totalRhoLiq * 1.0e-7;
    double[] muEq = new double[nc];
    for (int ic = 0; ic < nc; ic++) {
      double[] plus = rhoLiq.clone();
      double[] minus = rhoLiq.clone();
      plus[ic] += eps;
      minus[ic] -= eps;
      double fPlus = (Double) free.invoke(cdft, plus, aij, bComp, temperature, d1, d2, nc);
      double fMinus = (Double) free.invoke(cdft, minus, aij, bComp, temperature, d1, d2, nc);
      muEq[ic] = (fPlus - fMinus) / (2.0 * eps);
    }
    double pBulk = (Double) press.invoke(cdft, rhoLiq, aij, bComp, temperature, d1, d2, nc);
    System.out.println("eps=" + eps);
    System.out.println("muEq=" + arr(muEq));
    System.out.println("pBulk=" + pBulk);

    for (double width : new double[] { 0.3 * dMol, 0.5 * dMol, 1.0 * dMol, 2.0 * dMol }) {
      double value = (Double) sigma.invoke(cdft, rhoLiq, rhoVap, temperature, aij, bComp, d1,
          d2, dMol, width, muEq, pBulk, nc);
      System.out.println("sigma_at_width_" + width / dMol + "_dMol=" + value);
    }
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
