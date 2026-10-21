// A trace of `GTSurfaceTensionFullGT`'s own sequence: the class's private work system and
// bulk densities, then `runcase`'s outer loop replayed with its public statics so each
// refinement's Newton answer, its `sigmaCalc`, its width and its new grid are visible.
//
//     javac -proc:none -cp neqsim-f0c7436.jar FullGtTraceProbe.java
//     java -cp .:neqsim-f0c7436.jar FullGtTraceProbe > captures/full_gt_trace.tsv

import java.lang.reflect.Field;
import neqsim.physicalproperties.interfaceproperties.surfacetension.GTSurfaceTensionFullGT;
import neqsim.physicalproperties.interfaceproperties.surfacetension.GTSurfaceTensionUtils;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class FullGtTraceProbe {

  static int ncomp;
  static double t;
  static double[] mueq;

  public static void main(String[] args) throws Exception {
    String first = args.length > 0 ? args[0] : "methane";
    String second = args.length > 1 ? args[1] : "n-butane";
    double z1 = args.length > 2 ? Double.parseDouble(args[2]) : 0.5;
    double z2 = args.length > 3 ? Double.parseDouble(args[3]) : 0.5;
    double temperature = args.length > 4 ? Double.parseDouble(args[4]) : 300.0;
    double pressure = args.length > 5 ? Double.parseDouble(args[5]) : 20.0;
    SystemInterface system = new SystemSrkEos(temperature, pressure);
    system.addComponent(first, z1);
    system.addComponent(second, z2);
    system.createDatabase(true);
    system.setMixingRule(2);
    new ThermodynamicOperations(system).TPflash();
    system.initPhysicalProperties();

    GTSurfaceTensionFullGT full = new GTSurfaceTensionFullGT(system, 0, 1);
    SystemInterface sys = (SystemInterface) get(full, "sys");
    ncomp = (Integer) get(full, "ncomp");
    t = (Double) get(full, "t");
    double[] ci = (double[]) get(full, "ci");
    mueq = (double[]) get(full, "mueq");
    double[] rho_ph1 = (double[]) get(full, "rho_ph1");
    double[] rho_ph2 = (double[]) get(full, "rho_ph2");

    System.out.println("t=" + t);
    System.out.println("ncomp=" + ncomp);
    System.out.println("rho_ph1=" + arr(rho_ph1));
    System.out.println("rho_ph2=" + arr(rho_ph2));
    System.out.println("ci=" + arr(ci));
    System.out.println("mueq=" + arr(mueq));

    double[] mu = new double[ncomp];
    double[][] dmu = new double[ncomp][ncomp];
    double[] p0 = new double[1];
    GTSurfaceTensionUtils.mufun(sys, ncomp, t, rho_ph1, mu, dmu, p0);
    System.out.println("mufun_rho_ph1_mu=" + arr(mu));
    System.out.println("mufun_rho_ph1_p=" + p0[0]);
    for (int i = 0; i < ncomp; i++) {
      System.out.println("mufun_rho_ph1_dmu_" + i + "=" + arr(dmu[i]));
    }
    GTSurfaceTensionUtils.mufun(sys, ncomp, t, rho_ph2, mu, dmu, p0);
    System.out.println("mufun_rho_ph2_mu=" + arr(mu));
    System.out.println("mufun_rho_ph2_p=" + p0[0]);
    for (int i = 0; i < ncomp; i++) {
      System.out.println("mufun_rho_ph2_dmu_" + i + "=" + arr(dmu[i]));
    }

    double[][] cij = new double[ncomp][ncomp];
    for (int i = 0; i < ncomp; i++) {
      for (int j = i; j < ncomp; j++) {
        cij[i][j] = 1.e18 * Math.sqrt(ci[i] * ci[j]);
        cij[j][i] = cij[i][j];
      }
    }
    System.out.println("cij_00=" + cij[0][0] + " cij_01=" + cij[0][1] + " cij_11=" + cij[1][1]);

    GTSurfaceTensionFullGT.initmu(sys, ncomp, t, rho_ph1, rho_ph2, mueq, p0, 1e-6);
    System.out.println("mueq_filled=" + arr(mueq));

    int Nlevel = 3;
    int Ngrid = (1 << Nlevel) + 1;
    double[][] rhomat = new double[Ngrid][ncomp];
    double[] xgrid = GTSurfaceTensionFullGT.linspace(-3.0, 3.0, Ngrid);
    double L = 3.0;
    int N2 = 1 << (Nlevel - 1);
    for (int i = 0; i < N2; i++) {
      rhomat[i] = rho_ph1.clone();
    }
    int N1 = N2;
    N2 = 1 + (1 << Nlevel);
    for (int i = N1; i < N2; i++) {
      rhomat[i] = rho_ph2.clone();
    }

    int Nrefinements = 7;
    double StopTolerance = 0.01;
    double sigma_old = 9.9e99;
    double sigma = 0.0;
    double std;

    for (int i = 0; i < Nrefinements; i++) {
      System.out.println("refine=" + i + " Ngrid=" + Ngrid + " L=" + L);
      double sigma_newton =
          traceNewton(cij, L, 20, 0.5, true, rhomat, sys, t, mueq, i);
      double[][] drhodz = new double[rhomat.length - 1][ncomp];
      sigma = GTSurfaceTensionFullGT.sigmaCalc(xgrid[1] - xgrid[0], rhomat, cij, false,
          drhodz, ncomp);
      std = full_calc_std(xgrid, cij, drhodz, rho_ph1, rho_ph2);
      System.out.println("refine=" + i + " sigma_newton=" + sigma_newton + " sigma=" + sigma
          + " std=" + std);

      if (i > 0 && Math.abs(sigma - sigma_old) < StopTolerance * sigma) {
        System.out.println("stopped_at=" + i);
        break;
      } else {
        sigma_old = sigma;
      }

      double H = 2.0 * L / (Ngrid - 1.0);
      int Nhalf = (int) Math.round(Math.ceil(6.0 * std / H));
      int Nideal = 2 * Nhalf + 1;
      double Lnew = Nhalf * H;
      int NgridNew = 2 * Nideal - 1;
      double[] xgridNew = GTSurfaceTensionFullGT.linspace(-Lnew, Lnew, NgridNew);
      System.out.println("refine=" + i + " H=" + H + " Nhalf=" + Nhalf + " Lnew=" + Lnew
          + " NgridNew=" + NgridNew);

      int kk = 0;
      double[][] rhotmp = new double[NgridNew][ncomp];
      double alpha;
      for (int j = 0; j < rhotmp.length; j++) {
        if (xgridNew[j] < xgrid[0]) {
          rhotmp[j] = rho_ph1.clone();
        } else if (xgridNew[j] > xgrid[xgrid.length - 1]) {
          rhotmp[j] = rho_ph2.clone();
        } else {
          while (xgridNew[j] > xgrid[kk + 1]) {
            kk++;
          }
          alpha = (xgrid[kk + 1] - xgridNew[j]) / (xgrid[kk + 1] - xgrid[kk]);
          for (int k = 0; k < ncomp; k++) {
            rhotmp[j][k] = alpha * rhomat[kk][k] + (1 - alpha) * rhomat[kk + 1][k];
          }
        }
      }
      rhotmp[0] = rho_ph1.clone();
      rhotmp[rhotmp.length - 1] = rho_ph2.clone();
      L = Lnew;
      rhomat = rhotmp;
      Ngrid = NgridNew;
      xgrid = xgridNew;
    }
    System.out.println("sigma_final=" + sigma);
  }

  static double traceNewton(double[][] cij, double L, int N_Newton, double allowedRelChange,
      boolean highOrder, double[][] rhomat, SystemInterface sys, double t, double[] mueq,
      int refinement) {
    int Ngrid = rhomat.length;
    double H = 2.0 * L / (Ngrid - 1);
    double[][][] Jac = new double[Ngrid][ncomp][ncomp];
    double[][] dmu = new double[Ngrid][ncomp];
    double[][] rres = new double[Ngrid][ncomp];
    double[][] rrho_prev = new double[Ngrid][ncomp];
    double[][] drhodz = new double[Ngrid - 1][ncomp];
    double sigma = 0.0;
    double[] xgrid = GTSurfaceTensionFullGT.linspace(-L, L, Ngrid);

    for (int i = 0; i < Ngrid; i++) {
      GTSurfaceTensionFullGT.delta_mu(sys, ncomp, t, mueq, rhomat[i], dmu[i], Jac[i]);
      for (int j = 0; j < ncomp; j++) {
        rres[i][j] = dmu[i][j];
        for (int k = 0; k < ncomp; k++) {
          rres[i][j] -= Jac[i][j][k] * rhomat[i][k];
        }
      }
    }
    addCorrection(rres, dmu, Ngrid);

    for (int step = 0; step < N_Newton; step++) {
      rrho_prev = new double[Ngrid][ncomp];
      for (int i = 0; i < Ngrid; i++) {
        for (int j = 0; j < ncomp; j++) {
          rrho_prev[i][j] = rhomat[i][j];
        }
      }

      GTSurfaceTensionFullGT.directsolve(rres, Jac, cij, H, Ngrid, rhomat, ncomp);

      double maxrelchange = 0.0;
      for (int i = 1; i < Ngrid - 1; i++) {
        for (int j = 0; j < ncomp; j++) {
          if (Math.abs(rhomat[i][j] - rrho_prev[i][j]) > maxrelchange * Math.abs(rrho_prev[i][j])) {
            maxrelchange = Math.abs(rhomat[i][j] - rrho_prev[i][j]) / Math.abs(rrho_prev[i][j]);
          }
        }
      }
      double urel_Newton = Math.min(allowedRelChange / maxrelchange, 1.0);
      for (int i = 1; i < Ngrid - 1; i++) {
        for (int j = 0; j < ncomp; j++) {
          rhomat[i][j] = urel_Newton * rhomat[i][j] + (1.0 - urel_Newton) * rrho_prev[i][j];
        }
      }
      if (urel_Newton * maxrelchange < 0.00001) {
        System.out.println("refine=" + refinement + " newton_step=" + step
            + " broke maxrelchange=" + maxrelchange + " urel=" + urel_Newton
            + " sigma_returned=" + sigma);
        break;
      }
      sigma = GTSurfaceTensionFullGT.sigmaCalc(H, rhomat, cij, false, drhodz, ncomp);
      System.out.println("refine=" + refinement + " newton_step=" + step
          + " maxrelchange=" + maxrelchange + " urel=" + urel_Newton + " sigma=" + sigma
          + " rho_mid=" + arr(rhomat[Ngrid / 2]));

      double sum_ztmp = 0.0;
      double sum_tmp = 0.0;
      for (int i = 0; i < drhodz.length; i++) {
        double tmp = 0.0;
        for (int j = 0; j < ncomp; j++) {
          for (int k = 0; k < ncomp; k++) {
            tmp += drhodz[i][j] * cij[j][k] * drhodz[i][k];
          }
        }
        sum_ztmp += 0.5 * (xgrid[i + 1] + xgrid[i]) * tmp;
        sum_tmp += tmp;
      }
      double cg = sum_ztmp / sum_tmp;
      int icorr = (int) Math.round(cg / H);
      if (icorr > 0) {
        for (int i = 1; i < Ngrid - icorr; i++) {
          rhomat[i] = rhomat[i + icorr];
        }
        for (int i = Ngrid - icorr; i < Ngrid - 1; i++) {
          rhomat[i] = rhomat[Ngrid - 1].clone();
        }
      } else if (icorr < 0) {
        icorr = -icorr;
        for (int i = Ngrid - 2; i >= icorr; i--) {
          rhomat[i] = rhomat[i - icorr];
        }
        for (int i = 1; i < icorr; i++) {
          rhomat[i] = rhomat[0].clone();
        }
      }
      if (icorr != 0) {
        System.out.println("refine=" + refinement + " newton_step=" + step + " cg_shift=" + icorr);
      }

      for (int i = 1; i < Ngrid - 1; i++) {
        GTSurfaceTensionFullGT.delta_mu(sys, ncomp, t, mueq, rhomat[i], dmu[i], Jac[i]);
        for (int j = 0; j < ncomp; j++) {
          rres[i][j] = dmu[i][j];
          for (int k = 0; k < ncomp; k++) {
            rres[i][j] -= Jac[i][j][k] * rhomat[i][k];
          }
        }
      }
      addCorrection(rres, dmu, Ngrid);
    }
    return sigma;
  }

  static void addCorrection(double[][] rres, double[][] dmu, int Ngrid) {
    for (int i = 1; i < Ngrid - 1; i++) {
      for (int j = 0; j < ncomp; j++) {
        rres[i][j] += (dmu[i - 1][j] - 2.0 * dmu[i][j] + dmu[i + 1][j]) / 12.0;
      }
    }
  }

  static double full_calc_std(double[] z, double[][] C, double[][] drhodz, double[] a,
      double[] b) {
    double h1 = z[1] - z[0];
    int Ngrid = drhodz.length + 1;
    double zdum = 0.0;
    double z2dum = 0.0;
    double sumdum = 0.0;
    double dum;
    for (int j = 0; j < Ngrid - 1; j++) {
      dum = 0.0;
      for (int m = 0; m < ncomp; m++) {
        for (int n = 0; n < ncomp; n++) {
          dum = drhodz[j][m] * C[m][n] * drhodz[j][n];
        }
      }
      sumdum += dum * h1;
      zdum += dum * (z[j + 1] * z[j + 1] - z[j] * z[j]) / 2.0;
      z2dum += dum * (z[j + 1] * z[j + 1] * z[j + 1] - z[j] * z[j] * z[j]) / 3.0;
    }
    double mean = zdum / sumdum;
    return Math.sqrt(z2dum / sumdum - mean * mean);
  }

  static Object get(Object target, String name) throws Exception {
    Field field = target.getClass().getDeclaredField(name);
    field.setAccessible(true);
    return field.get(target);
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
