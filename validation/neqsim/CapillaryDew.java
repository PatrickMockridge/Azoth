// NeqSim's capillary (in-pore) dew point, against the bulk one.
//
//     javac -proc:none -cp neqsim-3.20.0.jar CapillaryDew.java
//     java -cp .:neqsim-3.20.0.jar CapillaryDew
//
// `ThermodynamicOperations.capillaryDewPointTemperatureFlash(r)` — and the two-argument form
// that also takes a contact angle — build a `CapillaryDewPointFlash` and run it. That is the
// public entry, and the class is the one NeqSim class in the saturation-op family carrying
// physics azoth does not have: a Kelvin shift on the dew point inside a pore.
//
// The bulk dew point is printed beside every capillary one because the two are the whole point
// of the correction: the Kelvin equation moves the dew point *upwards*, and how far is the
// only thing this class computes. `r` is printed in metres and in nanometres, because the
// interesting radii are nanometres and the constructor takes metres - the unit hazard this
// port has been bitten by three times.
//
// The surface tension is read by the class itself from the interphase properties, so it is
// printed here too: it is an input the caller never supplies, and a port that takes it as an
// argument needs to know what value the oracle used.

import java.util.Arrays;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class CapillaryDew {

  static SystemInterface build(String[] names, double[] moles, double temperatureK, double pressureBar) {
    SystemInterface fluid = new SystemPrEos(temperatureK, pressureBar);
    for (int i = 0; i < names.length; i++) {
      fluid.addComponent(names[i], moles[i]);
    }
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(1);
    return fluid;
  }

  /** One bulk/capillary pair at a radius, fresh fluid each time because the flash mutates it. */
  static void one(String label, String[] names, double[] moles, double temperatureK,
      double pressureBar, double poreRadiusM, double contactAngleRad) {
    // The bulk dew point, from the operation that has no curvature.
    SystemInterface bulk = build(names, moles, temperatureK, pressureBar);
    double bulkT = Double.NaN;
    try {
      new ThermodynamicOperations(bulk).dewPointTemperatureFlash();
      bulkT = bulk.getTemperature();
    } catch (Exception ex) {
      bulkT = Double.NaN;
    }
    double sigmaOfBulk = Double.NaN;
    try {
      sigmaOfBulk = bulk.getInterphaseProperties().getSurfaceTension(0, 1);
    } catch (Exception ex) {
      sigmaOfBulk = Double.NaN;
    }

    SystemInterface capillary = build(names, moles, temperatureK, pressureBar);
    String outcome;
    double capillaryT = Double.NaN;
    double sigma = Double.NaN;
    try {
      ThermodynamicOperations ops = new ThermodynamicOperations(capillary);
      ops.capillaryDewPointTemperatureFlash(poreRadiusM, contactAngleRad);
      capillaryT = capillary.getTemperature();
      try {
        sigma = capillary.getInterphaseProperties().getSurfaceTension(0, 1);
      } catch (Exception ex) {
        sigma = Double.NaN;
      }
      // The two phases' molar volumes, because the Kelvin exponent is `Vm_L dP/(R T)` and the
      // port takes that volume from the cubic's root as `z R T / P`. If the two disagree, the
      // volumes are where.
      double vm0 = Double.NaN;
      double vm1 = Double.NaN;
      try {
        vm0 = capillary.getPhase(0).getMolarVolume("m3/mol");
        vm1 = capillary.getPhase(1).getMolarVolume("m3/mol");
      } catch (Exception ex) {
        // left as NaN
      }
      outcome = String.format("T=%.10f K   Vm0=%.10g  Vm1=%.10g  zRT/P(0)=%.10g  zRT/P(1)=%.10g",
          capillaryT, vm0, vm1,
          capillary.getPhase(0).getZ() * 8.31446 * capillaryT / capillary.getPressure(),
          capillary.getPhase(1).getZ() * 8.31446 * capillaryT / capillary.getPressure());
    } catch (Exception ex) {
      outcome = ex.getClass().getSimpleName() + ": " + ex.getMessage();
    }

    System.out.printf("%n=== %s   T0=%.2f K  P0=%.2f bar%n", label, temperatureK, pressureBar);
    System.out.printf("bulk dew point        T=%.10f K   sigma=%.10g N/m%n", bulkT, sigmaOfBulk);
    System.out.printf("capillary r=%.3e m (%.2f nm)  theta=%.4f rad  -> %s%n", poreRadiusM,
        poreRadiusM * 1.0e9, contactAngleRad, outcome);
    if (Double.isFinite(capillaryT) && Double.isFinite(bulkT)) {
      System.out.printf("     shift = %.10f K   sigma=%s N/m%n", capillaryT - bulkT,
          Double.isFinite(sigma) ? String.format("%.10g", sigma) : "n/a");
    }
  }

  public static void main(String[] args) {
    String[] names = {"methane", "n-butane"};
    double[] moles = {0.5, 0.5};
    double t0 = 300.0;
    double p0 = 20.0;
    System.out.printf("mixture %s  z=%s%n", Arrays.toString(names), Arrays.toString(moles));
    for (double radius : new double[] {1.0e-8, 1.0e-7, 1.0e-6}) {
      one("capillary dew, perfectly wetting", names, moles, t0, p0, radius, 0.0);
    }
    // A non-zero contact angle, which the two-argument entry takes.
    one("capillary dew, theta = 60 deg", names, moles, t0, p0, 1.0e-7, Math.PI / 3.0);
  }
}
