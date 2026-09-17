// The gamma-phi flash of an SRK vapour over a `PhaseGENRTL` liquid, for
// `eos.ge_nrtl_flash`.
//
//     javac -proc:none -cp neqsim-3.20.0.jar NrtlGammaPhiFlash.java
//     java -cp .:neqsim-3.20.0.jar NrtlGammaPhiFlash
//
// **Why this drives the loop itself.** `SystemNRTL` does not opt into NeqSim's direct
// gamma-phi flash: `EosGeFlashModel.requiresDirectGammaPhiFlash` defaults to false,
// `SystemEosGE` does not override it, and in 3.20.0 the only system that returns true
// is `SystemVanLaarActivitySRK`. Measured, `SystemNRTL`'s own `TPflash` therefore runs
// the ordinary EOS successive substitution and returns a single SRK phase - at 21
// states from 300 K to 350 K and z = 0.1/0.5/0.9, one phase every time and the GE
// liquid never activated. So a comparison against `SystemNRTL`'s flash would measure
// the wrong thing.
//
// What *is* NeqSim's, and is what this prints, is the model and the update rule:
// `PhaseSrkEos` over `PhaseGENRTL`, the K-value `phi_i^L / phi_i^V` from
// `TPflash.sucsSubsGammaPhi`, and NeqSim's own `RachfordRice`. The loop below is that
// method transcribed, with the two model hooks left at their defaults - `presdiff` is
// one because both phases are at the feed's pressure, `relaxGammaPhiKValue` is the
// identity, and `constrainGammaPhiKValue` does nothing.

import neqsim.thermo.system.SystemNRTL;
import neqsim.thermo.system.EosGeFlashModel;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.component.ComponentGEInterface;
import neqsim.thermo.phase.PhaseGEInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermodynamicoperations.flashops.RachfordRice;

public class NrtlGammaPhiFlash {

  /** How close the logarithmic K-values must come to stop. */
  private static final double TOLERANCE = 1.0e-14;

  private static final int MAX_ITERATIONS = 500;

  static void one(String first, String second, double temperatureK, double pressureBar, double zFirst)
      throws Exception {
    SystemNRTL system = new SystemNRTL(temperatureK, pressureBar);
    system.addComponent(first, zFirst);
    system.addComponent(second, 1.0 - zFirst);
    system.createDatabase(true);
    system.setMixingRule("classic");
    system.init(0);

    PhaseInterface vapour = system.getPhase(0);
    PhaseInterface liquid = system.getPhase(1);
    EosGeFlashModel model = (EosGeFlashModel) system;
    RachfordRice rachfordRice = new RachfordRice();

    int n = vapour.getNumberOfComponents();
    double[] k = new double[n];
    // Wilson K-values, the same opening estimate the ordinary flash uses.
    for (int i = 0; i < n; i++) {
      k[i] = Math.exp(Math.log(vapour.getComponent(i).getPC() / pressureBar)
          + 5.373 * (1.0 + vapour.getComponent(i).getAcentricFactor())
          * (1.0 - vapour.getComponent(i).getTC() / temperatureK));
      vapour.getComponent(i).setK(k[i]);
      liquid.getComponent(i).setK(k[i]);
    }

    double beta = 0.5;
    double residual = Double.NaN;
    int iteration = 0;
    for (iteration = 1; iteration <= MAX_ITERATIONS; iteration++) {
      // `ComponentGE.getGamma` is a cached field that only `getExcessGibbsEnergy`
      // fills, so the liquid's coefficients are only meaningful after this call.
      // `GeGamma.java` documents the same requirement.
      ((PhaseGEInterface) liquid).getExcessGibbsEnergy(liquid, n, temperatureK, pressureBar,
          PhaseType.LIQUID);

      residual = 0.0;
      for (int i = 0; i < n; i++) {
        ComponentInterface inVapour = vapour.getComponent(i);
        ComponentInterface inLiquid = liquid.getComponent(i);
        double kOld = inVapour.getK();
        inLiquid.fugcoef(liquid);
        double phiVapour = model.getGammaPhiVapourFugacityCoefficient(inVapour, vapour);
        double targetK = inLiquid.getFugacityCoefficient() / phiVapour;
        double kNew = model.relaxGammaPhiKValue(kOld, targetK);
        inVapour.setK(kNew);
        inLiquid.setK(kNew);
        k[i] = kNew;
        residual += Math.pow(Math.log(kNew / kOld), 2.0);
      }
      residual = Math.sqrt(residual / n);

      beta = rachfordRice.calcBeta(k, system.getzvector());
      system.setBeta(beta);
      system.calc_x_y();
      system.init(1);

      if (residual <= TOLERANCE) {
        break;
      }
    }

    // One consistent converged state: the last K the loop set, with the compositions
    // it implies.
    for (int i = 0; i < n; i++) {
      vapour.getComponent(i).setK(k[i]);
      liquid.getComponent(i).setK(k[i]);
    }
    system.setBeta(beta);
    system.calc_x_y();
    system.init(1);
    ((PhaseGEInterface) liquid).getExcessGibbsEnergy(liquid, n, temperatureK, pressureBar,
        PhaseType.LIQUID);

    System.out.println(first + "/" + second + "  T=" + temperatureK + " P=" + pressureBar
        + "  z=" + zFirst);
    System.out.println("  iterations=" + iteration + " residual=" + residual + " beta=" + beta);
    for (int i = 0; i < n; i++) {
      ComponentInterface inVapour = vapour.getComponent(i);
      ComponentInterface inLiquid = liquid.getComponent(i);
      inLiquid.fugcoef(liquid);
      double phiVapour = model.getGammaPhiVapourFugacityCoefficient(inVapour, vapour);
      System.out.println("  " + inVapour.getName()
          + "  x=" + inLiquid.getx()
          + "  y=" + inVapour.getx()
          + "  K=" + inVapour.getK()
          + "  gamma=" + ((ComponentGEInterface) inLiquid).getGamma()
          + "  lnPhiLiquid=" + Math.log(inLiquid.getFugacityCoefficient())
          + "  lnPhiVapour=" + Math.log(phiVapour));
    }
    System.out.println("  zVapour=" + vapour.getZ());
  }

  public static void main(String[] args) throws Exception {
    one("methanol", "water", 350.0, 1.0, 0.5);
    one("methanol", "water", 353.0, 1.0, 0.4);
    one("ethanol", "water", 360.0, 1.0, 0.3);
  }
}
