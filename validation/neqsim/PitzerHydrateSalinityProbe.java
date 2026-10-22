// PitzerHydrateSalinityProbe - the salt suppression curve of the Pitzer hydrate saturation.
//
// The same fluid construction as PitzerHydrateProbe: a one-kilogram water basis, excess CO2,
// one chloride salt at a weight percent of the solution, the classic mixing rule, multiphase
// on, and the PHREEQC CO2-chloride parameters with the upstream test's zero zeta map. The
// upstream `PitzerHydrateSalinityValidationTest` walks a CSV of experimental points; this
// probe walks the concentration grid itself, so the suppression is a curve in the capture.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar PitzerHydrateSalinityProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar PitzerHydrateSalinityProbe > captures/pitzer_hydrate_salinity_probe.tsv

import java.util.LinkedHashMap;
import java.util.Map;
import neqsim.thermo.phase.PhaseHydrate;
import neqsim.thermo.system.SystemPitzer;
import neqsim.thermodynamicoperations.flashops.saturationops.PitzerHydrateFlash;

public class PitzerHydrateSalinityProbe {

  private static final double WATER_MOLES = 1.0 / 0.01801528;
  private static final String[] CATIONS = {"Na+", "K+", "Ca++", "Mg++"};
  private static final double[] SALT_MOLAR_MASS = {0.05844, 0.0745513, 0.11098, 0.095211};

  static Map<String, Double> zeroZeta() {
    Map<String, Double> values = new LinkedHashMap<String, Double>();
    for (String cation : CATIONS) {
      values.put(cation, 0.0);
    }
    return values;
  }

  static SystemPitzer brine(int salt, double weightPercent, double pressure) {
    SystemPitzer fluid = new SystemPitzer(280.0, pressure);
    fluid.addComponent("CO2", 10.0);
    fluid.addComponent("water", WATER_MOLES);
    double waterPart = 100.0 - weightPercent;
    double chloride = 0.0;
    for (int i = 0; i < CATIONS.length; i++) {
      double moles = (i == salt ? weightPercent : 0.0) / waterPart / SALT_MOLAR_MASS[i];
      fluid.addComponent(CATIONS[i], moles);
      chloride += moles * (i < 2 ? 1.0 : 2.0);
    }
    fluid.addComponent("Cl-", chloride);
    fluid.setMixingRule("classic");
    fluid.setMultiPhaseCheck(true);
    if (weightPercent > 0.0) {
      fluid.applyPhreeqcCo2ChlorideParameters(zeroZeta(), 0.0);
    }
    return fluid;
  }

  static void row(String label, int salt, double weightPercent, double pressure) {
    try {
      SystemPitzer fluid = brine(salt, weightPercent, pressure);
      PitzerHydrateFlash flash = new PitzerHydrateFlash(fluid, false);
      flash.run();
      System.out.println("  " + label + " status=solved temperature_k=" + fluid.getTemperature()
          + " structure=" + ((PhaseHydrate) fluid.getPhases()[4]).getStableHydrateStructure()
          + " residual=" + flash.getResidual() + " water_activity=" + flash.getWaterActivity());
    } catch (RuntimeException failure) {
      System.out.println("  " + label + " status=refused reason=" + failure.getClass().getSimpleName()
          + ":" + failure.getMessage());
    }
  }

  public static void main(String[] args) {
    System.out.println("pitzer_hydrate_salinity");
    double[] grid = {0.0, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0};
    System.out.println(" nacl_30bara");
    for (double percent : grid) {
      row("nacl_" + percent + "wt", 0, percent, 30.0);
    }
    System.out.println(" nacl_60bara");
    for (double percent : grid) {
      row("nacl_" + percent + "wt", 0, percent, 60.0);
    }
    System.out.println(" mgcl2_30bara");
    for (double percent : grid) {
      row("mgcl2_" + percent + "wt", 3, percent, 30.0);
    }
    System.out.println(" cacl2_30bara");
    for (double percent : new double[] {1.0, 5.0, 10.0, 15.0}) {
      row("cacl2_" + percent + "wt", 2, percent, 30.0);
    }
  }
}
