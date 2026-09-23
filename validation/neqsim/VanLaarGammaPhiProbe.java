// The direct gamma-phi flash, for `eos.ge_flash`.
//
// `TPflash.runInternal` dispatches to a **direct gamma-phi K-value loop** when the system's
// `EosGeFlashModel` answers `requiresDirectGammaPhiFlash()`, and `SystemVanLaarActivitySRK`
// is the one system in `src/main` that overrides it to `true`
// (`SystemVanLaarActivitySRK.java:386`, selected at `TPflash.java:226`). So this is the
// fluid that route is for, and it is the only one: the K-value the loop iterates is
// `phi_i^L / phi_i^V` with an EoS gas over a GE liquid, which is the same shape
// `eos.ge_nrtl_flash` has and the same shape a generalised one has.
//
// The fluid and the state are `SystemEosGEOperationsTest.createTwoPhaseSystem`'s - CO2 10 /
// water 0.70 / nitric acid 0.15 / sulfuric acid 0.15 at 273.15 K and 1 bar - which is where
// the two-phase state below comes from.
//
// The bubble and dew pressures are printed too, because that test's own assertions are on
// them: the bubble pressure is the **sum of the liquid model's three partial pressures** to
// 2%, and each gas mole fraction is `p_i / sum(p_i)` to 2e-3, which is what makes that route
// a bubble-point calculation and not only a flash.
//
//     javac -proc:none -cp neqsim-f0c7436.jar VanLaarGammaPhiProbe.java
//     java -cp .:neqsim-f0c7436.jar VanLaarGammaPhiProbe > captures/van_laar_gamma_phi_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemVanLaarActivitySRK;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.thermo.util.empiric.NitricSulfuricAcidVaporPressure;

public class VanLaarGammaPhiProbe {

  static final double PASCALS_PER_BAR = 1.0e5;

  public static void main(String[] args) throws Exception {
    twoPhase();
    bubbleAndDew();
  }

  /// The flash at the test's own state: the split, both phases, the K-values.
  static void twoPhase() {
    System.out.println("block=two-phase-flash");
    SystemInterface system = build(273.15, 1.0,
        new String[] { "CO2", "water", "nitric acid", "sulfuric acid" },
        new double[] { 10.0, 0.70, 0.15, 0.15 });
    ThermodynamicOperations flash = new ThermodynamicOperations(system);
    flash.TPflash();
    system.init(1);

    System.out.println("phases=" + system.getNumberOfPhases());
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      PhaseInterface phase = system.getPhase(p);
      System.out.println("  phase[" + p + "]_type=" + phase.getPhaseTypeName()
          + " beta=" + phase.getBeta() + " density=" + phase.getPhysicalProperties().getDensity()
          + " z=" + phase.getZ());
      for (int i = 0; i < phase.getNumberOfComponents(); i++) {
        ComponentInterface component = phase.getComponent(i);
        System.out.println("    x[" + p + "][" + component.getName() + "]=" + component.getx()
            + " K=" + component.getK() + " ln_phi=" + component.getLogFugacityCoefficient());
      }
    }
    System.out.println();
  }

  /// The bubble and dew pressures, and the gas composition at the bubble pressure.
  static void bubbleAndDew() throws Exception {
    System.out.println("block=bubble-and-dew");
    double temperature = 273.15;
    double[] liquid = NitricSulfuricAcidVaporPressure.moleFractionsFromMassFractions(60.0, 20.0,
        20.0);
    double[] partial = new double[] {
        NitricSulfuricAcidVaporPressure.partialPressureWater(liquid[0], liquid[1], liquid[2],
            temperature) / PASCALS_PER_BAR,
        NitricSulfuricAcidVaporPressure.partialPressureNitricAcid(liquid[0], liquid[1], liquid[2],
            temperature) / PASCALS_PER_BAR,
        NitricSulfuricAcidVaporPressure.partialPressureSulfuricAcid(liquid[0], liquid[1],
            liquid[2], temperature) / PASCALS_PER_BAR};
    System.out.println("mass_fractions=60 20 20");
    for (int i = 0; i < 3; i++) {
      System.out.println("  mole_fraction[" + i + "]=" + liquid[i]);
    }
    for (int i = 0; i < 3; i++) {
      System.out.println("  partial_pressure_bar[" + i + "]=" + partial[i]);
    }
    double expected = partial[0] + partial[1] + partial[2];
    System.out.println("expected_bubble_pressure_bar=" + expected);
    for (int i = 0; i < 3; i++) {
      System.out.println("  expected_gas_x[" + i + "]=" + partial[i] / expected);
    }

    SystemInterface bubble = build(temperature, 1.0,
        new String[] { "water", "nitric acid", "sulfuric acid" },
        new double[] { liquid[0], liquid[1], liquid[2] });
    new ThermodynamicOperations(bubble).bubblePointPressureFlash(false);
    System.out.println("bubble_pressure_bar=" + bubble.getPressure());
    for (int i = 0; i < bubble.getPhase(0).getNumberOfComponents(); i++) {
      System.out.println("  gas_x[" + bubble.getPhase(0).getComponent(i).getName() + "]="
          + bubble.getPhase(0).getComponent(i).getx());
    }

    double pureWater = NitricSulfuricAcidVaporPressure.pureVaporPressureWater(temperature)
        / PASCALS_PER_BAR;
    System.out.println("pure_water_pressure_bar=" + pureWater);
    SystemInterface dew = build(temperature, 1.0, new String[] { "water" }, new double[] { 1.0 });
    new ThermodynamicOperations(dew).dewPointPressureFlash();
    System.out.println("dew_pressure_bar=" + dew.getPressure());
    System.out.println();
  }

  /// `SystemVanLaarActivitySRK` with the acid components the caller asks for.
  static SystemInterface build(double temperature, double pressure, String[] names,
      double[] moles) {
    SystemInterface system = new SystemVanLaarActivitySRK(temperature, pressure);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.createDatabase(true);
    system.setMixingRule("classic");
    return system;
  }
}
