// The fixed-role EoS-gas / EoS-oil / GE-aqueous flash, for `eos.hybrid_eos_ge_flash`.
//
// `SystemPitzer` extends `SystemEosGE`, and with `setMultiPhaseCheck(true)` its `TPflash`
// takes the hybrid route: slot 0 is an EoS gas, slot 1 a GE aqueous phase and slot 2 an EoS
// oil, and the flash moves components between an EoS phase and a GE one. The roles are
// **fixed**, not decided by a stability test, which is what makes the answer pinnable: the
// test's acceptance contract is structural.
//
// The fluid is `SystemHybridEosGeFlashTest.createGasOilAqueousSystem`'s - methane 5 /
// n-heptane 2 / water 55.5 / Na+ 1 / Cl- 1 at 313.15 K and 50 bar - which is the issue-2862
// synthetic hydrocarbon/brine system that route exists for.
//
//     javac -proc:none -cp neqsim-f0c7436.jar HybridEosGeProbe.java
//     java -cp .:neqsim-f0c7436.jar HybridEosGeProbe > captures/hybrid_eos_ge_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPitzer;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class HybridEosGeProbe {

  /** The test's own material-balance and log-fugacity tolerances, printed as residuals. */
  static final double MATERIAL_PHASE_FRACTION = 1.0e-10;

  public static void main(String[] args) {
    System.out.println("block=topology");
    SystemInterface system = new SystemPitzer(313.15, 50.0);
    system.addComponent("methane", 5.0);
    system.addComponent("n-heptane", 2.0);
    system.addComponent("water", 55.5);
    system.addComponent("Na+", 1.0);
    system.addComponent("Cl-", 1.0);
    system.setMixingRule("classic");
    system.setMultiPhaseCheck(true);

    new ThermodynamicOperations(system).TPflash();

    System.out.println("phases=" + system.getNumberOfPhases());
    System.out.println("eos_gas_slot=" + ((SystemPitzer) system).getEosGasPhaseSlot());
    System.out.println("eos_oil_slot=" + ((SystemPitzer) system).getEosOilPhaseSlot());
    System.out.println("ge_liquid_slot=" + ((SystemPitzer) system).getGeLiquidPhaseSlot());

    int components = system.getPhase(0).getNumberOfComponents();
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      PhaseInterface phase = system.getPhase(p);
      System.out.println("  phase[" + p + "]_type=" + phase.getPhaseTypeName()
          + " beta=" + system.getBeta(p) + " z=" + phase.getZ());
      for (int i = 0; i < components; i++) {
        ComponentInterface component = phase.getComponent(i);
        System.out.println("    x[" + p + "][" + component.getName() + "]=" + component.getx()
            + " phi=" + component.getFugacityCoefficient()
            + " ln_f=" + logFugacity(phase, component)
            + " K=" + component.getK());
      }
    }

    System.out.println();
    System.out.println("block=acceptance");
    double betaSum = 0.0;
    double worstBalance = 0.0;
    double worstLogFugacity = 0.0;
    double worstIonOutsideAqueous = 0.0;
    int aqueousIndex = -1;
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      if ("aqueous".equals(system.getPhase(p).getPhaseTypeName())) {
        aqueousIndex = p;
      }
    }
    System.out.println("aqueous_index=" + aqueousIndex);
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      betaSum += system.getBeta(p);
    }
    for (int i = 0; i < components; i++) {
      double feedFraction = 0.0;
      for (int p = 0; p < system.getNumberOfPhases(); p++) {
        feedFraction += system.getBeta(p) * system.getPhase(p).getComponent(i).getx();
      }
      double z = system.getPhase(0).getComponent(i).getz();
      worstBalance = Math.max(worstBalance, Math.abs(z - feedFraction));
    }
    for (int i = 0; i < components; i++) {
      ComponentInterface reference = system.getPhase(0).getComponent(i);
      // The contract's own exclusion: a zero feed fraction or an ion takes no part in the
      // equal-fugacity assertion, because an ion's coefficient is a model constant rather
      // than an equilibrium quantity.
      if (reference.getz() <= 1.0e-30 || reference.getIonicCharge() != 0.0
          || reference.isIsIon()) {
        continue;
      }
      double first = Double.NaN;
      for (int p = 0; p < system.getNumberOfPhases(); p++) {
        if (system.getBeta(p) <= MATERIAL_PHASE_FRACTION) {
          continue;
        }
        double value = logFugacity(system.getPhase(p), system.getPhase(p).getComponent(i));
        if (Double.isNaN(first)) {
          first = value;
        } else {
          worstLogFugacity = Math.max(worstLogFugacity, Math.abs(first - value));
        }
      }
    }
    for (int i = 0; i < components; i++) {
      ComponentInterface reference = system.getPhase(0).getComponent(i);
      if (reference.getIonicCharge() == 0.0 && !reference.isIsIon()) {
        continue;
      }
      for (int p = 0; p < system.getNumberOfPhases(); p++) {
        if (p == aqueousIndex) {
          continue;
        }
        worstIonOutsideAqueous = Math.max(worstIonOutsideAqueous,
            system.getPhase(p).getComponent(i).getx());
      }
    }
    System.out.println("beta_sum=" + betaSum);
    System.out.println("worst_material_balance=" + worstBalance);
    System.out.println("worst_log_fugacity=" + worstLogFugacity);
    System.out.println("worst_ion_fraction_outside_aqueous=" + worstIonOutsideAqueous);
    System.out.println();
  }

  /// `ln(x_i phi_i P)`, the quantity the test's equilibrium assertion compares across phases.
  static double logFugacity(PhaseInterface phase, ComponentInterface component) {
    double value = component.getx() * component.getFugacityCoefficient() * phase.getPressure();
    return value > 0.0 ? Math.log(value) : Double.NaN;
  }
}
