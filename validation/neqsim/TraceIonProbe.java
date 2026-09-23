// The trace-ion short circuit, for `reactions.reactive_tp_flash`'s driver.
//
// `ReactiveMultiphaseTPflash.run` initialises a VLE split and then **accepts it unreactive**
// when every ionic species is a trace: `ionZTotal`, the sum of the ions' overall mole
// fractions, below `1e-8`. The comment says why - "the full RAND solve can spend thousands of
// iterations polishing near-zero ion amounts" - and what this probe measures is what the
// shortcut costs and saves: no solve at all, `totalIterations = 0`, and the VLE split as the
// answer.
//
// **Two fluids, and the second is the control.** Trace salt (`1e-12` mol) fires the circuit;
// a mole of salt does not, and the same fluid then runs the whole reactive stack. Both need
// `NR > 0`, or the class's earlier `NR = 0` short circuit returns first and the branch this
// probe is about is never reached - which is why the fluid is carbonated water rather than
// methane/water.
//
//     javac -proc:none -cp neqsim-f0c7436.jar TraceIonProbe.java
//     java -cp .:neqsim-f0c7436.jar TraceIonProbe > captures/trace_ion_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.flashops.reactiveflash.FormulaMatrix;
import neqsim.thermodynamicoperations.flashops.reactiveflash.ReactiveMultiphaseTPflash;

public class TraceIonProbe {

  public static void main(String[] args) {
    one("trace-salt", 1.0e-12);
    one("molar-salt", 1.0);
  }

  /// Carbonated water with a trace of sodium chloride: `NR = 1` from the CO2-water reaction,
  /// two phases at 300 K and 50 bar, and the ions below the shortcut's threshold.
  static SystemInterface carbonated(double salt) {
    SystemSrkEos system = new SystemSrkEos(300.0, 50.0);
    system.addComponent("CO2", 0.5);
    system.addComponent("water", 10.0);
    system.addComponent("Na+", salt);
    system.addComponent("Cl-", salt);
    system.chemicalReactionInit();
    system.setMixingRule("classic");
    system.setMultiPhaseCheck(true);
    return system;
  }

  static void one(String label, double salt) {
    SystemInterface system = carbonated(salt);
    System.out.println("fluid=" + label);
    system.init(0);
    system.init(1);
    System.out.println("is_chemical_system=" + system.isChemicalSystem());

    FormulaMatrix matrix = new FormulaMatrix(system);
    String[] names = matrix.getComponentNames();
    System.out.println("reactive_components=" + String.join(" ", names));
    System.out.println("has_ionic_species=" + matrix.hasIonicSpecies());
    System.out.println("independent_reactions=" + matrix.getNumberOfIndependentReactions());
    double ionZ = 0.0;
    for (int i = 0; i < names.length; i++) {
      if (matrix.isIon(i)) {
        ionZ += Math.abs(system.getPhase(0).getComponent(i).getz());
      }
    }
    System.out.println("ion_z_total=" + ionZ);
    System.out.println("is_trace_ion_case=" + (matrix.hasIonicSpecies() && ionZ < 1.0e-8));

    ReactiveMultiphaseTPflash flash = new ReactiveMultiphaseTPflash(system);
    flash.run();
    System.out.println("converged=" + flash.isConverged());
    System.out.println("total_iterations=" + flash.getTotalIterations());
    System.out.println("phases=" + system.getNumberOfPhases());
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      PhaseInterface held = system.getPhase(phase);
      StringBuilder x = new StringBuilder("  phase_x[" + phase + "]=");
      for (int i = 0; i < names.length; i++) {
        x.append(held.getComponent(i).getx()).append(i + 1 < names.length ? " " : "");
      }
      System.out.println(x);
      System.out.println("  phase_beta[" + phase + "]=" + system.getBeta(phase));
    }
    System.out.println("total_moles=" + system.getNumberOfMoles());
    System.out.println();
  }
}
