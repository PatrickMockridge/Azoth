// **Does the ion rows' filler data decide the answer?**
//
// The two dispatches S6c is about - `TPflash`'s sequential chemistry and
// `chemicalequilibrium/ChemicalEquilibrium` - fire only on a chemical system, and every fluid
// `chemicalReactionInit` makes chemical is one it has added ions to. On a `SystemSrkEos` fluid
// those ions are cubic components like any other, so every activity coefficient the chemistry
// reads is a *difference of log fugacity coefficients from a cubic built over them* - and the
// component table's critical columns for an ion are filler (sodium at `Tc = 717.3 K`,
// `Pc = 290.89 bar`, `omega = 0.344`, for a species with no critical point).
//
// The question this probe answers is whether that filler is *load-bearing*: perturb a reaction
// species' critical constants and see whether the converged chemistry moves. Ten percent is
// the perturbation, and three quantities are printed per run - the brine's species amounts, the
// `HCO3-` reaction heat, and the gas phase's composition - so that a quantity which is inert
// and one which is not are told apart rather than averaged.
//
//     javac -proc:none -cp neqsim-f0c7436.jar IonFillerSensitivityProbe.java
//     java -cp .:neqsim-f0c7436.jar IonFillerSensitivityProbe > captures/ion_filler_sensitivity_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class IonFillerSensitivityProbe {

  public static void main(String[] args) {
    one("unperturbed", 1.0);
    one("tc-x1.1", 1.1);
    one("tc-x0.9", 0.9);
    one("tc-x2.0", 2.0);
  }

  static void perturb(ComponentInterface component, double factor) {
    if (component.getIonicCharge() != 0 || component.isIsIon()) {
      component.setTC(component.getTC() * factor);
      component.setPC(component.getPC() * factor);
      component.setAcentricFactor(component.getAcentricFactor() * factor);
    }
  }

  static double first_ion_tc_of(neqsim.thermo.phase.PhaseInterface phase) {
    for (int i = 0; i < phase.getNumberOfComponents(); i++) {
      ComponentInterface component = phase.getComponent(i);
      if (component.getIonicCharge() != 0 || component.isIsIon()) {
        return component.getTC();
      }
    }
    return Double.NaN;
  }

  static double first_ion_tc(SystemInterface system) {
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      ComponentInterface component = system.getPhase(0).getComponent(i);
      if (component.getIonicCharge() != 0 || component.isIsIon()) {
        return component.getTC();
      }
    }
    return Double.NaN;
  }

  static void one(String label, double factor) {
    SystemInterface system = new SystemSrkEos(298.15, 1.01325);
    system.addComponent("CO2", 0.01);
    system.addComponent("water", 10.0);
    system.chemicalReactionInit();
    system.setMixingRule("classic");
    system.setMultiPhaseCheck(true);

    // The perturbation, applied to every component the reaction machinery added rather than to
    // the feed: those are the rows whose critical constants are filler.
    // **Every component of every phase, and the system's own array** - a NeqSim phase holds its
    // own component objects, so perturbing phase 0's reaches the gas and not the brine the
    // chemistry acts on. That was the first version of this probe, and it measured nothing.
    if (factor != 1.0) {
      for (int i = 0; i < system.getNumberOfComponents(); i++) {
        perturb(system.getComponent(i), factor);
      }
      for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
        for (int i = 0; i < system.getPhase(phase).getNumberOfComponents(); i++) {
          perturb(system.getPhase(phase).getComponent(i), factor);
        }
      }
    }

    System.out.println("fluid=" + label);
    // **Whether the perturbation survives the initialisation**, which is the difference between
    // an inert parameter and an experiment that did nothing: `init(0)` is where a system
    // re-reads the database.
    System.out.println("perturbed_tc_before_init=" + first_ion_tc(system) + " factor=" + factor);
    system.init(0);
    // After `init(0)` the phases are rebuilt, so the perturbation is re-applied to them - this
    // is what makes the experiment about the parameters rather than about the order of calls.
    if (factor != 1.0) {
      for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
        for (int i = 0; i < system.getPhase(phase).getNumberOfComponents(); i++) {
          perturb(system.getPhase(phase).getComponent(i), factor);
        }
      }
    }
    system.init(1);
    System.out.println("perturbed_tc_after_init=" + first_ion_tc(system));
    // **Which phase's components were perturbed**, and what the chemistry actually reads: the
    // activity coefficient of each component in the brine. If those are identical across a
    // factor of two the ion rows' parameters are inert; if they move, the chemistry rests on
    // them.
    System.out.println("aqueous_perturbed_tc=" + first_ion_tc_of(system.getPhase(1)));
    StringBuilder names = new StringBuilder();
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      names.append(system.getPhase(0).getComponent(i).getName()).append(" ");
    }
    System.out.println("components=" + names.toString().trim());
    System.out.println("aqueous_phase_class=" + system.getPhase(1).getClass().getSimpleName());

    ThermodynamicOperations operations = new ThermodynamicOperations(system);
    operations.TPflash();
    operations.chemicalEquilibrium();

    System.out.println("phases=" + system.getNumberOfPhases());
    for (int phase = 0; phase < system.getNumberOfPhases(); phase++) {
      StringBuilder x = new StringBuilder("  phase_x[" + phase + "]=");
      for (int i = 0; i < system.getPhase(phase).getNumberOfComponents(); i++) {
        x.append(system.getPhase(phase).getComponent(i).getName())
            .append(":")
            .append(system.getPhase(phase).getComponent(i).getx())
            .append(" ");
      }
      System.out.println(x.toString().trim());
    }
    for (int i = 0; i < system.getPhase(1).getNumberOfComponents(); i++) {
      ComponentInterface component = system.getPhase(1).getComponent(i);
      System.out.println("  log_activity[" + component.getName() + "]="
          + system.getPhase(1).getLogActivityCoefficient(i, 1));
    }
    System.out.println("delta_reaction_heat="
        + system.getChemicalReactionOperations().getDeltaReactionHeat());
    System.out.println("hco3_moles_overall="
        + system.getPhase(0).getComponent("HCO3-").getNumberOfmoles());
    System.out.println();
  }
}
