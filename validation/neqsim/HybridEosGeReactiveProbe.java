// The *reactive* hybrid flash, for `eos.hybrid_eos_ge_flash`'s reactive coupling.
//
// `TPHybridEosGeFlash.run` wraps the fixed-topology solve in a second loop when
// `system.isChemicalSystem()`: each pass re-equilibrates the aqueous phase with
// `ChemicalReactionOperations.solveChemEq`, projects the resulting mole changes onto the
// stoichiometric conservation null space, writes the adjusted species inventory back into all
// three roles, and only then solves the fractions again. The convergence test is then a pair -
// the composition deviation and the fraction residual - and the loop needs at least three
// passes.
//
// The fluid is `SystemHybridEosGeFlashTest.createReactiveScaleSystem`'s: methane 5 / CO2 0.05 /
// n-heptane 2 / water 55.5 / Ca++ 6e-4 / Cl- 2e-4 / HCO3- 1e-3 at 313.15 K and 50 bar, with and
// without the oil. Its assertions are that the aqueous phase retains bicarbonate and carbonate,
// that the carbonate moved from its `1e-3` seed by more than a tenth of a percent, that every
// element and the net charge are conserved, and that the calcite scale potential is finite,
// positive and repeatable.
//
//     javac -proc:none -cp neqsim-f0c7436.jar HybridEosGeReactiveProbe.java
//     java -cp .:neqsim-f0c7436.jar HybridEosGeReactiveProbe > captures/hybrid_eos_ge_reactive_probe.tsv

import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPitzer;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class HybridEosGeReactiveProbe {

  public static void main(String[] args) throws Exception {
    one("reactive-gas-aqueous", false);
    one("reactive-gas-oil-aqueous", true);
  }

  static void one(String label, boolean includeOil) throws Exception {
    System.out.println("fluid=" + label);
    SystemInterface system = new SystemPitzer(313.15, 50.0);
    system.addComponent("methane", 5.0);
    system.addComponent("CO2", 0.05);
    if (includeOil) {
      system.addComponent("n-heptane", 2.0);
    }
    system.addComponent("water", 55.5);
    system.addComponent("Ca++", 6.0e-4);
    system.addComponent("Cl-", 2.0e-4);
    system.addComponent("HCO3-", 1.0e-3);
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule("classic");
    system.setMultiPhaseCheck(true);

    System.out.println("is_chemical_system=" + system.isChemicalSystem());
    // The conservation matrix the projection is built from, and Commons Math's own
    // pseudo-inverse of it: the port has no SVD, so this is the quantity its own
    // alternative has to reproduce.
    double[][] conservation = system.getChemicalReactionOperations().getAmatrix();
    ComponentInterface[] reactive = system.getChemicalReactionOperations().getComponents();
    StringBuilder names = new StringBuilder();
    for (ComponentInterface component : reactive) {
      names.append(component.getName()).append(" ");
    }
    System.out.println("reactive_components=" + names.toString().trim());
    System.out.println("conservation_rows=" + conservation.length + " cols=" + conservation[0].length);
    for (double[] row : conservation) {
      StringBuilder line = new StringBuilder("  a_row=");
      for (double value : row) {
        line.append(value).append(" ");
      }
      System.out.println(line.toString().trim());
    }
    org.ejml.simple.SimpleMatrix matrix = new org.ejml.simple.SimpleMatrix(conservation);
    org.ejml.simple.SimpleMatrix pinv = matrix.pseudoInverse();
    System.out.println("pseudo_inverse_rows=" + pinv.numRows() + " cols=" + pinv.numCols());
    for (int r = 0; r < pinv.numRows(); r++) {
      StringBuilder line = new StringBuilder("  pinv_row=");
      for (int c = 0; c < pinv.numCols(); c++) {
        line.append(pinv.get(r, c)).append(" ");
      }
      System.out.println(line.toString().trim());
    }
    double[] before = conserved(system, true);
    double initialCarbonate = system.getPhase(0).getComponent("CO3--").getNumberOfmoles();
    System.out.println("initial_carbonate_moles=" + initialCarbonate);

    ThermodynamicOperations operations = new ThermodynamicOperations(system);
    operations.TPflash();

    System.out.println("phases=" + system.getNumberOfPhases());
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      PhaseInterface phase = system.getPhase(p);
      System.out.println("  phase[" + p + "]_type=" + phase.getPhaseTypeName()
          + " beta=" + system.getBeta(p) + " molecules=" + phase.getNumberOfMolesInPhase());
      for (int i = 0; i < phase.getNumberOfComponents(); i++) {
        ComponentInterface component = phase.getComponent(i);
        System.out.println("    x[" + p + "][" + component.getName() + "]=" + component.getx()
            + " moles=" + component.getNumberOfMolesInPhase());
      }
    }

    System.out.println("  aqueous_hco3_moles="
        + find(system, "aqueous", "HCO3-").getNumberOfMolesInPhase());
    System.out.println("  aqueous_co3_moles="
        + find(system, "aqueous", "CO3--").getNumberOfMolesInPhase());
    System.out.println("  aqueous_h3o_moles="
        + find(system, "aqueous", "H3O+").getNumberOfMolesInPhase());
    System.out.println("  aqueous_oh_moles="
        + find(system, "aqueous", "OH-").getNumberOfMolesInPhase());

    double[] after = conserved(system, false);
    double worstElement = 0.0;
    for (int i = 0; i < before.length - 1; i++) {
      worstElement = Math.max(worstElement, Math.abs(before[i] - after[i]));
    }
    System.out.println("worst_element_residual=" + worstElement);
    System.out.println("net_charge_residual=" + Math.abs(before[before.length - 1] - after[after.length - 1]));

    double first = operations.getRelativeScalePotential("CaCO3");
    System.out.println("relative_scale_potential_caco3=" + first);
    operations.TPflash();
    System.out.println("repeated_scale_potential_caco3=" + operations.getRelativeScalePotential("CaCO3"));
    System.out.println();
  }

  /// One component of one phase, by the phase's type name.
  static ComponentInterface find(SystemInterface system, String type, String name) {
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      PhaseInterface phase = system.getPhase(p);
      if (type.equalsIgnoreCase(phase.getPhaseTypeName())) {
        return phase.getComponent(name);
      }
    }
    throw new IllegalStateException("no " + type + " phase");
  }

  /// The conservation matrix applied to the component amounts, then the net charge.
  ///
  /// The same quantity `SystemHybridEosGeFlashTest.getReactiveConservedQuantities` computes,
  /// so what is printed here is what that test asserts.
  static double[] conserved(SystemInterface system, boolean useFeedTotals) {
    ChemicalReactionOperations reactions = system.getChemicalReactionOperations();
    ComponentInterface[] reactiveComponents = reactions.getComponents();
    double[][] conservation = reactions.getAmatrix();
    double[] quantities = new double[conservation.length];
    int components = system.getPhase(0).getNumberOfComponents();

    for (int reactiveIndex = 0; reactiveIndex < reactiveComponents.length; reactiveIndex++) {
      int componentIndex = reactiveComponents[reactiveIndex].getComponentNumber();
      double moles = useFeedTotals ? system.getPhase(0).getComponent(componentIndex).getNumberOfmoles()
          : summed(system, componentIndex);
      for (int row = 0; row < conservation.length - 1; row++) {
        quantities[row] += conservation[row][reactiveIndex] * moles;
      }
    }
    for (int componentIndex = 0; componentIndex < components; componentIndex++) {
      double moles = useFeedTotals ? system.getPhase(0).getComponent(componentIndex).getNumberOfmoles()
          : summed(system, componentIndex);
      quantities[quantities.length - 1] +=
          system.getPhase(0).getComponent(componentIndex).getIonicCharge() * moles;
    }
    return quantities;
  }

  static double summed(SystemInterface system, int componentIndex) {
    double total = 0.0;
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      total += system.getPhase(p).getComponent(componentIndex).getNumberOfMolesInPhase();
    }
    return total;
  }
}
