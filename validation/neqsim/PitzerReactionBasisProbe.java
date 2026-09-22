// The solute-molality concentration basis of `ChemicalEquilibrium`, for
// `reactions.chemical_equilibrium`.
//
// `ChemicalEquilibrium.getLogReactionActivity` has three branches, and the two captured
// fluids take only the first. On the mole-fraction basis the activity term is
// `ln n_i - ln n_t + ln gamma_i`; on `SOLUTE_MOLALITY` - which is `SystemPitzer`'s, and
// the only system that selects it - a **solvent** component keeps that form while a
// **solute** takes `ln n_i - ln w_solvent`, where `w_solvent` is the solvent's own mass in
// kg (`calculateSolventWeight`: `sum(n_j * M_j)` over the components whose reference state
// is `solvent`). The third branch is a *phase* and not a basis: `PhaseDeshmukhMather`
// converts its stored mole-fraction coefficient to the molality scale through
// `getMolality` and `getSolventMolarMass`.
//
// The fluid is the one NeqSim's own `PitzerHydrogenSulfideEquilibriumTest` builds, because
// that is the only place the branch is driven: `SystemPitzer` at 1.01325 bar, water and one
// acid, `setMixingRule("classic")`, then `solveChemEq(1, 0)` and `solveChemEq(1, 1)`.
//
//     javac -proc:none -cp neqsim-f0c7436.jar PitzerReactionBasisProbe.java
//     java -cp .:neqsim-f0c7436.jar PitzerReactionBasisProbe > captures/pitzer_reaction_basis_probe.tsv

import java.lang.reflect.Method;
import java.util.Map;
import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.thermo.ThermodynamicConstantsInterface;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemPitzer;

public class PitzerReactionBasisProbe {

  public static void main(String[] args) {
    one("co2-water-pitzer", 298.15, "CO2", 0.01, 10.0);
    one("h2s-water-pitzer", 298.15, "H2S", 0.01, 10.0);
  }

  static void one(String label, double temperature, String solute, double soluteMoles,
      double waterMoles) {
    System.out.println("fluid=" + label);
    SystemPitzer system = new SystemPitzer(temperature, 1.01325);
    system.addComponent("water", waterMoles);
    system.addComponent(solute, soluteMoles);
    system.setMultiPhaseCheck(false);
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);

    System.out.println("concentration_basis=" + system.getChemicalReactionConcentrationBasis());
    System.out.println("phases=" + system.getNumberOfPhases());
    System.out.println("temperature_K=" + system.getTemperature());
    System.out.println("pressure_bara=" + system.getPressure());

    ChemicalReactionOperations operations = system.getChemicalReactionOperations();
    int reactive = reactivePhaseIndex(operations);
    System.out.println("reactive_phase=" + reactive);
    int readPhase = reactive < 0 ? 0 : reactive;
    operations.setReactiveComponents(readPhase);
    ComponentInterface[] components = operations.getComponents();
    PhaseInterface phase = system.getPhase(readPhase);

    System.out.println("phase_moles=" + phase.getNumberOfMolesInPhase());
    for (int i = 0; i < components.length; i++) {
      System.out.println("  component[" + i + "]=" + components[i].getName()
          + " reference_state=" + components[i].getReferenceStateType()
          + " molar_mass=" + components[i].getMolarMass()
          + " moles=" + components[i].getNumberOfMolesInPhase()
          + " x=" + components[i].getx()
          + " log_activity_coefficient="
          + phase.getLogActivityCoefficient(components[i].getComponentNumber(), 0));
    }

    // `calculateSolventWeight`, recomputed here because the field is private: the solvent
    // component's moles times its molar mass, in kg.
    double solventWeight = 0.0;
    for (int i = 0; i < phase.getNumberOfComponents(); i++) {
      ComponentInterface held = phase.getComponent(i);
      if ("solvent".equalsIgnoreCase(held.getReferenceStateType())) {
        solventWeight += held.getNumberOfMolesInPhase() * held.getMolarMass();
      }
    }
    System.out.println("solvent_weight_kg=" + solventWeight);

    // What fills `logactivityVec`, which the solver reads once in its constructor.
    int solventIndex = -1;
    for (int i = 0; i < phase.getNumberOfComponents(); i++) {
      if ("solvent".equalsIgnoreCase(phase.getComponent(i).getReferenceStateType())) {
        solventIndex = i;
      }
    }
    for (int i = 0; i < components.length; i++) {
      int componentNumber = components[i].getComponentNumber();
      double logActivity = 0.0;
      if (components[i].calcActivity() && solventIndex >= 0) {
        logActivity = phase.getLogActivityCoefficient(componentNumber, solventIndex);
      }
      System.out.println("  log_activity[" + components[i].getName() + "]=" + logActivity);
    }

    printVector("moles_before", phaseMoles(components, phase));

    // --- the direct solve, which is the id's own answer ------------------------
    double[][] a = operations.getAmatrix();
    double[] b = operations.calcBVector();
    printMatrix("A", a);
    printVector("b", b);
    double[] chemRef = new double[components.length];
    for (int i = 0; i < components.length; i++) {
      chemRef[i] = components[i].getReferencePotential()
          / (ThermodynamicConstantsInterface.R * phase.getTemperature());
    }
    printVector("chem_ref_reduced", chemRef);
    neqsim.chemicalreactions.chemicalequilibrium.ChemicalEquilibrium direct =
        new neqsim.chemicalreactions.chemicalequilibrium.ChemicalEquilibrium(a, b, system,
            components, readPhase);
    boolean directConverged = direct.solve();
    System.out.println("direct_converged=" + directConverged);
    System.out.println("direct_iterations=" + direct.getLastIterationCount());
    System.out.println("direct_error=" + direct.getLastError());
    printVector("direct_moles", direct.getMoles());

    boolean first = operations.solveChemEq(readPhase, 0);
    boolean second = operations.solveChemEq(readPhase, 1);
    system.init(1);
    PhaseInterface after = system.getPhase(readPhase);
    System.out.println("solve_type0=" + first);
    System.out.println("solve_type1=" + second);
    printVector("moles_after", phaseMoles(components, after));
    System.out.println("max_reaction_log_residual="
        + operations.getMaximumAbsoluteReactionLogResidual());
    System.out.println("max_element_residual="
        + operations.getMaximumAbsoluteElementBalanceResidual());
    System.out.println("reactive_phase_charge_moles=" + operations.getReactivePhaseChargeMoles());
    Map<String, Double> residuals = operations.getReactionLogResiduals();
    for (Map.Entry<String, Double> entry : residuals.entrySet()) {
      System.out.println("  residual[" + entry.getKey() + "]=" + entry.getValue());
    }
    System.out.println();
  }

  /// `getReactivePhaseIndex()` is private; this calls it rather than copying it.
  static int reactivePhaseIndex(ChemicalReactionOperations operations) {
    try {
      Method method = ChemicalReactionOperations.class.getDeclaredMethod("getReactivePhaseIndex");
      method.setAccessible(true);
      return (Integer) method.invoke(operations);
    } catch (ReflectiveOperationException ex) {
      throw new RuntimeException(ex);
    }
  }

  static double[] phaseMoles(ComponentInterface[] components, PhaseInterface phase) {
    double[] out = new double[components.length];
    for (int i = 0; i < components.length; i++) {
      out[i] = phase.getComponent(components[i].getComponentNumber()).getNumberOfMolesInPhase();
    }
    return out;
  }

  static void printMatrix(String name, double[][] m) {
    System.out.println(name + "=" + m.length + "x" + (m.length == 0 ? 0 : m[0].length));
    for (int r = 0; r < m.length; r++) {
      StringBuilder row = new StringBuilder("  " + name + "[" + r + "]=");
      for (int c = 0; c < m[r].length; c++) {
        row.append(m[r][c]).append(c + 1 < m[r].length ? " " : "");
      }
      System.out.println(row);
    }
  }

  static void printVector(String name, double[] v) {
    StringBuilder out = new StringBuilder(name + "=");
    for (int i = 0; i < v.length; i++) {
      out.append(v[i]).append(i + 1 < v.length ? " " : "");
    }
    System.out.println(out);
  }
}
