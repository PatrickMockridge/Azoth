// The LP initial estimate `LinearProgrammingChemicalEquilibrium`, for the seed machinery
// of `reactions.reactive_phase_equilibrium`.
//
// `ChemicalReactionOperations.solveChemEq` seeds its Newton solve from a linear program:
// minimise `sum(mu_i n_i / R T)` subject to `A n = b` and `n >= 0`, where `A` is the
// element matrix with the electroneutrality row last. **The seed is not a suggestion** -
// `updateMoles` writes it back into the phase before the Newton solve starts - so the
// starting composition is part of the answer.
//
// **The LP is infeasible on the flashed states.** `seed_null=true` below is the finding
// this probe exists for: on both captured aqueous phases `generateInitialEstimates` throws
// `NoFeasibleSolutionException`, `solveChemEq` clears `newMoles`, and the Newton solve
// starts from the phase as it is. It is feasible on the unflashed fluids, which is why one
// of them is here, and there the converged composition is the same with and without it.
//
// **Two paths are printed for each state, and they are not "with and without the seed".**
// `TPflash` runs the reaction machinery itself, so `newMoles` survives from the flash. In
// the production path (`initCalc` alive) a null LP answer assigns `newMoles = null` and
// `updateMoles` is skipped. Nulling `initCalc` by reflection skips that assignment instead,
// and then `updateMoles` **writes the flash's stale vector into the phase** - which is what
// `without_initcalc_moles` measures, and it is a defect rather than a second experiment.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ReactionLpSeedProbe.java
//     java -cp .:neqsim-f0c7436.jar ReactionLpSeedProbe > captures/reaction_lp_seed_probe.tsv

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import neqsim.chemicalreactions.ChemicalReactionOperations;
import neqsim.chemicalreactions.chemicalequilibrium.LinearProgrammingChemicalEquilibrium;
import neqsim.thermo.ThermodynamicConstantsInterface;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemSrkEos;

public class ReactionLpSeedProbe {

  public static void main(String[] args) {
    one("co2-water", 298.15, new String[] { "CO2", "water" }, new double[] { 0.01, 10.0 }, false);
    one("co2-water-aqueous", 298.15, new String[] { "CO2", "water" }, new double[] { 0.01, 10.0 },
        true);
    one("co2-h2s-water-aqueous", 298.15, new String[] { "CO2", "H2S", "water" },
        new double[] { 0.01, 0.01, 10.0 }, true);
  }

  static SystemSrkEos build(double temperature, String[] names, double[] moles) {
    SystemSrkEos system = new SystemSrkEos(temperature, 1.01325);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.chemicalReactionInit();
    system.createDatabase(true);
    system.setMixingRule(2);
    system.init(0);
    return system;
  }

  static SystemSrkEos prepared(double temperature, String[] names, double[] moles, boolean flash) {
    SystemSrkEos system = build(temperature, names, moles);
    if (flash) {
      system.setMultiPhaseCheck(true);
      new neqsim.thermodynamicoperations.ThermodynamicOperations(system).TPflash();
      system.init(3);
    }
    return system;
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

  static Object readField(ChemicalReactionOperations operations, String name) {
    try {
      Field field = ChemicalReactionOperations.class.getDeclaredField(name);
      field.setAccessible(true);
      return field.get(operations);
    } catch (ReflectiveOperationException ex) {
      throw new RuntimeException(ex);
    }
  }

  static void one(String label, double temperature, String[] names, double[] moles, boolean flash) {
    System.out.println("fluid=" + label);
    System.out.println("flashed=" + flash);

    SystemSrkEos system = prepared(temperature, names, moles, flash);
    ChemicalReactionOperations operations = system.getChemicalReactionOperations();
    int reactive = reactivePhaseIndex(operations);
    int readPhase = reactive < 0 ? 0 : reactive;
    System.out.println("reactive_phase=" + reactive);
    System.out.println("phase_moles=" + system.getPhase(readPhase).getNumberOfMolesInPhase());
    // What the flash left behind, which is the value the stale path writes back.
    System.out.println("newMoles_left_by_flash_null="
        + (readField(operations, "newMoles") == null));

    operations.setReactiveComponents(readPhase);
    ComponentInterface[] components = operations.getComponents();
    System.out.println("components=" + components.length);
    for (int i = 0; i < components.length; i++) {
      System.out.println("  component[" + i + "]=" + components[i].getName());
    }

    double[] b = operations.calcBVector();
    printVector("b", b);
    printVector("n0", phaseMoles(system, components, readPhase));

    // The constructor recomputes the reference potentials and the matrix, so what is
    // printed here is the state `generateInitialEstimates` runs on.
    LinearProgrammingChemicalEquilibrium lp = new LinearProgrammingChemicalEquilibrium(
        operations.calcChemRefPot(readPhase).clone(), components, operations.getAllElements(),
        operations, readPhase);
    printMatrix("A", lp.getA());
    printVector("chem_ref", lp.getRefPot());

    double phaseTemperature = system.getPhase(readPhase).getTemperature();
    double[] objective = new double[components.length];
    for (int i = 0; i < components.length; i++) {
      objective[i] = lp.getRefPot()[i] / (ThermodynamicConstantsInterface.R * phaseTemperature);
    }
    printVector("objective", objective);

    double[] seed = lp.generateInitialEstimates(system, b, operations.calcInertMoles(readPhase),
        readPhase);
    System.out.println("seed_null=" + (seed == null));
    if (seed != null) {
      printVector("seed", seed);
      double value = 0.0;
      for (int i = 0; i < seed.length; i++) {
        value += objective[i] * seed[i];
      }
      System.out.println("seed_objective=" + value);
      StringBuilder active = new StringBuilder("seed_active=");
      for (int i = 0; i < seed.length; i++) {
        active.append(seed[i] > 1e-30 ? components[i].getName() : "-")
            .append(i + 1 < seed.length ? " " : "");
      }
      System.out.println(active);
    }

    solveAndPrint("with_lp", names, moles, flash, false);
    solveAndPrint("without_initcalc", names, moles, flash, true);
    System.out.println();
  }

  static void solveAndPrint(String key, String[] names, double[] moles, boolean flash,
      boolean nullTheSeed) {
    SystemSrkEos system = prepared(298.15, names, moles, flash);
    ChemicalReactionOperations operations = system.getChemicalReactionOperations();
    int reactive = reactivePhaseIndex(operations);
    int readPhase = reactive < 0 ? 0 : reactive;
    operations.setReactiveComponents(readPhase);
    if (nullTheSeed) {
      try {
        Field initCalc = ChemicalReactionOperations.class.getDeclaredField("initCalc");
        initCalc.setAccessible(true);
        initCalc.set(operations, null);
      } catch (ReflectiveOperationException ex) {
        throw new RuntimeException(ex);
      }
    }
    boolean converged = operations.solveChemEq(readPhase, 0);
    system.init(3);
    System.out.println(key + "_converged=" + converged);
    printVector(key + "_moles", phaseMoles(system, operations.getComponents(), readPhase));
    System.out.println(key + "_newMoles_used_null="
        + (readField(operations, "newMoles") == null));
  }

  static double[] phaseMoles(SystemSrkEos system, ComponentInterface[] components, int phase) {
    double[] out = new double[components.length];
    for (int i = 0; i < components.length; i++) {
      out[i] = system.getPhase(phase).getComponent(components[i].getComponentNumber())
          .getNumberOfMolesInPhase();
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
