// The effective-diffusion assembly's input, for `reactions.kinetics`.
//
// `PhysicalProperties.calcEffectiveDiffusionCoefficients` is eight lines: for each component,
// `D_eff_i = (1 - x_i) / sum_{j != i} x_j / D_ij` over a **binary** matrix the diffusivity
// model fills. This port takes the effective vector as an input because a first reading of
// the class concluded the matrix was unreachable - `getFickDiffusionCoefficient(i, j)`
// returned a diagonal array and the field behind it is package-private.
//
// **That reading asked the wrong method.** `getFickDiffusionCoefficient` and
// `getFickBinaryDiffusionCoefficient` are two different methods, and it is the *interface*
// that has to be asked:
//
//   DiffusivityInterface.calcDiffusionCoefficients(int, int)   returns the matrix itself
//   DiffusivityInterface.getFickBinaryDiffusionCoefficient(i,j)  reads one entry of it
//   PhysicalProperties.diffusivityCalc                          is a **public field**
//
// So the input is observable with no cast at all, and the assembly is portable. This probe
// prints the matrix, the vector the class computes from it and the vector recomputed from
// the matrix here - after re-running the assembly, so that both come from the same fill.
//
//     javac -proc:none -cp neqsim-f0c7436.jar EffectiveDiffusionProbe.java
//     java -cp .:neqsim-f0c7436.jar EffectiveDiffusionProbe > captures/effective_diffusion_probe.tsv

import neqsim.physicalproperties.methods.methodinterface.DiffusivityInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class EffectiveDiffusionProbe {

  public static void main(String[] args) {
    one("co2-water-298K", 298.15, new String[] { "CO2", "water" }, new double[] { 0.01, 10.0 });
    one("methane-nheptane-313K", 313.15, new String[] { "methane", "n-heptane" },
        new double[] { 0.5, 0.5 });
  }

  static void one(String label, double temperature, String[] names, double[] moles) {
    System.out.println("fluid=" + label);
    SystemInterface system = new SystemSrkEos(temperature, 1.01325);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.setMixingRule(2);
    system.init(0);
    ThermodynamicOperations flash = new ThermodynamicOperations(system);
    flash.TPflash();
    system.init(1);

    int phaseIndex = system.getNumberOfPhases() > 1 ? 1 : 0;
    PhaseInterface phase = system.getPhase(phaseIndex);
    int count = phase.getNumberOfComponents();
    System.out.println("phase_index=" + phaseIndex);
    System.out.println("phase_type=" + phase.getPhaseTypeName());
    System.out.println("components=" + count);

    DiffusivityInterface diffusivity = phase.getPhysicalProperties().diffusivityCalc;
    System.out.println("model=" + diffusivity.getClass().getName());

    // The matrix the model fills, from the interface's own returning method.
    double[][] binary = diffusivity.calcDiffusionCoefficients(0, 0);
    StringBuilder matrix = new StringBuilder("binary_diffusion=");
    for (int i = 0; i < count; i++) {
      for (int j = 0; j < count; j++) {
        matrix.append(binary[i][j]).append(" ");
      }
    }
    System.out.println(matrix.toString().trim());

    // **The same matrix read one entry at a time**, so the two getters are shown to be the
    // same array and not two spellings of one name.
    StringBuilder entries = new StringBuilder("binary_entries=");
    for (int i = 0; i < count; i++) {
      for (int j = 0; j < count; j++) {
        entries.append(diffusivity.getFickBinaryDiffusionCoefficient(i, j)).append(" ");
      }
    }
    System.out.println(entries.toString().trim());

    // **The assembly, re-run on that matrix.** `calcEffectiveDiffusionCoefficients` reads the
    // array the call above filled, so the vector printed here is the one this matrix gives.
    diffusivity.calcEffectiveDiffusionCoefficients();
    StringBuilder effective = new StringBuilder("effective_diffusion=");
    for (int i = 0; i < count; i++) {
      effective.append(diffusivity.getEffectiveDiffusionCoefficient(i)).append(" ");
    }
    System.out.println(effective.toString().trim());

    // And the same eight lines written out here, which is what a port has to reproduce.
    StringBuilder recomputed = new StringBuilder("effective_recomputed=");
    for (int i = 0; i < count; i++) {
      double sum = 0.0;
      for (int j = 0; j < count; j++) {
        if (i != j) {
          sum += phase.getComponent(j).getx() / binary[i][j];
        }
      }
      recomputed.append((1.0 - phase.getComponent(i).getx()) / sum).append(" ");
    }
    System.out.println(recomputed.toString().trim());

    for (int i = 0; i < count; i++) {
      System.out.println("  x[" + phase.getComponent(i).getName() + "]="
          + phase.getComponent(i).getx());
    }
    System.out.println();
  }
}
