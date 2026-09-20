// The Pitzer intermediates, printed one at a time so each can be ported and checked.
//
//     javac -proc:none -cp neqsim-3.20.0.jar PitzerArithmetic.java
//     java -cp .:neqsim-3.20.0.jar PitzerArithmetic
//
// `ComponentGePitzer`'s arithmetic is large and mostly private, and a phase-level
// comparison says only that the whole of it agrees or does not. This prints the pieces:
// the Debye-Huckel parameter, the ionic strength, the `g` functions, and the assembled
// `B`, `B'` and `C` - so a port can be wrong in one of them and be seen to be.
//
// `debyeHuckelAphi` is private, so it is reached by reflection. That is a probe's
// privilege and not a port's: the point is to obtain the number, and the alternative -
// inverting it out of an osmotic coefficient - would make the oracle depend on the thing
// it is checking.
//
// **`debyeHuckelAphi(T)` returns `3 * Aphi`, not `Aphi`.** Its own comment says so and
// every caller divides by three; the probe prints both, because the name and the value
// disagree and a port that took the name at face value would be out by a factor of three
// in every Debye-Huckel term - which is a plausible-looking number rather than a failure.

import java.lang.reflect.Method;

import neqsim.thermo.component.ComponentGePitzer;
import neqsim.thermo.phase.PhasePitzer;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPitzer;

public class PitzerArithmetic {

  private static Method privateMethod(String name, Class<?>... parameters) throws Exception {
    Method method = ComponentGePitzer.class.getDeclaredMethod(name, parameters);
    method.setAccessible(true);
    return method;
  }

  public static void main(String[] args) throws Exception {
    Method aphi = privateMethod("debyeHuckelAphi", double.class);
    // An **instance** method, so it needs a receiver: the parameter is a function of the
    // temperature alone, but the method that computes it is not static.
    ComponentGePitzer receiver = anyComponent();

    System.out.println("=== Debye-Huckel parameter against temperature ===");
    System.out.println("   T / K        Agamma = debyeHuckelAphi(T)      Aphi = Agamma / 3");
    for (double t : new double[] {273.15, 283.15, 298.15, 313.15, 323.15, 348.15, 373.15, 398.15}) {
      double agamma = (double) aphi.invoke(receiver, t);
      System.out.printf("  %8.4f    %22.15g    %22.15g%n", t, agamma, agamma / 3.0);
    }

    // The `g` and `g'` functions, which every binary term is built from. Pure functions of
    // `x = alpha * sqrt(I)`, so they are printed at the same `x` a brine produces.
    System.out.println();
    System.out.println("=== g(x) and g'(x) ===");
    System.out.println("   x               g(x)                    g'(x)");
    for (double x : new double[] {0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0}) {
      double g = 2.0 * (1.0 - (1.0 + x) * Math.exp(-x)) / (x * x);
      double gp = -2.0 * (1.0 - (1.0 + x + x * x / 2.0) * Math.exp(-x)) / (x * x);
      System.out.printf("  %6.2f    %22.15g    %22.15g%n", x, g, gp);
    }

    // And the same quantities read off a built phase, which is where the intermediates
    // meet the parameters the selection rule chose.
    for (String[] brine : new String[][] {
        {"water", "Na+", "Cl-"},
        {"water", "Na+", "Ca++", "Cl-"}}) {
      report(brine);
    }
  }

  /** Any built `ComponentGePitzer`, for reaching the instance methods by reflection. */
  private static ComponentGePitzer anyComponent() {
    SystemInterface system = new SystemPitzer(298.15, 1.0);
    system.addComponent("water", 1.0);
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);
    return (ComponentGePitzer) system.getPhase(1).getComponent(0);
  }

  private static void report(String[] names) {
    SystemInterface system = new SystemPitzer(298.15, 1.0);
    double[] z = new double[names.length];
    for (int i = 0; i < names.length; i++) {
      z[i] = names[i].equals("water") ? 0.9 : 0.1 / (names.length - 1);
      system.addComponent(names[i], z[i]);
    }
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);

    PhasePitzer phase = (PhasePitzer) system.getPhase(1);
    double ionicStrength = phase.getIonicStrength();
    double sqrtI = Math.sqrt(ionicStrength);

    System.out.printf("%n=== %s ===%n", String.join(" + ", names));
    System.out.printf("  I        = %.15g mol/kg%n", ionicStrength);
    System.out.printf("  sqrt(I)  = %.15g%n", sqrtI);
    System.out.printf("  dataset  = %s%n", phase.getParameterDatasetId());
    System.out.printf("  common-ion terms = %s%n", phase.isPhreeqcCommonIonTermsActive());

    double zsum = 0.0;
    for (int i = 0; i < phase.getNumberOfComponents(); i++) {
      double charge = phase.getComponent(i).getIonicCharge();
      if (Math.abs(charge) > 0.5) {
        zsum += phase.getComponent(i).getMolality(phase) * Math.abs(charge);
      }
    }
    System.out.printf("  Z = sum m|z| = %.15g%n", zsum);

    for (int i = 0; i < phase.getNumberOfComponents(); i++) {
      for (int j = 0; j < phase.getNumberOfComponents(); j++) {
        if (i == j) {
          continue;
        }
        double zi = phase.getComponent(i).getIonicCharge();
        double zj = phase.getComponent(j).getIonicCharge();
        if (zi * zj >= 0.0) {
          continue;
        }
        double alpha1 = Math.abs(zi) >= 1.5 && Math.abs(zj) >= 1.5 ? 1.4 : 2.0;
        double x1 = alpha1 * sqrtI;
        double g1 = 2.0 * (1.0 - (1.0 + x1) * Math.exp(-x1)) / (x1 * x1);
        double gp1 = -2.0 * (1.0 - (1.0 + x1 + x1 * x1 / 2.0) * Math.exp(-x1)) / (x1 * x1);
        double beta0 = phase.getBeta0ij(i, j, 298.15);
        double beta1 = phase.getBeta1ij(i, j, 298.15);
        double cphi = phase.getCphiij(i, j, 298.15);
        double beta2 = phase.getBeta2ij(i, j, 298.15);
        double bval = beta0 + beta1 * g1;
        double bprime = beta1 * gp1 / ionicStrength;
        double cval = cphi / (2.0 * Math.sqrt(Math.abs(zi * zj)));
        System.out.printf("  pair %-6s %-6s alpha1=%.2f g1=%.12g gp1=%.12g%n",
            phase.getComponent(i).getComponentName(), phase.getComponent(j).getComponentName(), alpha1, g1,
            gp1);
        System.out.printf("       b0=%.12g b1=%.12g b2=%.12g cphi=%.12g -> B=%.12g B'=%.12g C=%.12g%n",
            beta0, beta1, beta2, cphi, bval, bprime, cval);
      }
    }

    System.out.printf("  osmotic coefficient of water = %.15g%n", phase.getOsmoticCoefficientOfWater());
    for (int i = 0; i < phase.getNumberOfComponents(); i++) {
      ComponentGePitzer component = (ComponentGePitzer) phase.getComponent(i);
      double gamma = component.getGamma(phase, phase.getNumberOfComponents(), 298.15, 1.0, phase.getType());
      System.out.printf("  ln gamma(%-6s) = %.15g%n", component.getComponentName(), Math.log(gamma));
    }
  }
}
