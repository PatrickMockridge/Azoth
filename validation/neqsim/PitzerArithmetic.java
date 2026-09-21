// The Pitzer intermediates, printed one at a time so each can be ported and checked.
//
//     javac -proc:none -cp neqsim-f0c7436.jar PitzerArithmetic.java
//     java -cp .:neqsim-f0c7436.jar PitzerArithmetic
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
import neqsim.thermo.phase.PitzerElectrostaticMixing;
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

    // **The parameter getters, which is what the temperature forms are for.** The two
    // datasets use different forms, so each is printed across temperature: a catalogue
    // pair (the PHREEQC six coefficients) and a CSV pair (Silvester-Pitzer's two).
    parameters();
    activityCoefficients();
  }

  private static void activityCoefficients() {
    System.out.println();
    System.out.println("=== the ion activity coefficient, end to end ===");
    // **The third brine is the legacy path's.** The catalogue carries no `C0` row for
    // Na+/HCO3-, so that mixture falls back to `PitzerParameters.csv` and reads every
    // parameter through Silvester-Pitzer's two-coefficient form instead of the six. It is
    // also a single salt, so it is the one topology NeqSim initializes with an absent pair
    // rather than refusing - which is why the audit and the refusal are separate.
    String[][] brines = {
        {"water", "Na+", "Cl-"},
        {"water", "Na+", "Ca++", "Cl-"},
        {"water", "Na+", "HCO3-"}};
    double[][] compositions = {
        {0.88, 0.06, 0.06}, {0.88, 0.03, 0.03, 0.06}, {0.88, 0.06, 0.06}};
    for (int b = 0; b < brines.length; b++) {
      String[] names = brines[b];
      double[] z = compositions[b];
      SystemInterface system = new SystemPitzer(298.15, 1.0);
      for (int i = 0; i < names.length; i++) {
        system.addComponent(names[i], z[i]);
      }
      system.setMixingRule("classic");
      system.init(0);
      system.init(1);
      PhasePitzer phase = (PhasePitzer) system.getPhase(1);
      System.out.printf("  %s   (dataset=%s, common-ion=%s, unequal-charge-same-sign=%s, "
              + "non-2-2 beta2=%s)%n",
          String.join(" + ", names),
          phase.getParameterDatasetId().startsWith("usgs") ? "phreeqc" : "legacy",
          phase.isPhreeqcCommonIonTermsActive(),
          phase.hasUnequalChargeSameSignPair(), phase.isNonTwoTwoBeta2Active());
      // `getGamma` for water already *is* the water gamma - it dispatches to
      // `getWaterGamma` - so the ion loop below prints it and this adds the osmotic
      // coefficient, which `getGamma` does not return.
      System.out.printf("    I = %.15g   phi = %.15g%n", phase.getIonicStrength(),
          phase.getOsmoticCoefficientOfWater());
      for (int i = 0; i < phase.getNumberOfComponents(); i++) {
        ComponentGePitzer component = (ComponentGePitzer) phase.getComponent(i);
        double gamma = component.getGamma(phase, phase.getNumberOfComponents(), 298.15, 1.0, phase.getType());
        System.out.printf("    ln gamma(%-6s) = %22.15g   molality = %.12g   z = %+.0f%n",
            component.getComponentName(), Math.log(gamma), component.getMolality(phase),
            component.getIonicCharge());
      }
    }

    System.out.println();
    System.out.println("=== the neutral layer, on a CO2-bearing brine ===");
    SystemInterface system = new SystemPitzer(298.15, 1.0);
    system.addComponent("water", 0.86);
    system.addComponent("Na+", 0.03);
    system.addComponent("Cl-", 0.03);
    system.addComponent("CO2", 0.08);
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);
    PhasePitzer neutralPhase = (PhasePitzer) system.getPhase(1);
    System.out.printf("  dataset = %s%n",
        neutralPhase.getParameterDatasetId().startsWith("usgs") ? "phreeqc" : "legacy");
    System.out.printf("  neutral interactions active = %s%n",
        neutralPhase.hasNeutralPitzerInteractions());
    System.out.printf("  osmotic contribution = %.15g%n",
        neutralPhase.getNeutralPitzerOsmoticContribution(298.15));
    System.out.printf("  I = %.15g   phi = %.15g%n", neutralPhase.getIonicStrength(),
        neutralPhase.getOsmoticCoefficientOfWater());
    for (int i = 0; i < neutralPhase.getNumberOfComponents(); i++) {
      ComponentGePitzer component = (ComponentGePitzer) neutralPhase.getComponent(i);
      System.out.printf("    %-6s ln gamma = %22.15g   neutral contribution = %22.15g%n",
          component.getComponentName(),
          Math.log(component.getGamma(neutralPhase, neutralPhase.getNumberOfComponents(), 298.15, 1.0,
              neutralPhase.getType())),
          neutralPhase.getNeutralPitzerLogGammaContribution(i, 298.15));
    }

    System.out.println();
    System.out.println("=== PitzerElectrostaticMixing, the E_theta integral ===");
    System.out.printf("  %8s %8s %12s %10s %22s %22s%n", "z_j", "z_k", "I", "Aphi", "E_theta", "dE_theta/dI");
    double[] result = new double[2];
    double aphi = 0.392034451863750;
    for (double[] row : new double[][] {{1.0, 2.0, 6.0}, {1.0, 2.0, 0.5}, {2.0, 2.0, 6.0},
        {1.0, 3.0, 1.0}, {1.0, 2.0, 100.0}}) {
      PitzerElectrostaticMixing.calculate(row[0], row[1], row[2], aphi, result);
      System.out.printf("  %8.1f %8.1f %12.6g %10.6g %22.15g %22.15g%n", row[0], row[1], row[2], aphi,
          result[0], result[1]);
    }
  }

  private static void parameters() {
    System.out.println();
    System.out.println("=== the parameter getters across temperature ===");
    for (String label : new String[] {"catalogue", "legacy"}) {
      String[] names = label.equals("catalogue")
          ? new String[] {"water", "Na+", "Cl-"}
          : new String[] {"water", "Na+", "HCO3-"};
      SystemInterface system = new SystemPitzer(298.15, 1.0);
      system.addComponent(names[0], 0.9);
      system.addComponent(names[1], 0.05);
      system.addComponent(names[2], 0.05);
      system.setMixingRule("classic");
      system.init(0);
      system.init(1);
      PhasePitzer phase = (PhasePitzer) system.getPhase(1);
      System.out.printf("  %s: %s / %s   (dataset %s)%n", label, names[1], names[2],
          phase.getParameterDatasetId().startsWith("usgs") ? "phreeqc" : "legacy");
      int first = phase.getComponent(names[1]).getComponentNumber();
      int second = phase.getComponent(names[2]).getComponentNumber();
      System.out.printf("    %8s %18s %18s %18s %18s %18s%n", "T/K", "beta0", "beta1", "Cphi",
          "beta2", "theta");
      for (double temp : new double[] {273.15, 298.15, 298.15005, 323.15, 373.15}) {
        System.out.printf("    %8.4f %18.12g %18.12g %18.12g %18.12g %18.12g%n", temp,
            phase.getBeta0ij(first, second, temp), phase.getBeta1ij(first, second, temp),
            phase.getCphiij(first, second, temp), phase.getBeta2ij(first, second, temp),
            phase.getThetaij(first, second, temp));
      }
    }

    System.out.println();
    System.out.println("=== the alpha coefficients against charge ===");
    System.out.printf("  %6s %6s %10s %10s%n", "z1", "z2", "alpha1", "alpha2");
    SystemInterface system = new SystemPitzer(298.15, 1.0);
    system.addComponent("water", 0.99);
    system.addComponent("Na+", 0.01);
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);
    PhasePitzer phase = (PhasePitzer) system.getPhase(1);
    // The charge enters through the components, so this is read off the pair the phase
    // actually holds; a wider charge sweep would need wider components.
    System.out.printf("  %6.1f %6.1f %10.4f %10.4f%n", 1.0, -1.0,
        phase.getPitzerAlpha1(0, 1), phase.getPitzerAlpha2(0, 1));
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
