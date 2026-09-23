// The Pitzer phase's fugacity coefficients, for `eos.pitzer_phase`.
//
// `ComponentGePitzer.fugcoef(phase)` is where the coefficient is built, and it is **two
// methods rather than three branches in one**:
//
//   a neutral that is not water   gamma * H * (m / x) / P        `ComponentGePitzer`
//   anything else                 `super.fugcoef`, i.e. `ComponentGE`, which is
//                                 `gamma * P0 / P` for a solvent (water) and
//                                 `(gamma / gamma_inf) * H / P` otherwise (an ion)
//
// so the branch a component takes is `|z| < 0.5 && name != water` first, and only then the
// reference state. Pressure divides in **bar**, which is NeqSim's internal unit: the ion
// branch is `(gamma / gamma_inf) * 1e12 / P_bar`.
//
// **`H` has two arms and the second one is reachable.** `ComponentGE.getEffectiveHenryCoefficient`
// prefers the IAPWS pure-water table for a supported neutral solute in a water-bearing
// phase, and `ComponentGePitzer` overrides it with three gates in front: the IAPWS table
// only applies where the phase has **no ions**, the phase carries **no active neutral
// Pitzer interaction family**, and the solute is **not CO2 or H2S**. So every ion-bearing
// fluid in this capture takes the database correlation, and the IAPWS table is reached by
// the four **salt-free** fluids at the end. `IapwsHenryLaw` is a public class, so its
// support and range predicates and its own coefficient are printed beside the effective
// one - which is what makes the branch each component took readable off the capture rather
// than inferred.
//
// `gamma_inf` is `PhaseGE.getActivityCoefficientInfDilWater(k, waterIndex)`, which builds a
// **two-component reference phase** - the solute at `1e-10` mol in slot 0, the solvent at
// `10.0` mol in slot 1 (`Phase.initRefPhases`) - and returns the *solute's* gamma there.
//
// `getEffectiveHenryCoefficient` and `isHenryCoefficientCapped` are `protected`, so they are
// read through reflection up the class hierarchy - the alternative would be re-deriving them
// and checking a derivation against itself. `getActivityCoefficientInfDilWater` is read
// defensively: it re-inits the phase and throws a `NullPointerException` where a component
// array entry is unset, which is a fact about the class and not about this file.
//
// **No flash and no reaction initialisation**: the fluids are read at their constructed
// state, because `SystemPitzer`'s constructor already puts `PhasePitzer` in slot 1 and
// `init(1)` is what calls `fugcoef`. A `TPflash` here relabels and can collapse the phases,
// which would change what is being measured.
//
//     javac -proc:none -cp neqsim-f0c7436.jar PitzerFugacityProbe.java
//     java -cp .:neqsim-f0c7436.jar PitzerFugacityProbe > captures/pitzer_fugacity_probe.tsv

import java.lang.reflect.Method;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.component.IapwsHenryLaw;
import neqsim.thermo.phase.PhaseGE;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPitzer;

public class PitzerFugacityProbe {

  public static void main(String[] args) {
    one("co2-water-nacl", 298.15, 10.0, new String[] { "CO2", "water", "Na+", "Cl-" },
        new double[] { 0.5, 55.508, 1.0, 1.0 });
    one("methane-water-nacl", 298.15, 10.0, new String[] { "methane", "water", "Na+", "Cl-" },
        new double[] { 5.0, 55.508, 1.0, 1.0 });
    one("h2s-water", 298.15, 1.01325, new String[] { "H2S", "water" },
        new double[] { 0.01, 10.0 });
    one("water-methanol-nacl", 313.15, 5.0, new String[] { "water", "methanol", "Na+", "Cl-" },
        new double[] { 50.0, 5.0, 2.0, 2.0 });
    // The four above all carry ions, and `ComponentGePitzer.getEffectiveHenryCoefficient` sends
    // an ion-bearing phase to the database correlation. These four have **no salt**, which is
    // the only way the IAPWS table is reached at all.
    one("methane-water", 298.15, 1.01325, new String[] { "methane", "water" },
        new double[] { 0.01, 10.0 });
    one("nitrogen-water", 298.15, 1.01325, new String[] { "nitrogen", "water" },
        new double[] { 0.01, 10.0 });
    one("co2-water", 298.15, 1.01325, new String[] { "CO2", "water" },
        new double[] { 0.5, 10.0 });
    one("methane-water-333", 333.15, 1.01325, new String[] { "methane", "water" },
        new double[] { 0.01, 10.0 });
  }

  static void one(String label, double temperature, double pressure, String[] names,
      double[] moles) {
    System.out.println("fluid=" + label);
    SystemInterface system = new SystemPitzer(temperature, pressure);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.setMixingRule("classic");
    system.init(0);
    system.init(1);

    System.out.println("phases=" + system.getNumberOfPhases());
    System.out.println("pressure_pa=" + system.getPressure() * 1.0e5);
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      PhaseInterface phase = system.getPhase(p);
      if (!(phase instanceof PhaseGE)) {
        System.out.println("  phase[" + p + "]_type=" + phase.getPhaseTypeName()
            + " is_not_a_ge_phase=true");
        continue;
      }
      int waterIndex = -1;
      for (int i = 0; i < phase.getNumberOfComponents(); i++) {
        if ("water".equalsIgnoreCase(phase.getComponent(i).getName())) {
          waterIndex = i;
        }
      }
      System.out.println("  phase[" + p + "]_type=" + phase.getPhaseTypeName()
          + " density=" + phase.getPhysicalProperties().getDensity());
      for (int i = 0; i < phase.getNumberOfComponents(); i++) {
        ComponentInterface component = phase.getComponent(i);
        System.out.println("    component[" + component.getName() + "]_x=" + component.getx()
            + " charge=" + component.getIonicCharge()
            + " fugacity_coefficient=" + component.getFugacityCoefficient()
            + " ln_phi=" + Math.log(component.getFugacityCoefficient()));
        System.out.println("      gamma=" + gamma(component, phase)
            + " molality=" + molality(component, phase)
            + " antoine_p0=" + component.getAntoineVaporPressure(phase.getTemperature()));
        System.out.println("      activinf_water=" + activinf(phase, i, waterIndex)
            + " henry_effective=" + henry(component, phase)
            + " henry_capped_1e12=" + capped(component, 1.0e12)
            + " henry_capped_1=" + capped(component, 1.0)
            + " iapws_supported=" + IapwsHenryLaw.isSupportedSpecies(component.getName())
            + " iapws_usable=" + IapwsHenryLaw.isUsable(component.getName(),
                phase.getTemperature())
            + " iapws_henry_bar=" + iapwsHenry(component.getName(),
                phase.getTemperature()));
      }
    }
    System.out.println();
  }

  /// `ComponentGePitzer.getGamma(phase, n, T, P, phaseType)`.
  static Object gamma(ComponentInterface component, PhaseInterface phase) {
    // The fifth parameter is `PhaseType`, not an `int`, and there are two overloads: the
    // one carried here is the shorter, which is what `fugcoef` itself calls.
    return call(component, "getGamma",
        new Class<?>[] { PhaseInterface.class, int.class, double.class, double.class,
            phase.getType().getClass() },
        new Object[] { phase, phase.getNumberOfComponents(), phase.getTemperature(),
            phase.getPressure(), phase.getType() });
  }

  /// `getMolality`, which the Pitzer components override to `n_i / getSolventWeight()`.
  static Object molality(ComponentInterface component, PhaseInterface phase) {
    Object withPhase = call(component, "getMolality", new Class<?>[] { PhaseInterface.class },
        new Object[] { phase });
    if (!"no such method".equals(withPhase)) {
      return withPhase;
    }
    return call(component, "getMolality", new Class<?>[] {}, new Object[] {});
  }

  /// `getEffectiveHenryCoefficient(phase)`, which is `protected` on the Pitzer component.
  static Object henry(ComponentInterface component, PhaseInterface phase) {
    return call(component, "getEffectiveHenryCoefficient",
        new Class<?>[] { PhaseInterface.class }, new Object[] { phase });
  }

  /// `isHenryCoefficientCapped(double)`, also `protected`.
  static Object capped(ComponentInterface component, double value) {
    return call(component, "isHenryCoefficientCapped", new Class<?>[] { double.class },
        new Object[] { value });
  }

  /// The IAPWS mole-fraction Henry constant in bar, with the extrapolation allowed inside the
  /// liquid-water equation domain. Throws for the species the table does not carry - which is
  /// most of them, water included - so the throw is reported rather than propagated.
  static Object iapwsHenry(String name, double temperature) {
    try {
      return IapwsHenryLaw.getHenryCoefficientBarAllowExtrapolation(name, temperature);
    } catch (RuntimeException e) {
      return "threw_" + e.getClass().getSimpleName();
    }
  }

  /// `PhaseGE.getActivityCoefficientInfDilWater(k, solventIndex)`.
  ///
  /// **Read defensively: it re-inits the phase and throws a `NullPointerException` where a
  /// component array entry is unset**, which is a fact about the class and not about this
  /// file, so it is reported rather than propagated.
  static String activinf(PhaseInterface phase, int k, int solventIndex) {
    if (solventIndex < 0) {
      return "no_water";
    }
    try {
      return String.valueOf(
          ((PhaseGE) phase).getActivityCoefficientInfDilWater(k, solventIndex));
    } catch (RuntimeException e) {
      return "threw_" + e.getClass().getSimpleName();
    }
  }

  /// A method walked up the class hierarchy, by explicit signature.
  static Object call(ComponentInterface component, String method, Class<?>[] signature,
      Object[] arguments) {
    for (Class<?> type = component.getClass(); type != null; type = type.getSuperclass()) {
      Method candidate;
      try {
        candidate = type.getDeclaredMethod(method, signature);
      } catch (NoSuchMethodException ignored) {
        continue;
      }
      try {
        candidate.setAccessible(true);
        return candidate.invoke(component, arguments);
      } catch (ReflectiveOperationException ignored) {
        return "unreadable";
      }
    }
    return "no such method";
  }
}
