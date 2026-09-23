// The Pitzer phase's fugacity coefficients, for `eos.pitzer_phase`.
//
// `ComponentGePitzer.fugcoef` has three branches and the port carries none of them: this
// model is the **activity-coefficient surface only**, and its spec says so. The branches are
//
//   water (or a solvent)          gamma_i * P0_i / P
//   a neutral that is not water   gamma_i * H_m * (m_i / x_i) / P, with H_m Henry on the
//                                 molality scale - `getEffectiveHenryCoefficient` converts
//                                 the bar-scale coefficient to bar kg/mol
//   an ion                        activinf * H_capped / P, with `H_capped = 1e12` because
//                                 `isHenryCoefficientCapped` includes `isIsIon()`, and
//                                 `activinf = gamma / getActivityCoefficientInfDilWater`
//
// The second and third need `PhaseGE.getActivityCoefficientInfDilWater`, the one piece of
// the GE surface this library does not have. `getEffectiveHenryCoefficient` and
// `isHenryCoefficientCapped` are `protected`, so they are read through reflection up the
// class hierarchy - the alternative would be re-deriving them and checking a derivation
// against itself.
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
            + " henry_capped_1=" + capped(component, 1.0));
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
