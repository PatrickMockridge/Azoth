// `ReactiveMultiphasePHflash`, for `reactions.reactive_ph_flash`.
//
// The class looks for the temperature at which a reactive fluid's enthalpy matches a
// specification, and it does that by wrapping the reactive TP flash in an outer loop on T.
// **Its own docstring calls that loop a `1/T` Newton following Michelsen 1987; the code is a
// secant with a bisection fallback on T itself**, and this probe is what settles which one it
// is: the iteration count, the residual and the temperature it lands on are all printed.
//
// The state is the class's own test's: the water-gas shift at 600 K and 1 bar, brought to
// reactive equilibrium, its enthalpy recorded, then the temperature perturbed to 500 K and the
// PH flash asked to find its way back. The round trip is the oracle - a temperature that comes
// back is a temperature the loop found rather than a number it reported.
//
//     javac -proc:none -cp neqsim-f0c7436.jar ReactivePhFlashProbe.java
//     java -cp .:neqsim-f0c7436.jar ReactivePhFlashProbe > captures/reactive_ph_flash_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.flashops.reactiveflash.ReactiveMultiphasePHflash;
import neqsim.thermodynamicoperations.flashops.reactiveflash.ReactiveMultiphaseTPflash;

public class ReactivePhFlashProbe {

  public static void main(String[] args) {
    roundTrip("wgs-600K-from-500K", 600.0, 500.0, 1.0,
        new String[] { "CO", "water", "CO2", "hydrogen" },
        new double[] { 0.25, 0.25, 0.25, 0.25 });
    roundTrip("wgs-1000K-from-700K", 1000.0, 700.0, 1.0,
        new String[] { "CO", "water", "CO2", "hydrogen" },
        new double[] { 0.25, 0.25, 0.25, 0.25 });
    roundTrip("methane-water-co2-hydrogen-1000K-from-800K", 1000.0, 800.0, 1.0,
        new String[] { "methane", "water", "CO2", "hydrogen" },
        new double[] { 0.4, 0.2, 0.2, 0.2 });
  }

  /// One round trip: flash at `equilibriumTemperature`, record the enthalpy, perturb the
  /// temperature, and ask the PH flash for the temperature back.
  static void roundTrip(String label, double equilibriumTemperature, double perturbed,
      double pressure, String[] names, double[] moles) {
    System.out.println("fluid=" + label);
    System.out.println("flash_temperature_K=" + equilibriumTemperature);
    System.out.println("perturbed_temperature_K=" + perturbed);
    System.out.println("pressure_bara=" + pressure);

    SystemInterface system = build(equilibriumTemperature, pressure, names, moles);
    ReactiveMultiphaseTPflash tpFlash = new ReactiveMultiphaseTPflash(system);
    tpFlash.run();
    system.init(2);

    double enthalpySpec = system.getEnthalpy();
    System.out.println("specified_sensible_enthalpy_J=" + enthalpySpec);
    System.out.println("tp_flash_converged=" + tpFlash.isConverged());
    System.out.println("tp_flash_iterations=" + tpFlash.getTotalIterations());
    printComposition("tp_x", system);
    System.out.println("tp_cp_J_per_K=" + system.getCp());
    System.out.println("tp_entropy_J_per_K=" + system.getEntropy());

    // The perturbed state the loop starts from: the same composition at a different
    // temperature, which is what the class's own test does.
    system.setTemperature(perturbed);

    ReactiveMultiphasePHflash phFlash = new ReactiveMultiphasePHflash(system, enthalpySpec);
    printPrivate("formation_inventory_at_construction", phFlash,
        "getFormationEnthalpyInventory");
    printPrivate("thermochemical_enthalpy_spec", phFlash, "getThermochemicalEnthalpySpec");
    phFlash.run();

    System.out.println("converged=" + phFlash.isConverged());
    System.out.println("equilibrium_temperature_K=" + phFlash.getEquilibriumTemperature());
    System.out.println("outer_iterations=" + phFlash.getOuterIterations());
    System.out.println("total_inner_iterations=" + phFlash.getTotalInnerIterations());
    system.init(2);
    System.out.println("final_sensible_enthalpy_J=" + system.getEnthalpy());
    printPrivate("final_thermochemical_enthalpy", phFlash, "getThermochemicalEnthalpy");
    printPrivate("final_formation_inventory", phFlash, "getFormationEnthalpyInventory");
    System.out.println("final_entropy_J_per_K=" + system.getEntropy());
    System.out.println("final_cp_J_per_K=" + system.getCp());
    printComposition("final_x", system);
    System.out.println();
  }

  static SystemInterface build(double temperature, double pressure, String[] names, double[] moles) {
    SystemInterface system = new SystemSrkEos(temperature, pressure);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.setMixingRule("classic");
    return system;
  }

  static void printComposition(String label, SystemInterface system) {
    StringBuilder out = new StringBuilder(label + "=");
    for (int i = 0; i < system.getPhase(0).getNumberOfComponents(); i++) {
      ComponentInterface component = system.getPhase(0).getComponent(i);
      out.append(component.getName()).append(":").append(component.getx()).append(" ");
    }
    System.out.println(out.toString().trim());
  }

  /// A private reader the class never exposes, reached by reflection rather than recomputed:
  /// the thermochemical enthalpy is what the outer loop's residual is measured against, and
  /// the formation inventory is the term that separates it from the system's own.
  static void printPrivate(String label, Object target, String method) {
    try {
      java.lang.reflect.Method reader = target.getClass().getDeclaredMethod(method);
      reader.setAccessible(true);
      System.out.println(label + "=" + reader.invoke(target));
    } catch (ReflectiveOperationException ex) {
      System.out.println(label + "_failed=" + ex);
    }
  }
}
