// The wax phase's component models, run side by side on one fluid so a port can tell them apart.
//
// `PhaseWax` builds its component class from `waxComponentModelName` inside `addComponent`, so
// the model has to be chosen **before the wax phase exists**. `SystemThermo.setWaxModelType`
// says "If wax phase already exists, update it" and then calls `PhaseWax.setWaxComponentModel`,
// which sets the name and nothing else - the components stay `ComponentWax`, so a model set
// after the phase is built is a silent no-op and the phase answers with Pedersen's numbers
// under the other model's name. This probe sets it first, and prints the component class so a
// reader can see which one answered.
//
// **What it measured, and it is why there is nothing to port.** On this fluid the wax fraction
// is `0.146` at 275 K under Pedersen and **`1.0`** under Won, Wilson and Coutinho - the whole
// feed as wax, methane included, which is not a state. Two of the three do that by returning
// `NaN`, and the third by returning `8.0e+252`; the rows below say which layer went first.
// `ComponentWonWax.getWonParam` is `sqrt` of a negative at every temperature a wax exists at,
// so its activity coefficient and its coefficient are `NaN`; Wilson's is `NaN` the same way,
// and Coutinho's `ln gamma` is finite while `exp(thermal + heat-capacity + Poynting)` is not.
//
// **What each key is for.** `wax_volume` and `wax_param` are the Won model's two inputs and
// are printed for every model, because `ComponentWonWax` builds them from correlations of its
// own rather than reading the databank; `wax_gamma` is its activity coefficient, which is the
// whole of the difference between the two coefficients; `wax_phi` is each model's own answer.
// A port whose `wax_phi` is out while its `wax_gamma` agrees has the fusion term wrong, and one
// whose `wax_gamma` is out has the solubility parameter wrong.
//
// The fluid is `WaxProbe`'s - a characterisation this library cannot build, which is why the
// components are printed too.
//
// A capture, not a test: run it against the pinned jar and commit what it prints.
//   javac -proc:none -cp neqsim-f0c7436.jar WaxModelProbe.java
//   java -cp .:neqsim-f0c7436.jar WaxModelProbe > captures/wax_model_probe.tsv

import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.component.ComponentWonWax;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;
import neqsim.util.database.NeqSimDataBase;

public class WaxModelProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    for (String model : new String[] {"Pedersen", "Won", "Wilson", "Coutinho"}) {
      for (double temperatureK : new double[] {275.0, 261.0}) {
        report(model, temperatureK, 5.0);
      }
    }
  }

  /** The probe's fluid, with the wax model chosen before the phase is built. */
  private static SystemInterface build(String model) {
    NeqSimDataBase.setCreateTemporaryTables(true);
    SystemInterface system = new SystemSrkEos(298.0, 10.0);
    system.addComponent("methane", 6.78);
    system.addTBPfraction("C19", 10.13, 170.0 / 1000.0, 0.7814);
    system.addPlusFraction("C20", 10.62, 381.0 / 1000.0, 0.850871882888);
    system.getCharacterization().characterisePlusFraction();
    system.getWaxModel().addTBPWax();
    system.createDatabase(true);
    system.setMixingRule(2);
    // **The model, and it has to be here.** `addSolidComplexPhase` is what calls the wax
    // phase's `addComponent`, and that is where the class is chosen.
    system.setWaxModelType(model);
    NeqSimDataBase.setCreateTemporaryTables(true);
    system.addSolidComplexPhase("wax");
    system.setMultiphaseWaxCheck(true);
    system.setMultiPhaseCheck(true);
    NeqSimDataBase.setCreateTemporaryTables(false);
    system.init(0);
    system.init(1);
    return system;
  }

  private static void report(String model, double temperatureK, double pressureBara) {
    System.out.printf("# %s, T = %.15g K, P = %.15g bara%n", model, temperatureK, pressureBara);
    try {
      SystemInterface fluid = build(model);
      fluid.setTemperature(temperatureK);
      fluid.setPressure(pressureBara);
      new ThermodynamicOperations(fluid).TPflash();

      PhaseInterface wax = null;
      for (PhaseInterface phase : fluid.getPhases()) {
        if (phase != null && phase.getType() == PhaseType.WAX) {
          wax = phase;
        }
      }
      if (wax == null) {
        System.out.println("note = no wax phase in the system");
        System.out.println();
        return;
      }

      // **Which class answered**, which is the whole point of the file: a run that built
      // `ComponentWax` for both models would agree with itself on every key below.
      System.out.printf("# wax_component_class = %s%n",
          wax.getComponent(0).getClass().getSimpleName());
      row("wax_beta", wax.getBeta());
      row("phases", fluid.getNumberOfPhases());

      for (int i = 0; i < wax.getNumberOfComponents(); i++) {
        ComponentInterface component = wax.getComponent(i);
        String name = component.getName();
        boolean former = component.isWaxFormer();
        row("wax_x[" + name + "]", component.getx());
        row("wax_former[" + name + "]", former ? 1.0 : 0.0);
        if (!former) {
          continue;
        }
        row("wax_molar_mass_g[" + name + "]", component.getMolarMass() * 1000.0);
        row("wax_tc[" + name + "]", component.getTC());
        row("wax_heat_of_fusion[" + name + "]", component.getHeatOfFusion());
        row("wax_triple_point[" + name + "]", component.getTriplePointTemperature());
        if (component instanceof ComponentWonWax won) {
          row("wax_volume[" + name + "]", won.getWonVolume(wax));
          row("wax_param[" + name + "]", won.getWonParam(wax));
          row("wax_gamma[" + name + "]", won.getWonActivityCoefficient(wax));
        }
        if (component instanceof neqsim.thermo.component.ComponentWaxWilson wilson) {
          row("wax_gamma[" + name + "]", wilson.getWilsonActivityCoefficient(wax));
        }
        if (component instanceof neqsim.thermo.component.ComponentCoutinhoWax coutinho) {
          row("wax_ln_gamma[" + name + "]", coutinho.calcLnGammaUNIQUAC(wax));
        }
        row("wax_phi[" + name + "]", component.getFugacityCoefficient());
        // **Called directly**, because the flash does not evaluate the wax phase's own
        // coefficients: its field stays at its zero default and a model that returns NaN
        // would look like one that returns nothing.
        row("wax_phi_direct[" + name + "]", component.fugcoef(wax));
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
