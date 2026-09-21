// The two hydrate component models' layers, printed side by side so a port can tell them apart.
//
// `PhaseHydrate` picks its component class from the *system's model name*: the `CPA-SRK-EOS`
// family builds `ComponentHydrateGF` and everything else `ComponentHydratePVTsim`. That is not
// a property of the fluid, so the state is reached here by constructing the phase the way the
// dispatch does - `new PhaseHydrate("CPA-SRK-EOS")` - rather than by asking NeqSim for it,
// which would need a CPA system this port does not carry.
//
// **The state is a plain SRK fluid on purpose.** Both models read their guests' fugacities from
// the phase they are attached to, so the hydrate arithmetic is checked without also checking a
// mixing rule. The guests' reference fugacities are fed in exactly as
// `HydrateFormationTemperatureFlash` feeds them (`setRefFug(j, gas.getFugacity(j))`).
//
// **What each key is for.** `empty_vapour_pressure` and `molar_volume_hydrate` are the
// reference term the Guo-Finch route multiplies in and the PVTsim route does not, and they are
// printed per structure because that is the whole of the difference: the PVTsim reference is
// one number for both structures, so its comparison of the two is a comparison of the
// exponents, and the Guo-Finch one is not. `cki` and `yki` are the occupancy chain, which is
// where the two models' fitted pairs diverge, and `fugcoef` is each model's own answer -
// the lower of its two structures.
//
// A capture, not a test: run it against the pinned jar and commit what it prints.
//   javac -proc:none -cp neqsim-f0c7436.jar HydrateGFProbe.java
//   java -cp .:neqsim-f0c7436.jar HydrateGFProbe > captures/hydrate_gf_probe.tsv

import neqsim.thermo.component.ComponentHydrate;
import neqsim.thermo.phase.PhaseHydrate;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class HydrateGFProbe {

  /** NeqSim's own hydrate-equilibrium feed, as `HydrateProbe` uses it. */
  private static final String[] NAMES = {"methane", "ethane", "propane", "water"};

  private static final double[] MOLES = {79.0, 10.10, 2.050, 10.0};

  /** A `key = value` row, which is the shape `tools/neqsim_layer_diff.py` reads. */
  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    for (double pressureBara : new double[] {100.0, 50.0, 200.0}) {
      report(pressureBara);
    }
  }

  private static void report(double pressureBara) {
    System.out.printf("# methane/ethane/propane/water at P = %.15g bara, both models%n",
        pressureBara);
    try {
      SystemInterface fluid = new SystemSrkEos(273.15, pressureBara);
      for (int i = 0; i < NAMES.length; i++) {
        fluid.addComponent(NAMES[i], MOLES[i]);
      }
      fluid.setMixingRule(2);
      fluid.setMultiPhaseCheck(true);
      fluid.init(0);
      fluid.init(1);
      new ThermodynamicOperations(fluid).TPflash();

      // The guests' fugacities, from the phase the hydrate's own flash reads them off.
      double[] fugacities = new double[NAMES.length];
      for (int j = 0; j < NAMES.length; j++) {
        fugacities[j] = fluid.getPhase(0).getFugacity(j);
      }
      row("phases", fluid.getNumberOfPhases());
      for (int j = 0; j < NAMES.length; j++) {
        row("f[" + NAMES[j] + "]", fugacities[j]);
      }

      for (String model : new String[] {"pvtsim", "guo_finch"}) {
        PhaseHydrate hydrate = model.equals("guo_finch")
            ? new PhaseHydrate("CPA-SRK-EOS")
            : new PhaseHydrate();
        for (int i = 0; i < NAMES.length; i++) {
          hydrate.addComponent(NAMES[i], MOLES[i], MOLES[i], i);
        }
        hydrate.setTemperature(fluid.getTemperature());
        hydrate.setPressure(fluid.getPressure());
        // **The reference fluid phase, which `SystemThermo` gives the hydrate as it builds it**
        // (`SystemThermo.java:632` passes `phaseArray[0]`). The two models differ in what they
        // do with it - `ComponentHydratePVTsim` reads its water fugacity and `ComponentHydrateGF`
        // does not read it at all off the ice branch - but both refuse to run without it.
        hydrate.setSolidRefFluidPhase(fluid.getPhase(0));
        for (int i = 0; i < NAMES.length; i++) {
          if (hydrate.getComponent(i) instanceof ComponentHydrate component) {
            for (int j = 0; j < NAMES.length; j++) {
              component.setRefFug(j, fugacities[j]);
            }
          }
        }

        // **The component class is the model**, and printing it is what makes the rest of the
        // block readable: a run that silently built the same class twice would agree with
        // itself on every key below.
        System.out.printf("# model = %s%n", model);
        System.out.printf("# component_class = %s%n",
            hydrate.getComponent(0).getClass().getSimpleName());

        for (int structure = 0; structure < 2; structure++) {
          row("empty_vapour_pressure[" + structure + "]",
              ((ComponentHydrate) hydrate.getComponent(0))
                  .getEmptyHydrateStructureVapourPressure(structure, hydrate.getTemperature()));
          row("molar_volume_hydrate[" + structure + "]",
              ((ComponentHydrate) hydrate.getComponent(0))
                  .getMolarVolumeHydrate(structure, hydrate.getTemperature()));
          for (int cavity = 0; cavity < 2; cavity++) {
            row("cavprwat[" + structure + "," + cavity + "]",
                ((ComponentHydrate) hydrate.getComponent(0)).getCavprwat(structure, cavity));
          }
        }

        for (int i = 0; i < NAMES.length; i++) {
          if (!(hydrate.getComponent(i) instanceof ComponentHydrate component)) {
            continue;
          }
          for (int structure = 0; structure < 2; structure++) {
            for (int cavity = 0; cavity < 2; cavity++) {
              String key = NAMES[i] + "," + structure + "," + cavity;
              row("cki[" + key + "]", component.calcCKI(structure, cavity, hydrate));
              row("yki[" + key + "]", component.calcYKI(structure, cavity, hydrate));
            }
          }
        }
        row("fugcoef",
            ((ComponentHydrate) hydrate.getComponent("water"))
                .fugcoef(hydrate, hydrate.getNumberOfComponents(), hydrate.getTemperature(),
                    hydrate.getPressure()));
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
