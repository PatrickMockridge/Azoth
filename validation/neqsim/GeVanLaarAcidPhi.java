// What `PhaseGEVanLaarAcid` reports for a liquid, for `eos.ge_van_laar_acid_phase`.
//
//     javac -proc:none -cp neqsim-f0c7436.jar GeVanLaarAcidPhi.java
//     java -cp .:neqsim-f0c7436.jar GeVanLaarAcidPhi
//
// **Why a bare phase and not `SystemVanLaarActivitySRK`.** That class carries a great
// deal of tuning, and all of it is vapour-side or flash-iteration: tuned HNO3 SRK
// criticals, CO2-HNO3 and CO2-H2SO4 kij fits, initial K seeds, damped K relaxation, a
// nitric-acid carrier K clamp, trace-phase collapse thresholds. **None of it enters the
// liquid**: `ComponentGEVanLaarAcid.computeGamma` and `pureVaporPressureBar` read only
// `NitricSulfuricAcidVaporPressure` and the phase's temperature and acid composition.
// They are also both private, and the system exposes no gamma or P0 accessor.
//
// So the tuned system is passed over not because it would contaminate the liquid but
// because there is nothing to read from it. A bare `PhaseGEVanLaarAcid` is directly
// constructible, `addComponent` gives it `ComponentGEVanLaarAcid` instances, and
// `setAlpha`/`setDij`/`setDijT` throw deliberately because the model uses fixed
// correlation parameters. Driving it gives the same liquid the system would.
//
// Unlike `PhaseGEWilson` and `PhaseGEUniquac`, this component *does* publish gamma: its
// eight-argument `getGamma` calls `computeGamma` and stores the field, so
// `getExcessGibbsEnergy` works and `fugcoef` reads a real value.
//
// The P0 the phase uses is private, so the driver takes it from the same static utility
// the phase calls - which is also the oracle for `eos.nitric_sulfuric_acid_vapor_pressure`.
// What is NeqSim's here is gamma, P0 and the composition; what is the driver's is the one
// multiplication and the division by pressure, stated rather than hidden.

import neqsim.thermo.component.ComponentGEVanLaarAcid;
import neqsim.thermo.phase.PhaseGEInterface;
import neqsim.thermo.phase.PhaseGEVanLaarAcid;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.util.empiric.NitricSulfuricAcidVaporPressure;

public class GeVanLaarAcidPhi {

  static void one(String[] names, double[] x, double temperatureK, double pressureBar) {
    PhaseGEVanLaarAcid phase = new PhaseGEVanLaarAcid();
    phase.setTemperature(temperatureK);
    phase.setPressure(pressureBar);
    for (int i = 0; i < names.length; i++) {
      phase.addComponent(names[i], x[i], x[i], i);
    }
    for (int i = 0; i < names.length; i++) {
      phase.getComponent(i).setx(x[i]);
    }

    ((PhaseGEInterface) phase).getExcessGibbsEnergy(phase, names.length, temperatureK,
        pressureBar, PhaseType.LIQUID);

    System.out.println("T=" + temperatureK + " P=" + pressureBar
        + "  x=" + java.util.Arrays.toString(x));
    for (int i = 0; i < names.length; i++) {
      ComponentGEVanLaarAcid component = (ComponentGEVanLaarAcid) phase.getComponent(i);
      component.fugcoef(phase);
      // The same three functions `pureVaporPressureBar` chooses between, in pascal.
      double p0Pa = switch (component.getName().toLowerCase().trim()) {
        case "water", "h2o" -> NitricSulfuricAcidVaporPressure.pureVaporPressureWater(temperatureK);
        case "nitric acid", "hno3" ->
          NitricSulfuricAcidVaporPressure.pureVaporPressureNitricAcid(temperatureK);
        case "sulfuric acid", "sulphuric acid", "h2so4" ->
          NitricSulfuricAcidVaporPressure.pureVaporPressureSulfuricAcid(temperatureK);
        default -> Double.NaN;
      };
      System.out.println("  " + component.getName()
          + "  x=" + component.getx()
          + "  gamma=" + component.getGamma()
          + "  P0_Pa=" + p0Pa
          + "  phi=" + component.getFugacityCoefficient());
    }
  }

  public static void main(String[] args) {
    one(new String[] {"water", "nitric acid"}, new double[] {0.6, 0.4}, 250.0, 1.0);
    one(new String[] {"water", "nitric acid", "sulfuric acid"},
        new double[] {0.5, 0.3, 0.2}, 250.0, 1.0);
    one(new String[] {"water", "nitric acid", "sulfuric acid"},
        new double[] {0.7, 0.1, 0.2}, 273.15, 2.0);
  }
}
