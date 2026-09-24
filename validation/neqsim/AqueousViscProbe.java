// `eos.aqueous_viscosity`: the liquid viscosity NeqSim gives an **aqueous** phase.
//
// **The phase type is the whole point.** `getPhysicalProperties()` dispatches on it:
// `PhaseType.AQUEOUS` takes `WaterPhysicalProperties`, whose viscosity model is the liquid
// `Viscosity` class - the `polynom` one - where a gas or a hydrocarbon liquid takes
// `PFCTViscosityMethodHeavyOil`, which `eos.viscosity` ports. So the same water at the same
// state has two answers depending on which class the phase reaches, and this captures the
// aqueous one.
//
// # What the rows are for
//
// Water at five temperatures and two pressures, which is what pins the model *form* and the
// pressure correction separately; methanol and TEG, whose `LIQVISC` sets are their own; and
// a brine, whose ions are not.
//
// # The two things the capture carries that the numbers do not
//
// **`Gij` is zero everywhere.** `PhysicalPropertyMixingRule.initMixingRules` fills it with
// an inner loop that starts at `k = l` and breaks on `k == l`, so its body never runs - and
// the water/methanol row is printed with the whole matrix so that the zero is visible on a
// mixture whose `gijvisc` column is not empty. Grunberg-Nissan therefore reduces to
// `exp(sum_i w_i ln mu_i)` over mass fractions.
//
// **The ion rows are methanol's.** `na+` and `cl-` carry `-39.35, 4826.0, 0.1091, -1.13e-4`
// in NeqSim's own COMP.csv, which is methanol's set, and the brine row comes out *less*
// viscous than pure water - `8.4026e-4` against `8.5510e-4` - which is the wrong direction.
// That is why `eos.aqueous_viscosity` refuses a phase carrying an ion.
//
//     javac -proc:none -cp neqsim-f0c7436.jar AqueousViscProbe.java
//     java -cp .:neqsim-f0c7436.jar AqueousViscProbe > captures/aqueous_viscosity_probe.tsv

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class AqueousViscProbe {
  public static void main(String[] args) {
    row("water_300K_1bar", new String[] { "water" }, new double[] { 1.0 }, 300.0, 1.0);
    row("water_300K_5bar", new String[] { "water" }, new double[] { 1.0 }, 300.0, 5.0);
    row("water_280K_5bar", new String[] { "water" }, new double[] { 1.0 }, 280.0, 5.0);
    row("water_320K_5bar", new String[] { "water" }, new double[] { 1.0 }, 320.0, 5.0);
    row("water_350K_5bar", new String[] { "water" }, new double[] { 1.0 }, 350.0, 5.0);
    row("methanol_300K_5bar", new String[] { "methanol" }, new double[] { 1.0 }, 300.0, 5.0);
    row("teg_300K_5bar", new String[] { "TEG" }, new double[] { 1.0 }, 300.0, 5.0);
    row("water_methanol_300K_5bar", new String[] { "water", "methanol" },
        new double[] { 0.5, 0.5 }, 300.0, 5.0);
    row("brine_nacl_300K_5bar", new String[] { "water", "Na+", "Cl-" },
        new double[] { 0.98, 0.01, 0.01 }, 300.0, 5.0);
  }

  static void row(String label, String[] names, double[] z, double t, double p) {
    SystemInterface f = new SystemPrEos(t, p);
    for (int i = 0; i < names.length; i++) {
      f.addComponent(names[i], z[i]);
    }
    f.setMixingRule(2);
    f.setTotalFlowRate(1.0, "mol/sec");
    f.init(0);
    new ThermodynamicOperations(f).TPflash();
    f.initProperties();

    System.out.println(label);
    System.out.println("T_K=" + f.getPhase(0).getTemperature());
    System.out.println("P_bara=" + f.getPhase(0).getPressure("bara"));
    System.out.println("phase_type=" + f.getPhase(0).getType());
    System.out.println("viscosity_model="
        + f.getPhase(0).getPhysicalProperties().getViscosityModel().getClass().getSimpleName());
    System.out.println("viscosity_Pa_s=" + f.getPhase(0).getPhysicalProperties().getViscosity());

    var model = f.getPhase(0).getPhysicalProperties().getViscosityModel();
    if (model instanceof neqsim.physicalproperties.methods.liquidphysicalproperties.viscosity.Viscosity) {
      var v = (neqsim.physicalproperties.methods.liquidphysicalproperties.viscosity.Viscosity) model;
      v.calcPureComponentViscosity();
      for (int i = 0; i < f.getPhase(0).getNumberOfComponents(); i++) {
        String name = f.getPhase(0).getComponent(i).getName();
        System.out.println("pure_" + name + "_cP=" + v.getPureComponentViscosity(i));
        System.out.println("weight_" + name + "=" + f.getPhase(0).getWtFrac(i));
        System.out.println("liqviscmodel_" + name + "="
            + f.getPhase(0).getComponent(i).getLiquidViscosityModel());
        System.out.println("correction_" + name + "=" + v.getViscosityPressureCorrection(i));
      }
      // The whole matrix, so that its zero is in the capture rather than in a comment.
      int n = f.getPhase(0).getNumberOfComponents();
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          System.out.println("gij_" + i + "_" + j + "="
              + f.getPhase(0).getPhysicalProperties().getMixingRule().getViscosityGij(i, j));
        }
      }
    }
    System.out.println();
  }
}
