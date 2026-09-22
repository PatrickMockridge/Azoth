// How NeqSim labels a phase, printed a layer at a time so the rule can be ported.
//
// `PhaseEos.init` ends with three branches and one number:
//
//   if (getVolume() / getB() > 1.75)        GAS
//   else if (sumHydrocarbons > sumAqueous)  OIL
//   else                                    AQUEOUS
//
// where the two sums are over `isHydrocarbon() || isInert() || isIsTBPfraction()`, with water
// excluded from the first by name, and everything else counted aqueous. **Nothing states what
// units `getVolume()` and `getB()` are in**, so this prints both and their ratio rather than
// leaving the threshold's scale to be inferred - and it prints the sums so the second branch
// is checkable where the first does not take it.
//
// The fluids are chosen to land in each branch: a gas, an oil, an aqueous liquid, and a
// two-phase one whose two phases differ.
//
// A capture, not a test: run it against the pinned jar and commit what it prints.
//   javac -proc:none -cp neqsim-f0c7436.jar PhaseTypeProbe.java
//   java -cp .:neqsim-f0c7436.jar PhaseTypeProbe > captures/phase_type_probe.tsv

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class PhaseTypeProbe {

  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  public static void main(String[] args) {
    // (name, components, moles, T in K, P in bara)
    report("methane gas", new String[] {"methane"}, new double[] {1.0}, 300.0, 50.0);
    report("lng liquid", new String[] {"methane", "ethane", "propane"},
        new double[] {0.5, 0.3, 0.2}, 150.0, 50.0);
    report("meg and water", new String[] {"MEG", "water"}, new double[] {0.1, 1.0}, 300.0, 1.0);
    report("inhibitor feed", new String[] {"methane", "ethane", "propane", "i-butane", "MEG", "water"},
        new double[] {1.0, 0.10, 0.050, 0.0050, 0.1, 1.0}, 273.15, 100.0);
    report("gas and aqueous", new String[] {"methane", "water"}, new double[] {0.5, 0.5}, 300.0, 10.0);
    report("co2 liquid", new String[] {"CO2"}, new double[] {1.0}, 260.0, 60.0);
  }

  private static void report(String label, String[] names, double[] moles, double t, double pBara) {
    System.out.printf("# %s at %.15g K, %.15g bara%n", label, t, pBara);
    try {
      SystemInterface fluid = new SystemSrkEos(t, pBara);
      for (int i = 0; i < names.length; i++) {
        fluid.addComponent(names[i], moles[i]);
      }
      fluid.setMixingRule(2);
      fluid.setMultiPhaseCheck(true);
      fluid.init(0);
      fluid.init(1);
      new ThermodynamicOperations(fluid).TPflash();

      row("phases", fluid.getNumberOfPhases());
      for (int i = 0; i < fluid.getNumberOfPhases(); i++) {
        PhaseInterface phase = fluid.getPhase(i);
        row("phase[" + i + "].beta", phase.getBeta());
        row("phase[" + i + "].z", phase.getZ());
        // **Both of the ratio's halves**, because nothing says what units they are in.
        row("phase[" + i + "].volume", phase.getVolume());
        row("phase[" + i + "].molar_volume", phase.getMolarVolume());
        row("phase[" + i + "].moles", phase.getNumberOfMolesInPhase());
        row("phase[" + i + "].b", phase.getB());
        row("phase[" + i + "].volume_over_b", phase.getVolume() / phase.getB());

        double hydrocarbons = 0.0;
        double aqueous = 0.0;
        for (int j = 0; j < phase.getNumberOfComponents(); j++) {
          String name = phase.getComponent(j).getName();
          boolean isHydrocarbon = (phase.getComponent(j).isHydrocarbon()
              || phase.getComponent(j).isInert()
              || phase.getComponent(j).isIsTBPfraction())
              && !name.equals("water") && !name.equals("water_PC");
          if (isHydrocarbon) {
            hydrocarbons += phase.getComponent(j).getx();
          } else {
            aqueous += phase.getComponent(j).getx();
          }
        }
        row("phase[" + i + "].sum_hydrocarbons", hydrocarbons);
        row("phase[" + i + "].sum_aqueous", aqueous);
        System.out.printf("phase[%d].type = %s%n", i, phase.getType());
      }
    } catch (Exception error) {
      System.out.printf("# failed: %s: %s%n", error.getClass().getSimpleName(), error.getMessage());
    }
    System.out.println();
  }
}
