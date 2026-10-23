// TbpFractionComponentProbe - what NeqSim's `addTBPfraction` makes of a cut, at the pin.
//
// One state per (molar mass, specific gravity): the component's own getters after the call, and
// the intermediates that are not on the component - the model's own acentric factor and `m`, and
// the reference system's molar volume that `calcRacketZ` reads. The volume is printed because the
// Rackett compressibility is not a function of the cut's two numbers: it is the Peneloux shift of
// a one-component SRK phase at 15 C and one atmosphere, and a port that skipped that solve would
// state a plausible Z that is not this one.
//
//   ../../tools/gated.sh javac -proc:none -cp neqsim-f0c7436.jar TbpFractionComponentProbe.java
//   ../../tools/gated.sh java -cp .:neqsim-f0c7436.jar TbpFractionComponentProbe > captures/tbp_fraction_component_probe.tsv

import java.util.Locale;
import neqsim.thermo.ThermodynamicConstantsInterface;
import neqsim.thermo.characterization.TBPfractionModel;
import neqsim.thermo.component.ComponentInterface;
import neqsim.thermo.phase.PhaseType;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkEos;

public class TbpFractionComponentProbe {

  private static final String COLUMNS = "label|molar_mass|density|model_tc|model_pc|model_tb|"
      + "model_acentric_factor|model_m|reference_molar_volume|model_racket_z|component_tc|"
      + "component_pc|component_tb|component_acentric_factor|component_m|critical_volume|"
      + "racket_z|racket_z_cpa|parachor|critical_viscosity|triple_point_temperature|"
      + "heat_of_fusion|ideal_gas_enthalpy_of_formation|cp_a|cp_b|cp_c|cp_d|"
      + "lennard_jones_diameter|lennard_jones_energy|normal_liquid_density";

  private static final String UNITS = "string|kg/mol|kg/m3|K|Pa|K|1|1|m3/mol|1|K|Pa|K|1|1|"
      + "m3/mol|1|1|1|Pa*s|K|J/mol|J/mol|J/(mol*K)|J/(mol*K2)|J/(mol*K3)|J/(mol*K4)|angstrom|K|kg/m3";

  private static String number(double value) {
    return Double.isFinite(value) ? String.format(Locale.ROOT, "%.17g", value) : "-";
  }

  /** One cut through `addTBPfraction`, the class's own route. */
  private static ComponentInterface component(double molarMass, double density) {
    SystemInterface system = new SystemSrkEos(298.15, 1.0);
    system.addTBPfraction("CUT", 1.0, molarMass, density);
    system.setMixingRule(2);
    system.init(0);
    return system.getComponent(0);
  }

  /**
   * The reference system `addTBPfraction` flashes, and the model calls it makes.
   *
   * Reproduced line for line from the class: a one-component SRK system at 15 C and one
   * atmosphere, its component overwritten with the cut's constants and the fitted `m`, one
   * liquid phase, `init(1)`. Returns `{molar volume [m3/mol], Z, model Rackett Z, model m,
   * model acentric factor, model Tc, model Pc, model Tb}`.
   */
  private static double[] reference(double molarMass, double density) {
    TBPfractionModel model = new TBPfractionModel();
    TBPfractionModel.PedersenTBPModelSRK tbp = model.new PedersenTBPModelSRK();
    double molarMassGrams = molarMass * 1000.0;
    double tc = tbp.calcTC(molarMassGrams, density);
    double pc = tbp.calcPC(molarMassGrams, density);
    double m = tbp.calcm(molarMassGrams, density);
    double tb = tbp.calcTB(molarMassGrams, density);
    double acs = tbp.calcAcentricFactor(molarMassGrams, density);

    SystemInterface reference = new SystemSrkEos();
    reference.setTemperature(273.15 + 15.0);
    reference.setPressure(ThermodynamicConstantsInterface.referencePressure);
    reference.addComponent("default", 1.0, 273.15, 50.0, 0.1);
    reference.setMixingRule(1);
    reference.init(0);
    reference.setNumberOfPhases(1);
    reference.setPhaseType(0, PhaseType.LIQUID);
    ComponentInterface cut = reference.getPhase(0).getComponent(0);
    cut.setComponentName("CUT_PC");
    cut.setMolarMass(molarMass);
    cut.setAcentricFactor(acs);
    cut.setTC(tc);
    cut.setPC(pc);
    cut.setComponentType("TBPfraction");
    cut.setIsTBPfraction(true);
    cut.getAttractiveTerm().setm(m);
    acs = cut.getAcentricFactor();
    reference.setTemperature(273.15 + 15.0);
    reference.setPressure(ThermodynamicConstantsInterface.referencePressure);
    reference.init(1);

    double volume = reference.getPhase(0).getMolarVolume();
    double z = reference.getPhase(0).getZ();
    double racket = tbp.calcRacketZ(reference, molarMassGrams, density);
    return new double[] {volume, z, racket, m, acs, tc, pc, tb};
  }

  private static void state(String label, double molarMass, double density) {
    double[] ref = reference(molarMass, density);
    ComponentInterface c = component(molarMass, density);
    StringBuilder text = new StringBuilder(label).append('|');
    text.append(number(molarMass)).append('|').append(number(density * 1000.0)).append('|');
    text.append(number(ref[5])).append('|').append(number(ref[6] * 1.0e5)).append('|');
    text.append(number(ref[7])).append('|');
    text.append(number(ref[4])).append('|').append(number(ref[3])).append('|');
    text.append(number(ref[0])).append('|').append(number(ref[2])).append('|');
    text.append(number(c.getTC())).append('|').append(number(c.getPC() * 1.0e5)).append('|');
    text.append(number(c.getNormalBoilingPoint())).append('|');
    text.append(number(c.getAcentricFactor())).append('|');
    text.append(number(c.getAttractiveTerm() == null ? Double.NaN : c.getAttractiveTerm().getm()))
        .append('|');
    text.append(number(c.getCriticalVolume() / 1.0e6)).append('|');
    text.append(number(c.getRacketZ())).append('|').append(number(c.getRacketZCPA())).append('|');
    text.append(number(c.getParachorParameter())).append('|');
    text.append(number(c.getCriticalViscosity())).append('|');
    text.append(number(c.getTriplePointTemperature())).append('|');
    text.append(number(c.getHeatOfFusion())).append('|');
    text.append(number(c.getIdealGasEnthalpyOfFormation())).append('|');
    text.append(number(c.getCpA())).append('|').append(number(c.getCpB())).append('|');
    text.append(number(c.getCpC())).append('|').append(number(c.getCpD())).append('|');
    text.append(number(c.getLennardJonesMolecularDiameter())).append('|');
    text.append(number(c.getLennardJonesEnergyParameter())).append('|');
    text.append(number(c.getNormalLiquidDensity() * 1000.0));
    System.out.println(text);
  }

  public static void main(String[] args) {
    System.out.println("TBP fraction component at f0c7436c69");
    System.out.println("# columns = " + COLUMNS);
    System.out.println("# units = " + UNITS);
    System.out.println();
    state("c7", 0.095, 0.73);
    state("c10", 0.134, 0.81);
    state("c15", 0.206, 0.86);
    state("below_the_boiling_point_switch", 0.500, 0.95);
    state("above_it", 1.000, 1.00);
    state("above_the_heavy_switch", 1.500, 1.05);
  }
}
