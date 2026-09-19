import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

/**
 * Drives NeqSim's `VUflashSingleComp` on the states `eos.vu_flash_single_comp` carries cases for.
 *
 * <p>A pure component in the two-phase region: the pressure fixes the saturation line, and an
 * internal energy between the saturated liquid's and the vapour's fixes the split. The probe
 * finds the saturation temperature itself, prints the two saturated internal energies, and then
 * hands the halfway value back to `ThermodynamicOperations.VUflash` - which routes a pure
 * component to `VUflashSingleComp` - and prints the temperature and vapour fraction it returns.
 *
 * <p><b>The claim under test is that the `Vspec` it takes is never read.</b> So the probe runs
 * every state twice with the same pressure and internal energy and two wildly different volumes
 * - one the physical `z R T/P` of the saturation state, and one a cubic metre for a mole - and
 * prints both answers beside each other. If the volume were read anywhere in the flash the two
 * lines would differ; they are the measurement, and the printed pair is what makes this a test
 * rather than an assertion about the source.
 */
public class VuSingleCompProbe {

  private static void one(String name, double pressureBar) {
    SystemInterface sat = new SystemPrEos(300.0, pressureBar);
    sat.addComponent(name, 1.0);
    sat.setMixingRule("classic");
    sat.setAttractiveTerm(1);
    try {
      new ThermodynamicOperations(sat).bubblePointTemperatureFlash();
    } catch (Exception e) {
      System.out.printf("%-10s %6.1f bar  bubblePointTemperatureFlash THREW %s%n", name, pressureBar,
          e.getClass().getSimpleName());
      return;
    }
    double tsat = sat.getTemperature();
    sat.init(2);
    double uGas = sat.getPhase(0).getInternalEnergy("J/mol");
    double uLiq = sat.getPhase(1).getInternalEnergy("J/mol");
    double uSpec = 0.5 * (uLiq + uGas);
    System.out.printf("%-10s %6.1f bar  Tsat=%.12f  u_liq=%.9f  u_gas=%.9f  u_spec=%.9f%n", name,
        pressureBar, tsat, uLiq, uGas, uSpec);

    // Two volumes: the physical one for this state, and a cubic metre for a mole, which is
    // four orders of magnitude away. A flash that read it could not answer the same twice.
    double physical = sat.getPhase(0).getZ() * 8.3144621 * tsat / (pressureBar * 1.0e5);
    for (double volume : new double[] {physical, 1.0}) {
      SystemInterface flash = new SystemPrEos(298.0, pressureBar);
      flash.addComponent(name, 1.0);
      flash.setMixingRule("classic");
      flash.setAttractiveTerm(1);
      try {
        new ThermodynamicOperations(flash).VUflash(volume, uSpec);
        System.out.printf("   V=%-12.6g -> T=%.12f beta=%.12f P=%.6f bar%n", volume,
            flash.getTemperature(), flash.getBeta(), flash.getPressure());
      } catch (Exception e) {
        System.out.printf("   V=%-12.6g -> VUflash THREW %s: %s%n", volume,
            e.getClass().getSimpleName(), e.getMessage());
      }
    }
  }

  public static void main(String[] args) {
    one("propane", 10.0);
    one("propane", 20.0);
    one("n-butane", 2.0);
  }
}
