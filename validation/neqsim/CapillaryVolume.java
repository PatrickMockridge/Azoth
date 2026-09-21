// What `system.getPhase(1).getMolarVolume("m3/mol")` returns at the state
// `CapillaryDewPointFlash` evaluates its Kelvin exponent at.
//
//     javac -proc:none -cp neqsim-f0c7436.jar CapillaryVolume.java
//     java -cp .:neqsim-f0c7436.jar CapillaryVolume
//
// The class sets `setBeta(1, 1e-15)` and calls `getMolarVolume` on that phase inside its
// iteration. azoth takes the same exponent's volume from the cubic's own lower root as
// `z_liquid R T/P`. The two differ by about 6 per cent, which is the whole of the 5.3 per cent
// shortfall in the dew-point shift the case records; this asks whether the difference is the
// liquid's volume under a different definition, or a number that is not a liquid volume at all.
//
// NeqSim holds bar throughout.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemPrEos;

public class CapillaryVolume {

  static final double R = 8.31446261815324;

  public static void main(String[] args) {
    double t = 347.6658507316;
    double pBar = 20.0;
    SystemInterface fluid = new SystemPrEos(t, pBar);
    fluid.addComponent("methane", 0.5);
    fluid.addComponent("n-butane", 0.5);
    fluid.setMixingRule("classic");
    fluid.setAttractiveTerm(1);
    fluid.init(0);
    fluid.setBeta(0, 1.0 - 1e-15);
    fluid.setBeta(1, 1e-15);
    fluid.init(1);
    fluid.setNumberOfPhases(2);
    // The incipient liquid's composition, which is what the class iterates on.
    fluid.getPhases()[1].getComponent(0).setx(0.0560038706);
    fluid.getPhases()[1].getComponent(1).setx(0.9439961294);
    fluid.init(1);
    // The class reads the volume only after `getSurfaceTensionSafe` has run, which calls this.
    fluid.initPhysicalProperties();
    fluid.init(1);

    for (int ph = 0; ph < fluid.getNumberOfPhases(); ph++) {
      double v = fluid.getPhase(ph).getMolarVolume("m3/mol");
      System.out.printf("phase %d type=%-4s beta=%.3e moles=%.6e V=%.12e  Z(from V)=%.10f  "
              + "Z(phase)=%.10f  x=(%.6f, %.6f)%n",
          ph, fluid.getPhase(ph).getType(), fluid.getPhase(ph).getBeta(),
          fluid.getPhase(ph).getNumberOfMolesInPhase(), v,
          pBar * 1e5 * v / (R * fluid.getPhase(ph).getTemperature()), fluid.getPhase(ph).getZ(),
          fluid.getPhase(ph).getComponent(0).getx(), fluid.getPhase(ph).getComponent(1).getx());
    }
    System.out.printf("%nThe liquid root at that composition and state is Z = 0.076308064834, "
        + "so `z R T/P` = %.12e%n", 0.076308064834 * R * t / (pBar * 1e5));
  }
}
