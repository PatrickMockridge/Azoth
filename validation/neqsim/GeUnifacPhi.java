// What `PhaseGEUnifac` reports for a liquid, for `eos.ge_unifac_phase`.
//
//     javac -proc:none -cp neqsim-3.20.0.jar GeUnifacPhi.java
//     java -cp .:neqsim-3.20.0.jar GeUnifacPhi
//
// The phase's fugacity coefficient is `gamma_i * P0_i / P`, the same composition
// `GeNrtlPhi.java` prints for the NRTL phase - see `ComponentGE.fugcoef` - so this is the
// ported UNIFAC activity coefficient over the ported Antoine vapour pressure.
//
// **The group array has to be synchronised first.** `ComponentGEUnifac`'s constructor
// fills `unifacGroups` but not `unifacGroupsArray`, and `getNumberOfUNIFACgroups()` reads
// the one while `getUnifacGroup(i)` reads the other, so a freshly built phase computes on
// an empty array. `PhaseGEUnifac.checkGroups()` is the public method that repairs it: it
// adds every subgroup the phase's *other* components carry, with count zero, then sorts
// and reindexes. Without it the phase silently returns the combinatorial term alone.
//
// The evaluation before the read is `GeGamma.java`'s requirement: a component's `gamma`
// is a cached field that only `getExcessGibbsEnergy` fills, and `fugcoef` reads it.

import neqsim.thermo.system.SystemUNIFAC;
import neqsim.thermo.component.ComponentGEInterface;
import neqsim.thermo.phase.PhaseGEInterface;
import neqsim.thermo.phase.PhaseGEUnifac;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;

public class GeUnifacPhi {

  static void one(String first, String second, double temperatureK, double pressureBar, double xFirst) {
    SystemUNIFAC system = new SystemUNIFAC(temperatureK, pressureBar);
    system.addComponent(first, xFirst);
    system.addComponent(second, 1.0 - xFirst);
    // `ComponentGEUnifac`'s constructor fills each component's group *list* and leaves
    // its group *array* empty; `getNumberOfUNIFACgroups()` reads the list and
    // `getUnifacGroup(i)` the array. Every step that follows reads the array - the first
    // `init(0)`, `setMixingRule`'s `checkGroups()`, and every later `init(0)` - so the
    // array is synchronised before each of them.
    //
    // `checkGroups()` cannot do this repair itself: it rebuilds the list from the array,
    // so called on an empty array it empties the list too. The order is forced.
    system.createDatabase(true);
    synchronise(system);
    system.init(0);
    synchronise(system);
    system.setMixingRule("classic");
    synchronise(system);
    system.init(0);
    synchronise(system);

    PhaseInterface ge = null;
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      if (system.getPhase(p) instanceof PhaseGEUnifac) {
        ge = system.getPhase(p);
      }
    }
    if (ge == null) {
      System.out.println("no UNIFAC phase was built");
      return;
    }
    ((PhaseGEUnifac) ge).checkGroups();

    ((PhaseGEInterface) ge).getExcessGibbsEnergy(ge, ge.getNumberOfComponents(), temperatureK,
        pressureBar, PhaseType.LIQUID);

    System.out.println(first + "/" + second + "  T=" + temperatureK + " P=" + pressureBar
        + "  x=" + xFirst + "  (" + ge.getClass().getSimpleName() + ")");
    for (int i = 0; i < ge.getNumberOfComponents(); i++) {
      ge.getComponent(i).fugcoef(ge);
      double gamma = ((ComponentGEInterface) ge.getComponent(i)).getGamma();
      double p0 = ge.getComponent(i).getAntoineVaporPressure(temperatureK);
      double phi = ge.getComponent(i).getFugacityCoefficient();
      System.out.println("  " + ge.getComponent(i).getName()
          + "  x=" + ge.getComponent(i).getx()
          + "  gamma=" + gamma
          + "  P0_bar=" + p0
          + "  phi=" + phi
          + "  gamma*P0/P=" + (gamma * p0 / pressureBar));
    }
  }

  /**
   * Fill every component's group *array* from its group *list*.
   *
   * `ComponentGEUnifac`'s constructor fills the list and leaves the array empty;
   * `getNumberOfUNIFACgroups()` reads the list and `getUnifacGroup(i)` the array. Both
   * writers are public, so this is the whole of the repair.
   */
  static void synchronise(neqsim.thermo.system.SystemInterface system) {
    for (int p = 0; p < system.getNumberOfPhases(); p++) {
      if (!(system.getPhase(p) instanceof PhaseGEUnifac)) {
        continue;
      }
      PhaseInterface ge = system.getPhase(p);
      for (int i = 0; i < ge.getNumberOfComponents(); i++) {
        neqsim.thermo.component.ComponentGEUnifac component =
            (neqsim.thermo.component.ComponentGEUnifac) ge.getComponent(i);
        component.setUnifacGroups(component.getUnifacGroups2());
      }
    }
  }

  public static void main(String[] args) {
    one("methanol", "water", 298.15, 1.0, 0.5);
    one("methanol", "benzene", 330.0, 1.0, 0.3);
    one("water", "toluene", 340.0, 1.0, 0.7);
  }
}
