// NeqSim's UNIFAC activity coefficients, reached through public API only.
//
//     javac -proc:none -cp neqsim-f0c7436.jar UnifacGamma.java
//     java -cp .:neqsim-f0c7436.jar UnifacGamma
//
// `eos.unifac_activity_coefficients`'s spec says NeqSim's classic UNIFAC "fails
// at runtime, so there is no differential oracle". It does fail - but the failure is a
// desynchronised array inside NeqSim rather than a dead end, and both halves of it are
// reachable through public API:
//
//   1. `ComponentGEUnifac`'s constructor fills `unifacGroups` (an ArrayList) and never
//      `unifacGroupsArray`, which only `addUNIFACgroup`/`setUnifacGroups` write. The
//      count getter reads the list and the indexer reads the array, so
//      `PhaseGEUnifac.checkGroups` - reached from `setMixingRule` - throws on a
//      length-zero array before anything is computed. `setUnifacGroups` is public, so
//      one call puts the two back in step.
//
//   2. `ComponentGEUnifacPSRK` and `ComponentGEUnifacUMRPRU` never read groups at all:
//      the same constructor returns early for a subclass. There is nothing to
//      synchronise there, and this driver has not yet been shown to work for them.
//
// Every `getActivityCoefficient` below is preceded by an evaluation, for the reason
// `GeGamma.java` sets out: the accessor returns a *cached* field that only a GE
// evaluation fills.

import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemUNIFAC;
import neqsim.thermo.system.SystemUNIFACpsrk;
import neqsim.thermo.component.ComponentGEInterface;
import neqsim.thermo.component.ComponentGEUnifac;
import neqsim.thermo.phase.PhaseGEInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseType;

public class UnifacGamma {
  /**
   * Supply the decomposition to a component whose own constructor never read one.
   *
   * `ComponentGEUnifacPSRK` and `ComponentGEUnifacUMRPRU` inherit a constructor that
   * returns early for a subclass, so their components arrive with no groups at all -
   * there is no desynchronised list to put back in step, only an empty one. The
   * decomposition is read from **NeqSim's own `unifaccomp` table**, through the same
   * JDBC handle its components use, so nothing here supplies a number NeqSim does not
   * already carry.
   */
  static void injectGroups(PhaseInterface ge, String table, int subgroups) {
    for (int i = 0; i < ge.getNumberOfComponents(); i++) {
      ComponentGEUnifac c = (ComponentGEUnifac) ge.getComponent(i);
      try (neqsim.util.database.NeqSimDataBase database = new neqsim.util.database.NeqSimDataBase()) {
        java.sql.ResultSet row =
            database.getResultSet("SELECT * FROM " + table + " WHERE Name='" + c.getName() + "'");
        row.next();
        for (int p = 1; p < subgroups; p++) {
          int count = Integer.parseInt(row.getString("sub" + Integer.toString(p)));
          if (count > 0) {
            c.addUNIFACgroup(p, count);
          }
        }
        row.close();
      } catch (Exception ex) {
        System.out.println("    could not read " + c.getName() + ": " + ex);
      }
      System.out.println("    " + c.getName() + " groups=" + c.getNumberOfUNIFACgroups());
    }
  }

  /** Put every classic UNIFAC component's group list and group array back in step. */
  static void syncGroups(PhaseInterface ge) {
    for (int i = 0; i < ge.getNumberOfComponents(); i++) {
      ComponentGEUnifac c = (ComponentGEUnifac) ge.getComponent(i);
      // The constructor filled the list; the array is still empty.
      c.setUnifacGroups(c.getUnifacGroups2());
      StringBuilder g = new StringBuilder();
      for (int k = 0; k < c.getNumberOfUNIFACgroups(); k++) {
        if (k > 0) {
          g.append(", ");
        }
        g.append("sub").append(c.getUnifacGroup(k).getSubGroup())
            .append(" n=").append(c.getUnifacGroup(k).getN())
            .append(" R=").append(c.getUnifacGroup(k).getR())
            .append(" Q=").append(c.getUnifacGroup(k).getQ());
      }
      System.out.println("    " + c.getName() + " groups=" + c.getNumberOfUNIFACgroups()
          + " array=" + c.getUnifacGroups().length + "  [" + g + "]"
          + "  R_i=" + c.getR() + " Q_i=" + c.getQ());
    }
  }

  /** The GE phase of a system, or `null` when it built none. */
  static PhaseInterface findGePhase(SystemInterface fluid) {
    PhaseInterface ge = null;
    for (int p = 0; p < fluid.getNumberOfPhases(); p++) {
      if (fluid.getPhase(p) instanceof PhaseGEInterface) {
        ge = fluid.getPhase(p);
      }
    }
    return ge;
  }

  static void one(String label, SystemInterface fluid, String[] names, double t) {
    System.out.println(label);
    try {
      for (String name : names) {
        fluid.addComponent(name, 0.5);
      }
      fluid.createDatabase(true);

      // The sync has to happen twice, and the reason is that NeqSim does two things at
      // once in each of the calls around it. `setMixingRule` calls `checkGroups`, which
      // walks the desynchronised array; `init(0)` calls `getExcessGibbsEnergy`, which
      // walks it again - and `init` is also what gives the phase the moles a UNIFAC
      // evaluation weights by, so it cannot simply be dropped or moved.
      PhaseInterface ge = findGePhase(fluid);
      if (ge == null) {
        System.out.println("    no GE phase was built");
        return;
      }
      if (ge.getComponent(0) instanceof ComponentGEUnifac
          && ((ComponentGEUnifac) ge.getComponent(0)).getNumberOfUNIFACgroups() == 0) {
        injectGroups(ge, "unifaccomp", 139);
      } else {
        syncGroups(ge);
      }
      fluid.setMixingRule("classic");
      fluid.init(0);
      ge = findGePhase(fluid);
      syncGroups(ge);
      ((PhaseGEInterface) ge).getExcessGibbsEnergy(ge, ge.getNumberOfComponents(), t, 1.0,
          PhaseType.LIQUID);

      StringBuilder gamma = new StringBuilder();
      for (int i = 0; i < ge.getNumberOfComponents(); i++) {
        if (i > 0) {
          gamma.append(", ");
        }
        gamma.append(((ComponentGEInterface) ge.getComponent(i)).getGamma());
      }
      System.out.println("    molesInPhase " + ge.getNumberOfMolesInPhase());
      System.out.println("    gamma [" + gamma + "]");
    } catch (Throwable failure) {
      System.out.println("    THREW " + failure);
      for (int i = 0; i < Math.min(8, failure.getStackTrace().length); i++) {
        System.out.println("      at " + failure.getStackTrace()[i]);
      }
    }
  }

  /** What NeqSim's own PSRK tables say for the methanol/water main-group pair. */
  static void psrkTables() {
    System.out.println("NeqSim's own tables, main groups 6 and 7:");
    try (neqsim.util.database.NeqSimDataBase database = new neqsim.util.database.NeqSimDataBase()) {
      for (String table : new String[] {"unifacinterparam", "unifacinterparamb", "unifacinterparamc"}) {
        for (String main : new String[] {"6", "7"}) {
          java.sql.ResultSet row =
              database.getResultSet("SELECT * FROM " + table + " WHERE MainGroup=" + main);
          row.next();
          System.out.println("  " + table + "[" + main + "]  n6=" + row.getString("n6")
              + "  n7=" + row.getString("n7") + "  n15=" + row.getString("n15")
              + "  n16=" + row.getString("n16"));
          row.close();
        }
      }
    } catch (Exception ex) {
      System.out.println("  could not read: " + ex);
    }
  }

  public static void main(String[] args) {
    psrkTables();
    String[] methanolWater = {"methanol", "water"};
    String[] acetoneHexane = {"acetone", "n-hexane"};
    one("classic UNIFAC, methanol/water at 298.15 K",
        new SystemUNIFAC(298.15, 1.0), methanolWater, 298.15);
    one("classic UNIFAC, acetone/n-hexane at 298.15 K",
        new SystemUNIFAC(298.15, 1.0), acetoneHexane, 298.15);
    one("UNIFAC-PSRK, methanol/water at 298.15 K",
        new SystemUNIFACpsrk(298.15, 1.0), methanolWater, 298.15);
  }
}
