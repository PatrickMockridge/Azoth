// The Furst electrolyte layers, printed one at a time so each can be ported and checked.
//
//     javac -proc:none -cp neqsim-f0c7436.jar FurstProbe.java
//     java -cp .:neqsim-f0c7436.jar FurstProbe
//
// `PhaseModifiedFurstElectrolyteEos` is `PhaseSrkEos` plus three additive Helmholtz terms:
//
//     getF() = super.getF() + FSR2()*sr2On + FLR()*lrOn + FBorn()*bornOn
//
// The short-range `W` term, the MSA long-range term and the Born solvation term. Everything
// each of them needs is computed in `volInit()`: the solvent dielectric constant and its
// derivatives, the packing fraction, the shielding parameter, `XLR`, `alphaLR2` and `bornX`.
// This prints all of it, in the order the phase builds it, so a port can be checked layer by
// layer rather than only at `ln phi`.
//
// **The rows are `key = value` and the blocks are `#`-opened**, which is the shape
// `tools/neqsim_layer_diff.py` reads. The other electrolyte probes print aligned tables and
// so cannot feed it; this one can, and `eos.furst_electrolyte_phase` gets a layer diff
// because of it.
//
// Three things the source does not make obvious, which the rows below are here to measure:
//
//   1. **The solvent dielectric constant sums over `ionicCharge == 0` only.** The ions are
//      excluded from the average, so adding salt does not move the dielectric constant
//      directly - it moves it through the mole fractions of what is left.
//   2. **`calcSolventDiElectricConstantdT` is the molar-average formula whatever rule is
//      selected.** `setDielectricMixingRule` switches the value and not its derivative, so
//      the two disagree under VOLUME_AVERAGE and LOOYENGA.
//   3. **The Born term is absent from `dFdV`, `dFdTdV`, `dFdVdV` and `dFdVdVdV`.** It is in
//      `getF`, `dFdT` and `dFdTdT` only - which is consistent, because `FBorn` carries no
//      volume, but it is the kind of omission a port has to confirm rather than assume.

import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.phase.PhaseModifiedFurstElectrolyteEos;
import neqsim.thermo.system.SystemFurstElectrolyteEos;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermodynamicoperations.ThermodynamicOperations;

public class FurstProbe {

  /** One `key = value` row. */
  private static void row(String key, double value) {
    System.out.printf("%s = %.15g%n", key, value);
  }

  /** One package-private field, by name. */
  private static double field(Object target, String name) {
    try {
      java.lang.reflect.Field f = target.getClass().getDeclaredField(name);
      f.setAccessible(true);
      return f.getDouble(target);
    } catch (Throwable e) {
      return Double.NaN;
    }
  }

  private static void report(String label, String[] names, double[] moles, double tC, double pBara) {
    System.out.printf("# %s%n", label);
    SystemInterface system = new SystemFurstElectrolyteEos(298.15, 10.01325);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.setMixingRule(4);
    system.setTemperature(tC, "C");
    system.setPressure(pBara, "bara");
    try {
      ThermodynamicOperations ops = new ThermodynamicOperations(system);
      ops.TPflash();
      system.initProperties();
      row("phases", system.getNumberOfPhases());
      for (int p = 0; p < system.getNumberOfPhases(); p++) {
        PhaseInterface phase = system.getPhase(p);
        System.out.printf("# %s phase %d%n", label, p);
        reportPhase(phase);
      }
    } catch (Throwable error) {
      System.out.printf("# %s failed%n", label);
      System.out.println("error = " + error.getClass().getSimpleName() + ": " + error.getMessage());
      for (StackTraceElement frame : error.getStackTrace()) {
        if (frame.getClassName().startsWith("neqsim.")) {
          System.out.println("#   at " + frame);
        }
      }
    }
  }

  private static void reportPhase(PhaseInterface phase) {
    double temperature = phase.getTemperature();
    row("T", temperature);
    row("P", phase.getPressure());
    row("n_total", phase.getNumberOfMolesInPhase());
    row("V", phase.getMolarVolume());
    row("Z", phase.getZ());

    int n = phase.getNumberOfComponents();
    // **The covolumes, which is the one thing the port cannot get from a formula.** An ion's
    // `b` is not `0.08664 R Tc/Pc`: `ComponentModifiedFurstElectrolyteEos` overwrites it
    // with `(p0 d^3 + p1) 1e5` from the fitted parameters and sets `a = 1e-35`, and azoth's
    // `Component` cannot express either. Reading both sides is what fixes the conversion.
    row("A_phase", ((neqsim.thermo.phase.PhaseEos) phase).getA());
    row("B_phase", phase.getB());
    for (int i = 0; i < n; i++) {
      neqsim.thermo.component.ComponentEosInterface component =
          (neqsim.thermo.component.ComponentEosInterface) phase.getComponent(i);
      row("b[" + i + "]", component.getb());
      row("a[" + i + "]", component.geta());
      row("aT[" + i + "]", component.getaT());
      // **The Huron-Vidal rule's own per-component quantity.** `ader_i` is
      // `a_i/(b_i R T) - ln(gamma_i)/lambda`, and `alphaMix = sum_i x_i ader_i` is what
      // `calcA` multiplies by `n B R T`. Printed beside `A` because a difference in the
      // sum is a difference in one of these, and every way of getting one wrong looks the
      // same in `A` alone.
      row("ader[" + i + "]", component.getAder());
      row("alpha[" + i + "]", component.getAttractiveTerm().alpha(temperature));
      // The component's own constants, which need not be the table's: `TC`/`PC` are the
      // table's in degrees C and bar, and `getTC()` is kelvin.
      row("omega[" + i + "]", component.getAcentricFactor());
      row("tc[" + i + "]", component.getTC());
      {
        java.lang.Object term = component.getAttractiveTerm();
        java.lang.Object params = null;
        for (java.lang.Class<?> k = term.getClass(); k != null; k = k.getSuperclass()) {
          try {
            java.lang.reflect.Field f = k.getDeclaredField("parameters");
            f.setAccessible(true);
            params = f.get(term);
            break;
          } catch (Throwable ignored) {
            // keep walking up
          }
        }
        System.out.printf("params[%d] = %s  term = %s%n", i, java.util.Arrays.toString((double[]) params),
            term.getClass().getSimpleName());
      }
      row("pc[" + i + "]", component.getPC());
    }
    for (int i = 0; i < n; i++) {
      row("x[" + i + "]", phase.getComponent(i).getx());
      row("charge[" + i + "]", phase.getComponent(i).getIonicCharge());
      row("lj[" + i + "]", phase.getComponent(i).getLennardJonesMolecularDiameter());
      row("eps_i[" + i + "]", phase.getComponent(i).getDielectricConstant(temperature));
    }

    if (!(phase instanceof PhaseModifiedFurstElectrolyteEos)) {
      System.out.println("note = not a Furst phase: " + phase.getClass().getSimpleName());
      return;
    }
    PhaseModifiedFurstElectrolyteEos furst = (PhaseModifiedFurstElectrolyteEos) phase;

    // What `volInit` built, in the order it built it.
    row("eps", furst.getSolventDiElectricConstant());
    row("eps_dT", furst.getSolventDiElectricConstantdT());
    row("eps_dTdT", furst.getSolventDiElectricConstantdTdT());
    row("eps_phase", furst.calcDiElectricConstant(temperature));
    row("eps_phase_dT", furst.calcDiElectricConstantdT(temperature));
    row("eps_phase_dTdT", furst.calcDiElectricConstantdTdT(temperature));
    row("eps_phase_dV", furst.calcDiElectricConstantdV(temperature));
    row("eps_phase_dVdV", furst.calcDiElectricConstantdVdV(temperature));
    row("eps_phase_dTdV", furst.calcDiElectricConstantdTdV(temperature));
    row("packing", furst.getEps());
    row("packing_V", furst.calcEpsV());
    row("packing_VV", furst.calcEpsVV());
    row("packing_ionic", furst.getEpsIonic());
    row("packing_ionic_V", furst.calcEpsIonicdV());
    row("gamma", furst.getShieldingParameter());
    row("gamma_dT", furst.calcShieldingParameterdT());
    row("alphaLR2", furst.getAlphaLR2());
    row("XLR", furst.getXLR());
    row("XLR_dT", furst.calcXLRdT());
    row("bornX", furst.calcBornX());

    // The short-range parameter table and the three terms, with the temperature
    // derivatives that say each is live rather than merely present.
    row("W", furst.getW());
    row("WT", furst.getWT());
    // `WTT` is a package-private field with no accessor, so it is read by reflection -
    // the same reach `PitzerArithmetic` makes for `debyeHuckelAphi`.
    row("WTT", field(furst, "WTT"));
    row("FSR2", furst.FSR2());
    row("FSR2_dT", furst.dFSR2dT());
    row("FSR2_dTdT", furst.dFSR2dTdT());
    row("FSR2_dV", furst.dFSR2dV());
    row("FSR2_dVdV", furst.dFSR2dVdV());
    row("FSR2_dTdV", furst.dFSR2dTdV());
    row("FLR", furst.FLR());
    row("FLR_dT", furst.dFLRdT());
    row("FLR_dTdT", furst.dFLRdTdT());
    row("FLR_dV", furst.dFLRdV());
    row("FLR_dVdV", furst.dFLRdVdV());
    row("FLR_dTdV", furst.dFLRdTdV());
    row("FBorn", furst.FBorn());
    row("FBorn_dT", furst.dFBorndT());
    row("FBorn_dTdT", furst.dFBorndTdT());

    // The Helmholtz energy and the derivatives the cubic's root and fugacity read.
    row("F", furst.getF());
    row("F_srk", furst.getF() - furst.FSR2() - furst.FLR() - furst.FBorn());
    row("dFdT", furst.dFdT());
    row("dFdV", furst.dFdV());
    row("dFdTdT", furst.dFdTdT());
    row("dFdVdV", furst.dFdVdV());
    row("dFdTdV", furst.dFdTdV());

    for (int i = 0; i < n; i++) {
      row("lnPhi[" + i + "]", Math.log(phase.getComponent(i).getFugacityCoefficient()));
    }

    // **The composition derivatives, which are what `ln phi` is actually built from.**
    // `ComponentEos.fugcoef` is `exp(dFdN - ln(PV/RT))`, so a port that reproduces `F` and
    // not `dFdN` has reproduced the energy and not the fugacity. These are the three
    // electrolyte contributions separately, so the model's `ln phi` can be checked layer by
    // layer rather than only at the end.
    for (int i = 0; i < n; i++) {
      neqsim.thermo.component.ComponentModifiedFurstElectrolyteEos c =
          (neqsim.thermo.component.ComponentModifiedFurstElectrolyteEos) phase.getComponent(i);
      row("dFdN[" + i + "]", c.dFdN(phase, n, temperature, phase.getPressure()));
      row("dFSR2dN[" + i + "]", c.dFSR2dN(phase, n, temperature, phase.getPressure()));
      row("dFLRdN[" + i + "]", c.dFLRdN(phase, n, temperature, phase.getPressure()));
      row("dFBorndN[" + i + "]", c.dFBorndN(phase, n, temperature, phase.getPressure()));
      row("xlnN[" + i + "]", c.getx() * phase.getNumberOfMolesInPhase());
    }
  }

  /** The same rows for the Mod2004 phase, whose class shares its method names. */
  private static void reportMod2004(String label, String[] names, double[] moles, double tC,
      double pBara) {
    System.out.printf("# %s%n", label);
    neqsim.thermo.system.SystemInterface system =
        new neqsim.thermo.system.SystemFurstElectrolyteEosMod2004(298.15, 10.01325);
    for (int i = 0; i < names.length; i++) {
      system.addComponent(names[i], moles[i]);
    }
    system.setMixingRule(4);
    system.setTemperature(tC, "C");
    system.setPressure(pBara, "bara");
    try {
      ThermodynamicOperations ops = new ThermodynamicOperations(system);
      ops.TPflash();
      system.initProperties();
      for (int p = 0; p < system.getNumberOfPhases(); p++) {
        PhaseInterface phase = system.getPhase(p);
        System.out.printf("# %s phase %d%n", label, p);
        row("T", phase.getTemperature());
        row("n_total", phase.getNumberOfMolesInPhase());
        row("V", phase.getMolarVolume());
        row("Z", phase.getZ());
        int n = phase.getNumberOfComponents();
        for (int i = 0; i < n; i++) {
          row("x[" + i + "]", phase.getComponent(i).getx());
          row("lnPhi[" + i + "]", Math.log(phase.getComponent(i).getFugacityCoefficient()));
        }
        // The five quantities Mod2004 zeroes, read reflectively so one method serves both
        // classes - and printed so the difference from the base model is a measurement.
        for (String name : new String[] {"getSolventDiElectricConstantdT", "getShieldingParameter",
            "getXLR", "getSolventDiElectricConstant", "getDielectricConstant"}) {
          row(name, (Double) phase.getClass().getMethod(name).invoke(phase));
        }
      }
    } catch (Throwable error) {
      System.out.printf("# %s failed: %s%n", label, error);
    }
  }

  public static void main(String[] args) {
    // `SystemFurstElectrolyteEosTest`'s own mixture, at its own state.
    report("the shipped test: methane water Na+ Cl-",
        new String[] {"methane", "water", "Na+", "Cl-"},
        new double[] {0.1, 1.0, 0.001, 0.001}, 25.0, 10.01325);

    // The same at four times the salt, so the electrostatic terms move and the SRK part
    // barely does.
    report("the same at four times the salt",
        new String[] {"methane", "water", "Na+", "Cl-"},
        new double[] {0.1, 1.0, 0.004, 0.004}, 25.0, 10.01325);

    // A mixed solvent: the dielectric mixing rule is what separates this from the above.
    report("a mixed solvent: methanol joins the water",
        new String[] {"methane", "water", "methanol", "Na+", "Cl-"},
        new double[] {0.1, 0.6, 0.4, 0.001, 0.001}, 25.0, 10.01325);

    // **The same composition at ten times the moles.** A Helmholtz energy's volume
    // derivative is a pressure and is intensive: `dFdV` must not move with the phase's
    // size. `FSR2V` is written against `(V n)^2` and `epsdV` against `V n`, which is the
    // *total* volume - so if the `1e-5` beside them is the only scaling, every volume
    // derivative here carries one factor of `n` too many or too few, and this is the
    // measurement that says which. `lnPhi` is the control: it does not depend on them.
    report("the same at ten times the moles",
        new String[] {"methane", "water", "Na+", "Cl-"},
        new double[] {1.0, 10.0, 0.01, 0.01}, 25.0, 10.01325);

    reportMod2004("MOD2004 shipped test", new String[] {"methane", "water", "Na+", "Cl-"},
        new double[] {0.1, 1.0, 0.001, 0.001}, 25.0, 10.01325);

    // A gas-rich state, where the aqueous phase barely exists.
    report("at 60 C and 40 bara",
        new String[] {"methane", "water", "Na+", "Cl-"},
        new double[] {0.1, 1.0, 0.001, 0.001}, 60.0, 40.0);
  }
}
