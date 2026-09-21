import neqsim.thermo.phase.PhaseCPAInterface;
import neqsim.thermo.phase.PhaseInterface;
import neqsim.thermo.system.SystemInterface;
import neqsim.thermo.system.SystemSrkCPA;

/**
 * What becomes of a CPA association energy whose scheme carries one charge sign.
 *
 * <p>
 * {@code CPAMixingRuleHandler}'s bond test was {@code charge[i] * charge[j] < 0}, so the
 * {@code 1A} scheme - one site of one sign - and the {@code 2A} - two of the same - could not
 * bond with themselves, and the fitted association energy the databank carries for them reached
 * no calculation. {@code #3832} replaced the sign test with an all-ones matrix for those two
 * schemes and left it for {@code 2B} and {@code 4C}, whose sites really are donor and acceptor.
 * <b>Cross-solvation keeps the sign test</b>, so both numbers here matter: the self term and
 * the mixture's.
 *
 * <p>
 * <b>The databank's own numbers make it checkable.</b> Five {@code 1A} components carry an
 * association energy of 40323 or 41917 J/mol - acetic, formic, hydrochloric, sulfuric and
 * nitric acid - and four {@code 2A} components carry 5000 J/mol: H2S, SF6, R12 and R134a. Those
 * are fitted values rather than defaulted zeros, and this prints what becomes of them.
 *
 * <p>
 * <b>The second defect is in the table rather than the code, and this still shows it.</b> Seven
 * of those rows name a {@code 1A} scheme and a site count of zero, so they have no sites before
 * any bond test is reached; the last state below is one, and its {@code hcpa} is zero for a
 * different reason. {@code hcpa} is the association's Helmholtz contribution, so zero means the
 * energy was read by nothing.
 */
/**
 * Whether NeqSim's CPA bonds anything for a component whose scheme carries one charge sign.
 *
 * <p>
 * {@code CPAMixingRuleHandler}'s bond test is {@code charge[i] * charge[j] < 0}: two sites bond
 * when their charges have opposite signs. The {@code 1A} scheme is one site of one sign and the
 * {@code 2A} scheme is two of the same, so neither can bond with itself - a {@code 2A}
 * component in a mixture of {@code 2A} components has an all-zero interaction matrix, and the
 * fitted association energy the databank carries for it is read by no calculation.
 *
 * <p>
 * <b>The databank's own numbers make the claim checkable.</b> Five {@code 1A} components carry
 * an association energy of 40323 or 41917 J/mol - acetic, formic, hydrochloric, sulfuric and
 * nitric acid - and four {@code 2A} components carry 5000 J/mol: H2S, SF6, R12 and R134a. Those
 * are not defaulted zeros; they are fitted values, and this prints what becomes of them.
 *
 * <p>
 * Three states, because the claim is not "a 2A component never associates". Water is the
 * control that association is visible at all; H2S with water is the control that a 2A
 * component's sites still bond with *another* scheme's, which is real physics; and H2S alone is
 * the finding. A probe with only the third state could not tell "the scheme is inert" from "the
 * probe cannot see association".
 *
 * <pre>
 * javac -proc:none -cp neqsim-f0c7436.jar CpaSchemeProbe.java
 * java -cp .:neqsim-f0c7436.jar CpaSchemeProbe
 * </pre>
 */
public final class CpaSchemeProbe {

  private CpaSchemeProbe() {}

  private static void probe(String label, String[][] components, double[] fractions) {
    SystemInterface system = new SystemSrkCPA(300.0, 100.0);
    for (int i = 0; i < components.length; i++) {
      system.addComponent(components[i][0], fractions[i]);
    }
    system.setMixingRule(10);
    system.init(0);
    system.init(1);

    PhaseInterface phase = system.getPhase(0);
    System.out.printf("%n# %s at 300 K, 100 bara%n", label);
    for (int i = 0; i < components.length; i++) {
      System.out.printf("   %-14s scheme=%-3s sites=%d%n", components[i][0], components[i][1],
          phase.getComponent(i).getNumberOfAssociationSites());
    }
    System.out.printf("   phase sites              = %d%n",
        ((PhaseCPAInterface) phase).getTotalNumberOfAccociationSites());
    System.out.printf("   hcpa  the association    = %.15g%n",
        ((PhaseCPAInterface) phase).getHcpatot());
    System.out.printf("   (hcpa = 0 means every association term was computed over nothing)%n");
  }

  public static void main(String[] args) {
    System.out.println("# azoth CpaSchemeProbe - NeqSim master's CPA bond test.");
    System.out.println("# A site bonds when charge[i]*charge[j] < 0, except within the 1A and 2A");
    System.out.println("# schemes, whose equivalent sites bond with each other; cross-solvation");
    System.out.println("# keeps the sign test. hcpa is the association's Helmholtz contribution:");
    System.out.println("# zero means the fitted energy was read by nothing.");

    // The control: 4C is two of each sign, so water associates with itself.
    probe("water alone (4C)", new String[][] {{"water", "4C"}}, new double[] {1.0});
    // The 2A component with the strongest fitted energy in the table.
    probe("H2S alone (2A, eps = 5000 J/mol)",
        new String[][] {{"H2S", "2A"}}, new double[] {1.0});
    // And the same 2A sites, where a 4C partner gives them something to bond to.
    probe("H2S and water (2A with 4C)",
        new String[][] {{"H2S", "2A"}, {"water", "4C"}}, new double[] {0.5, 0.5});
    // A `1A` row that carries one site, so "1A cannot bond with itself" is measured apart
    // from "this row has no sites at all" - which is what acetic acid below turns out to be.
    probe("asphaltene alone (1A with one site, eps = 3500 J/mol)",
        new String[][] {{"asphaltene", "1A"}}, new double[] {1.0});
    // A 1A component with the largest fitted energy in the table - and, in the table, no
    // site count at all, so the energy has nothing to act on before the scheme is reached.
    probe("acetic acid alone (1A, eps = 40323 J/mol)",
        new String[][] {{"acetic acid", "1A"}}, new double[] {1.0});
  }
}
