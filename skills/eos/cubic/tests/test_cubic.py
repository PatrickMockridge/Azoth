"""The Peng-Robinson coefficient chain `SKILL.md` documents must run."""

import azoth


def test_cubic_chain_runs() -> None:
    kappa = azoth.eos.pr_kappa(0.0115)
    assert isinstance(kappa.kappa, float)
    ab = azoth.eos.pr_alpha_ab(kappa.kappa, 2.0, 1.0)
    assert isinstance(ab.a_reduced, float) and isinstance(ab.b_reduced, float)
    z = azoth.eos.pr_z_factor(ab.a_reduced, ab.b_reduced)
    assert isinstance(z.z_min, float) and isinstance(z.z_max, float)
