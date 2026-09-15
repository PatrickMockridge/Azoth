"""The Peng-Robinson coefficient chain: kappa -> a/b -> Z."""

import azoth

kappa = azoth.eos.pr_kappa(0.0115)
ab = azoth.eos.pr_alpha_ab(kappa.kappa, 2.0, 1.0)
z = azoth.eos.pr_z_factor(ab.a_reduced, ab.b_reduced)

print("kappa:", kappa.kappa)
print("a_reduced, b_reduced:", ab.a_reduced, ab.b_reduced)
print("z roots:", z.z_min, z.z_max)
