"""Author a keycard and read a component through it."""

import azoth

card = azoth.keycard.use(
    {
        "schema_version": 2,
        "keyholder": {"name": "example"},
        "components": {"methane": {"omega": {"value": 0.1111, "unit": "dimensionless"}}},
    }
)

base = azoth.eos.component("methane")
overridden = azoth.eos.component("methane", card=card)
print("databank omega:", base.omega)
print("card omega:    ", overridden.omega)
print("Tc kept from databank:", overridden.Tc)
