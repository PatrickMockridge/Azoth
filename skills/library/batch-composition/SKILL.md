# Build a mixture or look up a component

Every equation-of-state calculation takes a `Mixture` and a composition, and the
substances are looked up by name from the databank. This skill is how a caller turns
names into the values a calculation takes.

## When to Use

- When a flash or phase-boundary calculation needs its components and `kij`.
- When a caller must look a substance up by name, or override one via a keycard.

## Inputs

- Substance names, or the `Tc`/`Pc`/`omega` values directly.
- A keycard, where the caller's values differ from the databank's.

## Outputs

- A `Component` (`Tc`, `Pc`, `omega`), a `Mixture`, and the `kij` pairs.

## How a calculation runs

`component(name)` looks a substance up and returns what a cubic reads. `from_names`
builds a `Mixture` and pulls the published `kij` for each pair. The databank is 286
substances read from `data/components/components.csv`, generated from NeqSim and
attributed in `NOTICE`. A name the databank does not have is refused, not approximated.

## Python usage pattern

```python
import azoth

methane = azoth.eos.component("methane")
methane.Tc  # 190.56 K

fluid = azoth.eos.from_names(["methane", "n-butane"])
len(fluid)  # 2, with the databank's kij pulled in

# A caller's own values: build the components directly, or override via a keycard.
```

## The keycard

A keycard wins over the databank by name, parameter by parameter: overriding `omega`
keeps the shipped `Tc` and `Pc`. A new name must be complete — a partial component is
refused rather than completed from a similar substance, which would be inventing data.

## Validation checklist

- [ ] The name is in the databank or the keycard; otherwise it is refused, not
      approximated.
- [ ] A keycard override names the parameters it changes, and nothing is silently lost.

## Common mistakes

| Symptom | Cause | Fix |
|---|---|---|
| `PropertyUnavailableError` | A name not in the databank or card | Check `azoth.eos.available()` |
| A plausible wrong fluid | A similar substance substituted | Supply the substance's own values |
| The wrong `kij` | A pair omitted and treated as zero | `from_names` pulls the databank pair |

## Limitations

The databank is NeqSim's `COMP.csv` filtered to what a cubic can describe: ions are
excluded because a cubic has no notion of one. `COMP_EXT.csv` (heavy fluids) is not
vendored. Every pure-component constant is still an argument, so the databank is a
default a caller can replace.

## Related Azoth functionality

- `azoth.eos.component`, `azoth.eos.from_names`, `azoth.eos.mixture`, `azoth.eos.available`
  and `azoth.eos.kij_for` — `python/src/azoth/eos/components.py`.
- The `eos` calc pages under `docs/src/eos/`.

## References

- `NOTICE` — attribution for the vendored databank.
