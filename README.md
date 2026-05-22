# SMART

**S**implified **M**ass **A**blation and **R**elocation **T**reatment — a fast,
reduced-physics model of high-field-side / low-field-side pellet injection
into a tokamak plasma, packaged as an [IMAS](https://imas.iter.org/) actor and
a standalone Fortran executable.

SMART couples a Kuteev/Parks ablation model [1,2] to a Strauss-style
post-ablation mass relocation along the magnetic field [3], following the
formulation of Polevoi & Shimada [4]. It reads `equilibrium` and
`core_profiles` IDSs (and optionally a `pellets` IDS) and returns an updated
`core_profiles` IDS suitable for use inside an integrated-modelling workflow.

## Getting Started

### Dependencies

SMART is written in Fortran 90 / Fortran 77 and Python and links against the
IMAS Access Layer. To build and run you need:

| Dependency | Purpose | Where to obtain |
| ---------- | ------- | --------------- |
| A Fortran compiler (`ifx`, `ifort`, or `gfortran`) | Build SMART | vendor / distribution |
| IMAS Access Layer (Fortran + Python bindings) | IDS I/O | https://github.com/iterorganization/IMAS-AL-Core (and language bindings) |
| MDSplus or HDF5 | IDS backend | https://www.mdsplus.org / https://www.hdfgroup.org |
| XMLlib | XML parsing in Fortran | https://github.com/iterorganization (XMLlib) |
| Fundamental-Constants | Physical constants module | https://github.com/iterorganization (Fundamental-Constants) |
| iWrap (optional) | Generates the IMAS Python actor from `libsmart.a` | https://github.com/iterorganization (iWrap) |
| `pkg-config` | Discovers compile/link flags | distribution package |
| `xmllint` (optional) | XSD validation of input files | `libxml2-utils` |

On the ITER SDCC environment all of the above are available as
[environment-modules](https://lmod.readthedocs.io/) — see
[`config_sdcc.sh`](config_sdcc.sh) for the loadout we use, and the
[`ci-sdcc/`](ci-sdcc/) directory for the CI scripts.

### Building

```sh
git clone https://github.com/iterorganization/smart.git
cd smart

# On ITER SDCC:
source config_sdcc.sh

# Elsewhere: set FC and ensure pkg-config finds al-fortran, xmllib,
#            fundamental-constants. Then:
make            # builds libsmart.a + standalone executable
make actor      # builds the Python IMAS actor (requires iWrap)
make validate   # checks the input XML against the XSD
```

### Input Parameters (`input/smart.xml`)

#### `STEPUP` — transport step parameters
| Name | Unit | Description |
| ---- | ---- | ----------- |
| TAU  | s    | Time step   |
| dtau | s    | Time interval between consecutive pellet injections |
| QNB  | 10<sup>19</sup> s<sup>-1</sup> | Boundary influx of neutrals (D+T+H) for transport equation in CX approximation |
| F01B | 10<sup>19</sup> m<sup>-3</sup> | Edge neutral density (puffing, H component) |
| F02B | 10<sup>19</sup> m<sup>-3</sup> | Edge neutral density (puffing, D component) |
| F03B | 10<sup>19</sup> m<sup>-3</sup> | Edge neutral density (puffing, T component) |
| GN2E | - | Coefficient for convective term in heat transport (0 / 1.5 / 2.5) — electrons |
| GN2I | - | Coefficient for convective term in heat transport (0 / 1.5 / 2.5) — ions |
| sw_stdout | - | Toggle for debug information to stdout (off when 0, on when ≠ 0) |

#### `SMART` — pellet parameters
| Name | Unit | Description |
| ---- | ---- | ----------- |
| YAM  | AMU  | Pellet mass in atomic units (2.5 for 50:50 DT) |
| YVP  | km/s | Pellet velocity |
| YVOL | mm<sup>3</sup>  | Pellet size |
| YCOS0| -    | Cosine of angle to normal of the last closed surface |
| YEFF | -    | Efficiency of loss in the guide tube (1 = 100% intact pellet) |
| YDL  | -    | Ratio of distance along the trajectory to that mapped on the mid-plane |
| yswitch | - | Switch between distinct-pellet and permanent-pellet models |
| sw_smart | - | Toggle for SMART (off when 0, on when ≠ 0) |

#### `ECH2a` — ECH heating term
| Name | Unit | Description |
| ---- | ---- | ----------- |
| ROCEC | - | Normalised EC location w.r.t. toroidal flux coordinate = sqrt(phi/(π·b0)) |
| ROCDR | - | Normalised EC width w.r.t. toroidal flux coordinate = sqrt(phi/(π·b0)) |
| QECR | MW | ECRH power |
| YEFFec | MA/MW | I<sub>EC</sub> / Q<sub>EC</sub> |
| sw_ech2a | - | Toggle for ECH2a (off when 0, on when ≠ 0) |

### Running

#### Standalone Fortran wrapper

```sh
# Edit input/standalone.xml and/or input/smart.xml as needed
./smart
```

#### IMAS actor via Python

```sh
# Edit input/scenario.yaml and/or input/smart.xml as needed
python ./run_smart.py
```

The Python wrapper expects an `equilibrium` IDS and a `core_profiles` IDS in
the IMAS database referenced by `scenario.yaml`. By default it reads the ITER
public scenarios.

## References

1. Kuteev, B., *Pellet ablation model*, Nucl. Fusion **35** (1995) 431
2. Parks, P. B., *Cloud size*, Phys. Plasmas **7** (2000) 1968
3. Strauss, H. R., *Mass relocation model*, Phys. Plasmas **5** (1998) 2676
4. **Polevoi, A. R. and Shimada, M.**, *Simplified mass ablation and relocation
   treatment for pellet injection optimization*, Plasma Phys. Control. Fusion
   **43** (2001) 1525 — *primary reference*
5. Pereverzev, G. V. and Yushmanov, P. N., *ASTRA automated system for transport
   analysis in a tokamak*, Max-Planck IPP Report **5/98**
6. Polevoi, A. R. et al., *Reassessment of steady-state operation in ITER with
   NBI and EC heating and current drive*, Nucl. Fusion **60** (2020) 096024
7. Polevoi, A. R. et al., *Impact of pellet injection on ionization source from
   gas puffing and NBI*, 50th EPS Conference on Plasma Physics, P4-092

## Authors

* A. R. Polevoi (ITER Organization)
* M. Hosokawa (ITER Organization)
* M. Schneider (ITER Organization)
* S. D. Pinches (ITER Organization)

## How to cite

If you use SMART in published work, please cite reference [4] above and this
repository (see [`CITATION.cff`](CITATION.cff) for a machine-readable form).

## Contributing

Bug reports, questions and pull requests are welcome — see
[`CONTRIBUTING.md`](CONTRIBUTING.md) for the workflow.

## License

SMART is distributed under the GNU Lesser General Public License, version 3.0
or any later version. See [`COPYING.LESSER`](COPYING.LESSER) (LGPL-3.0) and
[`COPYING`](COPYING) (GPL-3.0, which the LGPL extends).

> The views and opinions expressed herein do not necessarily reflect those of
> the ITER Organization.
