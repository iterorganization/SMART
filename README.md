# SMART

IMAS actor for Simplified Mass Ablation and Relocation Treatment of high-field-side/lowfield-side pellet injection (SMART)

## Getting Started

### Installing

* Download the program
```
$ git clone ssh://git@git.iter.org/fuel/smart.git
```
* Setup environment and build the actor
```
$ cd smart
$ source config_sdcc.sh
$ make
$ make actor
```

### Input Parameters (smart.xml)

#### STEPUP
| Name | Unit | Description |
| ---- | ---- | ----------- |
| TAU  | s    | Time step   |
| dtau | s    | Time interval between consequential pellet injection   |
| QNB  | 10<sup>19</sup> s<sup>-1</sup> | Boundary influx of neutrals (D+T+H) for transport equation in CX approximation |
| F01B | 10<sup>19</sup> m<sup>-3</sup> | Edge neutral density (Puffing, H component) |
| F02B | 10<sup>19</sup> m<sup>-3</sup> | Edge neutral density (Puffing, D component) |
| F03B | 10<sup>19</sup> m<sup>-3</sup> | Edge neutral density (Puffing, T component) |
| GN2E | - | Coefficient for convective term in heat transport (0/1.5/2.5) for electron |
| GN2I | - | Coefficient for convective term in heat transport (0/1.5/2.5) for ion |
| GN2I | - | Coefficient for convective term in heat transport (0/1.5/2.5) for ion |
| sw_stdout | - | Toggle for debug information to stdout (Off when =0, On when !=0) |

#### SMART
| Name | Unit | Description |
| ---- | ---- | ----------- |
| YAM  | AMU  | Pellet mass in Atomic Units (2.5 for 50:50 DT) |
| YVP  | km/s | Vpel  |
| YVOL | mm<sup>3</sup>  | Pellet size |
| YCOS0| -    | Angle to normal to LCS  |
| YEFF | -    | Efficiency of loss in the tube 1 corresponds to 100% intact pellet |
| YDL  | -    | Ratio of distance along the trajectory to that mapped on the mid-plane |
| yswitch | - | Switch between distinct pellets and permanent pellet model |
| sw_smart | - | Toggle for SMART (Off when =0, On when !=0) |

#### ECH2a
| Name | Unit | Description |
| ---- | ---- | ----------- |
| ROCEC | - | Normalized EC location w.r.t. toroidal flux coordinate = sqrt(phi/(pi\*b0)) |
| ROCDR | - | Normalized EC width w.r.t. toroidal flux coordinate = sqrt(phi/(pi\*b0)) |
| QECR | MW | ECRH Power |
| YEFFec | MA/MW | IEC/QEC |
| sw_ech2a | - | Toggle for ECH2A (Off when =0, On when !=0) |

### Executing program

#### The wrapper code `smart` 
* Edit `input/standalone.xml` and/or `input/smart.xml` accordingly
* Run the command
```
$ ./smart
```
#### The IMAS actor using `run_smart.py`
* Edit `input/scenario.yaml` and/or  `input/smart.xml` accordingly
* Run the command
```
$ python ./run_smart.py
```

## Reference

[1] Pellet ablation model by B.Kuteev, NF, 35 (1995) 431

[2] Cloud size by P.B.Parks, Phys. of Plasmas, 7 (2000) 1968

[3] Mass relocation model by H.R.Strauss, Phys. of Plasmas, 5 (1998) 2676

[4] Simplified mass ablation and relocation treatment for pellet injection optimization by A.R. Polevoi and M. Shimada, PPCF, 43 (2001) 1525

[5] ASTRA automated system for transport analysis in a tokamak by G.V. Pereverzev and P.N. Yushmanov, Max-Planck IPP Report vol 5/98

[6] Reassessment of steady-state operation in ITER with NBI and EC heating and current drive by A.R. Polevoi et al, NF, 60 (2020) 096024

[7] Impact of pellet injection on ionization source from gas puffing and NBI by A.R. Polevoi et al, 50th EPS Conference on Plasma Physics, P4-092

## Authors

* A.R. Polevoi (PMA/SCD/SID, ITER Organization)
* M. Hosokawa (PMA/SCD/SID, ITER Organization)
