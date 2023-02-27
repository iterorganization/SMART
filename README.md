# SMART

IMAS actor for Simplified Mass Ablation and Relocation Treatment of high-field-side/lowfield-side pellet injection (SMART)

## Description

[1] Pellet ablation model by B.Kuteev, NF, 35 (1995) 431

[2] Cloud size by P.B.Parks, Phys. of Plasmas, 7 (2000) 1968

[3] Mass relocation model by H.R.Strauss, Phys. of Plasmas, 5 (1998) 2676

[4] Simplified mass ablation and relocation treatment for pellet injection optimization by A.R. Polevoi and M. Shimada, PPCF, 43 (2001) 1525

## Getting Started

### Installing

* Download the program
```
$ git clone ssh://git@git.iter.org/...
```
* Setup environment and build the actor
```
$ cd smart
$ ml IMAS iWrap XMLlib
$ make
```

### Input Parameters (smart.xml)

| Name | Unit | Description |
| ---- | ---- | ----------- |
| YAM  | AMU  | Pellet mass in Atomic Units (2.5 for 50:50 DT) |
| YVP  | km/s | Vpel  |
| YVOL | mm3  | Pellet size |
| YCOS0| -    | Angle to normal to LCS  |
| YEFF | -    | Efficiency of loss in the tube 1 corresponds to 100% intact pellet |
| YDL  | -    | Ratio of distance along the trajectory to that mapped on the mid-plane |
| yswitch | - | Switch between distinct pellets and permanent pellet model |

### Executing program

#### By using `standalone.exe` 
* Edit `input/standalone.xml` and/or `input/smart.xml` accordingly
* Run the command
```
$ ./standalone.exe
```
#### By using `run_smart.py`
* Edit `input/scenario.yaml` and/or  `input/smart.xml` accordingly
* Run the command
```
$ python ./run_smart.py
```

## Authors

* A.R. Polevoi (PMA/SCD/SCOD, ITER Organization)
* M. Hosokawa (PMA/SCD/SCOD, ITER Organization)
