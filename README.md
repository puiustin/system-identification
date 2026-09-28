# System Identification Labs

MATLAB coursework for a System Identification course (_Identificarea Sistemelor_), covering data acquisition, stochastic process analysis, parametric model estimation (ARX/OE/ARMAX/BJ), model structure/order selection, and a final case study identifying the physical parameters of a DC motor.

Each lab folder is self-contained: it includes the lab handout (PDF), the assignment writeup (`.docx`/`.txt`), the helper/library functions used across problems, and the scripts solving each problem.

## Structure

| Folder    | Topic                                                                                                                                                             |
| --------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `LAB I`   | IDDATA/IDSS objects, vectorizing data, generating and importing acquired data (MIMO).                                                                             |
| `LAB II`  | Stochastic processes: noise filtering/shaping, spectral factorization (`spefac`), spectrum analysis (`d_spektr`).                                                 |
| `LAB III` | Parametric identification with ARX and OE models; filtered vs. unfiltered estimation, model comparison.                                                           |
| `LAB IV`  | Least-Squares (LS) and Instrumental Variable (IV) estimation, model validation indexes, structure/order selection via the generalized Akaike criterion (GAIC).    |
| `LAB V`   | ARMAX and Box-Jenkins (BJ) model estimation and validation, GAIC-based order selection for 2/3/4-parameter structures.                                            |
| `Tema XI` | Case study: off-line identification of a DC engine's physical parameters (gain, time constant) from input/output data, including Hz-to-Hs conversion and scaling. |

`LAB IV` has no PDF handout/`.docx` writeup in the repo; its assignment (`tema.m`) and scripts (`ISLAB_5A`–`5D`) are otherwise structured the same as the other labs.

## Requirements

- MATLAB (developed/tested with a version compatible with the `IDDATA`/`IDSS` object API — the older System Identification Toolbox object model).
- MATLAB System Identification Toolbox.

## Usage

Open the desired `LAB *` folder in MATLAB (so its helper functions are on the path) and run the top-level problem scripts, e.g.:

```matlab
cd 'LAB III'
Problema4_1_ARX
```

Most problem scripts (`Problema*.m`, `ISLAB_*.m`, `tema.m`) are runnable directly and will generate/load data, estimate models, and plot validation results. Refer to each lab's PDF handout and assignment document for the underlying theory and requirements.
