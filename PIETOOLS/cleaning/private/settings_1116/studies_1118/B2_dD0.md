## B2 dD 0 against control

Rows in both files: 124. Summed over them: variables 1462393 -> 1037045 (0.71x), solves 332 -> 332, seconds 654 -> 527.

### Gained (both gates): 0

### Lost (both gates): 6

| row | status ctl / B2 dD 0 | result | variables | solves | seconds | gates |
|---|---|---|---|---|---|---|
| PDE:29 | certified / unknown | 1.000004 / 12337.791 | 1274 / 994 | 1 / 1 | 0.77 / 1.43 | H 1/0 E 1/0 eta 3.9e-09/2.9e-09 |
| PDE:30 | certified / inaccurate | 0.6499404 / 1.1252811 | 1394 / 1102 | 1 / 1 | 1.57 / 3.03 | H 1/0 E 1/0 eta 1.8e-10/3.1e-12 |
| PDE:31 | certified / unknown | 0.00025642962 / 0.00049278362 | 4387 / 3363 | 1 / 1 | 1.18 / 2.50 | H 1/0 E 1/0 eta 3.2e-11/2.6e-11 |
| PDE:51 | certified / unknown | 1.2699242 / 4012.6832 | 2416 / 1952 | 1 / 1 | 0.82 / 1.45 | H 1/0 E 1/0 eta 2.8e-11/1.1e-08 |
| PDE:52 | certified / unknown | 1.0963826 / 13014.204 | 2523 / 2051 | 1 / 1 | 1.32 / 1.17 | H 1/0 E 1/0 eta 4.9e-12/9.4e-09 |
| PDE:7/stability_dual | certified / infeasible | 1 / 0 | 1368 / 768 | 1 / 1 | 0.59 / 0.29 | H 1/0 E 1/0 eta 4e-12/1 |

### Gates disagree: 9

| row | status ctl / B2 dD 0 | result | variables | solves | seconds | gates |
|---|---|---|---|---|---|---|
| DDE:DDE12 | unknown / certified | 4.9839982 / 4.9972143 | 3837 / 2429 | 1 / 1 | 0.46 / 0.58 | H 0/1 E 1/1 eta 7.3e-10/2.7e-09 |
| DDE:DDE12/Hinf_gain_dual | certified / unknown | 4.9801891 / 4.9852763 | 3837 / 2429 | 1 / 1 | 0.47 / 0.48 | H 1/0 E 1/1 eta 3.1e-09/2.1e-09 |
| DDE:DDE25 | certified / unknown | 3.6746856 / 3.6749151 | 1746 / 1042 | 1 / 1 | 0.58 / 0.78 | H 1/0 E 1/1 eta 6.8e-09/2.1e-07 |
| DDE:DDE26/H2_norm_dual | unknown / certified | 3.8306734 / 3.830769 | 54996 / 31572 | 1 / 1 | 5.73 / 4.19 | H 0/1 E 1/1 eta 6.9e-09/2.8e-08 |
| NDSDDF:NDS11 | certified / unknown | 1 / 0 | 1712 / 1428 | 1 / 1 | 0.32 / 0.27 | H 1/0 E 1/1 eta 1.9e-06/1.1e-06 |
| PDE:18 | certified / unknown | 1 / 0 | 1417 / 977 | 1 / 1 | 0.48 / 0.40 | H 1/0 E 1/1 eta 2.3e-09/9.1e-09 |
| PDE:22 | unknown / certified | 0.20847412 / 0.20847418 | 929 / 577 | 1 / 1 | 0.78 / 0.71 | H 0/1 E 1/1 eta 7.4e-09/6.6e-09 |
| PDE:27 | unknown / certified | 1.0000223 / 1.0000268 | 1180 / 624 | 1 / 1 | 0.66 / 0.88 | H 0/1 E 1/1 eta 2.1e-12/2.3e-12 |
| PDE:49 | certified / unknown | 0.1016557 / 0.10522219 | 1563 / 899 | 1 / 1 | 0.90 / 1.22 | H 1/0 E 1/1 eta 6e-10/8.7e-10 |

### Bisection thresholds that differ: 3

| row | status ctl / B2 dD 0 | result | variables | solves | seconds | change |
|---|---|---|---|---|---|---|
| DDE:DDE01 | converged / converged | 2.0180859 / 2.0200781 | 6196 / 15028 | 11 / 11 | 9.23 / 11.48 | +0.099% |
| NDSDDF:NDS08 | converged / converged | 0.027388574 / 0.027496191 | 24784 / 20368 | 11 / 11 | 32.71 / 18.73 | +0.393% |
| PDE:5 | converged / converged | 9.8696044 / 9.8599661 | 4944 / 1128 | 11 / 11 | 7.29 / 2.40 | -0.098% |

### Other differences: 29

| row | status ctl / B2 dD 0 | result | variables | solves | seconds | gates |
|---|---|---|---|---|---|---|
| DDE:DDE09 | certified / certified | 25.912755 / 25.942216 | 1040 / 668 | 1 / 1 | 0.69 / 0.73 | H 1/1 E 1/1 eta 2.1e-09/1.7e-08 |
| DDE:DDE09/Hinf_gain_dual | certified / certified | 25.909445 / 25.919043 | 1040 / 668 | 1 / 1 | 0.36 / 0.39 | H 1/1 E 1/1 eta 2.9e-09/1.6e-08 |
| DDE:DDE10 | certified / certified | 0.2364431 / 0.2364434 | 3837 / 2429 | 1 / 1 | 0.50 / 0.44 | H 1/1 E 1/1 eta 8.7e-10/1.3e-09 |
| DDE:DDE10/Hinf_gain_dual | certified / certified | 0.23645056 / 0.23646638 | 3837 / 2429 | 1 / 1 | 0.42 / 0.28 | H 1/1 E 1/1 eta 1.5e-10/1.1e-09 |
| DDE:DDE11 | certified / certified | 2.8772529 / 2.8796899 | 1246 / 834 | 1 / 1 | 0.48 / 0.50 | H 1/1 E 1/1 eta 7.8e-10/4e-09 |
| DDE:DDE11/Hinf_gain_dual | certified / certified | 2.8779832 / 2.8815455 | 1246 / 834 | 1 / 1 | 0.37 / 0.50 | H 1/1 E 1/1 eta 9.2e-11/3.2e-09 |
| DDE:DDE13 | certified / certified | 1.3653282 / 1.4469767 | 4983 / 2703 | 1 / 1 | 0.79 / 0.87 | H 1/1 E 1/1 eta 2.2e-11/2.1e-11 |
| DDE:DDE14 | certified / certified | 0.39528506 / 0.39528477 | 4983 / 2703 | 1 / 1 | 0.45 / 0.44 | H 1/1 E 1/1 eta 1e-09/1.6e-10 |
| DDE:DDE15 | certified / certified | 0.61043363 / 0.61043362 | 5173 / 2837 | 1 / 1 | 0.42 / 0.55 | H 1/1 E 1/1 eta 3.6e-10/2.2e-10 |
| DDE:DDE16 | certified / certified | 0.96293737 / 0.96293752 | 86188 / 44176 | 1 / 1 | 8.66 / 5.49 | H 1/1 E 1/1 eta 2e-09/1.9e-09 |
| DDE:DDE17/Hinf_control | certified / certified | 9.3513439 / 9.3513437 | 29146 / 15526 | 1 / 1 | 1.83 / 1.40 | H 1/1 E 1/1 eta 5.8e-10/8.1e-10 |
| DDE:DDE18 | certified / certified | 0.37791699 / 0.41518046 | 13672 / 7912 | 1 / 1 | 1.73 / 1.77 | H 1/1 E 1/1 eta 3.8e-12/3.9e-11 |
| DDE:DDE19 | inaccurate / inaccurate | 0.0013984255 / 0.0012126863 | 2242 / 1242 | 1 / 1 | 1.36 / 1.43 | H 0/0 E 0/0 eta 3e-12/2.7e-12 |
| DDE:DDE20 | certified / certified | 0.23570228 / 0.23570226 | 1626 / 958 | 1 / 1 | 0.74 / 0.77 | H 1/1 E 1/1 eta 1e-09/7.6e-11 |
| DDE:DDE21 | certified / certified | 2.3270003 / 2.327002 | 4983 / 2703 | 1 / 1 | 0.67 / 0.66 | H 1/1 E 1/1 eta 1.1e-09/7.8e-10 |
| DDE:DDE23 | certified / certified | 1.8080605 / 1.8080609 | 10446 / 5526 | 1 / 1 | 1.14 / 0.79 | H 1/1 E 1/1 eta 2.8e-09/3.5e-09 |
| DDE:DDE24 | certified / certified | 1.0439637 / 1.0440831 | 17909 / 9349 | 1 / 1 | 1.32 / 1.48 | H 1/1 E 1/1 eta 9.8e-11/1.2e-10 |
| DDE:DDE25/H2_norm_dual | certified / certified | 3.6723928 / 3.6723991 | 1746 / 1042 | 1 / 1 | 0.45 / 0.53 | H 1/1 E 1/1 eta 1.1e-08/5.1e-09 |
| DDE:DDE26 | certified / certified | 3.8308551 / 3.8308627 | 54996 / 31572 | 1 / 1 | 5.52 / 3.85 | H 1/1 E 1/1 eta 4.1e-08/1.1e-07 |
| NDSDDF:DDF03 | unknown / inaccurate | 0.00052024274 / 0.00066921933 | 7869 / 4157 | 1 / 1 | 0.88 / 0.67 | H 0/0 E 0/0 eta 5.4e-12/8.9e-12 |
| PDE:20 | certified / certified | 0.50081458 / 0.50285584 | 929 / 577 | 1 / 1 | 1.52 / 1.84 | H 1/1 E 1/1 eta 5.1e-11/1.9e-09 |
| PDE:21 | certified / certified | 0.84187746 / 0.8418811 | 3409 / 2081 | 1 / 1 | 1.03 / 1.37 | H 1/1 E 1/1 eta 9.1e-09/4.7e-09 |
| PDE:21/Hinf_gain_dual | certified / certified | 0.84188347 / 0.84188459 | 3409 / 3361 | 1 / 1 | 0.84 / 1.68 | H 1/1 E 1/1 eta 3.1e-11/2.2e-10 |
| PDE:23 | certified / certified | 23.719575 / 23.742847 | 929 / 577 | 1 / 1 | 0.59 / 0.68 | H 1/1 E 1/1 eta 5.1e-10/1.5e-11 |
| PDE:24 | certified / certified | 8.2137035 / 8.2153465 | 929 / 917 | 1 / 1 | 0.55 / 0.79 | H 1/1 E 1/1 eta 2.7e-09/9.4e-09 |
| PDE:25 | certified / certified | 0.81029155 / 0.81037008 | 3409 / 2081 | 1 / 1 | 0.74 / 0.95 | H 1/1 E 1/1 eta 1.1e-09/1.2e-09 |
| PDE:26 | certified / certified | 2.1472936 / 2.1474248 | 7893 / 4845 | 1 / 1 | 1.06 / 1.02 | H 1/1 E 1/1 eta 5e-09/5.4e-09 |
| PDE:28 | certified / certified | 1.0000533 / 1.0000573 | 1180 / 624 | 1 / 1 | 0.49 / 0.56 | H 1/1 E 1/1 eta 2e-12/2.1e-12 |
| PDE:32 | unknown / unknown | 0.5510909 / 5.5521662 | 7249 / 3681 | 1 / 1 | 3.43 / 4.05 | H 0/0 E 0/0 eta 1e-09/1.2e-08 |

