## C3 vs light against control

Rows in both files: 124. Summed over them: variables 1462393 -> 684275 (0.47x), solves 332 -> 332, seconds 654 -> 351.

### Gained (both gates): 0

### Lost (both gates): 9

| row | status ctl / C3 vs light | result | variables | solves | seconds | gates |
|---|---|---|---|---|---|---|
| DDE:DDE01/stability_dual | certified / infeasible | 1 / 0 | 5996 / 2204 | 1 / 1 | 1.08 / 0.50 | H 1/0 E 1/0 eta 1.5e-11/1 |
| PDE:21/Hinf_gain_dual | certified / unknown | 0.84188347 / 0.8859611 | 3409 / 2081 | 1 / 1 | 0.84 / 0.76 | H 1/0 E 1/0 eta 3.1e-11/2.3e-09 |
| PDE:24 | certified / inaccurate | 8.2137035 / 8.8517602 | 929 / 577 | 1 / 1 | 0.55 / 0.34 | H 1/0 E 1/0 eta 2.7e-09/2.3e-08 |
| PDE:29 | certified / unknown | 1.000004 / 43707.218 | 1274 / 450 | 1 / 1 | 0.77 / 0.85 | H 1/0 E 1/0 eta 3.9e-09/3.4e-09 |
| PDE:30 | certified / infeasible | 0.6499404 / 0.00066628506 | 1394 / 526 | 1 / 1 | 1.57 / 1.70 | H 1/0 E 1/0 eta 1.8e-10/1 |
| PDE:31 | certified / unknown | 0.00025642962 / 0.00062185726 | 4387 / 1443 | 1 / 1 | 1.18 / 1.23 | H 1/0 E 1/0 eta 3.2e-11/7.1e-11 |
| PDE:51 | certified / unknown | 1.2699242 / 3692.216 | 2416 / 732 | 1 / 1 | 0.82 / 0.89 | H 1/0 E 1/0 eta 2.8e-11/1.4e-08 |
| PDE:52 | certified / unknown | 1.0963826 / 11368.181 | 2523 / 791 | 1 / 1 | 1.32 / 0.62 | H 1/0 E 1/0 eta 4.9e-12/5.6e-09 |
| PDE:7/stability_dual | certified / infeasible | 1 / 0 | 1368 / 468 | 1 / 1 | 0.59 / 0.33 | H 1/0 E 1/0 eta 4e-12/1 |

### Gates disagree: 9

| row | status ctl / C3 vs light | result | variables | solves | seconds | gates |
|---|---|---|---|---|---|---|
| DDE:DDE12 | unknown / certified | 4.9839982 / 4.9972143 | 3837 / 2429 | 1 / 1 | 0.46 / 0.48 | H 0/1 E 1/1 eta 7.3e-10/2.7e-09 |
| DDE:DDE12/Hinf_gain_dual | certified / unknown | 4.9801891 / 4.9852763 | 3837 / 2429 | 1 / 1 | 0.47 / 0.37 | H 1/0 E 1/1 eta 3.1e-09/2.1e-09 |
| DDE:DDE25 | certified / unknown | 3.6746856 / 3.6749151 | 1746 / 1042 | 1 / 1 | 0.58 / 0.65 | H 1/0 E 1/1 eta 6.8e-09/2.1e-07 |
| DDE:DDE26/H2_norm_dual | unknown / certified | 3.8306734 / 3.830769 | 54996 / 31572 | 1 / 1 | 5.73 / 3.94 | H 0/1 E 1/1 eta 6.9e-09/2.8e-08 |
| PDE:18 | certified / unknown | 1 / 0 | 1417 / 977 | 1 / 1 | 0.48 / 0.50 | H 1/0 E 1/1 eta 2.3e-09/9.1e-09 |
| PDE:22 | unknown / certified | 0.20847412 / 0.20847418 | 929 / 577 | 1 / 1 | 0.78 / 0.61 | H 0/1 E 1/1 eta 7.4e-09/6.6e-09 |
| PDE:27 | unknown / certified | 1.0000223 / 1.0000431 | 1180 / 400 | 1 / 1 | 0.66 / 0.75 | H 0/1 E 1/1 eta 2.1e-12/4e-12 |
| PDE:28 | certified / unknown | 1.0000533 / 1.000072 | 1180 / 400 | 1 / 1 | 0.49 / 0.42 | H 1/0 E 1/1 eta 2e-12/1.6e-12 |
| PDE:49 | certified / unknown | 0.1016557 / 0.10522219 | 1563 / 899 | 1 / 1 | 0.90 / 0.92 | H 1/0 E 1/1 eta 6e-10/8.7e-10 |

### Bisection thresholds that differ: 7

| row | status ctl / C3 vs light | result | variables | solves | seconds | change |
|---|---|---|---|---|---|---|
| DDE:DDE01 | converged / converged | 2.0180859 / 2.0160938 | 6196 / 6892 | 11 / 11 | 9.23 / 8.90 | -0.099% |
| DDE:DDE03 | converged / converged | 1.7161724 / 1.6708775 | 1712 / 1898 | 11 / 11 | 5.50 / 4.92 | -2.639% |
| DDE:DDE04 | converged / converged | 6.1665521 / 6.0580497 | 6196 / 3388 | 11 / 11 | 7.17 / 5.40 | -1.760% |
| DDE:DDE05 | converged / converged | 1.5707 / 1.5660983 | 5233 / 847 | 11 / 11 | 9.50 / 4.07 | -0.293% |
| DDE:DDE07 | converged / converged | 1.3721821 / 1.3478511 | 6196 / 3388 | 11 / 11 | 12.37 / 7.82 | -1.773% |
| DDE:DDE08 | converged / converged | 1.42229 / 1.3435913 | 6160 / 1246 | 12 / 12 | 8.39 / 8.19 | -5.533% |
| NDSDDF:NDS08 | converged / converged | 0.027388574 / 0.027415479 | 24784 / 13552 | 11 / 11 | 32.71 / 14.22 | +0.098% |

### Other differences: 25

| row | status ctl / C3 vs light | result | variables | solves | seconds | gates |
|---|---|---|---|---|---|---|
| DDE:DDE09 | certified / certified | 25.912755 / 25.942216 | 1040 / 668 | 1 / 1 | 0.69 / 0.77 | H 1/1 E 1/1 eta 2.1e-09/1.7e-08 |
| DDE:DDE09/Hinf_gain_dual | certified / certified | 25.909445 / 25.919043 | 1040 / 668 | 1 / 1 | 0.36 / 0.39 | H 1/1 E 1/1 eta 2.9e-09/1.6e-08 |
| DDE:DDE10 | certified / certified | 0.2364431 / 0.2364434 | 3837 / 2429 | 1 / 1 | 0.50 / 0.63 | H 1/1 E 1/1 eta 8.7e-10/1.3e-09 |
| DDE:DDE10/Hinf_gain_dual | certified / certified | 0.23645056 / 0.23646638 | 3837 / 2429 | 1 / 1 | 0.42 / 0.53 | H 1/1 E 1/1 eta 1.5e-10/1.1e-09 |
| DDE:DDE11 | certified / certified | 2.8772529 / 2.8796899 | 1246 / 834 | 1 / 1 | 0.48 / 0.45 | H 1/1 E 1/1 eta 7.8e-10/4e-09 |
| DDE:DDE11/Hinf_gain_dual | certified / certified | 2.8779832 / 2.8815455 | 1246 / 834 | 1 / 1 | 0.37 / 0.46 | H 1/1 E 1/1 eta 9.2e-11/3.2e-09 |
| DDE:DDE13 | certified / certified | 1.3653282 / 1.4591196 | 4983 / 1775 | 1 / 1 | 0.79 / 1.06 | H 1/1 E 1/1 eta 2.2e-11/2.6e-12 |
| DDE:DDE14 | certified / certified | 0.39528506 / 0.3952848 | 4983 / 1775 | 1 / 1 | 0.45 / 0.59 | H 1/1 E 1/1 eta 1e-09/2.6e-10 |
| DDE:DDE15 | certified / certified | 0.61043363 / 0.6104337 | 5173 / 1877 | 1 / 1 | 0.42 / 0.64 | H 1/1 E 1/1 eta 3.6e-10/8.2e-10 |
| DDE:DDE16 | certified / certified | 0.96293737 / 0.96293831 | 86188 / 27760 | 1 / 1 | 8.66 / 4.22 | H 1/1 E 1/1 eta 2e-09/2.4e-09 |
| DDE:DDE18 | certified / certified | 0.37791699 / 0.41809405 | 13672 / 5464 | 1 / 1 | 1.73 / 1.25 | H 1/1 E 1/1 eta 3.8e-12/7.4e-11 |
| DDE:DDE19 | inaccurate / inaccurate | 0.0013984255 / 0.0013963761 | 2242 / 634 | 1 / 1 | 1.36 / 0.86 | H 0/0 E 0/0 eta 3e-12/3e-12 |
| DDE:DDE20 | certified / certified | 0.23570228 / 0.23570226 | 1626 / 670 | 1 / 1 | 0.74 / 0.64 | H 1/1 E 1/1 eta 1e-09/7.6e-11 |
| DDE:DDE21 | certified / certified | 2.3270003 / 2.3270476 | 4983 / 1775 | 1 / 1 | 0.67 / 0.72 | H 1/1 E 1/1 eta 1.1e-09/5.5e-09 |
| DDE:DDE23 | certified / certified | 1.8080605 / 1.808062 | 10446 / 3558 | 1 / 1 | 1.14 / 0.59 | H 1/1 E 1/1 eta 2.8e-09/1.9e-09 |
| DDE:DDE24 | certified / certified | 1.0439637 / 1.0443303 | 17909 / 5957 | 1 / 1 | 1.32 / 0.94 | H 1/1 E 1/1 eta 9.8e-11/1.2e-10 |
| DDE:DDE25/H2_norm_dual | certified / certified | 3.6723928 / 3.6723991 | 1746 / 1042 | 1 / 1 | 0.45 / 0.26 | H 1/1 E 1/1 eta 1.1e-08/5.1e-09 |
| DDE:DDE26 | certified / certified | 3.8308551 / 3.8308627 | 54996 / 31572 | 1 / 1 | 5.52 / 3.46 | H 1/1 E 1/1 eta 4.1e-08/1.1e-07 |
| NDSDDF:DDF03 | unknown / unknown | 0.00052024274 / 0.00068329739 | 7869 / 1981 | 1 / 1 | 0.88 / 0.75 | H 0/0 E 0/0 eta 5.4e-12/7.3e-12 |
| PDE:20 | certified / certified | 0.50081458 / 0.50285584 | 929 / 577 | 1 / 1 | 1.52 / 1.30 | H 1/1 E 1/1 eta 5.1e-11/1.9e-09 |
| PDE:21 | certified / certified | 0.84187746 / 0.8418811 | 3409 / 2081 | 1 / 1 | 1.03 / 1.14 | H 1/1 E 1/1 eta 9.1e-09/4.7e-09 |
| PDE:23 | certified / certified | 23.719575 / 23.742847 | 929 / 577 | 1 / 1 | 0.59 / 0.40 | H 1/1 E 1/1 eta 5.1e-10/1.5e-11 |
| PDE:25 | certified / certified | 0.81029155 / 0.81037008 | 3409 / 2081 | 1 / 1 | 0.74 / 0.65 | H 1/1 E 1/1 eta 1.1e-09/1.2e-09 |
| PDE:26 | certified / certified | 2.1472936 / 2.1474248 | 7893 / 4845 | 1 / 1 | 1.06 / 0.58 | H 1/1 E 1/1 eta 5e-09/5.4e-09 |
| PDE:32 | unknown / unknown | 0.5510909 / 4030.3693 | 7249 / 1633 | 1 / 1 | 3.43 / 2.59 | H 0/0 E 0/0 eta 1e-09/3.3e-09 |

