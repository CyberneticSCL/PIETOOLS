## C2 vs light against control

Rows in both files: 124. Summed over them: variables 1462393 -> 882751 (0.60x), solves 332 -> 332, seconds 654 -> 385.

### Gained (both gates): 0

### Lost (both gates): 4

| row | status ctl / C2 vs light | result | variables | solves | seconds | gates |
|---|---|---|---|---|---|---|
| DDE:DDE01/stability_dual | certified / inaccurate | 1 / 0.5 | 5996 / 3532 | 1 / 1 | 1.08 / 0.89 | H 1/0 E 1/0 eta 1.5e-11/7e-12 |
| PDE:30 | certified / inaccurate | 0.6499404 / 0.65068774 | 1394 / 866 | 1 / 1 | 1.57 / 1.28 | H 1/0 E 1/0 eta 1.8e-10/9.5e-11 |
| PDE:31 | certified / unknown | 0.00025642962 / 0.0002810559 | 4387 / 2563 | 1 / 1 | 1.18 / 1.05 | H 1/0 E 1/0 eta 3.2e-11/2.7e-11 |
| PDE:7/stability_dual | certified / unknown | 1 / 0 | 1368 / 780 | 1 / 1 | 0.59 / 0.47 | H 1/0 E 1/0 eta 4e-12/9.5e-11 |

### Gates disagree: 3

| row | status ctl / C2 vs light | result | variables | solves | seconds | gates |
|---|---|---|---|---|---|---|
| NDSDDF:NDS11 | certified / unknown | 1 / 0 | 1712 / 2278 | 1 / 1 | 0.32 / 0.79 | H 1/0 E 1/1 eta 1.9e-06/1.2e-06 |
| PDE:27 | unknown / certified | 1.0000223 / 1.0000233 | 1180 / 700 | 1 / 1 | 0.66 / 0.52 | H 0/1 E 1/1 eta 2.1e-12/2.2e-12 |
| PDE:68 | certified / unknown | 1 / 0 | 12564 / 6714 | 1 / 1 | 1.92 / 2.53 | H 1/0 E 1/1 eta 4.2e-09/8.1e-09 |

### Bisection thresholds that differ: 7

| row | status ctl / C2 vs light | result | variables | solves | seconds | change |
|---|---|---|---|---|---|---|
| DDE:DDE01 | converged / converged | 2.0180859 / 2.0160938 | 6196 / 8348 | 11 / 11 | 9.23 / 9.67 | -0.099% |
| DDE:DDE03 | converged / converged | 1.7161724 / 1.6708775 | 1712 / 2278 | 11 / 11 | 5.50 / 7.18 | -2.639% |
| DDE:DDE04 | converged / converged | 6.1665521 / 6.0580497 | 6196 / 3436 | 11 / 11 | 7.17 / 7.48 | -1.760% |
| DDE:DDE05 | converged / converged | 1.5707 / 1.5660983 | 5233 / 859 | 11 / 11 | 9.50 / 5.75 | -0.293% |
| DDE:DDE07 | converged / converged | 1.3721821 / 1.3478511 | 6196 / 3436 | 11 / 11 | 12.37 / 7.86 | -1.773% |
| DDE:DDE08 | converged / converged | 1.42229 / 1.3435913 | 6160 / 1258 | 12 / 12 | 8.39 / 8.69 | -5.533% |
| NDSDDF:NDS08 | converged / converged | 0.027388574 / 0.027523096 | 24784 / 13744 | 11 / 11 | 32.71 / 14.81 | +0.491% |

### Other differences: 15

| row | status ctl / C2 vs light | result | variables | solves | seconds | gates |
|---|---|---|---|---|---|---|
| DDE:DDE13 | certified / certified | 1.3653282 / 1.4303312 | 4983 / 3015 | 1 / 1 | 0.79 / 0.97 | H 1/1 E 1/1 eta 2.2e-11/1.8e-11 |
| DDE:DDE14 | certified / certified | 0.39528506 / 0.39528477 | 4983 / 3015 | 1 / 1 | 0.45 / 0.68 | H 1/1 E 1/1 eta 1e-09/1.6e-10 |
| DDE:DDE15 | certified / certified | 0.61043363 / 0.61043373 | 5173 / 3157 | 1 / 1 | 0.42 / 0.59 | H 1/1 E 1/1 eta 3.6e-10/4.3e-10 |
| DDE:DDE16 | certified / certified | 0.96293737 / 0.96293763 | 86188 / 49900 | 1 / 1 | 8.66 / 6.60 | H 1/1 E 1/1 eta 2e-09/2.4e-09 |
| DDE:DDE17/Hinf_control | certified / certified | 9.3513439 / 9.3513437 | 29146 / 17386 | 1 / 1 | 1.83 / 1.17 | H 1/1 E 1/1 eta 5.8e-10/7.9e-10 |
| DDE:DDE18 | certified / certified | 0.37791699 / 0.3922283 | 13672 / 8704 | 1 / 1 | 1.73 / 1.49 | H 1/1 E 1/1 eta 3.8e-12/2.8e-11 |
| DDE:DDE19 | inaccurate / inaccurate | 0.0013984255 / 0.0013589698 | 2242 / 994 | 1 / 1 | 1.36 / 0.96 | H 0/0 E 0/0 eta 3e-12/2.9e-12 |
| DDE:DDE21 | certified / certified | 2.3270003 / 2.3270005 | 4983 / 3015 | 1 / 1 | 0.67 / 0.57 | H 1/1 E 1/1 eta 1.1e-09/1.1e-09 |
| DDE:DDE23 | certified / certified | 1.8080605 / 1.8080604 | 10446 / 6198 | 1 / 1 | 1.14 / 0.79 | H 1/1 E 1/1 eta 2.8e-09/7.5e-10 |
| DDE:DDE24 | certified / certified | 1.0439637 / 1.0440079 | 17909 / 10517 | 1 / 1 | 1.32 / 1.08 | H 1/1 E 1/1 eta 9.8e-11/2.2e-10 |
| NDSDDF:DDF03 | unknown / unknown | 0.00052024274 / 0.00055921051 | 7869 / 3261 | 1 / 1 | 0.88 / 0.70 | H 0/0 E 0/0 eta 5.4e-12/6.6e-12 |
| PDE:28 | certified / certified | 1.0000533 / 1.0000321 | 1180 / 700 | 1 / 1 | 0.49 / 0.39 | H 1/1 E 1/1 eta 2e-12/1.3e-12 |
| PDE:32 | unknown / unknown | 0.5510909 / 0.59829485 | 7249 / 2833 | 1 / 1 | 3.43 / 1.74 | H 0/0 E 0/0 eta 1e-09/6.3e-09 |
| PDE:51 | certified / certified | 1.2699242 / 1.2707136 | 2416 / 1344 | 1 / 1 | 0.82 / 0.85 | H 1/1 E 1/1 eta 2.8e-11/2.3e-11 |
| PDE:52 | certified / certified | 1.0963826 / 1.0974921 | 2523 / 1423 | 1 / 1 | 1.32 / 0.81 | H 1/1 E 1/1 eta 4.9e-12/3.9e-12 |

