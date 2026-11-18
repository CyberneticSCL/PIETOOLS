## Summary per tier (11/17 pass vs 10/28 sopvar reference)
| tier | rows | same | diff | new | certified | infeasible | unknown | error | skipped | seconds |
|---|---|---|---|---|---|---|---|---|---|---|
| heavy | 158 | 15 | 123 | 20 | 109 | 4 | 15 | 26 | 0 | 3515 |
| light | 158 | 16 | 122 | 20 | 110 | 8 | 15 | 21 | 0 | 2270 |

## Rows that differ from the reference
| library | id | tier | reference | 11/17 | difference |
|---|---|---|---|---|---|
| DDE | DDE01/stability_dual | heavy | certified 1 6100/572 0.48s | certified 1 6100/559 0.46s | m 572->559 |
| DDE | DDE01/stability_dual | light | certified 1 5996/556 0.86s | certified 1 5996/545 0.89s | m 556->545 |
| DDE | DDE01 | heavy | converged 1.9085156 6100/580 6.42s | converged 1.8587109 6100/567 5.94s | m 580->567; value 1.9085156->1.8587109 (2.6e-02 rel) |
| DDE | DDE01 | light | converged 2.0160938 5996/568 7.49s | converged 2.0101172 3532/499 6.08s | ndv 5996->3532; m 568->499; value 2.0160938->2.0101172 (3.0e-03 rel) |
| DDE | DDE02/stability_dual | heavy | certified 1 23312/2096 0.80s | certified 1 23312/2023 0.64s | m 2096->2023 |
| DDE | DDE02/stability_dual | light | certified 1 6800/1616 0.48s | certified 1 6800/1573 0.38s | m 1616->1573 |
| DDE | DDE02 | heavy | converged 2.0839844 23312/2112 32.61s | converged 2.0820312 23312/2039 30.68s | m 2112->2039; value 2.0839844->2.0820312 (9.4e-04 rel) |
| DDE | DDE02 | light | converged 2.0839844 6800/1624 14.18s | converged 2.0839844 6800/1581 14.53s | m 1624->1581 |
| DDE | DDE03/stability_dual | heavy | certified 1 1667/167 0.28s | certified 1 1667/166 0.27s | m 167->166 |
| DDE | DDE03/stability_dual | light | certified 1 545/132 0.31s | certified 1 545/131 0.19s | m 132->131 |
| DDE | DDE03 | heavy | converged 1.7161724 1667/170 4.04s | converged 1.7161724 1667/169 2.92s | m 170->169 |
| DDE | DDE03 | light | converged 1.7144948 1637/165 4.37s | converged 1.7144948 1637/164 3.31s | m 165->164 |
| DDE | DDE04/stability_dual | heavy | certified 1 6100/568 0.41s | certified 1 6100/555 0.24s | m 568->555 |
| DDE | DDE04/stability_dual | light | certified 1 1852/444 0.27s | certified 1 1852/436 0.24s | m 444->436 |
| DDE | DDE04 | heavy | converged 6.17258 13380/700 12.13s | converged 6.17258 13380/683 10.81s | m 700->683 |
| DDE | DDE04 | light | converged 6.1665521 3532/508 4.64s | converged 6.1665521 3532/499 4.51s | m 508->499 |
| DDE | DDE06/stability_dual | heavy | certified 1 5828/522 0.34s | certified 1 5828/510 0.24s | m 522->510 |
| DDE | DDE06/stability_dual | light | certified 1 1700/402 0.25s | certified 1 1700/395 0.20s | m 402->395 |
| DDE | DDE06 | heavy | converged 3 12996/644 11.44s | converged 3 12996/628 10.62s | m 644->628 |
| DDE | DDE06 | light | converged 3 5732/518 7.11s | converged 3 5732/508 7.32s | m 518->508 |
| DDE | DDE07/stability_dual | heavy | certified 1 6100/566 0.32s | certified 1 6100/553 0.23s | m 566->553 |
| DDE | DDE07/stability_dual | light | certified 1 1852/438 0.26s | certified 1 1852/430 0.22s | m 438->430 |
| DDE | DDE07 | heavy | converged 1.3713711 13380/692 9.85s | converged 1.3713711 13380/675 8.24s | m 692->675 |
| DDE | DDE07 | light | converged 1.3713711 5996/558 5.57s | converged 1.3713711 5996/547 5.24s | m 558->547 |
| DDE | DDE08/stability_dual | heavy | certified 1 1975/229 0.32s | certified 1 1975/223 0.18s | m 229->223 |
| DDE | DDE08/stability_dual | light | certified 1 733/186 0.23s | certified 1 733/180 0.15s | m 186->180 |
| DDE | DDE08 | heavy | converged 1.42229 3963/268 6.17s | converged 1.42229 3963/262 4.75s | m 268->262 |
| DDE | DDE08 | light | converged 1.4185425 1237/205 4.72s | converged 1.4185425 1237/199 4.29s | m 205->199 |
| DDE | DDE10/Hinf_gain_dual | heavy | certified 0.23644362 6517/608 0.46s | certified 0.23644372 6517/595 0.37s | m 608->595 |
| DDE | DDE10/Hinf_gain_dual | light | certified 0.23645073 3837/532 0.37s | certified 0.23645055 3837/523 0.32s | m 532->523 |
| DDE | DDE10 | heavy | certified 0.23644325 6517/620 0.48s | certified 0.23644327 6517/607 0.47s | m 620->607 |
| DDE | DDE10 | light | certified 0.23644309 3837/544 0.48s | certified 0.23644309 3837/535 0.43s | m 544->535 |
| DDE | DDE11/Hinf_gain_dual | heavy | unknown 2.8763708 2004/205 0.44s | certified 2.876346 2004/203 0.63s | status unknown->certified; m 205->203; value 2.8763708->2.876346 (8.6e-06 rel) |
| DDE | DDE11/Hinf_gain_dual | light | certified 2.8779855 1246/181 0.35s | certified 2.8779834 1246/179 0.33s | m 181->179 |
| DDE | DDE11 | heavy | unknown 2.8758942 2004/208 0.77s | certified 2.8758939 2004/206 0.57s | status unknown->certified; m 208->206 |
| DDE | DDE11 | light | certified 2.8772544 1246/184 0.50s | certified 2.8772595 1246/182 0.36s | m 184->182; value 2.8772544->2.8772595 (1.8e-06 rel) |
| DDE | DDE12/Hinf_gain_dual | heavy | certified 4.9793913 6517/606 0.61s | certified 4.9794268 6517/593 0.44s | m 606->593; value 4.9793913->4.9794268 (7.1e-06 rel) |
| DDE | DDE12/Hinf_gain_dual | light | certified 4.9796898 3837/528 0.41s | certified 4.9796696 3837/519 0.40s | m 528->519; value 4.9796898->4.9796696 (4.1e-06 rel) |
| DDE | DDE12 | heavy | certified 4.9793843 6517/612 0.58s | certified 4.9794672 6517/599 0.41s | m 612->599; value 4.9793843->4.9794672 (1.7e-05 rel) |
| DDE | DDE12 | light | certified 4.9840028 3837/534 0.41s | certified 4.9835374 3837/525 0.45s | m 534->525; value 4.9840028->4.9835374 (9.3e-05 rel) |
| DDE | DDE13 | heavy | certified 1.3531041 7623/360 0.93s | certified 1.3527283 7623/357 0.81s | m 360->357; value 1.3531041->1.3527283 (2.8e-04 rel) |
| DDE | DDE13 | light | certified 1.3649751 4983/302 0.79s | certified 1.3657495 4983/299 0.70s | m 302->299; value 1.3649751->1.3657495 (5.7e-04 rel) |
| DDE | DDE14 | heavy | certified 0.39528489 7623/360 0.50s | certified 0.39528489 7623/357 0.32s | m 360->357 |
| DDE | DDE14 | light | certified 0.39528506 4983/302 0.41s | certified 0.39528506 4983/299 0.33s | m 302->299 |
| DDE | DDE15 | heavy | certified 0.61043361 7861/386 0.53s | certified 0.61043361 7861/382 0.39s | m 386->382 |
| DDE | DDE15 | light | certified 0.61043363 5173/326 0.46s | certified 0.61043363 5173/322 0.36s | m 326->322 |
| DDE | DDE16 | heavy | certified 0.96293716 135913/5615 10.26s | certified 0.96293716 135913/5602 9.46s | m 5615->5602 |
| DDE | DDE16 | light | certified 0.96293737 86188/4589 5.45s | certified 0.96293737 86188/4576 5.60s | m 4589->4576 |
| DDE | DDE17/Hinf_control | heavy | certified 9.3513437 45101/1987 1.91s | certified 9.351344 45101/1970 1.44s | m 1987->1970 |
| DDE | DDE17/Hinf_control | light | certified 9.3513439 29146/1647 1.21s | certified 9.3513439 29146/1630 1.04s | m 1647->1630 |
| DDE | DDE17 | heavy | converged 0.14833984 84599/4522 89.30s | converged 0.14833984 84599/4287 76.19s | m 4522->4287 |
| DDE | DDE17 | light | converged 0.14833984 38149/3682 39.71s | converged 0.14833984 38149/3492 38.05s | m 3682->3492 |
| DDE | DDE18 | heavy | unknown 0.37579925 20179/1126 1.27s | certified 0.3750075 20179/1090 1.42s | status unknown->certified; m 1126->1090; value 0.37579925->0.3750075 (2.1e-03 rel) |
| DDE | DDE18 | light | certified 0.3791484 13672/973 1.01s | certified 0.37801951 13672/937 0.96s | m 973->937; value 0.3791484->0.37801951 (3.0e-03 rel) |
| DDE | DDE19 | heavy | error NaN NaN/NaN 0.12s | unknown 0.00139538 3111/141 0.86s | status error->unknown; ndv NaN->3111; m NaN->141; value NaN->0.00139538 |
| DDE | DDE19 | light | error NaN NaN/NaN 0.11s | unknown 0.0014027311 2242/124 0.77s | status error->unknown; ndv NaN->2242; m NaN->124; value NaN->0.0014027311 |
| DDE | DDE20 | heavy | certified 0.23570229 2373/142 0.56s | certified 0.23570229 2373/139 0.35s | m 142->139 |
| DDE | DDE20 | light | certified 0.23570228 1626/124 0.51s | certified 0.23570228 1626/121 0.35s | m 124->121 |
| DDE | DDE21 | heavy | certified 2.3270004 7623/360 0.59s | certified 2.3270004 7623/357 0.39s | m 360->357 |
| DDE | DDE21 | light | certified 2.3270003 4983/302 0.43s | certified 2.3270003 4983/299 0.35s | m 302->299 |
| DDE | DDE22 | heavy | certified 1.3500949 27965/1206 1.03s | certified 1.3500949 27965/1202 0.79s | m 1206->1202 |
| DDE | DDE22 | light | certified 1.3500949 17909/994 0.70s | certified 1.3500949 17909/990 0.61s | m 994->990 |
| DDE | DDE23 | heavy | certified 1.8080604 16197/720 0.74s | certified 1.8080605 16197/715 0.54s | m 720->715 |
| DDE | DDE23 | light | certified 1.8080604 10446/597 0.60s | certified 1.8080605 10446/592 0.47s | m 597->592 |
| DDE | DDE24 | heavy | certified 1.0439678 27965/1206 1.35s | certified 1.0439681 27965/1202 0.93s | m 1206->1202 |
| DDE | DDE24 | light | certified 1.0439636 17909/994 1.22s | certified 1.0439644 17909/990 0.69s | m 994->990 |
| DDE | DDE26/H2_norm_dual | heavy | certified 3.8306405 99540/7045 6.21s | certified 3.8306126 99540/6754 6.18s | m 7045->6754; value 3.8306405->3.8306126 (7.3e-06 rel) |
| DDE | DDE26/H2_norm_dual | light | certified 3.8306886 54996/5899 3.31s | certified 3.830674 54996/5698 4.18s | m 5899->5698; value 3.8306886->3.830674 (3.8e-06 rel) |
| DDE | DDE26 | heavy | certified 3.8305602 99540/7098 5.96s | certified 3.8305426 99540/6807 5.69s | m 7098->6807; value 3.8305602->3.8305426 (4.6e-06 rel) |
| DDE | DDE26 | light | unknown 3.8308423 54996/5952 3.34s | unknown 3.8308533 54996/5751 3.04s | m 5952->5751; value 3.8308423->3.8308533 (2.9e-06 rel) |
| NDSDDF | DDF01/stability_dual | heavy | certified 1 5564/496 0.27s | certified 1 5564/476 0.14s | m 496->476 |
| NDSDDF | DDF01/stability_dual | light | certified 1 1552/376 0.17s | certified 1 1552/362 0.15s | m 376->362 |
| NDSDDF | DDF01 | heavy | certified 1 5564/496 0.33s | certified 1 5564/476 0.21s | m 496->476 |
| NDSDDF | DDF01 | light | certified 1 1552/376 0.27s | certified 1 1552/362 0.23s | m 376->362 |
| NDSDDF | DDF03 | heavy | error NaN NaN/NaN 0.13s | inaccurate 0.00065446149 11145/441 1.06s | status error->inaccurate; ndv NaN->11145; m NaN->441; value NaN->0.00065446149 |
| NDSDDF | DDF03 | light | error NaN NaN/NaN 0.10s | inaccurate 0.00054836065 7869/381 1.03s | status error->inaccurate; ndv NaN->7869; m NaN->381; value NaN->0.00054836065 |
| NDSDDF | NDS01/stability_dual | heavy | certified 1 6100/584 0.32s | certified 1 6100/571 0.25s | m 584->571 |
| NDSDDF | NDS01/stability_dual | light | certified 1 1852/460 0.25s | certified 1 1852/452 0.20s | m 460->452 |
| NDSDDF | NDS01 | heavy | converged 2.2233267 6100/596 5.35s | converged 2.2233267 6100/583 3.97s | m 596->583 |
| NDSDDF | NDS01 | light | converged 2.2233267 1852/468 3.97s | converged 2.2233267 1852/460 3.38s | m 468->460 |
| NDSDDF | NDS02/stability_dual | heavy | certified 1 6100/584 0.34s | certified 1 6100/555 0.26s | m 584->555 |
| NDSDDF | NDS02/stability_dual | light | certified 1 1852/460 0.24s | certified 1 1852/436 0.27s | m 460->436 |
| NDSDDF | NDS02 | heavy | certified 1 6100/596 0.35s | certified 1 6100/567 0.26s | m 596->567 |
| NDSDDF | NDS02 | light | certified 1 1852/468 0.24s | certified 1 1852/444 0.22s | m 468->444 |
| NDSDDF | NDS03/stability_dual | heavy | certified 1 23312/2152 0.72s | certified 1 23312/2079 0.66s | m 2152->2079 |
| NDSDDF | NDS03/stability_dual | light | certified 1 6800/1672 0.39s | certified 1 6800/1629 0.41s | m 1672->1629 |
| NDSDDF | NDS03 | heavy | certified 1 23312/2176 0.77s | certified 1 23312/2103 0.68s | m 2176->2103 |
| NDSDDF | NDS03 | light | certified 1 6800/1688 0.42s | certified 1 6800/1645 0.47s | m 1688->1645 |
| NDSDDF | NDS04/stability_dual | heavy | certified 1 23312/2152 0.72s | certified 1 23312/2079 0.64s | m 2152->2079 |
| NDSDDF | NDS04/stability_dual | light | certified 1 6800/1672 0.41s | certified 1 6800/1629 0.43s | m 1672->1629 |
| NDSDDF | NDS04 | heavy | hi_pass 0.69 23312/2176 1.42s | hi_pass 0.69 23312/2103 1.24s | m 2176->2103 |
| NDSDDF | NDS04 | light | hi_pass 0.69 6800/1688 0.79s | hi_pass 0.69 6800/1645 0.90s | m 1688->1645 |
| NDSDDF | NDS05/stability_dual | heavy | certified 1 23312/2120 0.96s | certified 1 23312/2047 0.83s | m 2120->2047 |
| NDSDDF | NDS05/stability_dual | light | certified 1 6800/1640 0.44s | certified 1 6800/1597 0.39s | m 1640->1597 |
| NDSDDF | NDS05 | heavy | converged 1.5797836 51984/2640 41.97s | converged 1.5797836 51984/2511 40.45s | m 2640->2511 |
| NDSDDF | NDS05 | light | converged 1.5797836 22928/2136 19.13s | converged 1.5797836 22928/2043 21.06s | m 2136->2043 |
| NDSDDF | NDS07/stability_dual | heavy | certified 1 5828/530 0.31s | certified 1 5828/518 0.22s | m 530->518 |
| NDSDDF | NDS07/stability_dual | light | certified 1 1700/410 0.25s | certified 1 1700/403 0.21s | m 410->403 |
| NDSDDF | NDS07 | heavy | converged 0.44453613 12996/660 8.37s | converged 0.44453613 12996/628 7.35s | m 660->628 |
| NDSDDF | NDS07 | light | converged 0.44453613 5732/534 5.45s | converged 0.44453613 5732/508 4.98s | m 534->508 |
| NDSDDF | NDS08/stability_dual | heavy | certified 1 24400/2336 1.00s | certified 1 24400/2258 0.89s | m 2336->2258 |
| NDSDDF | NDS08/stability_dual | light | certified 1 7408/1840 0.46s | certified 1 7408/1792 0.46s | m 1840->1792 |
| NDSDDF | NDS08 | heavy | converged 0.027334766 24400/2384 20.28s | converged 0.027334766 24400/2306 17.99s | m 2384->2306 |
| NDSDDF | NDS08 | light | converged 0.027307861 7408/1872 9.74s | converged 0.027415479 7408/1824 9.28s | m 1872->1824; value 0.027307861->0.027415479 (3.9e-03 rel) |
| NDSDDF | NDS09/stability_dual | heavy | certified 1 38125/3650 1.74s | unknown 0 38125/3520 2.61s | status certified->unknown; m 3650->3520; value 1->0 (1.0e+00 rel) |
| NDSDDF | NDS09/stability_dual | light | unknown 0 11575/2875 0.86s | unknown 0 11575/2795 0.76s | m 2875->2795 |
| NDSDDF | NDS09 | heavy | converged 0.5994668 83625/4475 52.84s | converged 0.5994668 83625/4305 51.16s | m 4475->4305 |
| NDSDDF | NDS09 | light | converged 0.5994668 37475/3650 24.60s | converged 0.5994668 37475/3540 23.93s | m 3650->3540 |
| NDSDDF | NDS10/stability_dual | heavy | certified 1 6100/584 0.31s | certified 1 6100/571 0.21s | m 584->571 |
| NDSDDF | NDS10/stability_dual | light | certified 1 1852/460 0.26s | certified 1 1852/452 0.17s | m 460->452 |
| NDSDDF | NDS10 | heavy | certified 1 6100/596 0.30s | certified 1 6100/583 0.23s | m 596->583 |
| NDSDDF | NDS10 | light | certified 1 1852/468 0.23s | certified 1 1852/460 0.18s | m 468->460 |
| NDSDDF | NDS11/stability_dual | heavy | unknown 0 1667/171 0.28s | unknown 0 1667/166 0.18s | m 171->166 |
| NDSDDF | NDS11/stability_dual | light | certified 1 545/136 0.22s | unknown 0 545/131 0.18s | status certified->unknown; m 136->131; value 1->0 (1.0e+00 rel) |
| NDSDDF | NDS11 | heavy | unknown 0 1667/178 0.30s | unknown 0 1667/169 0.22s | m 178->169 |
| NDSDDF | NDS11 | light | unknown 0 545/143 0.27s | unknown 0 545/133 0.20s | m 143->133 |
| NDSDDF | NDS12/stability_dual | heavy | certified 1 5828/530 0.32s | certified 1 5828/518 0.20s | m 530->518 |
| NDSDDF | NDS12/stability_dual | light | certified 1 1700/410 0.24s | certified 1 1700/403 0.21s | m 410->403 |
| NDSDDF | NDS12 | heavy | converged 0.44436562 5828/544 4.89s | converged 0.44436562 5828/516 3.96s | m 544->516 |
| NDSDDF | NDS12 | light | converged 0.44436562 1700/422 3.33s | converged 0.44436562 1700/399 2.79s | m 422->399 |
| PDE | 10 | heavy | certified 1 1525/168 0.52s | certified 1 1525/161 0.76s | m 168->161 |
| PDE | 10 | light | unknown 0 1499/165 0.77s | infeasible 0 1499/158 1.53s | status unknown->infeasible; m 165->158 |
| PDE | 11/stability_dual | heavy | certified 1 5564/584 0.30s | certified 1 5564/544 0.27s | m 584->544 |
| PDE | 11/stability_dual | light | certified 1 1552/512 0.23s | certified 1 1552/458 0.21s | m 512->458 |
| PDE | 11 | heavy | converged 2.9296875 5564/584 7.64s | converged 2.9296875 5564/544 9.22s | m 584->544 |
| PDE | 11 | light | converged 2.9296875 3120/544 5.39s | converged 2.9296875 3120/496 6.40s | m 544->496 |
| PDE | 12 | heavy | hi_pass 37.5 12519/1295 0.89s | hi_pass 37.5 12519/1236 1.16s | m 1295->1236 |
| PDE | 12 | light | hi_pass 37.5 3492/1115 0.61s | hi_pass 37.5 3492/1071 0.83s | m 1115->1071 |
| PDE | 13.2/stability_dual | heavy | certified 1 13052/1092 0.88s | certified 1 13052/1075 1.03s | m 1092->1075 |
| PDE | 13.2/stability_dual | light | certified 1 3280/908 0.46s | certified 1 3280/895 0.59s | m 908->895 |
| PDE | 13.2 | heavy | certified 1 13052/1092 1.22s | unknown 0 13052/1075 2.44s | status certified->unknown; m 1092->1075; value 1->0 (1.0e+00 rel) |
| PDE | 13.2 | light | certified 1 3280/908 1.21s | certified 1 3280/895 1.95s | m 908->895 |
| PDE | 13/stability_dual | heavy | certified 1 5564/540 0.33s | certified 1 5564/522 0.33s | m 540->522 |
| PDE | 13/stability_dual | light | certified 1 1552/420 0.23s | certified 1 1552/408 0.20s | m 420->408 |
| PDE | 13 | heavy | certified 1 5564/540 0.41s | certified 1 5564/522 0.54s | m 540->522 |
| PDE | 13 | light | certified 1 1552/420 0.33s | certified 1 1552/408 0.36s | m 420->408 |
| PDE | 14 | heavy | certified 1 22256/1932 0.87s | certified 1 22256/1852 1.01s | m 1932->1852 |
| PDE | 14 | light | certified 1 6208/1476 0.47s | certified 1 6208/1432 0.62s | m 1476->1432 |
| PDE | 15 | heavy | infeasible 0 50480/2483 4.55s | infeasible 0 43892/2351 11.93s | ndv 50480->43892; m 2483->2351 |
| PDE | 15 | light | infeasible 0 21888/2003 2.26s | infeasible 0 19188/1915 8.45s | ndv 21888->19188; m 2003->1915 |
| PDE | 16 | heavy | infeasible 0 122992/6864 11.44s | infeasible 0 117685/6779 19.65s | ndv 122992->117685; m 6864->6779 |
| PDE | 16 | light | infeasible 0 83904/6512 7.03s | infeasible 0 80337/6441 14.21s | ndv 83904->80337; m 6512->6441 |
| PDE | 17 | heavy | certified 1 5564/488 0.30s | certified 1 5564/481 0.61s | m 488->481 |
| PDE | 17 | light | certified 1 1552/368 0.28s | certified 1 1552/365 0.56s | m 368->365 |
| PDE | 18 | heavy | certified 1 5564/530 0.35s | certified 1 4889/513 0.48s | ndv 5564->4889; m 530->513 |
| PDE | 18 | light | certified 1 1552/436 0.29s | certified 1 1417/429 0.43s | ndv 1552->1417; m 436->429 |
| PDE | 19 | heavy | certified 1 12519/1100 0.41s | certified 1 11494/1068 1.24s | ndv 12519->11494; m 1100->1068 |
| PDE | 19 | light | certified 1 3492/848 0.31s | certified 1 3287/836 0.85s | ndv 3492->3287; m 848->836 |
| PDE | 2.2 | heavy | certified 1 5564/500 0.28s | certified 1 5564/488 0.38s | m 500->488 |
| PDE | 2.2 | light | certified 1 1552/404 0.35s | certified 1 1552/398 0.34s | m 404->398 |
| PDE | 2.3 | heavy | certified 1 5564/500 0.31s | certified 1 5564/488 0.61s | m 500->488 |
| PDE | 2.3 | light | certified 1 1552/404 0.46s | certified 1 1552/398 0.60s | m 404->398 |
| PDE | 2.4 | heavy | certified 1 5564/500 0.29s | certified 1 5564/488 0.31s | m 500->488 |
| PDE | 2.4 | light | certified 1 1552/404 0.25s | certified 1 1552/398 0.29s | m 404->398 |
| PDE | 21/Hinf_gain_dual | heavy | certified 0.84187959 5965/528 0.39s | certified 0.84187953 5965/521 0.86s | m 528->521 |
| PDE | 21/Hinf_gain_dual | light | certified 0.84188348 3409/462 0.37s | certified 0.84188347 3409/457 0.73s | m 462->457 |
| PDE | 21 | heavy | certified 0.84187729 5965/526 0.56s | certified 0.84187729 5965/519 1.06s | m 526->519 |
| PDE | 21 | light | certified 0.84187746 3409/460 0.55s | certified 0.84187746 3409/455 0.98s | m 460->455 |
| PDE | 22 | heavy | certified 0.20833385 1596/157 0.39s | certified 0.20833385 1596/155 0.73s | m 157->155 |
| PDE | 22 | light | certified 0.20847421 929/140 0.41s | certified 0.20847416 929/138 0.47s | m 140->138 |
| PDE | 23 | heavy | certified 23.706053 1596/179 0.49s | certified 23.706102 1596/176 0.64s | m 179->176; value 23.706053->23.706102 (2.1e-06 rel) |
| PDE | 23 | light | certified 23.709405 929/170 0.42s | certified 23.718499 929/167 0.51s | m 170->167; value 23.709405->23.718499 (3.8e-04 rel) |
| PDE | 24 | heavy | certified 8.2140553 1596/167 0.45s | certified 8.213836 1596/160 0.39s | m 167->160; value 8.2140553->8.213836 (2.7e-05 rel) |
| PDE | 24 | light | certified 8.2137125 929/156 0.38s | certified 8.2138923 929/146 0.49s | m 156->146; value 8.2137125->8.2138923 (2.2e-05 rel) |
| PDE | 25 | heavy | certified 0.81021203 5965/624 0.54s | certified 0.81021204 5965/584 0.69s | m 624->584 |
| PDE | 25 | light | certified 0.81029224 3409/582 0.49s | certified 0.81029249 3409/534 0.73s | m 582->534 |
| PDE | 26 | heavy | certified 2.1456961 13728/1519 0.79s | certified 2.1456961 13728/1377 0.91s | m 1519->1377 |
| PDE | 26 | light | certified 2.147294 7893/1414 0.73s | certified 2.1472865 7893/1274 0.93s | m 1414->1274; value 2.147294->2.1472865 (3.5e-06 rel) |
| PDE | 28 | heavy | certified 1.0000414 1827/88 0.65s | certified 1.0000281 1827/88 0.37s | value 1.0000414->1.0000281 (1.3e-05 rel) |
| PDE | 29 | heavy | certified 1.000004 1945/107 0.65s | certified 1.000004 1945/106 0.68s | m 107->106 |
| PDE | 29 | light | certified 1.000004 1274/93 0.55s | certified 1.000004 1274/92 0.72s | m 93->92 |
| PDE | 2 | heavy | certified 1 5564/500 0.70s | certified 1 5564/488 0.89s | m 500->488 |
| PDE | 2 | light | certified 1 1552/404 0.62s | certified 1 1552/398 0.81s | m 404->398 |
| PDE | 30 | heavy | inaccurate 0.52324381 2907/135 1.47s | inaccurate 0.52327256 2907/134 2.68s | m 135->134; value 0.52324381->0.52327256 (5.5e-05 rel) |
| PDE | 30 | light | unknown 0.64991793 1394/105 1.09s | unknown 0.64975777 1394/104 1.84s | m 105->104; value 0.64991793->0.64975777 (2.5e-04 rel) |
| PDE | 31 | heavy | certified 0.00024497011 6875/310 0.91s | unknown 0.00024517701 6875/309 1.35s | status certified->unknown; m 310->309; value 0.00024497011->0.00024517701 (8.4e-04 rel) |
| PDE | 31 | light | certified 0.00025719111 4387/260 0.79s | certified 0.00025708184 4387/259 1.21s | m 260->259; value 0.00025719111->0.00025708184 (4.2e-04 rel) |
| PDE | 32 | heavy | unknown 0.54215546 10413/402 1.32s | certified 0.5430436 10413/400 3.64s | status unknown->certified; m 402->400; value 0.54215546->0.5430436 (1.6e-03 rel) |
| PDE | 32 | light | unknown 0.55109001 7249/348 1.41s | unknown 0.55109092 7249/346 2.91s | m 348->346; value 0.55109001->0.55109092 (1.7e-06 rel) |
| PDE | 33 | heavy | certified 1 751849/5202 27.65s | certified 1 751849/5101 32.76s | m 5202->5101 |
| PDE | 33 | light | certified 1 103744/2106 4.10s | certified 1 103744/2041 4.77s | m 2106->2041 |
| PDE | 34 | heavy | certified 1 751849/5287 42.47s | certified 1 751849/5120 57.11s | m 5287->5120 |
| PDE | 34 | light | certified 1 187984/3075 13.96s | certified 1 187984/2968 16.47s | m 3075->2968 |
| PDE | 35 | heavy | certified 1 751849/5297 33.05s | certified 1 751849/5125 44.18s | m 5297->5125 |
| PDE | 35 | light | certified 1 103744/2139 3.91s | certified 1 103744/2052 3.22s | m 2139->2052 |
| PDE | 36 | heavy | skipped NaN NaN/NaN 0.02s | error NaN NaN/20460 13.67s | status skipped->error; m NaN->20460 |
| PDE | 36 | light | infeasible 0 3006976/20888 1652.97s | infeasible 0 751936/11867 221.02s | ndv 3006976->751936; m 20888->11867 |
| PDE | 37 | heavy | certified 1 751849/5284 44.33s | certified 1 751849/5118 49.89s | m 5284->5118 |
| PDE | 37 | light | certified 1 103744/2138 3.73s | certified 1 103744/2043 4.34s | m 2138->2043 |
| PDE | 38 | heavy | skipped NaN NaN/NaN 0.03s | error NaN NaN/20456 7.08s | status skipped->error; m NaN->20456 |
| PDE | 38 | light | infeasible 0 3006976/20888 1638.45s | infeasible 0 751936/11867 216.47s | ndv 3006976->751936; m 20888->11867 |
| PDE | 39 | heavy | skipped NaN NaN/NaN 0.02s | error NaN NaN/46075 26.65s | status skipped->error; m NaN->46075 |
| PDE | 39 | light | certified 1 933696/19143 250.21s | certified 1 933696/18411 315.48s | m 19143->18411 |
| PDE | 3 | heavy | unknown 0 12519/1125 0.66s | certified 1 12519/1065 0.39s | status unknown->certified; m 1125->1065; value 0->1 (1.0e+00 rel) |
| PDE | 3 | light | unknown 0 3492/909 0.38s | certified 1 3492/822 0.37s | status unknown->certified; m 909->822; value 0->1 (1.0e+00 rel) |
| PDE | 4.2 | heavy | certified 1 5564/492 0.29s | certified 1 5564/472 0.19s | m 492->472 |
| PDE | 4.2 | light | certified 1 1552/390 0.22s | certified 1 1552/370 0.23s | m 390->370 |
| PDE | 4.3 | heavy | certified 1 5564/500 0.31s | certified 1 5564/472 0.18s | m 500->472 |
| PDE | 4.3 | light | certified 1 1552/404 0.23s | certified 1 1552/370 0.17s | m 404->370 |
| PDE | 4.4 | heavy | certified 1 5564/492 0.28s | certified 1 5564/472 0.20s | m 492->472 |
| PDE | 4.4 | light | certified 1 1552/390 0.25s | certified 1 1552/370 0.20s | m 390->370 |
| PDE | 4.5 | heavy | infeasible 0 12620/612 0.81s | infeasible 0 12620/580 0.78s | m 612->580 |
| PDE | 4.5 | light | infeasible 0 5472/492 0.50s | infeasible 0 5472/466 0.47s | m 492->466 |
| PDE | 40 | heavy | unknown 0 755625/5308 60.16s | unknown 0 755625/5162 52.97s | m 5308->5162 |
| PDE | 40 | light | unknown 0 105206/2154 5.26s | unknown 0 105206/2079 3.76s | m 2154->2079 |
| PDE | 41 | heavy | unknown 0 12620/667 1.95s | unknown 0 12620/631 1.94s | m 667->631 |
| PDE | 41 | light | unknown 0 5472/547 1.42s | unknown 0 5472/515 1.68s | m 547->515 |
| PDE | 42 | heavy | error NaN NaN/NaN 4.77s | certified 1 1341893/7397 109.14s | status error->certified; ndv NaN->1341893; m NaN->7397; value NaN->1 |
| PDE | 42 | light | error NaN NaN/NaN 3.47s | certified 1 336926/4397 18.87s | status error->certified; ndv NaN->336926; m NaN->4397; value NaN->1 |
| PDE | 43 | heavy | error NaN NaN/NaN 2.09s | error NaN NaN/9386 109.93s | m NaN->9386 |
| PDE | 43 | light | error NaN NaN/NaN 1.64s | infeasible 0 508228/4757 45.30s | status error->infeasible; ndv NaN->508228; m NaN->4757; value NaN->0 |
| PDE | 44 | heavy | lo_fail NaN 51228/2715 7.02s | lo_fail NaN 47851/2582 13.40s | ndv 51228->47851; m 2715->2582 |
| PDE | 44 | light | lo_fail NaN 22412/2231 2.43s | lo_fail NaN 21023/2135 8.36s | ndv 22412->21023; m 2231->2135 |
| PDE | 45 | heavy | unknown 54.832297 1398292/9732 291.32s | unknown 64.243854 1398292/9476 251.15s | m 9732->9476; value 54.832297->64.243854 (1.7e-01 rel) |
| PDE | 45 | light | unknown 2593.8227 350746/5134 59.90s | unknown 3329.3592 350746/5018 37.32s | m 5134->5018; value 2593.8227->3329.3592 (2.8e-01 rel) |
| PDE | 46 | heavy | unknown 8.5889181 1398292/9922 277.64s | unknown 7.1623276 1398292/9509 178.05s | m 9922->9509; value 8.5889181->7.1623276 (1.7e-01 rel) |
| PDE | 46 | light | unknown 962.67199 350746/5222 55.37s | unknown 455.37038 350746/5036 43.43s | m 5222->5036; value 962.67199->455.37038 (5.3e-01 rel) |
| PDE | 47 | heavy | skipped NaN NaN/NaN 0.04s | unknown 348.20998 1398292/9444 271.29s | status skipped->unknown; ndv NaN->1398292; m NaN->9444; value NaN->348.20998 |
| PDE | 47 | light | unknown 22601.94 350746/5062 62.69s | unknown 355.67186 350746/5001 38.50s | m 5062->5001; value 22601.94->355.67186 (9.8e-01 rel) |
| PDE | 48 | heavy | unknown 0.0046891088 767456/7186 79.01s | unknown 0.0045330352 767456/5849 86.11s | m 7186->5849; value 0.0046891088->0.0045330352 (3.3e-02 rel) |
| PDE | 48 | light | unknown 0.005050555 195815/5298 16.58s | unknown 0.0061794314 195815/3764 13.58s | m 5298->3764; value 0.005050555->0.0061794314 (2.2e-01 rel) |
| PDE | 4 | heavy | certified 1 5564/500 0.29s | certified 1 5564/472 0.20s | m 500->472 |
| PDE | 4 | light | certified 1 1552/404 0.22s | certified 1 1552/370 0.18s | m 404->370 |
| PDE | 50/H2_norm_dual | heavy | error NaN NaN/NaN 4.78s | error NaN NaN/31501 9.25s | m NaN->31501 |
| PDE | 50/H2_norm_dual | light | error NaN NaN/NaN 5.26s | error NaN NaN/30166 10.88s | m NaN->30166 |
| PDE | 50 | heavy | error NaN NaN/NaN 4.90s | error NaN NaN/31673 13.59s | m NaN->31673 |
| PDE | 50 | light | error NaN NaN/NaN 4.67s | error NaN NaN/30372 16.01s | m NaN->30372 |
| PDE | 5 | heavy | converged 9.8696044 3155/174 5.28s | converged 9.8696044 3155/167 4.15s | m 174->167 |
| PDE | 5 | light | converged 9.8599661 780/136 3.25s | converged 9.8599661 780/126 2.68s | m 136->126 |
| PDE | 6/stability_dual | heavy | certified 1 1391/140 0.23s | certified 1 1391/138 0.13s | m 140->138 |
| PDE | 6/stability_dual | light | certified 1 388/117 0.20s | certified 1 388/114 0.15s | m 117->114 |
| PDE | 6 | heavy | converged 2.4674011 3155/168 4.93s | converged 2.4674011 3155/166 3.61s | m 168->166 |
| PDE | 6 | light | converged 2.4674011 1368/138 3.60s | converged 2.4674011 1368/136 2.94s | m 138->136 |
| PDE | 7/stability_dual | heavy | certified 1 1391/154 0.28s | certified 1 1391/152 0.28s | m 154->152 |
| PDE | 7/stability_dual | light | certified 1 1368/151 0.37s | certified 1 1368/149 0.29s | m 151->149 |
| PDE | 7 | heavy | converged 4.6508984 3155/180 3.80s | converged 4.6508984 3155/178 3.23s | m 180->178 |
| PDE | 7 | light | converged 4.6508984 1368/155 3.38s | converged 4.6508984 1368/152 2.54s | m 155->152 |
| PDE | 8/stability_dual | heavy | infeasible 0 14172/904 1.28s | infeasible 0 14172/854 1.89s | m 904->854 |
| PDE | 8/stability_dual | light | infeasible 0 6548/768 0.96s | infeasible 0 6548/716 0.99s | m 768->716 |
| PDE | 8 | heavy | lo_fail NaN 14172/904 1.99s | lo_fail NaN 14172/854 3.59s | m 904->854 |
| PDE | 8 | light | lo_fail NaN 6548/758 1.80s | lo_fail NaN 6548/714 2.08s | m 758->714 |
| PDE | 9/stability_dual | heavy | certified 1 1525/163 0.38s | certified 1 1525/156 0.55s | m 163->156 |
| PDE | 9/stability_dual | light | certified 1 463/130 0.33s | certified 1 463/123 0.40s | m 130->123 |
| PDE | 9 | heavy | converged -0.00024414062 1525/168 6.43s | converged -0.00024414062 1525/161 9.60s | m 168->161 |
| PDE | 9 | light | converged -0.00024414062 463/136 5.49s | converged -0.00024414062 463/129 9.11s | m 136->129 |

## Rows new in the 11/17 pass
| library | id | tier | status | result | ndv | m | seconds | note |
|---|---|---|---|---|---|---|---|---|
| PDE | 53 | heavy | error | NaN | NaN | 55274 | 148.72 | sop:PIETOOLS_PIE2PDEstability / lb:size: LB_SIZE ndv=9382545 m=55274 above the gate (ndv 2e+06; m 45000): not solved @ lpisolve:16 |
| PDE | 53 | light | certified | 1 | 907729 | 26123 | 267.91 | sop:PIETOOLS_PIE2PDEstability / PIETOOLS_stability_sop/Q fr=1.02 pinf=0 dinf=0 numerr=0 iter=18 cpu=231.2s m=26123 nx=907729 |
| PDE | 54 | heavy | error | NaN | NaN | 122663 | 299.30 | sop:PIETOOLS_PIE2PDEstability / lb:size: LB_SIZE ndv=29219345 m=122663 above the gate (ndv 2e+06; m 45000): not solved @ lpisolve:16 |
| PDE | 54 | light | error | NaN | NaN | 66295 | 129.76 | sop:PIETOOLS_PIE2PDEstability / lb:size: LB_SIZE ndv=7000529 m=66295 above the gate (ndv 2e+06; m 45000): not solved @ lpisolve:16 |
| PDE | 55 | heavy | error | NaN | NaN | NaN | 1305.48 | sop:PIETOOLS_PIE2PDEstability / MATLAB:nomem: Out of memory. @ ChangeDecVar:113 |
| PDE | 55 | light | error | NaN | NaN | 333851 | 581.02 | sop:PIETOOLS_PIE2PDEstability / lb:size: LB_SIZE ndv=28002116 m=333851 above the gate (ndv 2e+06; m 45000): not solved @ lpisolve:16 |
| PDE | 56 | heavy | error | NaN | NaN | NaN | 0.15 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman' for input argument |
| PDE | 56 | light | error | NaN | NaN | NaN | 0.06 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman' for input argument |
| PDE | 57 | heavy | error | NaN | NaN | NaN | 0.04 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman' for input  |
| PDE | 57 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman' for input  |
| PDE | 58 | heavy | error | NaN | NaN | NaN | 0.03 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Riemann_Fridman' fo |
| PDE | 58 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Riemann_Fridman' fo |
| PDE | 59 | heavy | error | NaN | NaN | NaN | 0.03 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Wave_Delayed_Boundary_Damping_Xu' for input ar |
| PDE | 59 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Wave_Delayed_Boundary_Damping_Xu' for input ar |
| PDE | 60 | heavy | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Heat_Delayed_Robin_Boundary_Nicaise' for input |
| PDE | 60 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Heat_Delayed_Robin_Boundary_Nicaise' for input |
| PDE | 61 | heavy | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Heat_Laplace_Beltrami_Edge' for input argument |
| PDE | 61 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Heat_Laplace_Beltrami_Edge' for input argument |
| PDE | 62 | heavy | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Road_Field_Berestycki' for input arguments of  |
| PDE | 62 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Road_Field_Berestycki' for input arguments of  |
| PDE | 63 | heavy | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff' for input ar |
| PDE | 63 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff' for input ar |
| PDE | 64 | heavy | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Two_Layer_Heat_Rectangles' for input arguments |
| PDE | 64 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Two_Layer_Heat_Rectangles' for input arguments |
| PDE | 65 | heavy | error | NaN | NaN | NaN | 0.02 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair' for input arguments of  |
| PDE | 65 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair' for input arguments of  |
| PDE | 66 | heavy | error | NaN | NaN | NaN | 0.02 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Age_Structured_Diffusion_Kang' for input argum |
| PDE | 66 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Age_Structured_Diffusion_Kang' for input argum |
| PDE | 67 | heavy | error | NaN | NaN | NaN | 0.00 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Elastic_Thermoelastic_Transmission_Marzocchi'  |
| PDE | 67 | light | error | NaN | NaN | NaN | 0.00 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Elastic_Thermoelastic_Transmission_Marzocchi'  |
| PDE | 68 | heavy | error | NaN | NaN | NaN | 0.00 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Wave_Heat_Transmission_1D_Batty' for input arg |
| PDE | 68 | light | error | NaN | NaN | NaN | 0.00 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Wave_Heat_Transmission_1D_Batty' for input arg |
| PDE | 69 | heavy | error | NaN | NaN | NaN | 0.02 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty' for input arg |
| PDE | 69 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty' for input arg |
| PDE | 70 | heavy | error | NaN | NaN | NaN | 0.02 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_2D_Wave_Acoustic_Boundary_Beale' for input arg |
| PDE | 70 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_2D_Wave_Acoustic_Boundary_Beale' for input arg |
| PDE | 71 | heavy | error | NaN | NaN | NaN | 0.02 | sop:PIETOOLS_Hinf_gain / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_2D_Heat_Eq_Exact_Gain_Lidstrom' for input arguments o |
| PDE | 71 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_Hinf_gain / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_2D_Heat_Eq_Exact_Gain_Lidstrom' for input arguments o |
| PDE | 72 | heavy | error | NaN | NaN | NaN | 0.02 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Transport_Edge_Diffusion' for input arguments  |
| PDE | 72 | light | error | NaN | NaN | NaN | 0.01 | sop:PIETOOLS_PIE2PDEstability / MATLAB:UndefinedFunction: Undefined function 'PIETOOLS_PDE_Ex_Transport_Edge_Diffusion' for input arguments  |
