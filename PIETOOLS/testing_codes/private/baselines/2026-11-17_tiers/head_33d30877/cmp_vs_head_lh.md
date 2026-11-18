## Summary per tier (11/17 pass vs 10/28 sopvar reference)
| tier | rows | same | diff | new | certified | infeasible | unknown | error | skipped | seconds |
|---|---|---|---|---|---|---|---|---|---|---|
| heavy | 158 | 54 | 104 | 0 | 128 | 4 | 11 | 9 | 2 | 3071 |
| light | 158 | 50 | 108 | 0 | 117 | 5 | 22 | 10 | 0 | 2519 |

## Rows that differ from the reference
| library | id | tier | reference | 11/17 | difference |
|---|---|---|---|---|---|
| DDE | DDE01 | heavy | converged 1.8587109 6100/567 5.94s | converged 2.0180859 6196/582 8.81s | ndv 6100->6196; m 567->582; value 1.8587109->2.0180859 (8.6e-02 rel) |
| DDE | DDE01 | light | converged 2.0101172 3532/499 6.08s | converged 2.0180859 15580/792 13.38s | ndv 3532->15580; m 499->792; value 2.0101172->2.0180859 (4.0e-03 rel) |
| DDE | DDE02 | heavy | converged 2.0820312 23312/2039 30.68s | converged 2.0839844 23540/2081 39.70s | ndv 23312->23540; m 2039->2081; value 2.0820312->2.0839844 (9.4e-04 rel) |
| DDE | DDE02 | light | converged 2.0839844 6800/1581 14.53s | converged 2.0839844 12820/1797 37.19s | ndv 6800->12820; m 1581->1797 |
| DDE | DDE03 | heavy | converged 1.7161724 1667/169 2.92s | converged 1.7161724 1712/179 4.70s | ndv 1667->1712; m 169->179 |
| DDE | DDE03 | light | converged 1.7144948 1637/164 3.31s | converged 1.6708775 4154/231 8.28s | ndv 1637->4154; m 164->231; value 1.7144948->1.6708775 (2.5e-02 rel) |
| DDE | DDE04 | heavy | converged 6.17258 13380/683 10.81s | converged 6.1665521 6196/582 7.12s | ndv 13380->6196; m 683->582; value 6.17258->6.1665521 (9.8e-04 rel) |
| DDE | DDE04 | light | converged 6.1665521 3532/499 4.51s | converged 6.0580497 3436/508 9.21s | ndv 3532->3436; m 499->508; value 6.1665521->6.0580497 (1.8e-02 rel) |
| DDE | DDE05 | heavy | converged 1.5707 3345/175 4.88s | converged 1.5707 5233/219 9.62s | ndv 3345->5233; m 175->219 |
| DDE | DDE05 | light | converged 1.5707 1499/142 4.69s | converged 1.5660983 859/128 5.81s | ndv 1499->859; m 142->128; value 1.5707->1.5660983 (2.9e-03 rel) |
| DDE | DDE06 | heavy | converged 3 12996/628 10.62s | converged 3 20349/799 14.27s | ndv 12996->20349; m 628->799 |
| DDE | DDE06 | light | converged 3 5732/508 7.32s | converged 3 15077/727 13.61s | ndv 5732->15077; m 508->727 |
| DDE | DDE07 | heavy | converged 1.3713711 13380/675 8.24s | converged 1.3721821 6196/582 11.28s | ndv 13380->6196; m 675->582; value 1.3713711->1.3721821 (5.9e-04 rel) |
| DDE | DDE07 | light | converged 1.3713711 5996/547 5.24s | converged 1.3478511 3436/500 8.17s | ndv 5996->3436; m 547->500; value 1.3713711->1.3478511 (1.7e-02 rel) |
| DDE | DDE08 | heavy | converged 1.42229 3963/262 4.75s | converged 1.42229 6160/336 8.07s | ndv 3963->6160; m 262->336 |
| DDE | DDE08 | light | converged 1.4185425 1237/199 4.29s | converged 1.3435913 1258/218 9.60s | ndv 1237->1258; m 199->218; value 1.4185425->1.3435913 (5.3e-02 rel) |
| DDE | DDE09 | light | certified 25.912728 1040/147 0.56s | certified 25.912755 1040/147 0.81s | value 25.912728->25.912755 (1.0e-06 rel) |
| DDE | DDE10/Hinf_gain_dual | heavy | certified 0.23644372 6517/595 0.37s | certified 0.2364442 6517/595 0.60s | value 0.23644372->0.2364442 (2.0e-06 rel) |
| DDE | DDE11/Hinf_gain_dual | heavy | certified 2.876346 2004/203 0.63s | unknown 2.8763565 2004/203 0.72s | status certified->unknown; value 2.876346->2.8763565 (3.7e-06 rel) |
| DDE | DDE11 | heavy | certified 2.8758939 2004/206 0.57s | certified 2.8758982 2004/206 1.33s | value 2.8758939->2.8758982 (1.5e-06 rel) |
| DDE | DDE11 | light | certified 2.8772595 1246/182 0.36s | certified 2.8772529 1246/182 0.56s | value 2.8772595->2.8772529 (2.3e-06 rel) |
| DDE | DDE12/Hinf_gain_dual | heavy | certified 4.9794268 6517/593 0.44s | certified 4.979455 6517/593 0.50s | value 4.9794268->4.979455 (5.7e-06 rel) |
| DDE | DDE12/Hinf_gain_dual | light | certified 4.9796696 3837/519 0.40s | certified 4.9801891 3837/519 0.44s | value 4.9796696->4.9801891 (1.0e-04 rel) |
| DDE | DDE12 | heavy | certified 4.9794672 6517/599 0.41s | certified 4.9793907 6517/599 0.60s | value 4.9794672->4.9793907 (1.5e-05 rel) |
| DDE | DDE12 | light | certified 4.9835374 3837/525 0.45s | unknown 4.9839982 3837/525 0.52s | status certified->unknown; value 4.9835374->4.9839982 (9.2e-05 rel) |
| DDE | DDE13 | heavy | certified 1.3527283 7623/357 0.81s | certified 1.3530096 7623/357 0.97s | value 1.3527283->1.3530096 (2.1e-04 rel) |
| DDE | DDE13 | light | certified 1.3657495 4983/299 0.70s | certified 1.3653282 4983/299 0.87s | value 1.3657495->1.3653282 (3.1e-04 rel) |
| DDE | DDE17 | heavy | converged 0.14833984 84599/4287 76.19s | converged 0.14833984 132300/5547 94.24s | ndv 84599->132300; m 4287->5547 |
| DDE | DDE17 | light | converged 0.14833984 38149/3492 38.05s | converged 0.14833984 98650/5027 81.34s | ndv 38149->98650; m 3492->5027 |
| DDE | DDE18 | heavy | certified 0.3750075 20179/1090 1.42s | certified 0.37550673 20179/1090 1.23s | value 0.3750075->0.37550673 (1.3e-03 rel) |
| DDE | DDE18 | light | certified 0.37801951 13672/937 0.96s | certified 0.37791699 13672/937 1.37s | value 0.37801951->0.37791699 (2.7e-04 rel) |
| DDE | DDE19 | heavy | unknown 0.00139538 3111/141 0.86s | unknown 0.0013951273 3111/141 0.82s | value 0.00139538->0.0013951273 (1.8e-04 rel) |
| DDE | DDE19 | light | unknown 0.0014027311 2242/124 0.77s | inaccurate 0.0013984255 2242/124 1.15s | status unknown->inaccurate; value 0.0014027311->0.0013984255 (3.1e-03 rel) |
| DDE | DDE25 | heavy | unknown 3.6742593 3060/233 0.43s | certified 3.674256 3060/233 0.43s | status unknown->certified |
| DDE | DDE25 | light | certified 3.6746237 1746/198 0.42s | certified 3.6746856 1746/198 0.62s | value 3.6746237->3.6746856 (1.7e-05 rel) |
| DDE | DDE26/H2_norm_dual | heavy | certified 3.8306126 99540/6754 6.18s | certified 3.8306529 99540/6754 4.90s | value 3.8306126->3.8306529 (1.1e-05 rel) |
| DDE | DDE26/H2_norm_dual | light | certified 3.830674 54996/5698 4.18s | unknown 3.8306734 54996/5698 4.79s | status certified->unknown |
| DDE | DDE26 | heavy | certified 3.8305426 99540/6807 5.69s | certified 3.83056 99540/6807 5.44s | value 3.8305426->3.83056 (4.5e-06 rel) |
| DDE | DDE26 | light | unknown 3.8308533 54996/5751 3.04s | certified 3.8308551 54996/5751 4.56s | status unknown->certified |
| NDSDDF | DDF01 | heavy | certified 1 5564/476 0.21s | certified 1 5584/471 0.23s | ndv 5564->5584; m 476->471 |
| NDSDDF | DDF01 | light | certified 1 1552/362 0.23s | certified 1 2984/401 0.35s | ndv 1552->2984; m 362->401 |
| NDSDDF | DDF03 | heavy | inaccurate 0.00065446149 11145/441 1.06s | inaccurate 0.00065442059 11145/441 0.99s | value 0.00065446149->0.00065442059 (6.2e-05 rel) |
| NDSDDF | DDF03 | light | inaccurate 0.00054836065 7869/381 1.03s | unknown 0.00052024274 7869/381 1.14s | status inaccurate->unknown; value 0.00054836065->0.00052024274 (5.1e-02 rel) |
| NDSDDF | NDS01 | heavy | converged 2.2233267 6100/583 3.97s | converged 2.2233267 6196/590 5.13s | ndv 6100->6196; m 583->590 |
| NDSDDF | NDS01 | light | converged 2.2233267 1852/460 3.38s | converged 2.2233267 3436/516 6.70s | ndv 1852->3436; m 460->516 |
| NDSDDF | NDS02 | heavy | certified 1 6100/567 0.26s | certified 1 6196/582 0.26s | ndv 6100->6196; m 567->582 |
| NDSDDF | NDS02 | light | certified 1 1852/444 0.22s | certified 1 3436/508 0.37s | ndv 1852->3436; m 444->508 |
| NDSDDF | NDS03 | heavy | certified 1 23312/2103 0.68s | certified 1 23540/2113 0.56s | ndv 23312->23540; m 2103->2113 |
| NDSDDF | NDS03 | light | certified 1 6800/1645 0.47s | certified 1 12820/1829 0.63s | ndv 6800->12820; m 1645->1829 |
| NDSDDF | NDS04 | heavy | hi_pass 0.69 23312/2103 1.24s | hi_pass 0.69 23540/2113 1.18s | ndv 23312->23540; m 2103->2113 |
| NDSDDF | NDS04 | light | hi_pass 0.69 6800/1645 0.90s | hi_pass 0.69 12820/1829 1.31s | ndv 6800->12820; m 1645->1829 |
| NDSDDF | NDS05 | heavy | converged 1.5797836 51984/2511 40.45s | converged 1.5784509 23540/2097 34.40s | ndv 51984->23540; m 2511->2097; value 1.5797836->1.5784509 (8.4e-04 rel) |
| NDSDDF | NDS05 | light | converged 1.5797836 22928/2043 21.06s | converged 1.5784509 31924/2361 39.24s | ndv 22928->31924; m 2043->2361; value 1.5797836->1.5784509 (8.4e-04 rel) |
| NDSDDF | NDS07 | heavy | converged 0.44453613 12996/628 7.35s | converged 0.44453613 20349/799 9.22s | ndv 12996->20349; m 628->799 |
| NDSDDF | NDS07 | light | converged 0.44453613 5732/508 4.98s | converged 0.44453613 15077/727 10.96s | ndv 5732->15077; m 508->727 |
| NDSDDF | NDS08 | heavy | converged 0.027334766 24400/2306 17.99s | converged 0.027388574 24784/2348 19.07s | ndv 24400->24784; m 2306->2348; value 0.027334766->0.027388574 (2.0e-03 rel) |
| NDSDDF | NDS08 | light | converged 0.027415479 7408/1824 9.28s | converged 0.027523096 13744/2056 18.69s | ndv 7408->13744; m 1824->2056; value 0.027415479->0.027523096 (3.9e-03 rel) |
| NDSDDF | NDS09/stability_dual | light | unknown 0 11575/2795 0.76s | certified 1 11575/2795 1.57s | status unknown->certified; value 0->1 (1.0e+00 rel) |
| NDSDDF | NDS09 | heavy | converged 0.5994668 83625/4305 51.16s | converged 0.5994668 130825/5425 54.93s | ndv 83625->130825; m 4305->5425 |
| NDSDDF | NDS09 | light | converged 0.5994668 37475/3540 23.93s | converged 0.5994668 97375/4970 60.53s | ndv 37475->97375; m 3540->4970 |
| NDSDDF | NDS10 | heavy | certified 1 6100/583 0.23s | certified 1 6196/590 0.28s | ndv 6100->6196; m 583->590 |
| NDSDDF | NDS10 | light | certified 1 1852/460 0.18s | certified 1 3436/516 0.46s | ndv 1852->3436; m 460->516 |
| NDSDDF | NDS11/stability_dual | heavy | unknown 0 1667/166 0.18s | certified 1 1667/166 0.21s | status unknown->certified; value 0->1 (1.0e+00 rel) |
| NDSDDF | NDS11 | heavy | unknown 0 1667/169 0.22s | certified 1 1712/179 0.27s | status unknown->certified; ndv 1667->1712; m 169->179; value 0->1 (1.0e+00 rel) |
| NDSDDF | NDS11 | light | unknown 0 545/133 0.20s | unknown 0 2278/193 0.92s | ndv 545->2278; m 133->193 |
| NDSDDF | NDS12 | heavy | converged 0.44436562 5828/516 3.96s | converged 0.44436562 5885/523 4.93s | ndv 5828->5885; m 516->523 |
| NDSDDF | NDS12 | light | converged 0.44436562 1700/399 2.79s | converged 0.44436562 3205/451 7.38s | ndv 1700->3205; m 399->451 |
| PDE | 10 | heavy | certified 1 1525/161 0.76s | certified 1 2093/171 0.84s | ndv 1525->2093; m 161->171 |
| PDE | 10 | light | infeasible 0 1499/158 1.53s | certified 1 1403/152 0.61s | status infeasible->certified; ndv 1499->1403; m 158->152; value 0->1 (1.0e+00 rel) |
| PDE | 11 | heavy | converged 2.9296875 5564/544 9.22s | converged 2.9296875 5584/488 8.29s | ndv 5564->5584; m 544->488 |
| PDE | 11 | light | converged 2.9296875 3120/496 6.40s | converged 2.9296875 2984/418 5.65s | ndv 3120->2984; m 496->418 |
| PDE | 12 | heavy | hi_pass 37.5 12519/1236 1.16s | hi_pass 37.5 12564/1086 0.82s | ndv 12519->12564; m 1236->1086 |
| PDE | 12 | light | hi_pass 37.5 3492/1071 0.83s | hi_pass 37.5 12546/1086 0.74s | ndv 3492->12546; m 1071->1086 |
| PDE | 13.2 | heavy | unknown 0 13052/1075 2.44s | certified 1 10037/712 1.25s | status unknown->certified; ndv 13052->10037; m 1075->712; value 0->1 (1.0e+00 rel) |
| PDE | 13.2 | light | certified 1 3280/895 1.95s | certified 1 7437/636 1.28s | ndv 3280->7437; m 895->636 |
| PDE | 13 | heavy | certified 1 5564/522 0.54s | certified 1 5584/485 0.30s | ndv 5564->5584; m 522->485 |
| PDE | 13 | light | certified 1 1552/408 0.36s | certified 1 2984/415 0.30s | ndv 1552->2984; m 408->415 |
| PDE | 14 | heavy | certified 1 22256/1852 1.01s | certified 1 22336/1858 0.69s | ndv 22256->22336; m 1852->1858 |
| PDE | 14 | light | certified 1 6208/1432 0.62s | certified 1 11936/1582 0.54s | ndv 6208->11936; m 1432->1582 |
| PDE | 16 | heavy | infeasible 0 117685/6779 19.65s | infeasible 0 119514/4139 11.96s | ndv 117685->119514; m 6779->4139 |
| PDE | 16 | light | infeasible 0 80337/6441 14.21s | infeasible 0 99882/3847 11.84s | ndv 80337->99882; m 6441->3847 |
| PDE | 17 | heavy | certified 1 5564/481 0.61s | certified 1 5584/507 0.53s | ndv 5564->5584; m 481->507 |
| PDE | 17 | light | certified 1 1552/365 0.56s | certified 1 2984/437 0.51s | ndv 1552->2984; m 365->437 |
| PDE | 19 | heavy | certified 1 11494/1068 1.24s | certified 1 11875/1154 1.10s | ndv 11494->11875; m 1068->1154 |
| PDE | 19 | light | certified 1 3287/836 0.85s | certified 1 6449/1000 1.05s | ndv 3287->6449; m 836->1000 |
| PDE | 1 | heavy | certified 1 1391/120 2.64s | certified 1 1396/118 1.76s | ndv 1391->1396; m 120->118 |
| PDE | 1 | light | certified 1 388/90 4.14s | certified 1 746/100 2.35s | ndv 388->746; m 90->100 |
| PDE | 2.2 | heavy | certified 1 5564/488 0.38s | certified 1 5584/475 0.26s | ndv 5564->5584; m 488->475 |
| PDE | 2.2 | light | certified 1 1552/398 0.34s | certified 1 2984/405 0.35s | ndv 1552->2984; m 398->405 |
| PDE | 2.3 | heavy | certified 1 5564/488 0.61s | certified 1 5584/475 0.24s | ndv 5564->5584; m 488->475 |
| PDE | 2.3 | light | certified 1 1552/398 0.60s | certified 1 2984/405 0.49s | ndv 1552->2984; m 398->405 |
| PDE | 2.4 | heavy | certified 1 5564/488 0.31s | certified 1 5584/475 0.22s | ndv 5564->5584; m 488->475 |
| PDE | 2.4 | light | certified 1 1552/398 0.29s | certified 1 2984/405 0.28s | ndv 1552->2984; m 398->405 |
| PDE | 20/Hinf_gain_dual | light | certified 0.50000223 929/125 0.72s | unknown 0.50000222 929/125 0.65s | status certified->unknown |
| PDE | 20 | heavy | certified 0.50066378 1596/141 1.53s | certified 0.50066131 1596/141 0.91s | value 0.50066378->0.50066131 (4.9e-06 rel) |
| PDE | 22 | light | certified 0.20847416 929/138 0.47s | unknown 0.20847412 929/138 0.43s | status certified->unknown |
| PDE | 23 | light | certified 23.718499 929/167 0.51s | certified 23.719575 929/167 0.40s | value 23.718499->23.719575 (4.5e-05 rel) |
| PDE | 24 | heavy | certified 8.213836 1596/160 0.39s | certified 8.2140658 1596/160 0.25s | value 8.213836->8.2140658 (2.8e-05 rel) |
| PDE | 24 | light | certified 8.2138923 929/146 0.49s | certified 8.2137035 929/146 0.30s | value 8.2138923->8.2137035 (2.3e-05 rel) |
| PDE | 25 | light | certified 0.81029249 3409/534 0.73s | certified 0.81029155 3409/534 0.53s | value 0.81029249->0.81029155 (1.2e-06 rel) |
| PDE | 26 | light | certified 2.1472865 7893/1274 0.93s | certified 2.1472936 7893/1274 0.82s | value 2.1472865->2.1472936 (3.3e-06 rel) |
| PDE | 27 | heavy | certified 1.0000244 1827/85 0.69s | certified 1.000022 1827/85 0.50s | value 1.0000244->1.000022 (2.4e-06 rel) |
| PDE | 27 | light | unknown 1.0000249 1180/71 0.75s | unknown 1.0000223 1180/71 0.63s | value 1.0000249->1.0000223 (2.6e-06 rel) |
| PDE | 28 | heavy | certified 1.0000281 1827/88 0.37s | certified 1.0000415 1827/88 0.29s | value 1.0000281->1.0000415 (1.3e-05 rel) |
| PDE | 28 | light | certified 1.000055 1180/75 0.48s | certified 1.0000533 1180/75 0.32s | value 1.000055->1.0000533 (1.7e-06 rel) |
| PDE | 2 | heavy | certified 1 5564/488 0.89s | certified 1 5584/475 0.58s | ndv 5564->5584; m 488->475 |
| PDE | 2 | light | certified 1 1552/398 0.81s | certified 1 2984/405 0.75s | ndv 1552->2984; m 398->405 |
| PDE | 30 | heavy | inaccurate 0.52327256 2907/134 2.68s | unknown 0.52622977 2907/134 1.68s | status inaccurate->unknown; value 0.52327256->0.52622977 (5.7e-03 rel) |
| PDE | 30 | light | unknown 0.64975777 1394/104 1.84s | certified 0.6499404 1394/104 1.38s | status unknown->certified; value 0.64975777->0.6499404 (2.8e-04 rel) |
| PDE | 31 | heavy | unknown 0.00024517701 6875/309 1.35s | certified 0.00024503739 6875/309 1.02s | status unknown->certified; value 0.00024517701->0.00024503739 (5.7e-04 rel) |
| PDE | 31 | light | certified 0.00025708184 4387/259 1.21s | certified 0.00025642962 4387/259 1.07s | value 0.00025708184->0.00025642962 (2.5e-03 rel) |
| PDE | 32 | heavy | certified 0.5430436 10413/400 3.64s | certified 0.54324228 7341/346 1.35s | ndv 10413->7341; m 400->346; value 0.5430436->0.54324228 (3.7e-04 rel) |
| PDE | 33 | heavy | certified 1 751849/5101 32.76s | certified 1 210361/7123 10.91s | ndv 751849->210361; m 5101->7123 |
| PDE | 33 | light | certified 1 103744/2041 4.77s | certified 1 52660/2552 3.07s | ndv 103744->52660; m 2041->2552 |
| PDE | 34 | heavy | certified 1 751849/5120 57.11s | certified 1 210361/7454 10.40s | ndv 751849->210361; m 5120->7454 |
| PDE | 34 | light | certified 1 187984/2968 16.47s | certified 1 72100/3281 3.40s | ndv 187984->72100; m 2968->3281 |
| PDE | 35 | heavy | certified 1 751849/5125 44.18s | certified 1 210361/7454 10.78s | ndv 751849->210361; m 5125->7454 |
| PDE | 35 | light | certified 1 103744/2052 3.22s | certified 1 72100/3281 2.91s | ndv 103744->72100; m 2052->3281 |
| PDE | 36 | heavy | error NaN NaN/20460 13.67s | skipped NaN NaN/NaN 0.01s | status error->skipped; m 20460->NaN |
| PDE | 36 | light | infeasible 0 751936/11867 221.02s | unknown 0 1775200/25947 501.30s | status infeasible->unknown; ndv 751936->1775200; m 11867->25947 |
| PDE | 37 | heavy | certified 1 751849/5118 49.89s | certified 1 210361/7385 9.25s | ndv 751849->210361; m 5118->7385 |
| PDE | 37 | light | certified 1 103744/2043 4.34s | certified 1 72100/3265 2.86s | ndv 103744->72100; m 2043->3265 |
| PDE | 38 | heavy | error NaN NaN/20456 7.08s | certified 1 744164/30999 260.62s | status error->certified; ndv NaN->744164; m 20456->30999; value NaN->1 |
| PDE | 38 | light | infeasible 0 751936/11867 216.47s | unknown 0 245200/13167 34.33s | status infeasible->unknown; ndv 751936->245200; m 11867->13167 |
| PDE | 39 | heavy | error NaN NaN/46075 26.65s | skipped NaN NaN/NaN 0.01s | status error->skipped; m 46075->NaN |
| PDE | 39 | light | certified 1 933696/18411 315.48s | unknown 0 648900/30721 260.49s | status certified->unknown; ndv 933696->648900; m 18411->30721; value 1->0 (1.0e+00 rel) |
| PDE | 3 | heavy | certified 1 12519/1065 0.39s | certified 1 12564/1053 0.31s | ndv 12519->12564; m 1065->1053 |
| PDE | 3 | light | certified 1 3492/822 0.37s | certified 1 6714/897 0.31s | ndv 3492->6714; m 822->897 |
| PDE | 4.2 | heavy | certified 1 5564/472 0.19s | certified 1 5584/467 0.23s | ndv 5564->5584; m 472->467 |
| PDE | 4.2 | light | certified 1 1552/370 0.23s | certified 1 2984/397 0.20s | ndv 1552->2984; m 370->397 |
| PDE | 4.3 | heavy | certified 1 5564/472 0.18s | certified 1 5584/467 0.22s | ndv 5564->5584; m 472->467 |
| PDE | 4.3 | light | certified 1 1552/370 0.17s | certified 1 2984/397 0.20s | ndv 1552->2984; m 370->397 |
| PDE | 4.4 | heavy | certified 1 5564/472 0.20s | certified 1 5584/467 0.22s | ndv 5564->5584; m 472->467 |
| PDE | 4.4 | light | certified 1 1552/370 0.20s | certified 1 2984/397 0.20s | ndv 1552->2984; m 370->397 |
| PDE | 4.5 | heavy | infeasible 0 12620/580 0.78s | infeasible 0 19776/735 1.03s | ndv 12620->19776; m 580->735 |
| PDE | 4.5 | light | infeasible 0 5472/466 0.47s | infeasible 0 14584/665 0.73s | ndv 5472->14584; m 466->665 |
| PDE | 40 | heavy | unknown 0 755625/5162 52.97s | certified 1 212694/7522 10.50s | status unknown->certified; ndv 755625->212694; m 5162->7522; value 0->1 (1.0e+00 rel) |
| PDE | 40 | light | unknown 0 105206/2079 3.76s | certified 1 73463/3356 4.43s | status unknown->certified; ndv 105206->73463; m 2079->3356; value 0->1 (1.0e+00 rel) |
| PDE | 41 | heavy | unknown 0 12620/631 1.94s | unknown 0 19776/797 2.18s | ndv 12620->19776; m 631->797 |
| PDE | 41 | light | unknown 0 5472/515 1.68s | unknown 0 14584/727 3.09s | ndv 5472->14584; m 515->727 |
| PDE | 42 | heavy | certified 1 1341893/7397 109.14s | certified 1 249298/10329 18.91s | ndv 1341893->249298; m 7397->10329 |
| PDE | 42 | light | certified 1 336926/4397 18.87s | certified 1 69915/4283 13.89s | ndv 336926->69915; m 4397->4283 |
| PDE | 43 | heavy | error NaN NaN/9386 109.93s | certified 1 215188/8686 13.37s | status error->certified; ndv NaN->215188; m 9386->8686; value NaN->1 |
| PDE | 43 | light | infeasible 0 508228/4757 45.30s | certified 1 374827/6576 27.67s | status infeasible->certified; ndv 508228->374827; m 4757->6576; value 0->1 (1.0e+00 rel) |
| PDE | 44 | heavy | lo_fail NaN 47851/2582 13.40s | lo_fail NaN 93828/3648 14.79s | ndv 47851->93828; m 2582->3648 |
| PDE | 44 | light | lo_fail NaN 21023/2135 8.36s | lo_fail NaN 93780/3642 22.30s | ndv 21023->93780; m 2135->3642 |
| PDE | 45 | heavy | unknown 64.243854 1398292/9476 251.15s | unknown 50.415781 1399419/9517 253.84s | ndv 1398292->1399419; m 9476->9517; value 64.243854->50.415781 (2.2e-01 rel) |
| PDE | 45 | light | unknown 3329.3592 350746/5018 37.32s | unknown 4375.5723 351307/5027 41.07s | ndv 350746->351307; m 5018->5027; value 3329.3592->4375.5723 (3.1e-01 rel) |
| PDE | 46 | heavy | unknown 7.1623276 1398292/9509 178.05s | unknown 8.5071045 1399419/9550 265.53s | ndv 1398292->1399419; m 9509->9550; value 7.1623276->8.5071045 (1.9e-01 rel) |
| PDE | 46 | light | unknown 455.37038 350746/5036 43.43s | unknown 812.63584 351307/5045 28.13s | ndv 350746->351307; m 5036->5045; value 455.37038->812.63584 (7.8e-01 rel) |
| PDE | 47 | heavy | unknown 348.20998 1398292/9444 271.29s | unknown 306.34762 1399419/9485 259.62s | ndv 1398292->1399419; m 9444->9485; value 348.20998->306.34762 (1.2e-01 rel) |
| PDE | 47 | light | unknown 355.67186 350746/5001 38.50s | unknown 10931.453 351307/5010 28.80s | ndv 350746->351307; m 5001->5010; value 355.67186->10931.453 (3.0e+01 rel) |
| PDE | 48 | heavy | unknown 0.0045330352 767456/5849 86.11s | certified 0.0039726464 768159/5493 96.98s | status unknown->certified; ndv 767456->768159; m 5849->5493; value 0.0045330352->0.0039726464 (1.2e-01 rel) |
| PDE | 48 | light | unknown 0.0061794314 195815/3764 13.58s | unknown 0.0068670362 195866/3152 10.75s | ndv 195815->195866; m 3764->3152; value 0.0061794314->0.0068670362 (1.1e-01 rel) |
| PDE | 49/H2_norm_dual | heavy | certified 0.10405054 2818/200 0.83s | certified 0.10401646 2818/200 1.06s | value 0.10405054->0.10401646 (3.3e-04 rel) |
| PDE | 49/H2_norm_dual | light | certified 0.10585578 1563/169 0.77s | unknown 0.10585986 1563/169 0.71s | status certified->unknown; value 0.10585578->0.10585986 (3.9e-05 rel) |
| PDE | 4 | heavy | certified 1 5564/472 0.20s | certified 1 5584/467 0.25s | ndv 5564->5584; m 472->467 |
| PDE | 4 | light | certified 1 1552/370 0.18s | certified 1 2984/397 0.25s | ndv 1552->2984; m 370->397 |
| PDE | 50/H2_norm_dual | heavy | error NaN NaN/31501 9.25s | error NaN NaN/56369 20.30s | m 31501->56369 |
| PDE | 50/H2_norm_dual | light | error NaN NaN/30166 10.88s | error NaN NaN/55535 231.35s | m 30166->55535 |
| PDE | 50 | heavy | error NaN NaN/31673 13.59s | error NaN NaN/57377 20.48s | m 31673->57377 |
| PDE | 50 | light | error NaN NaN/30372 16.01s | error NaN NaN/56623 139.63s | m 30372->56623 |
| PDE | 51 | light | certified 1.2699296 2416/137 0.68s | certified 1.2699242 2416/137 0.92s | value 1.2699296->1.2699242 (4.3e-06 rel) |
| PDE | 52 | heavy | certified 1.0900775 3986/176 0.80s | certified 1.0900946 3986/176 0.97s | value 1.0900775->1.0900946 (1.6e-05 rel) |
| PDE | 52 | light | certified 1.0962356 2523/148 0.70s | certified 1.0963826 2523/148 1.40s | value 1.0962356->1.0963826 (1.3e-04 rel) |
| PDE | 53 | heavy | error NaN NaN/55274 148.72s | error NaN NaN/NaN 0.81s | m 55274->NaN |
| PDE | 53 | light | certified 1 907729/26123 267.91s | error NaN NaN/NaN 1.22s | status certified->error; ndv 907729->NaN; m 26123->NaN; value 1->NaN |
| PDE | 54 | heavy | error NaN NaN/122663 299.30s | error NaN NaN/NaN 0.62s | m 122663->NaN |
| PDE | 54 | light | error NaN NaN/66295 129.76s | error NaN NaN/NaN 0.94s | m 66295->NaN |
| PDE | 55 | light | error NaN NaN/333851 581.02s | error NaN NaN/NaN 1.20s | m 333851->NaN |
| PDE | 57 | heavy | error NaN NaN/NaN 0.04s | certified 1 248497/10311 22.37s | status error->certified; ndv NaN->248497; m NaN->10311; value NaN->1 |
| PDE | 57 | light | error NaN NaN/NaN 0.01s | certified 1 433725/7871 39.13s | status error->certified; ndv NaN->433725; m NaN->7871; value NaN->1 |
| PDE | 58 | heavy | error NaN NaN/NaN 0.03s | certified 1 302957/10064 22.07s | status error->certified; ndv NaN->302957; m NaN->10064; value NaN->1 |
| PDE | 58 | light | error NaN NaN/NaN 0.01s | infeasible 0 563813/8126 50.93s | status error->infeasible; ndv NaN->563813; m NaN->8126; value NaN->0 |
| PDE | 59 | heavy | error NaN NaN/NaN 0.03s | certified 1 10581/753 1.96s | status error->certified; ndv NaN->10581; m NaN->753; value NaN->1 |
| PDE | 59 | light | error NaN NaN/NaN 0.01s | certified 1 5214/454 3.18s | status error->certified; ndv NaN->5214; m NaN->454; value NaN->1 |
| PDE | 5 | heavy | converged 9.8696044 3155/167 4.15s | converged 9.8696044 4944/191 5.24s | ndv 3155->4944; m 167->191 |
| PDE | 5 | light | converged 9.8599661 780/126 2.68s | converged 9.8696044 3646/173 4.80s | ndv 780->3646; m 126->173; value 9.8599661->9.8696044 (9.8e-04 rel) |
| PDE | 60 | heavy | error NaN NaN/NaN 0.01s | certified 1 4700/322 1.20s | status error->certified; ndv NaN->4700; m NaN->322; value NaN->1 |
| PDE | 60 | light | error NaN NaN/NaN 0.01s | certified 1 2496/201 1.58s | status error->certified; ndv NaN->2496; m NaN->201; value NaN->1 |
| PDE | 61 | heavy | error NaN NaN/NaN 0.01s | certified 1 841042/19855 204.11s | status error->certified; ndv NaN->841042; m NaN->19855; value NaN->1 |
| PDE | 61 | light | error NaN NaN/NaN 0.01s | error NaN NaN/19151 88.85s | m NaN->19151 |
| PDE | 62 | heavy | error NaN NaN/NaN 0.01s | certified 1 241236/9193 22.87s | status error->certified; ndv NaN->241236; m NaN->9193; value NaN->1 |
| PDE | 62 | light | error NaN NaN/NaN 0.01s | certified 1 86078/4132 6.63s | status error->certified; ndv NaN->86078; m NaN->4132; value NaN->1 |
| PDE | 63 | heavy | error NaN NaN/NaN 0.01s | certified 1 1110293/12861 116.50s | status error->certified; ndv NaN->1110293; m NaN->12861; value NaN->1 |
| PDE | 63 | light | error NaN NaN/NaN 0.01s | certified 1 165182/4910 10.70s | status error->certified; ndv NaN->165182; m NaN->4910; value NaN->1 |
| PDE | 64 | heavy | error NaN NaN/NaN 0.01s | certified 1 841444/32439 198.19s | status error->certified; ndv NaN->841444; m NaN->32439; value NaN->1 |
| PDE | 64 | light | error NaN NaN/NaN 0.01s | certified 1 288400/13915 40.40s | status error->certified; ndv NaN->288400; m NaN->13915; value NaN->1 |
| PDE | 65 | heavy | error NaN NaN/NaN 0.02s | error NaN NaN/148155 35.04s | m NaN->148155 |
| PDE | 65 | light | error NaN NaN/NaN 0.01s | error NaN NaN/67500 21.78s | m NaN->67500 |
| PDE | 66 | heavy | error NaN NaN/NaN 0.02s | certified 1 210361/7284 13.56s | status error->certified; ndv NaN->210361; m NaN->7284; value NaN->1 |
| PDE | 66 | light | error NaN NaN/NaN 0.01s | certified 1 56548/2775 3.95s | status error->certified; ndv NaN->56548; m NaN->2775; value NaN->1 |
| PDE | 67 | heavy | error NaN NaN/NaN 0.00s | certified 1 34900/3339 6.73s | status error->certified; ndv NaN->34900; m NaN->3339; value NaN->1 |
| PDE | 67 | light | error NaN NaN/NaN 0.00s | certified 1 18650/2909 8.86s | status error->certified; ndv NaN->18650; m NaN->2909; value NaN->1 |
| PDE | 68 | heavy | error NaN NaN/NaN 0.00s | certified 1 12564/1208 1.70s | status error->certified; ndv NaN->12564; m NaN->1208; value NaN->1 |
| PDE | 68 | light | error NaN NaN/NaN 0.00s | unknown 0 6714/1052 2.43s | status error->unknown; ndv NaN->6714; m NaN->1052; value NaN->0 |
| PDE | 69 | heavy | error NaN NaN/NaN 0.02s | error NaN NaN/72333 10.76s | m NaN->72333 |
| PDE | 69 | light | error NaN NaN/NaN 0.01s | error NaN NaN/60678 275.15s | m NaN->60678 |
| PDE | 6 | heavy | converged 2.4674011 3155/166 3.61s | converged 2.4674011 4944/190 4.94s | ndv 3155->4944; m 166->190 |
| PDE | 6 | light | converged 2.4674011 1368/136 2.94s | converged 2.4674011 3646/172 4.01s | ndv 1368->3646; m 136->172 |
| PDE | 70 | heavy | error NaN NaN/NaN 0.02s | error NaN NaN/63789 347.34s | m NaN->63789 |
| PDE | 70 | light | error NaN NaN/NaN 0.01s | error NaN NaN/30624 77.68s | m NaN->30624 |
| PDE | 71 | heavy | error NaN NaN/NaN 0.02s | inaccurate 3.1839994 1411226/9853 299.04s | status error->inaccurate; ndv NaN->1411226; m NaN->9853; value NaN->3.1839994 |
| PDE | 71 | light | error NaN NaN/NaN 0.01s | unknown 169.91135 354659/5228 40.63s | status error->unknown; ndv NaN->354659; m NaN->5228; value NaN->169.91135 |
| PDE | 72 | heavy | error NaN NaN/NaN 0.02s | unknown 0 267968/8755 24.48s | status error->unknown; ndv NaN->267968; m NaN->8755; value NaN->0 |
| PDE | 72 | light | error NaN NaN/NaN 0.01s | unknown 0 676410/7989 29.08s | status error->unknown; ndv NaN->676410; m NaN->7989; value NaN->0 |
| PDE | 7 | heavy | converged 4.6508984 3155/178 3.23s | converged 4.6508984 8508/226 4.76s | ndv 3155->8508; m 178->226 |
| PDE | 7 | light | converged 4.6508984 1368/152 2.54s | converged 4.6508984 6562/208 3.34s | ndv 1368->6562; m 152->208 |
| PDE | 8 | heavy | lo_fail NaN 14172/854 3.59s | lo_fail NaN 26480/1118 3.51s | ndv 14172->26480; m 854->1118 |
| PDE | 8 | light | lo_fail NaN 6548/714 2.08s | lo_fail NaN 16616/936 2.42s | ndv 6548->16616; m 714->936 |
| PDE | 9 | heavy | converged -0.00024414062 1525/161 9.60s | converged -0.00024414062 2093/171 11.13s | ndv 1525->2093; m 161->171 |
| PDE | 9 | light | converged -0.00024414062 463/129 9.11s | converged -0.00024414062 1403/152 8.26s | ndv 463->1403; m 129->152 |

## Rows new in the 11/17 pass
| library | id | tier | status | result | ndv | m | seconds | note |
|---|---|---|---|---|---|---|---|---|
