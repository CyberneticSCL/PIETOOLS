| row | form | tier | st (from lpisettings_sop(tier)) | variables start -> final | factor | reason |
|---|---|---|---|---|---|---|
| PDE:1 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 180 -> 9 | 20.00 | converged |
| DDE:DDE03 | P | stripped | `st.lpi.ineq.dD = 0;` | 440 -> 440 | 1.00 | converged |
| DDE:DDE10 | P | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0;` | 1669 -> 1329 | 1.26 | converged |
| PDE:33 | Q | light | `st.lpi.ineq.dD = 0; st.lpi.ineq.psatz{2} = 'none'; st.sop.deg_Q{2} = 0;` | 409 -> 89 | 4.60 | converged |
| DDE:DDE01 | P | stripped | `st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0;` | 1484 -> 1196 | 1.24 | converged |
| DDE:DDE01/stability_dual | Pd | stripped | `st.lpi.pos.len_com{1} = 0;` | 1356 -> 1068 | 1.27 | converged |
| DDE:DDE02 | P | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 3244 -> 2380 | 1.36 | converged |
| DDE:DDE02/stability_dual | Pd | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 3244 -> 1324 | 2.45 | converged |
| DDE:DDE03/stability_dual | Pd | stripped | `st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 280 -> 236 | 1.19 | converged |
| DDE:DDE04 | P | stripped | `st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0;` | 908 -> 748 | 1.21 | converged |
| DDE:DDE04/stability_dual | Pd | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 908 -> 440 | 2.06 | converged |
| DDE:DDE05 | P | stripped | `st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0;` | 227 -> 187 | 1.21 | converged |
| DDE:DDE05/stability_dual | Pd | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.pos.len_com{1} = 1;` | 227 -> 155 | 1.46 | converged |
| DDE:DDE06 | P | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0;` | 811 -> 187 | 4.34 | converged |
| DDE:DDE06/stability_dual | Pd | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 1; st.lpi.ineq.dD = 0;` | 811 -> 187 | 4.34 | converged |
| DDE:DDE07 | P | stripped | `st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0;` | 1356 -> 748 | 1.81 | converged |
| DDE:DDE07/stability_dual | Pd | stripped | `st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 908 -> 748 | 1.21 | converged |
| DDE:DDE08 | P | veryheavy | `st.lpi.pos.len_com{1} = 2; st.lpi.pos.len_sos{1} = [2;3;4]; st.lpi.ineq.dD = 0; st.sop.dw_feas = [0 2 2];` | 6452 -> 2864 | 2.25 | converged |
| DDE:DDE08/stability_dual | Pd | light | `st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 404 -> 352 | 1.15 | converged |
| DDE:DDE09 | P | stripped | `st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0;` | 468 -> 428 | 1.09 | converged |
| DDE:DDE09/Hinf_gain_dual | Pd | light | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 708 -> 651 | 1.09 | converged |
| DDE:DDE10/Hinf_gain_dual | Pd | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.pos.len_com{1} = 0;` | 1669 -> 1285 | 1.30 | converged |
| DDE:DDE11 | P | stripped | `st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 603 -> 559 | 1.08 | converged |
| DDE:DDE11/Hinf_gain_dual | Pd | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.pos.len_com{1} = 3; st.lpi.pos.len_sos{1} = [2;2;1]; st.sop.dw_P = 0;` | 963 -> 732 | 1.32 | converged |
| DDE:DDE12 | P | stripped | `st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 1669 -> 1509 | 1.11 | converged |
| DDE:DDE12/Hinf_gain_dual | Pd | stripped | `st.lpi.ineq.dD = 0;` | 1669 -> 1669 | 1.00 | converged |
| DDE:DDE13 |  | veryheavy | `st.lpi.pos.len_inc{1} = 2; st.lpi.pos.len_com{1} = 3; st.lpi.pos.len_sos{1} = [2;3;4]; st.lpi.ineq.dD = 0; st.lpi.var.deg{1} = [0 6 6];` | 32143 -> 11759 | 2.73 | converged |
| DDE:DDE14 |  | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.pos.len_com{1} = 1; st.lpi.var.deg{1} = [0 2 2]; st.sop.dw_P = 0;` | 1775 -> 811 | 2.19 | converged |
| DDE:DDE15 |  | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.var.deg{1} = [0 2 2]; st.sop.dw_P = 0;` | 1877 -> 1425 | 1.32 | converged |
| DDE:DDE16 |  | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0; st.lpi.var.deg{1} = [0 2 2];` | 27760 -> 23449 | 1.18 | converged |
| DDE:DDE17 | P | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0;` | 5928 -> 1373 | 4.32 | converged |
| DDE:DDE17/Hinf_control |  | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0; st.lpi.var.deg{1} = [0 2 2];` | 10086 -> 7901 | 1.28 | converged |
| DDE:DDE18 |  | veryheavy | `st.lpi.pos.len_com{1} = 2; st.lpi.pos.len_sos{1} = [2;3;4]; st.lpi.var.deg{1} = [2 6 6]; st.sop.dw_P = 1;` | 78328 -> 51580 | 1.52 | converged |
| DDE:DDE20 |  | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0; st.lpi.var.deg{1} = [0 2 2]; st.sop.dw_P = 0;` | 670 -> 217 | 3.09 | converged |
| DDE:DDE21 |  | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0; st.lpi.var.deg{1} = [0 2 2];` | 1775 -> 1543 | 1.15 | converged |
| DDE:DDE22 |  | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0; st.lpi.var.deg{1} = [0 2 2]; st.sop.dw_P = 0;` | 5957 -> 2717 | 2.19 | converged |
| DDE:DDE23 |  | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0; st.lpi.var.deg{1} = [0 2 2];` | 3558 -> 3057 | 1.16 | converged |
| DDE:DDE24 |  | stripped | `st.lpi.ineq.dD = 0; st.lpi.var.deg{1} = [1 2 2];` | 5957 -> 5953 | 1.00 | converged |
| DDE:DDE25 | o | stripped | `st.lpi.pos.len_sos{1} = [0;2;2]; st.lpi.ineq.dD = 0;` | 1746 -> 865 | 2.02 | converged |
| DDE:DDE25/H2_norm_dual | c | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [0;1;2]; st.lpi.ineq.dD = 0;` | 1042 -> 820 | 1.27 | converged |
| DDE:DDE26 | o | stripped | `st.lpi.pos.len_sos{1} = [0;2;1]; st.lpi.ineq.dD = 0;` | 31572 -> 24672 | 1.28 | converged |
| DDE:DDE26/H2_norm_dual | c | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [0;1;2]; st.lpi.ineq.dD = 0;` | 31572 -> 24240 | 1.30 | converged |
| NDSDDF:DDF01 | P | light | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 720 -> 36 | 20.00 | converged |
| NDSDDF:DDF01/stability_dual | Qd | extreme | `st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 544 -> 144 | 3.78 | converged |
| NDSDDF:NDS01 | P | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 908 -> 680 | 1.34 | converged |
| NDSDDF:NDS01/stability_dual | Pd | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 908 -> 440 | 2.06 | converged |
| NDSDDF:NDS02 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 908 -> 88 | 10.32 | converged |
| NDSDDF:NDS02/stability_dual | Pd | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 908 -> 88 | 10.32 | converged |
| NDSDDF:NDS03 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 3244 -> 236 | 13.75 | converged |
| NDSDDF:NDS03/stability_dual | Pd | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 3244 -> 236 | 13.75 | converged |
| NDSDDF:NDS04 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 3244 -> 236 | 13.75 | converged |
| NDSDDF:NDS04/stability_dual | Pd | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 3244 -> 236 | 13.75 | converged |
| NDSDDF:NDS05 | P | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.ineq.dD = 0;` | 3244 -> 2636 | 1.23 | converged |
| NDSDDF:NDS05/stability_dual | Pd | light | `st.lpi.pos.len_com{1} = 1; st.lpi.ineq.dD = 0;` | 3244 -> 2636 | 1.23 | converged |
| NDSDDF:NDS07 | P | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 811 -> 595 | 1.36 | converged |
| NDSDDF:NDS07/stability_dual | Pd | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 811 -> 331 | 2.45 | converged |
| NDSDDF:NDS08 | P | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0;` | 3632 -> 944 | 3.85 | converged |
| NDSDDF:NDS08/stability_dual | Pd | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 3632 -> 1584 | 2.29 | converged |
| NDSDDF:NDS09 | P | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0;` | 5675 -> 2475 | 2.29 | converged |
| NDSDDF:NDS09/stability_dual | Pd | stripped | `st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0;` | 5675 -> 4675 | 1.21 | converged |
| NDSDDF:NDS10 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 908 -> 88 | 10.32 | converged |
| NDSDDF:NDS10/stability_dual | Pd | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 908 -> 88 | 10.32 | converged |
| NDSDDF:NDS11 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [3;2;2]; st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 280 -> 57 | 4.91 | converged |
| NDSDDF:NDS11/stability_dual | Pd | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.ineq.dD = 0;` | 280 -> 68 | 4.12 | converged |
| NDSDDF:NDS12 | P | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 811 -> 595 | 1.36 | converged |
| NDSDDF:NDS12/stability_dual | Pd | light | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 811 -> 331 | 2.45 | converged |
| PDE:1/stability_dual | Qd | extreme | `st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 136 -> 36 | 3.78 | converged |
| PDE:2 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.ineq.dD = 0;` | 720 -> 96 | 7.50 | converged |
| PDE:2.2 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 720 -> 36 | 20.00 | converged |
| PDE:2.3 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.ineq.dD = 0;` | 720 -> 96 | 7.50 | converged |
| PDE:2.4 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 720 -> 36 | 20.00 | converged |
| PDE:3 | P | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;1]; st.lpi.ineq.dD = 0;` | 1620 -> 801 | 2.02 | converged |
| PDE:4 | P | stripped | `st.lpi.pos.len_inc{1} = 1; st.lpi.pos.len_com{1} = 1;` | 720 -> 464 | 1.55 | converged |
| PDE:4.2 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 720 -> 36 | 20.00 | converged |
| PDE:4.3 | P | stripped | `st.lpi.ineq.dD = 0;` | 720 -> 720 | 1.00 | converged |
| PDE:4.4 | P | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0;` | 720 -> 420 | 1.71 | converged |
| PDE:5 | P | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [0;2;2];` | 568 -> 484 | 1.17 | converged |
| PDE:6 | P | light | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [0;2;2]; st.lpi.ineq.dD = 0;` | 500 -> 416 | 1.20 | converged |
| PDE:6/stability_dual | Qd | stripped | `st.lpi.ineq.dD = 0;` | 268 -> 268 | 1.00 | converged |
| PDE:7 | P | light | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [0;2;2]; st.lpi.ineq.dD = 0;` | 500 -> 416 | 1.20 | converged |
| PDE:7/stability_dual | Qd | light | `st.lpi.ineq.cap_inc = 'auto';` | 1368 -> 988 | 1.38 | converged |
| PDE:9 | P | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [2;2;1]; st.lpi.ineq.dD = 0;` | 227 -> 131 | 1.73 | converged |
| PDE:9/stability_dual | Pd | light | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 227 -> 110 | 2.06 | converged |
| PDE:10 | Q | stripped | `st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto'; st.sop.deg_Q{1} = [2 1 1];` | 1303 -> 385 | 3.38 | converged |
| PDE:11 | P | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2];` | 2272 -> 1972 | 1.15 | converged |
| PDE:11/stability_dual | Qd | stripped | `` | 1552 -> 1552 | 1.00 | converged |
| PDE:12 | P | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 2772 -> 2097 | 1.32 | converged |
| PDE:13 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 720 -> 36 | 20.00 | converged |
| PDE:13.2 | Q | stripped | `st.lpi.ineq.dD = 0; st.sop.deg_Q{1} = [0 0 0];` | 4821 -> 3841 | 1.26 | converged |
| PDE:13.2/stability_dual | Qd | stripped | `` | 3280 -> 3280 | 1.00 | converged |
| PDE:13/stability_dual | Pd | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 720 -> 36 | 20.00 | converged |
| PDE:14 | P | stripped | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 2880 -> 144 | 20.00 | converged |
| PDE:17 | P | light | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 720 -> 36 | 20.00 | converged |
| PDE:18 | Pd | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;1]; st.lpi.ineq.dD = 0;` | 720 -> 356 | 2.02 | converged |
| PDE:19 | P | light | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [1;2;2]; st.lpi.ineq.dD = 0;` | 1620 -> 81 | 20.00 | converged |
| PDE:20 | P | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;2;1];` | 697 -> 633 | 1.10 | converged |
| PDE:20/Hinf_gain_dual | Pd | stripped | `st.lpi.ineq.dD = 0;` | 397 -> 397 | 1.00 | converged |
| PDE:21 | P | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 1401 -> 1197 | 1.17 | converged |
| PDE:21/Hinf_gain_dual | Pd | stripped | `st.lpi.pos.len_sos{1} = [2;1;2]; st.lpi.ineq.dD = 0;` | 1401 -> 1257 | 1.11 | converged |
| PDE:22 | P | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0;` | 397 -> 313 | 1.27 | converged |
| PDE:23 | P | stripped | `st.lpi.pos.patterns{1} = 'separable';` | 697 -> 633 | 1.10 | converged |
| PDE:24 | Q | light | `st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 917 -> 833 | 1.10 | converged |
| PDE:25 | Q | stripped | `st.lpi.ineq.dD = 0;` | 2081 -> 2081 | 1.00 | converged |
| PDE:26 | Q | stripped | `st.lpi.ineq.dD = 0;` | 4845 -> 4845 | 1.00 | converged |
| PDE:27 |  | extreme | `st.lpi.pos.len_sos{1} = [1;0;0]; st.lpi.var.deg{1} = [0 2 0]; st.sop.dw_P = 0;` | 304 -> 75 | 4.05 | converged |
| PDE:28 |  | light | `st.lpi.pos.len_sos{1} = [2;1;1]; st.lpi.ineq.dD = 0; st.lpi.var.deg{1} = [1 2 2]; st.sop.dw_P = [0 2 1];` | 624 -> 175 | 3.57 | build budget |
| PDE:29 |  | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [0;2;2]; st.lpi.ineq.cap_inc = 'auto'; st.lpi.var.deg{1} = [0 2 2];` | 770 -> 604 | 1.27 | converged |
| PDE:30 |  | light | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [3;2;1]; st.lpi.ineq.cap_inc = 'auto'; st.lpi.var.deg{1} = [0 2 2];` | 1394 -> 579 | 2.41 | converged |
| PDE:31 |  | veryheavy | `st.lpi.pos.len_inc{1} = 2; st.lpi.pos.len_com{1} = 3; st.lpi.pos.len_sos{1} = [4;4;3]; st.lpi.ineq.dD = 0; st.lpi.var.deg{1} = [1 2 6]; st.sop.dw_P = 0;` | 30683 -> 8369 | 3.67 | converged |
| PDE:32 |  | veryheavy | `st.lpi.pos.len_com{1} = 2; st.lpi.pos.len_sos{1} = [3;4;4]; st.lpi.var.deg{1} = [3 6 6];` | 31537 -> 29227 | 1.08 | converged |
| PDE:49 | o | stripped | `st.lpi.pos.len_sos{1} = [0;2;1];` | 1563 -> 1318 | 1.19 | converged |
| PDE:49/H2_norm_dual | c | veryheavy | `st.lpi.pos.len_inc{1} = 3; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [0;4;3]; st.lpi.ineq.dD = 0;` | 11847 -> 2507 | 4.73 | converged |
| PDE:51 |  | heavy | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [3;2;1]; st.lpi.var.deg{1} = [2 2 2]; st.sop.dw_P = 0; st.sop.deg_Q{1} = [6 2 2];` | 3851 -> 1293 | 2.98 | build budget |
| PDE:52 |  | heavy | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [3;2;1]; st.lpi.var.deg{1} = [1 2 2]; st.sop.dw_P = 1; st.sop.deg_Q{1} = [6 2 2];` | 3986 -> 2471 | 1.61 | build budget |
| PDE:67 | P | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [0;2;1]; st.lpi.ineq.dD = 0;` | 4500 -> 900 | 5.00 | converged |
| PDE:68 | P | stripped | `st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{1} = [0;2;1]; st.lpi.ineq.dD = 0;` | 1620 -> 324 | 5.00 | converged |
| PDE:34 | Q | stripped | `st.sop.deg_Q{2} = [0 1 1];` | 72100 -> 72081 | 1.00 | converged |
| PDE:35 | P | stripped | `st.lpi.pos.len_sos{2} = [0 1;1 1;0 0]; st.lpi.ineq.dD = 0; st.lpi.ineq.psatz{2} = 'none';` | 5409 -> 4112 | 1.32 | build budget |
| PDE:37 | Q | light | `st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto'; st.sop.deg_Q{2} = [1 0 0];` | 16656 -> 11040 | 1.51 | converged |
| PDE:38 | Q | light | `st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 1 1];` | 380880 -> 380804 | 1.00 | converged |
| PDE:40 | Q | light | `st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0;` | 8518 -> 8518 | 1.00 | converged |
| PDE:42 | Q | stripped | `st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 0 1];` | 23331 -> 6197 | 3.76 | converged |
| PDE:43 | Q | light | `st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto'; st.sop.deg_Q{2} = [1 2 2];` | 196513 -> 132609 | 1.48 | converged |
| PDE:45 | P | heavy | `st.lpi.pos.len_sos{2} = [2 2;1 0;0 0]; st.sop.dw_P = 1;` | 799170 -> 91945 | 8.69 | build budget |
| PDE:46 | P | heavy | `st.lpi.pos.len_sos{2} = [2 2;1 1;0 0]; st.lpi.pos.patterns{2} = 'uniform'; st.sop.dw_P = 1;` | 799170 -> 200425 | 3.99 | build budget |
| PDE:47 | P | heavy | `st.lpi.pos.len_sos{2} = [2 2;1 1;0 0]; st.lpi.pos.patterns{2} = 'uniform'; st.lpi.ineq.dD = 0;` | 170690 -> 169458 | 1.01 | build budget |
| PDE:48 |  | heavy | `` | 768159 -> 768159 | 1.00 | build budget |
| PDE:57 | Q | light | `st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [2 1 1];` | 23121 -> 23121 | 1.00 | converged |
| PDE:58 | Q | light | `st.lpi.ineq.dD = 0; st.lpi.ineq.psatz{2} = 'none'; st.lpi.ineq.cap_inc = 'auto'; st.sop.deg_Q{2} = [1 2 2];` | 225705 -> 81604 | 2.77 | converged |
| PDE:59 | Q | light | `st.lpi.ineq.dD = 0; st.lpi.ineq.psatz{2} = 'none'; st.sop.deg_Q{2} = 0; st.sop.max_raise_Q = [2 2 0]; st.sop.raise_Q = 'deg';` | 604 -> 136 | 4.44 | converged |
| PDE:60 | Q | stripped | `st.lpi.ineq.dD = 0; st.lpi.ineq.psatz{2} = 'none'; st.lpi.ineq.cap_inc = 'auto'; st.sop.deg_Q{2} = [2 1 1];` | 1340 -> 552 | 2.43 | build budget |
| PDE:61 | Q | light | `st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto';` | 197415 -> 132999 | 1.48 | converged |
| PDE:62 | P | light | `st.lpi.pos.len_inc{1} = 0; st.lpi.pos.len_com{1} = 0; st.lpi.pos.len_sos{2} = [1 1;1 2;0 0];` | 84205 -> 83269 | 1.01 | build budget |
| PDE:63 | Q | light | `st.lpi.ineq.cap_inc = 'auto';` | 165182 -> 149054 | 1.11 | converged |
| PDE:64 | Q | extreme | `st.sop.deg_Q{2} = [0 1 1];` | 288400 -> 288324 | 1.00 | converged |
| PDE:66 | Q | light | `st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 1 1];` | 16740 -> 16721 | 1.00 | converged |
| PDE:71 | P | light | `st.lpi.pos.len_com{2} = [2 1]; st.lpi.pos.len_sos{2} = [1 1;1 1;0 0]; st.lpi.pos.patterns{2} = 'uniform'; st.lpi.ineq.wcross = 1;` | 203522 -> 94987 | 2.14 | build budget |
| PDE:4.5 | P | veryheavy | `` | NaN -> NaN |  | no start certified |
| PDE:41 | P | heavy | `` | NaN -> NaN |  | no start certified |
| PDE:72 | Q | light | `st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1; st.sop.deg_Q{2} = 2;` | 256948 -> 256948 | 1.00 | converged |
| PDE:36 | Q | light | `` | NaN -> NaN |  | no start certified |
| PDE:50 | oco | light | `st.lpi.ineq.dD = 0;` | 503633 -> 503633 | 1.00 | build budget |
| PDE:50/H2_norm_dual | cco | stripped | `` | NaN -> NaN |  | no start certified |
| PDE:56 | Q | light | `` | NaN -> NaN |  | no start certified |
| PDE:65 | Q | light | `` | NaN -> NaN |  | no start certified |
| PDE:69 | Q | light | `` | NaN -> NaN |  | no start certified |
| PDE:70 | Q | light | `` | NaN -> NaN |  | no start certified |
| PDE:44 | Q | light | `` | NaN -> NaN |  | no start certified |
| PDE:15 | Qd>Q | light | `` | NaN -> NaN |  | no start certified |
| PDE:16 | Q | light | `` | NaN -> NaN |  | no start certified |
| PDE:53 | P | stripped | `` | 668443 -> 580612 | 1.15 | build budget |
| PDE:54 | P | stripped | `` | 1835012 -> 1835012 | 1.00 | build budget |
| PDE:55 | Q | light | `` | NaN -> NaN |  | no start certified |
| PDE:39 | Q | light | `` | NaN -> NaN |  | no start certified |
| PDE:X1 | P | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [0;1;1];` | 280 -> 189 | 1.48 | build budget |
| PDE:X1/stability_dual | Pd | stripped | `st.lpi.pos.len_com{1} = 1; st.lpi.pos.len_sos{1} = [0;1;1];` | 280 -> 189 | 1.48 | build budget |

Rows with a result: 143; reduced: 128; summed variables 8174899 -> 6175524.
