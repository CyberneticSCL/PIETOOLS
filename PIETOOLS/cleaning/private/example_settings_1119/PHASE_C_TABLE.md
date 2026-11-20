Corrected 11/19 18:05 (maintainer: "Apply all four"; PHASE_C.md §4, §7, §10): PDE 15 certifies (Qd, 74,593 variables), PDE 4.5 certifies (P, 720), DDE19 1.963870 at 221, DDF03 0.420543 at 4,673; the rows below predate it.

DDE19 and DDF03 certify since 11/19 17:40 with their regulated outputs (PHASE_C.md §7: 1.963986 at 265 variables; 0.439636 at 4,673); their rows below predate it.

Library-path sizes of 1e6 or more were written to 6 significant digits (libcheck before 05:55): PDE 54 is 1,835,012, the PDE 50 dual about 1.66e6.

| Example | Case | Form | Tier and changes | Measured (harness) | Variables | Light as is | Library path |
|---|---|---|---|---|---|---|---|
| PIETOOLS_PDDE_Ex | 50 | oco | light: lpi.ineq.dD = 0 |  |  |  | H2_norm optimal 0.877971 at 503,633, H2_norm_dual infeasible 3.97632 at 1,660,060 |
| DDE library | DDE01 | stab P | light: lpi.ineq.dD = 0 | certified up to tau(1) = 2.04; certified (stability_dual) | 2,316; 2,316 (stability_dual) | 6,196; 5,996 (stability_dual) | stability certified, stability_dual certified |
| DDE library | DDE02 | stab P | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;2]; lpi.ineq.dD = 0 | certified up to tauN = 2.08398; certified (stability_dual) | 2,380; 2,380 (stability_dual) | 23,540; 6,800 (stability_dual) | stability certified, stability_dual certified |
| DDE library | DDE03 | stab P | stripped: lpi.ineq.dD = 0 | certified up to tau(1) = 1.71617; certified (stability_dual) | 440; 280 (stability_dual) | 1,712; 545 (stability_dual) | stability certified, stability_dual certified |
| DDE library | DDE04 | stab P | stripped: lpi.pos.len_sos{1} = [2;2;1]; lpi.ineq.dD = 0 | certified up to tau(1) = 6.16655; certified (stability_dual) | 1,324; 748 (stability_dual) | 6,196; 1,852 (stability_dual) | stability certified, stability_dual certified |
| DDE library | DDE05 | stab P | stripped: lpi.pos.len_sos{1} = [2;2;1]; lpi.ineq.dD = 0 | certified up to tau(1) = 1.5707; certified (stability_dual) | 331; 187 (stability_dual) | 5,233; 463 (stability_dual) | stability certified, stability_dual certified |
| DDE library | DDE06 | stab P | stripped: lpi.pos.len_inc{1} = 1; lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;1]; lpi.ineq.dD = 0 | certified up to b (Ai{1}) = 3; certified (stability_dual) | 331; 187 (stability_dual) | 20,349; 1,700 (stability_dual) | stability certified, stability_dual certified |
| DDE library | DDE07 | stab P | stripped: lpi.pos.len_sos{1} = [2;2;1]; lpi.ineq.dD = 0 | certified up to tauN = 1.37218; certified (stability_dual) | 1,324; 748 (stability_dual) | 6,196; 1,852 (stability_dual) | stability certified, stability_dual certified |
| DDE library | DDE08 | stab P | veryheavy: lpi.pos.len_com{1} = 2; lpi.pos.len_sos{1} = [2;3;4]; lpi.ineq.dD = 0; sop.dw_feas = [0 2 2] | certified up to tau(1) = 3.93315; certified (stability_dual) | 7,088; 2,864 (stability_dual) | 6,160; 733 (stability_dual) | stability certified, stability_dual certified |
| DDE library | DDE09 | gain P | light: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;2]; lpi.ineq.dD = 0 | gamma 25.9228; gamma 25.9169 (Hinf_gain_dual) | 651; 651 (Hinf_gain_dual) | 1,040; 1,040 (Hinf_gain_dual) | Hinf_gain certified, Hinf_gain_dual certified |
| DDE library | DDE10 | gain P | stripped: lpi.ineq.dD = 0 | gamma 0.236443; gamma 0.236539 (Hinf_gain_dual) | 1,669; 1,669 (Hinf_gain_dual) | 3,837; 3,837 (Hinf_gain_dual) | Hinf_gain certified, Hinf_gain_dual certified |
| DDE library | DDE11 | gain P | stripped: lpi.pos.len_inc{1} = 1; lpi.pos.len_com{1} = 3; lpi.pos.len_sos{1} = [2;2;1]; sop.dw_P = 0 | gamma 2.86672; gamma 2.86706 (Hinf_gain_dual) | 732; 732 (Hinf_gain_dual) | 1,246; 1,246 (Hinf_gain_dual) | Hinf_gain certified, Hinf_gain_dual certified |
| DDE library | DDE12 | gain P | stripped: lpi.ineq.dD = 0 | gamma 4.97845; gamma 4.9785 (Hinf_gain_dual) | 1,669; 1,669 (Hinf_gain_dual) | 3,837; 3,837 (Hinf_gain_dual) | Hinf_gain certified, Hinf_gain_dual certified |
| DDE library | DDE13 |  | veryheavy: lpi.pos.len_inc{1} = 2; lpi.pos.len_com{1} = 3; lpi.pos.len_sos{1} = [2;3;4]; lpi.ineq.dD = 0; lpi.var.deg{1} = [0 6 6] | closed-loop gain 1.34225 | 11,759 | 4,983 | Hinf_control certified |
| DDE library | DDE14 |  | stripped: lpi.pos.len_inc{1} = 1; lpi.pos.len_com{1} = 1; lpi.var.deg{1} = [0 2 2]; sop.dw_P = 0 | closed-loop gain 0.395305 | 811 | 4,983 | Hinf_control certified |
| DDE library | DDE15 |  | stripped: lpi.pos.len_com{1} = 0; lpi.var.deg{1} = [0 2 2]; sop.dw_P = 0 | closed-loop gain 0.610589 | 1,425 | 5,173 | Hinf_control certified |
| DDE library | DDE16 |  | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;2]; lpi.ineq.dD = 0; lpi.var.deg{1} = [0 2 2] | closed-loop gain 0.963509 | 23,449 | 86,188 | Hinf_control certified |
| DDE library | DDE17 | stab P | stripped: lpi.ineq.dD = 0 | certified up to h = 0.14834; closed-loop gain 9.35134 (Hinf_control) | 9,608; 10,086 (Hinf_control) | 132,300; 29,146 (Hinf_control) | Hinf_control certified, stability certified |
| DDE library | DDE18 |  | veryheavy: lpi.pos.len_com{1} = 2; lpi.pos.len_sos{1} = [2;3;4]; lpi.var.deg{1} = [2 6 6]; sop.dw_P = 1 | closed-loop gain 0.333523 | 51,580 | 13,672 | Hinf_control certified |
| DDE library | DDE19 |  | – | no certified setting |  |  |  |
| DDE library | DDE20 |  | stripped: lpi.pos.len_inc{1} = 1; lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;2;1]; lpi.ineq.dD = 0; lpi.var.deg{1} = [0 2 2]; sop.dw_P = 0 | closed-loop gain 0.235702 | 217 | 1,626 | Hinf_estimator certified |
| DDE library | DDE21 |  | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;2]; lpi.ineq.dD = 0; lpi.var.deg{1} = [0 2 2] | closed-loop gain 2.32816 | 1,543 | 4,983 | Hinf_estimator certified |
| DDE library | DDE22 |  | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;2]; lpi.ineq.dD = 0; lpi.var.deg{1} = [0 2 2]; sop.dw_P = 0 | closed-loop gain 1.35012 | 2,717 | 17,909 | Hinf_estimator certified |
| DDE library | DDE23 |  | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;2]; lpi.ineq.dD = 0; lpi.var.deg{1} = [0 2 2] | closed-loop gain 1.8081 | 3,057 | 10,446 | Hinf_estimator certified |
| DDE library | DDE24 |  | stripped: lpi.ineq.dD = 0; lpi.var.deg{1} = [1 2 2] | closed-loop gain 1.04437 | 5,953 | 17,909 | Hinf_estimator certified |
| DDE library | DDE25 |  | stripped: lpi.pos.len_sos{1} = [0;2;2]; lpi.ineq.dD = 0 | H2 norm 3.67514; H2 norm 3.6724 (H2_norm_dual) | 865; 865 (H2_norm_dual) | 1,746; 1,746 (H2_norm_dual) | H2_norm certified, H2_norm_dual certified |
| DDE library | DDE26 |  | stripped: lpi.pos.len_sos{1} = [0;2;1]; lpi.ineq.dD = 0 | H2 norm 3.83098; H2 norm 3.83136 (H2_norm_dual) | 24,672; 24,672 (H2_norm_dual) | 54,996; 54,996 (H2_norm_dual) | H2_norm certified, H2_norm_dual certified |
| NDSDDF library | DDF01 | stab P | extreme: lpi.ineq.dD = 0; lpi.ineq.cap_inc = 'auto' | certified; certified (stability_dual) | 48; 48 (stability_dual) | 5,584; 1,552 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | DDF03 |  | – | no certified setting |  |  |  |
| NDSDDF library | NDS01 | stab P | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;2]; lpi.ineq.dD = 0 | certified up to tau(1) = 2.22333; certified (stability_dual) | 680; 680 (stability_dual) | 6,196; 1,852 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | NDS02 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified; certified (stability_dual) | 88; 88 (stability_dual) | 6,196; 1,852 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | NDS03 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified; certified (stability_dual) | 236; 236 (stability_dual) | 23,540; 6,800 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | NDS04 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified up to tau(2) = 0.69; certified (stability_dual) | 236; 236 (stability_dual) | 23,540; 6,800 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | NDS05 | stab P | stripped: lpi.pos.len_inc{1} = 1; lpi.ineq.dD = 0 | certified up to h = 1.57978; certified (stability_dual) | 4,812; 2,636 (stability_dual) | 23,540; 6,800 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | NDS07 | stab P | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;2]; lpi.ineq.dD = 0 | certified up to h = 0.444536; certified (stability_dual) | 1,139; 595 (stability_dual) | 20,349; 1,700 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | NDS08 | stab P | light: lpi.ineq.dD = 0 | certified up to tau(1) = 0.0275231; certified (stability_dual) | 3,632; 3,632 (stability_dual) | 24,784; 7,408 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | NDS09 | stab P | stripped: lpi.pos.len_sos{1} = [2;2;1]; lpi.ineq.dD = 0 | certified up to tau(1) = 0.599467; certified (stability_dual) | 8,275; 4,675 (stability_dual) | 130,825; 11,575 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | NDS10 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified; certified (stability_dual) | 88; 88 (stability_dual) | 6,196; 1,852 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | NDS11 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.ineq.dD = 0 | certified; certified (stability_dual) | 68; 68 (stability_dual) | 1,712; 545 (stability_dual) | stability certified, stability_dual certified |
| NDSDDF library | NDS12 | stab P | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;2]; lpi.ineq.dD = 0 | certified up to h = 0.444366; certified (stability_dual) | 1,139; 595 (stability_dual) | 5,885; 1,700 (stability_dual) | stability certified, stability_dual certified |
| 2D_Delayed_Heat_Eq_Fridman | 56 | stab Q | – | no certified setting |  |  |  |
| 2D_Heat_Eq_Exact_Gain_Lidstrom | 71 | gain P | light: lpi.pos.len_com{2} = [2 1]; lpi.pos.len_sos{2} = [1 1;1 1;0 0]; lpi.pos.patterns{2} = 'uniform'; lpi.ineq.wcross = 1 | gamma 0.0459442 | 202,922 | 354,659 | Hinf_gain optimal 0.0459804 at 94,987 |
| 2D_Heat_Eq_with_ODE | 40 | stab Q | light: lpi.ineq.dD = 0; sop.deg_Q{2} = 0 | certified | 8,518 | 73,463 | stability optimal at 8,518 |
| 2D_KISS_Model | 34 | stab Q | stripped: sop.deg_Q{2} = [0 1 1] | certified | 72,081 | 72,100 | stability optimal at 72,081 |
| 2D_Linear_Turing_Pair | 65 | stab Q | – | no certified setting |  |  |  |
| 2D_NS_Antonelli | 39 | stab Q | – | no certified setting |  |  |  |
| 2D_Parabolic_Eq | 37 | stab Q | light: lpi.ineq.dD = 0; lpi.ineq.cap_inc = 'auto'; sop.deg_Q{2} = [1 0 0] | certified | 11,040 | 72,100 | stability optimal at 11,040 |
| 2D_ReactionDiffusion_DDDD | 46 | gain P | heavy: lpi.pos.len_sos{2} = [2 2;1 1;0 0]; lpi.pos.patterns{2} = 'uniform'; sop.dw_P = 1 | gamma 0.171259 | 200,425 | 351,307 | Hinf_gain optimal 0.171259 at 200,425 |
| 2D_ReactionDiffusion_DDDN | 45 | gain P | heavy: lpi.pos.len_sos{2} = [2 2;1 0;0 0]; sop.dw_P = 1 | gamma 0.111119 | 91,945 | 351,307 | Hinf_gain optimal 0.111119 at 91,945 |
| 2D_ReactionDiffusion_DNDN | 47 | gain P | heavy: lpi.pos.len_sos{2} = [2 2;1 1;0 0]; lpi.pos.patterns{2} = 'uniform'; lpi.ineq.dD = 0 | gamma 0.206727 | 169,458 | 351,307 | Hinf_gain optimal 0.206727 at 169,458 |
| 2D_ReactionDiffusion_DNDN_EdgeObserve | 48 |  | heavy | closed-loop gain 0.00397265 | 768,159 | 195,866 | Hinf_estimator optimal 0.00397265 at 768,159 |
| 2D_Reaction_Diffusion_Eq | 35 | stab P | stripped: lpi.pos.len_sos{2} = [0 1;1 1;0 0]; lpi.ineq.dD = 0; lpi.ineq.psatz{2} = 'none' | certified | 4,112 | 72,100 | stability optimal at 4,112 |
| 2D_Telegraph_Eq | 36 | stab Q | – | no certified setting |  |  |  |
| 2D_Transport_Eq | 33 | stab Q | light: lpi.ineq.dD = 0; lpi.ineq.psatz{2} = 'none'; sop.deg_Q{2} = 0 | certified | 89 | 52,660 | stability optimal at 89 |
| 2D_Wave_Acoustic_Boundary_Beale | 70 | stab Q | – | no certified setting |  |  |  |
| 2D_Wave_Eq | 38 | stab Q | light: lpi.ineq.dD = 0; sop.deg_Q{2} = [0 1 1] | certified | 380,804 | 245,200 | stability optimal at 380,804 |
| 3D_Advection_Diffusion_Eq | 54 | stab P | stripped |  |  |  | stability optimal at 1,835,010 |
| 3D_Coupled_Heat_Eq | 55 | stab Q | – | no certified setting |  |  |  |
| 3D_Reaction_Diffusion_Eq | 53 | stab P | stripped |  |  |  | stability optimal at 580,612 |
| Age_Structured_Diffusion_Kang | 66 | stab Q | light: lpi.ineq.dD = 0; sop.deg_Q{2} = [0 1 1] | certified | 16,721 | 56,548 | stability optimal at 16,721 |
| Coupled_ODE_Heat_Eq_w_Delay | 42 | stab Q | stripped: lpi.ineq.dD = 0; sop.deg_Q{2} = [0 0 1] | certified | 6,197 | 69,915 | stability optimal at 6,197 |
| Damped_Wave_Interior_Delay_Fridman | 57 | stab Q | light: lpi.ineq.dD = 0; sop.deg_Q{2} = [2 1 1] | certified | 23,121 | 433,725 | stability optimal at 23,121 |
| Damped_Wave_Interior_Delay_Riemann_Fridman | 58 | stab Q | light: lpi.ineq.dD = 0; lpi.ineq.psatz{2} = 'none'; lpi.ineq.cap_inc = 'auto'; sop.deg_Q{2} = [1 2 2] | certified | 81,604 | 563,813 | stability optimal at 81,604 |
| Diagne | 3 | stab P | stripped: lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;1]; lpi.ineq.dD = 0 | certified | 801 | 12,564 | stability optimal at 801 |
| Diffusive_Eq_with_Distributed_Disturbance | 25 | gain Q | stripped: lpi.ineq.dD = 0 | gamma 0.81037 | 2,081 | 3,409 | Hinf_gain optimal 0.81037 at 2,081 |
| Diffusive_Eq_with_Distributed_Disturbance_2 | 26 | gain Q | stripped: lpi.ineq.dD = 0 | gamma 2.14742 | 4,845 | 7,893 | Hinf_gain optimal 2.14742 at 4,845 |
| Elastic_Thermoelastic_Transmission_Marzocchi | 67 | stab P | stripped: lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [0;2;1]; lpi.ineq.dD = 0 | certified | 900 | 34,900 | stability optimal at 900 |
| Euler_Bernoulli_Beam_Eq | 13.0 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified; certified (stability_dual) | 36; 36 (stability_dual) | 5,584; 1,552 (stability_dual) | stability optimal at 36, stability_dual optimal at 36 |
| Euler_Bernoulli_Beam_Eq | 13.2 | stab Q | stripped: lpi.ineq.dD = 0; sop.deg_Q{1} = [0 0 0] | certified; certified (stability_dual) | 3,841; 5,072 (stability_dual) | 10,037; 3,280 (stability_dual) | stability optimal at 3,841, stability_dual optimal at 5,072 |
| Euler_Bernoulli_Beam_with_Control | 31 |  | veryheavy: lpi.pos.len_inc{1} = 2; lpi.pos.len_com{1} = 3; lpi.pos.len_sos{1} = [4;4;3]; lpi.ineq.dD = 0; lpi.var.deg{1} = [1 2 6]; sop.dw_P = 0 | closed-loop gain 0.000169037 | 8,369 | 4,387 | Hinf_control optimal 0.000169037 at 8,369 |
| Heat_Delayed_Robin_Boundary_Nicaise | 60 | stab Q | stripped: lpi.ineq.dD = 0; lpi.ineq.psatz{2} = 'none'; lpi.ineq.cap_inc = 'auto'; sop.deg_Q{2} = [2 1 1] | certified | 552 | 2,496 | stability optimal at 552 |
| Heat_Edge_Cascade_Variable_Coeff | 63 | stab Q | light: lpi.ineq.cap_inc = 'auto' | certified | 149,054 | 165,182 | stability optimal at 149,054 |
| Heat_Eq_w_Delayed_Boundary_Input | 41 | stab P | – | no certified setting |  |  |  |
| Heat_Eq_w_Interior_Delay | 43 | stab Q | light: lpi.ineq.dD = 0; lpi.ineq.cap_inc = 'auto'; sop.deg_Q{2} = [1 2 2] | certified | 132,609 | 374,827 | stability optimal at 132,609 |
| Heat_Eq_with_Distributed_Disturbance | 22 | gain P | stripped: lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [2;1;1]; lpi.ineq.dD = 0 | gamma 0.208434 | 313 | 929 | Hinf_gain optimal 0.208434 at 313 |
| Heat_Eq_with_ODE_Das | 8 |  | – | no certified setting |  |  |  |
| Heat_Eq_with_ODE_Tang | 10 | stab Q | stripped: lpi.ineq.dD = 0; lpi.ineq.cap_inc = 'auto'; sop.deg_Q{1} = [2 1 1] | certified | 385 | 2,093 | stability optimal at 385 |
| Heat_Eq_with_ODE_in_BC | 9 | stab P | light: lpi.ineq.dD = 0 | certified up to k = -0.000244141; certified (stability_dual) | 227; 227 (stability_dual) | 2,093; 463 (stability_dual) | stability optimal at 227, stability_dual optimal at 227 |
| Heat_Eq_with_Observed_Output | 28 |  | light: lpi.pos.len_sos{1} = [2;1;1]; lpi.ineq.dD = 0; lpi.var.deg{1} = [1 2 2]; sop.dw_P = [0 2 1] | closed-loop gain 1.00041 | 175 | 1,180 | Hinf_estimator optimal 1.00041 at 175 |
| Heat_Laplace_Beltrami_Edge | 61 | stab Q | light: lpi.ineq.dD = 0; lpi.ineq.cap_inc = 'auto' | certified | 132,999 |  | stability optimal at 132,999 |
| Lamare | 2.0 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.ineq.dD = 0 | certified | 96 | 5,584 | stability optimal at 96 |
| Lamare | 2.2 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified | 36 | 5,584 | stability optimal at 36 |
| Lamare | 2.3 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.ineq.dD = 0 | certified | 96 | 5,584 | stability optimal at 96 |
| Lamare | 2.4 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified | 36 | 5,584 | stability optimal at 36 |
| Parabolic_Eq_Gahlawat | 7 | stab Q | light: lpi.ineq.cap_inc = 'auto' | certified up to lam = 4.6509; certified (stability_dual) | 8,136; 988 (stability_dual) | 8,508; 1,368 (stability_dual) | stability optimal at 3,664, stability_dual optimal at 988 |
| Parabolic_Eq_with_Disturbance | 23 | gain P | stripped: lpi.pos.patterns{1} = 'separable' | gamma 23.7285 | 633 | 929 | Hinf_gain optimal 23.7285 at 633 |
| Reaction_Diffusion_Ahmadi | 11 | stab P | stripped: lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2] | certified up to R = 2.92969; certified (stability_dual) | 1,972; 1,972 (stability_dual) | 5,584; 1,552 (stability_dual) | stability optimal at 1,972, stability_dual optimal at 1,972 |
| Reaction_Diffusion_Ahmadi_2 | 12 | stab P | light: lpi.ineq.dD = 0 | certified up to R = 37.5 | 2,772 | 12,564 | stability optimal at 2,772 |
| Reaction_Diffusion_Boundary_ODE_Control | 30 |  | light: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [3;2;1]; lpi.ineq.cap_inc = 'auto'; lpi.var.deg{1} = [0 2 2] | closed-loop gain 0.626125 | 579 | 1,394 | Hinf_control optimal 0.626125 at 579 |
| Reaction_Diffusion_Eq | 5 | stab P | stripped: lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [0;2;2] | certified up to lam = 9.85997 | 484 | 4,944 | stability optimal at 484 |
| Reaction_Diffusion_Eq_with_Controlled_Input | 29 |  | stripped: lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [0;2;2]; lpi.ineq.cap_inc = 'auto'; lpi.var.deg{1} = [0 2 2] | closed-loop gain 1.00002 | 604 | 1,274 | Hinf_control optimal 1.00002 at 604 |
| Reaction_Diffusion_Eq_with_Disturbance | 24 | gain Q | light: lpi.ineq.dD = 0; lpi.ineq.cap_inc = 'auto' | gamma 8.2156 | 833 | 929 | Hinf_gain optimal 8.2156 at 833 |
| Reaction_Diffusion_Gahlawat | 6 | stab P | light: lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [0;2;2]; lpi.ineq.dD = 0 | certified up to lam = 2.46499; certified (stability_dual) | 416; 416 (stability_dual) | 4,944; 388 (stability_dual) | stability optimal at 416, stability_dual optimal at 416 |
| Reaction_Diffusion_Robin | X | stab P | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [0;1;1] | certified up to lam = 1.70543; certified (stability_dual) | 189; 189 (stability_dual) |  | stability optimal at 189, stability_dual optimal at 189 |
| Road_Field_Berestycki | 62 | stab P | light: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{2} = [1 1;1 2;0 0] | certified | 83,269 | 86,078 | stability optimal at 83,269 |
| Saba | 4.0 | stab P | stripped: lpi.pos.len_inc{1} = 1; lpi.pos.len_com{1} = 1 | certified | 464 | 5,584 | stability optimal at 464 |
| Saba | 4.2 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified | 36 | 5,584 | stability optimal at 36 |
| Saba | 4.3 | stab P | stripped: lpi.ineq.dD = 0 | certified | 720 | 5,584 | stability optimal at 720 |
| Saba | 4.4 | stab P | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;1;1]; lpi.ineq.dD = 0 | certified | 420 | 5,584 | stability optimal at 420 |
| Saba | 4.5 | stab P | – | no certified setting |  |  |  |
| Timoshenko_Beam_1 | 14 | stab P | stripped: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified | 144 | 22,336 | stability optimal at 144 |
| Timoshenko_Beam_2 | 15 |  | – | no certified setting |  |  |  |
| Timoshenko_Beam_3 | 16 | stab Q | – | no certified setting |  |  |  |
| Transport_Edge_Diffusion | 72 | stab Q | light: lpi.ineq.dD = 0; lpi.ineq.wcross = 1; sop.deg_Q{2} = 2 | certified | 256,948 | 676,410 | stability optimal at 256,948 |
| Transport_Eq | 1 | stab Q | extreme: lpi.ineq.dD = 0; lpi.ineq.cap_inc = 'auto' | certified; certified (stability_dual) | 134; 36 (stability_dual) | 1,396; 388 (stability_dual) | stability optimal at 134, stability_dual optimal at 36 |
| Transport_Eq_with_Disturbance | 20 | gain P | stripped: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [2;2;1] | gamma 0.500037; gamma 0.500195 (Hinf_gain_dual) | 633; 633 (Hinf_gain_dual) | 929; 929 (Hinf_gain_dual) | Hinf_gain optimal 0.500037 at 633, Hinf_gain_dual optimal 0.500195 at 633 |
| Transport_Eq_with_Disturbance2 | 49 |  | veryheavy | H2 norm 0.101576; H2 norm 0.101576 (H2_norm_dual) | 11,847; 11,847 (H2_norm_dual) | 1,563; 1,563 (H2_norm_dual) | H2_norm optimal 0.101576 at 11,847, H2_norm_dual optimal 0.101576 at 11,847 |
| Transport_Eq_with_Observed_Output | 27 |  | extreme: lpi.pos.len_sos{1} = [1;0;0]; lpi.var.deg{1} = [0 2 0]; sop.dw_P = 0 | closed-loop gain 1.00007 | 75 | 1,180 | Hinf_estimator optimal 1.00007 at 75 |
| Two_Layer_Heat_Rectangles | 64 | stab Q | extreme: sop.deg_Q{2} = [0 1 1] | certified | 288,324 | 288,400 | stability optimal at 288,324 |
| Unstable_ReactionDiffusion | 51 |  | heavy: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [3;2;1]; lpi.var.deg{1} = [2 2 2]; sop.dw_P = 0; sop.deg_Q{1} = [6 2 2] | closed-loop gain 1.26957 | 1,293 | 2,416 | H2_estimator optimal 1.26957 at 1,293 |
| Unstable_ReactionDiffusion_Control | 52 |  | heavy: lpi.pos.len_com{1} = 1; lpi.pos.len_sos{1} = [3;2;1]; lpi.var.deg{1} = [1 2 2]; sop.dw_P = 1; sop.deg_Q{1} = [6 2 2] | closed-loop gain 1.0908 | 2,471 | 2,523 | H2_control optimal 1.0908 at 2,471 |
| Wave_Delayed_Boundary_Damping_Xu | 59 | stab Q | light: lpi.ineq.dD = 0; lpi.ineq.psatz{2} = 'none'; sop.deg_Q{2} = 0; sop.max_raise_Q = [2 2 0]; sop.raise_Q = 'deg' | certified | 136 | 5,214 | stability optimal at 136 |
| Wave_Eq_Boundary_Damped | 17 | stab P | light: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified | 36 | 5,584 | stability optimal at 36 |
| Wave_Eq_Boundary_ODE_Control | 32 |  | veryheavy: lpi.pos.len_com{1} = 2; lpi.pos.len_sos{1} = [3;4;4]; lpi.var.deg{1} = [3 6 6] | closed-loop gain 0.518052 | 29,227 | 7,249 | Hinf_control optimal 0.518052 at 29,227 |
| Wave_Eq_Datko_Boundary_Damped_1 | 18 | stab P | stripped: lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;1]; lpi.ineq.dD = 0 | certified | 356 | 1,417 | stability_dual optimal at 356 |
| Wave_Eq_Datko_Boundary_Damped_2 | 19 | stab P | light: lpi.pos.len_inc{1} = 0; lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [1;2;2]; lpi.ineq.dD = 0 | certified | 81 | 11,875 | stability optimal at 81 |
| Wave_Eq_Tip_Damped | 21 | gain P | stripped: lpi.ineq.dD = 0 | gamma 0.841902; gamma 0.842084 (Hinf_gain_dual) | 1,401; 1,401 (Hinf_gain_dual) | 3,409; 3,409 (Hinf_gain_dual) | Hinf_gain optimal 0.841902 at 1,401, Hinf_gain_dual optimal 0.842084 at 1,401 |
| Wave_Eq_w_Boundary_Delay | 44 | stab Q | – | no certified setting |  |  |  |
| Wave_Eq_w_Boundary_Delay_2D | X |  | – | no certified setting |  |  |  |
| Wave_Heat_Transmission_1D_Batty | 68 | stab P | stripped: lpi.pos.len_com{1} = 0; lpi.pos.len_sos{1} = [0;2;1]; lpi.ineq.dD = 0 | certified | 324 | 12,564 | stability optimal at 324 |
| Wave_Heat_Transmission_2D_Batty | 69 | stab Q | – | no certified setting |  |  |  |
| PIETOOLS_PDE_Example_Diagne | X |  | – | no certified setting |  |  |  |
