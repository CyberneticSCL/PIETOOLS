# M1: term sets at the N-D default storage (light, MOSEK, dw pinned)
| system | parameter | storage | set | dw | status | cert | ndv | m | t (s) | numerr | feasratio | eta |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | none | 0 | infeasible | 0 | 9700 | 1078 | 3 | 0 | -1 | 2.1e-02 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | none | 1 | infeasible | 0 | 83428 | 2716 | 7 | 0 | -1 | 1.1e-02 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | box | 0 | infeasible | 0 | 21220 | 1870 | 1 | 0 | -1 | 2.6e-02 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | box | 1 | unknown | 0 | 141028 | 2716 | 7 | 2 | 1.77 | 2.3e-10 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | facesS | 0 | infeasible | 0 | 18916 | 1294 | 1 | 0 | -1 | 1.3e-02 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | facesS | 1 | unknown | 0 | 120292 | 2716 | 9 | 2 | 3.28 | 5.5e-12 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | faces0 | 0 | infeasible | 0 | 46564 | 1717 | 1 | 0 | -1 | 1.1e-02 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | faces0 | 1 | unknown | 0 | 415204 | 3715 | 16 | 2 | 0.0309 | 7.0e-12 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | corners | 0 | infeasible | 0 | 41956 | 1474 | 2 | 0 | -1 | 8.7e-03 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | def | corners | 1 | unknown | 0 | 235492 | 2716 | 10 | 2 | 2.83 | 1.2e-10 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | none | 0 | infeasible | 0 | 7540 | 1022 | 1 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | none | 1 | infeasible | 0 | 63988 | 2556 | 3 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | box | 0 | infeasible | 0 | 16360 | 1764 | 1 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | box | 1 | infeasible | 0 | 108088 | 2556 | 4 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | facesS | 0 | infeasible | 0 | 14596 | 1225 | 1 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | facesS | 1 | infeasible | 0 | 92212 | 2556 | 5 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | faces0 | 0 | infeasible | 0 | 35764 | 1621 | 1 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | faces0 | 1 | infeasible | 0 | 318004 | 3491 | 12 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | corners | 0 | infeasible | 0 | 32236 | 1393 | 1 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | corners | 1 | infeasible | 0 | 180412 | 2556 | 7 | 0 | -1 | 1.0e+00 |

# M3: storage variants (every basis operator at the light degrees; sepAB: direction split where the digit is 0), facesS
| system | parameter | storage | set | dw | status | cert | ndv | m | t (s) | numerr | feasratio | eta |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | sep11 | facesS | 0 | unknown | 0 | 47376 | 1882 | 5 | 2 | 1.32 | 4.3e-11 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | sep11 | facesS | 1 | unknown | 0 | 335376 | 3880 | 17 | 2 | 0.0749 | 1.9e-11 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | sep01 | facesS | 0 | unknown | 0 | 49801 | 1903 | 3 | 2 | -0.118 | 1.2e-10 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | sep01 | facesS | 1 | unknown | 0 | 337801 | 3895 | 17 | 2 | -0.0593 | 6.6e-11 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | sep10 | facesS | 0 | unknown | 0 | 48996 | 2012 | 3 | 2 | -0.254 | 5.0e-11 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | sep10 | facesS | 1 | unknown | 0 | 336996 | 3936 | 18 | 2 | -0.0104 | 8.4e-11 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | sep00 | facesS | 0 | unknown | 0 | 54180 | 2252 | 4 | 2 | 0.499 | 5.4e-10 |
| S17_transport_edge_diffusion | {"r2":4.934802200544679} | sep00 | facesS | 1 | unknown | 0 | 342180 | 4044 | 16 | 2 | -0.0433 | 2.9e-11 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep11 | facesS | 0 | infeasible | 0 | 36576 | 1775 | 2 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep11 | facesS | 1 | infeasible | 0 | 257076 | 3645 | 12 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep10 | facesS | 0 | infeasible | 0 | 38196 | 1825 | 2 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep10 | facesS | 1 | infeasible | 0 | 258696 | 3672 | 14 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep01 | facesS | 0 | infeasible | 0 | 39001 | 2039 | 2 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep01 | facesS | 1 | infeasible | 0 | 259501 | 3768 | 16 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep00 | facesS | 0 | infeasible | 0 | 43380 | 2329 | 4 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep00 | facesS | 1 | unknown | 0 | 263880 | 3920 | 19 | 2 | -1 | 1.8e-09 |

# M6: S18 with the storage degrees raised (mult 1, int 2 per variable, joint 4, every basis operator; PF_POS override), facesS
| system | parameter | storage | set | dw | status | cert | ndv | m | t (s) | numerr | feasratio | eta |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep10 | facesS | 0 | optimal | 1 | 790884 | 7408 | 90 | 0 | 1.24 | 1.7e-09 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep10 | facesS | 1 | optimal | 1 | 3045924 | 12006 | 740 | 0 | 0.48 | 9.7e-09 |
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep00 | facesS | 1 | optimal | 1 | 3125056 | 12550 | 789 | 0 | 0.995 | 1.7e-08 |

# M8: S18 at the raised degrees (mult 1, int 2, joint 4): the separable storage with every basis operator ('sep11', control for the split) at 0.5 tau*; the s2-split storage at 1.1 tau* (unstable) and 0.9 tau*
| system | parameter | storage | set | dw | status | cert | ndv | m | t (s) | numerr | feasratio | eta |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | sep11 | facesS | 0 | optimal | 1 | 767305 | 7280 | 64 | 0 | 0.987 | 3.7e-09 |
| S18_lib43_interior_delay | {"tau":1.1381986000000002,"sgn":-1} | sep10 | facesS | 0 | infeasible | 0 | 790884 | 7408 | 146 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":0.9312534,"sgn":-1} | sep10 | facesS | 0 | optimal | 1 | 790884 | 7408 | 138 | 0 | -0.363 | 6.2e-11 |

# M9: S18, the DEG-omitted default storage ('def') with lpi.pos raised to mult 1, int 2, joint 4 (PF_POS read into st.lpi.pos), facesS, dw 0
| system | parameter | storage | set | dw | status | cert | ndv | m | t (s) | numerr | feasratio | eta |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S18_lib43_interior_delay | {"tau":0.517363,"sgn":-1} | def | facesS | 0 | infeasible | 0 | 109921 | 3735 | 18 | 0 | -1 | 1.0e+00 |
| S18_lib43_interior_delay | {"tau":1.1381986000000002,"sgn":-1} | def | facesS | 0 | infeasible | 0 | 109921 | 3736 | 14 | 0 | -1 | 1.0e+00 |
