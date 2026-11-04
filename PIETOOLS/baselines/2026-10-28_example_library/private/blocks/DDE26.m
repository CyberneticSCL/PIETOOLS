% examples_DDE_library_PIETOOLS.m lines 420-438, one comment level removed (lb_blocks)
% % % Example 2b from [13]
H2_norm=1;
H2_norm_dual=1;
DDE.A0=[-1 1 2
                1 -3 2
                0 0 -1];
nx=size(DDE.A0,1);
DDE.Ai{1}=(1/5).*[-3 0 1
                              1 -2 0
                              0 2 -2];
DDE.Ai{2}=(1/5).*[-4 1 0
                              0 -2 1
                              2 1 -3];
DDE.B1=[1
                1
                1];
DDE.C1=[1 1 1];
DDE.tau=[pi/10 1];
% % % %%%% answer - 3.8307
