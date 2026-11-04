% examples_DDE_library_PIETOOLS.m lines 74-80, one comment level removed (lb_blocks)
 % % % Example from [1], page 192
 stability=1
 stability_dual=1
    DDE.A0=[-2 0;0 -.9];%
    DDE.Ai{1}=[-1 0; -1 -1];%
    DDE.tau(1) = 6;
 % %%%% [1] says stable range is tau < 6.17258
