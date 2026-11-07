% hw_repro.m -- reproduce the 19.7 GB case with output UNSUPPRESSED.
% The probe showed the BUILD is flat at ~2.1 GB for every n/settings/w
% combination, and n=2 heavy is only m=1064. So either the solve is the wall or
% bl_b_io1 does something the mirror does not. bl_run wraps the builder in
% evalc, which hid whatever the executive printed last; here it prints, so the
% stage is visible in the log when memory starts climbing.
cuadmm_path;
M0 = memory; fprintf('HR start mem=%.2f GB\n',M0.MemUsedMATLAB/2^30);
t0 = tic;
[sol,Mx] = cuadmm_private('bl_b_io1','Hinf_gain','heavy',2);
M1 = memory;
S = cuadmm_private('sdpshape',sol);
fprintf('HR done t=%.1f s mem=%.2f GB m=%d Ks=%s\n', ...
        toc(t0),M1.MemUsedMATLAB/2^30,S.m,mat2str(S.Ks));
fprintf('HRDONE\n');
