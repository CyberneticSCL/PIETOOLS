% PF_RUN_M4  Drop-a-face in the transport direction on S17 (r2 = 0.5 r2*):
% the faces of s1 only, of s2 only, the box factor of s1 only, of s2 only
% (per-direction offset 1), against facesS, at dw 1, on the storage STO
% (set below from M1/M3). SCRATCH, MMP, 11/16/2026.
STO = {'def'};
X = struct('name',{'faces_s1','faces_s2','box_s1','box_s2'}, ...
           'rows',{[1 0;2 0],[0 1;0 2],[3 0],[0 3]}, ...
           'offs',{[1 1],[1 1],1,1},'pd',{true,true,true,true});
C = pf_cfg({'facesS'},2,X);
pf_exec_m('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},C,1,STO,'light','m4_s17');
disp('## M4 DONE');
