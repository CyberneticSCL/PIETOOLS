function B = t1d_build(PIE,spec)                                            % CC, 09/22/2026
% PROVENANCE.  scratchpad/reach1d/r1build.m verbatim, renamed (B6: no two
% copies of one name on the path) and with lset -> t1d_set for the same
% reason.  bm_setup is the SHIPPED file, byte-identity asserted by T0.
%                                                                 % CC, 09/22/2026
% r1build(PIE,spec) -- assemble the direct-form 1-D stability LPI at one degree
% setting and extract its SeDuMi data.  spec is a stock name or an lset knob
% vector.  Same LPI as PIETOOLS_PDEstability (build_stab_st is build_stability2
% with a settings STRUCT), so the sizes are comparable to the banked ladder.
st = t1d_set(spec);                                                         % CC, 09/22/2026
t0 = tic;
[prog,H] = build_stab_st(PIE,st);
B.t_build = toc(t0);
[Atf,bf,Ns,Kf] = raw_data(prog);
t0 = tic;
P = bm_setup(Atf,bf,Ns,Kf,1);   % pre=1: the whitened metric both arms report in
B.t_pre = toc(t0);
B.prog=prog; B.H=H; B.Atf=Atf; B.bf=bf; B.Ns=Ns; B.Kf=Kf; B.P=P; B.st=st;
B.m = size(Atf,2);   B.Ntot = size(Atf,1);   B.sumNs2 = sum(Ns.^2);
B.rankA = P.rankA;
end
