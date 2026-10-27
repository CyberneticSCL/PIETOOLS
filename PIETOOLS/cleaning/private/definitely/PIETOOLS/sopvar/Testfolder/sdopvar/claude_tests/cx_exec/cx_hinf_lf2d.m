function [prog,Pm] = cx_hinf_lf2d(prog,Tm,PIE,S,addeppos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,PM] = CX_HINF_LF2D(PROG,TM,PIE,S,ADDEPPOS) declares the 2-D storage
% operator with poscopvar over TM's state spaces, as the 2-D executives do:
%
%   [prog,Pop] = poslpivar_2d(prog,Top.dim,LF_deg,LF_opts);
%   for j: if LF_use_psatz(j)~=0, Pop = Pop + poslpivar_2d(...psatz j...)
%   if ~all(eppos==0), Pop = Pop + opvar2d(blkdiag(eppos(k)*I),...)  [ADDEPPOS]
%
% (PIETOOLS_Hinf_gain_2D.m:195-211, Hinf_gain_dual_2D:164-180; the
% non-coercive executive has no eppos, _non_coercive:191-199). S is from
% cx_hinf_set2d; degrees/options via cx_hinf_2dpos. The domain is passed
% name-keyed, so the sorted registry order cannot permute it.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dom = struct();     dom.vars = reshape(PIE.vars(:,1).varname,1,[]);    dom.dom = PIE.dom;
[sp,dm] = cx_space_list(Tm,'out');
[deg,co] = cx_hinf_2dpos(S.LF_deg,S.LF_opts,sp);
[prog,Pm] = poscopvar(prog,dm,sp,dom,deg,co);
for j = 1:numel(S.LF_use_psatz)
    if S.LF_use_psatz(j)~=0
        [deg,co] = cx_hinf_2dpos(S.LF_deg_psatz{j},S.LF_opts_psatz{j},sp);
        [prog,P2] = poscopvar(prog,dm,sp,dom,deg,co);
        Pm = Pm + P2;
    end
end
if addeppos && ~all(S.eppos==0)
    n = PIE.T.dim(:,1);     e = S.eppos;
    Ip = blkdiag(e(1)*eye(n(1)),e(2)*eye(n(2)),e(3)*eye(n(3)),e(4)*eye(n(4)));
    % cx_hinf_op2d: a zero eppos entry on a present space leaves a zero
    % diagonal component, which opvar2d2copvar cannot convert.
    Pm = Pm + cx_hinf_op2d(opvar2d(Ip,PIE.T.dim,PIE.dom,PIE.vars));
end
end
