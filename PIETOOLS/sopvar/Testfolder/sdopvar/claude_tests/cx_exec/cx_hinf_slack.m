function prog = cx_hinf_slack(prog,Km,PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_HINF_SLACK(PROG,KM,PIE,ST) imposes KM <= 0 exactly as the 1-D
% H-infinity executives' default branch does:
%
%   [prog,De1op] = poslpivar(prog,Dop.dim,dd2,options2);
%   if override2~=1
%       [prog,De2op] = poslpivar(prog,Dop.dim,dd3,options3);  De = De1+De2;
%   end
%   prog = lpi_eq(prog,Deop+Dop,'symmetric');
%
% (PIETOOLS_Hinf_gain.m:189-197, _coercive:171-179, _dual:194-202,
% _dual_coercive:188-196, Hinf_control:195-205, Hinf_estimator:191-199).
%
% The slack is declared over KM's OWN spaces (w, z, state kept separate by
% the container concatenation), where the stock legacy concatenation merges
% all R parts into one R block and all L_2 parts into one L_2 block. One Gram
% over the separate spaces is a permutation of the stock Gram over the
% merged ones (I_n kron Z vs n copies of Z), so the cones coincide.
%
% The sosineq_on branch (lpi_ineq(prog,-Dop,opts)) has no container form:
% there is no lpi_ineq for containers. It is an error, not a silent switch.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if st.sosineq_on
    error('cx_hinf_slack:noLpiIneq',['settings.sosineq_on=1 selects lpi_ineq, '...
          'which has no container counterpart.'])
end
[sp,dm] = cx_space_list(Km,'out');
[deg,co] = cx_hinf_posdeg(st.dd2,st.options2,sp);
[prog,Nm] = poscopvar(prog,dm,sp,PIE.dom,deg,co);
if st.override2~=1
    [deg,co] = cx_hinf_posdeg(st.dd3,st.options3,sp);  % light: psatz = 1
    [prog,N2] = poscopvar(prog,dm,sp,PIE.dom,deg,co);
    Nm = Nm + N2;
end
prog = lpi_eq_cdopvar(prog,Nm + Km,'symmetric');       % Deop + Dop = 0
end
