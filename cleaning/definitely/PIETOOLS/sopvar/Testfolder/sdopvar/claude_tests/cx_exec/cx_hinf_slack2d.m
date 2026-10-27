function [prog,eq_opts] = cx_hinf_slack2d(prog,Km,PIE,S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,EQ_OPTS] = CX_HINF_SLACK2D(PROG,KM,PIE,S) imposes KM <= 0 as the
% 2-D H-infinity executives' default branch does:
%
%   Qop = clean_opvar(Qop,1e-12);
%   eq_opts = get_eq_opts_2D(Qop,eq_opts,ztol);
%   [prog,Qeop] = poslpivar_2d(prog,Qop.dim,eq_deg,eq_opts);
%   for j: psatz terms with exclude/sep OR-ed with eq_opts
%   prog = lpi_eq_2d(prog,Qeop+Qop,'symmetric');
%
% (PIETOOLS_Hinf_gain_2D.m:241-296, _non_coercive:222-277,
% Hinf_gain_dual_2D:204-259). get_eq_opts_2D -> cx_hinf_eqopts2d (which
% also stands in for clean_opvar in its tests); the KYP operator itself is
% NOT cleaned - there is no clean_opvar for containers - so coefficients
% below 1e-12 that the stock drops stay in KM. The lpi_ineq_2d branch
% (use_sosineq) has no container form: error. The checkdeg_lpi_eq_2d loop
% is behind toggle = 0 in every executive, so it is not mirrored.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if S.use_sosineq
    error('cx_hinf_slack2d:noLpiIneq','use_sosineq selects lpi_ineq_2d, which has no container counterpart.')
end
ztol = 1e-12;
dom = struct();     dom.vars = reshape(PIE.vars(:,1).varname,1,[]);    dom.dom = PIE.dom;
eq_opts = cx_hinf_eqopts2d(Km,S.eq_opts,ztol);
[sp,dm] = cx_space_list(Km,'out');
[deg,co] = cx_hinf_2dpos(S.eq_deg,eq_opts,sp);
[prog,Nm] = poscopvar(prog,dm,sp,dom,deg,co);
for j = 1:numel(S.eq_use_psatz)
    if S.eq_use_psatz(j)~=0
        o = S.eq_opts_psatz{j};
        o.exclude = o.exclude | eq_opts.exclude;    o.sep = o.sep | eq_opts.sep;
        [deg,co] = cx_hinf_2dpos(S.eq_deg_psatz{j},o,sp);
        [prog,N2] = poscopvar(prog,dm,sp,dom,deg,co);
        Nm = Nm + N2;
    end
end
prog = lpi_eq_cdopvar(prog,Nm + Km,'symmetric');       % Qeop + Qop = 0
end
