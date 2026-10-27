function prog = cx_Hinf_gain_2D_non_coercive(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_HINF_GAIN_2D_NON_COERCIVE(PIE,ST,GAM) builds, with copvar/cdopvar
% only, the LPI of 'PIETOOLS_Hinf_gain_2D_non_coercive(PIE,settings,gain)'
% (2-D primal KYP, R = T'*PT, R >= 0) at the FIXED gain GAM, unsolved.
% Structure-only in this comparison.
%
% Mirrors executives/2D/PIETOOLS_Hinf_gain_2D_non_coercive.m:
%   59-64    operators; both Tw and B1 nonzero: error
%   67-77    Tw ~= 0: forwarded to the coercive executive (cx_Hinf_gain_2D)
%   87-126   settings                                       -> cx_hinf_set2d
%   178-182  gain given: gam = gain (objective branch not used)
%   191-199  R = poslpivar_2d(LF) [+ psatz], NO eppos         -> cx_hinf_lf2d
%   202-204  PT = lpivar(Top.dim, get_lpivar_degs(R,T)), T'*PT - R = 0.
%            lpivar_cdopvar cannot take lpivar_2d's per-cell joint caps:
%            per-role maxima of R's degrees (cx_hinf_qdeg2d), a family
%            that CONTAINS the stock PT's, so this program is NESTED around
%            the stock one, not equal. lpi_eq_cdopvar without 'symmetric'.
%   214-220  identities and the 3 x 3 KYP operator
%   222-277  clean_opvar, get_eq_opts_2D, slack, lpi_eq_2d 'symmetric'
%                                                            -> cx_hinf_slack2d
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || ~isnumeric(gam) || gam==0
    error('cx_Hinf_gain_2D_non_coercive:gam','A fixed nonzero numeric gamma is required.')
end
PIE = initialize(PIE);
Top = PIE.T;    Twop = PIE.Tw;    Bop = PIE.B1;    Cop = PIE.C1;
if ~(Twop==0) && ~(Bop==0)
    error('cx_Hinf_gain_2D_non_coercive:TwB',['The PIE takes both the input w and '...
          'its derivative dw/dt; not supported (as the executive).'])
end
if ~(Twop==0),  prog = cx_Hinf_gain_2D(PIE,st,gam);     return,     end
S = cx_hinf_set2d(st,[1e-4;1e-6;1e-6;1e-6]);    % eppos is not used here

Tm = cx_hinf_op2d(Top);     Am = cx_hinf_op2d(PIE.A);
Bw = cx_hinf_op2d(Bop);     Cz = cx_hinf_op2d(Cop);     Dzw = cx_hinf_op2d(PIE.D11);
Iw = cx_hinf_op2d(mat2opvar(eye(size(Bop,2)),Bop.dim(:,2),PIE.vars,PIE.dom));
Iz = cx_hinf_op2d(mat2opvar(eye(size(Cop,1)),Cop.dim(:,1),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Rm] = cx_hinf_lf2d(prog,Tm,PIE,S,false);
dom = struct();     dom.vars = reshape(PIE.vars(:,1).varname,1,[]);    dom.dom = PIE.dom;
[sp,dm] = cx_space_list(Tm,'out');
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,dom,cx_hinf_qdeg2d(Rm));
prog = lpi_eq_cdopvar(prog,Tm'*Qm - Rm);        % NOT symmetric (line 204)

Km = [-(gam*Iw),   Dzw',        Bw'*Qm;
       Dzw,        -(gam*Iz),   Cz;
       Qm'*Bw,     Cz',         Am'*Qm + Qm'*Am];
prog = cx_hinf_slack2d(prog,Km,PIE,S);
end
