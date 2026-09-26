function prog = cx_Hinf_gain_dual_2D(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_HINF_GAIN_DUAL_2D(PIE,ST,GAM) builds, with copvar/cdopvar only,
% the LPI of 'PIETOOLS_Hinf_gain_dual_2D(PIE,settings,gain)' (2-D dual KYP,
% coercive, P >= eppos I) at the FIXED gain GAM, unsolved. Structure-only
% in this comparison.
%
% Mirrors executives/2D/PIETOOLS_Hinf_gain_dual_2D.m:
%   62-67    operators; Tw ~= 0: error, as the executive
%   77-133   settings; eppos default [1e-4;0;0;1e-6]         -> cx_hinf_set2d
%   151-155  gain given: gam = gain (objective branch not used)
%   164-180  P = poslpivar_2d(LF) [+ psatz] + eppos           -> cx_hinf_lf2d
%   189-196  identities and the 3 x 3 KYP operator, transcribed AS WRITTEN:
%            (C*P)*T and T'*(C*P)' in the (1,3)/(3,1) blocks, where the 1-D
%            dual coercive executive has C*P*T' and T*P*C'. The two agree
%            only when T = T* (true for io2, Dirichlet on every edge); the
%            stock form is kept so the SDPs can be compared.
%   204-259  clean_opvar, get_eq_opts_2D, slack, lpi_eq_2d 'symmetric'
%                                                            -> cx_hinf_slack2d
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || ~isnumeric(gam) || gam==0
    error('cx_Hinf_gain_dual_2D:gam','A fixed nonzero numeric gamma is required.')
end
PIE = initialize(PIE);
Top = PIE.T;    Bop = PIE.B1;    Cop = PIE.C1;
if ~(PIE.Tw==0)
    error('cx_Hinf_gain_dual_2D:Tw',['The PIE takes the derivative dw/dt as input; '...
          'dual H_infty gain analysis is not supported (as the executive).'])
end
S = cx_hinf_set2d(st,[1e-4;0;0;1e-6]);

Tm = cx_hinf_op2d(Top);     Am = cx_hinf_op2d(PIE.A);
Bw = cx_hinf_op2d(Bop);     Cz = cx_hinf_op2d(Cop);     Dzw = cx_hinf_op2d(PIE.D11);
Iw = cx_hinf_op2d(mat2opvar(eye(size(Bop,2)),Bop.dim(:,2),PIE.vars,PIE.dom));
Iz = cx_hinf_op2d(mat2opvar(eye(size(Cop,1)),Cop.dim(:,1),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pm] = cx_hinf_lf2d(prog,Tm,PIE,S,true);

CP = Cz*Pm;     AP = Am*Pm;
Km = [-(gam*Iz),    Dzw,         CP*Tm;
       Dzw',        -(gam*Iw),   Bw';
       Tm'*CP',     Bw,          Tm*AP' + AP*Tm'];
prog = cx_hinf_slack2d(prog,Km,PIE,S);
end
