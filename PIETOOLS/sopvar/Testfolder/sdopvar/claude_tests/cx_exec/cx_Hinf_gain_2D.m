function prog = cx_Hinf_gain_2D(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_HINF_GAIN_2D(PIE,ST,GAM) builds, with copvar/cdopvar only, the
% LPI of 'PIETOOLS_Hinf_gain_2D(PIE,settings,gain)' (2-D primal KYP,
% coercive, P >= eppos I) at the FIXED gain GAM, unsolved. Structure-only
% in this comparison: the program is assembled and its shape compared with
% the stock program, never solved.
%
% Mirrors executives/2D/PIETOOLS_Hinf_gain_2D.m:
%   57-62    operators; both Tw and B1 nonzero: error, as the executive
%   72-128   settings                                       -> cx_hinf_set2d
%   183-186  gain given: gam = gain (the objective branch 174-181 is not
%            used: a dpvar cannot multiply a container)
%   195-211  P = poslpivar_2d(LF) [+ psatz terms] + eppos    -> cx_hinf_lf2d
%   220-230  identities and the 3 x 3 KYP operator (vertcat/horzcat_legacy
%            merge w, z, state spaces; the container keeps them apart and
%            the slack is declared over the separate spaces, the same cone)
%   242-296  clean_opvar, get_eq_opts_2D, slack, lpi_eq_2d 'symmetric'
%                                                            -> cx_hinf_slack2d
% Operators are converted with cx_hinf_op2d, since opvar2d2copvar errors on
% an all-zero row or column (D11 = 0 here). Tw ~= 0: -gam*Iw is restated
% over the registry of Tw'*(P*B) (cx_on_registry) before the sum.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || ~isnumeric(gam) || gam==0
    error('cx_Hinf_gain_2D:gam','A fixed nonzero numeric gamma is required.')
end
PIE = initialize(PIE);
Top = PIE.T;    Twop = PIE.Tw;    Bop = PIE.B1;    Cop = PIE.C1;
if ~(Twop==0) && ~(Bop==0)
    error('cx_Hinf_gain_2D:TwB',['The PIE takes both the input w and its derivative '...
          'dw/dt; not supported (as the executive).'])
end
S = cx_hinf_set2d(st,[1e-4;1e-6;1e-6;1e-6]);

Tm = cx_hinf_op2d(Top);     Am = cx_hinf_op2d(PIE.A);
Bw = cx_hinf_op2d(Bop);     Cz = cx_hinf_op2d(Cop);     Dzw = cx_hinf_op2d(PIE.D11);
Iw = cx_hinf_op2d(mat2opvar(eye(size(Bop,2)),Bop.dim(:,2),PIE.vars,PIE.dom));
Iz = cx_hinf_op2d(mat2opvar(eye(size(Cop,1)),Cop.dim(:,1),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pm] = cx_hinf_lf2d(prog,Tm,PIE,S,true);

PB = Pm*Bw;     PA = Pm*Am;
K11 = -(gam*Iw);
K13 = PB'*Tm;               K31 = Tm'*PB;
if ~(Twop==0)
    Twm = cx_hinf_op2d(Twop);
    X   = Twm'*PB;
    K11 = cx_on_registry(K11,X) + X + PB'*Twm;
    K13 = K13 + Twm'*PA;    K31 = K31 + PA'*Twm;
end
Km = [K11,      Dzw',        K13;
      Dzw,      -(gam*Iz),   Cz;
      K31,      Cz',         PA'*Tm + Tm'*PA];
prog = cx_hinf_slack2d(prog,Km,PIE,S);
end
