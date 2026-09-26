function prog = cx_Hinf_gain_coercive(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_HINF_GAIN_COERCIVE(PIE,ST,GAM) builds, with copvar/cdopvar only,
% the LPI of 'PIETOOLS_Hinf_gain_coercive' (primal KYP, V = <Tx,PTx>) at the
% FIXED gain GAM, unsolved.
%
% Mirrors executives/PIETOOLS_Hinf_gain_coercive.m:
%   62-70    2-D PIE: forwarded (cx_Hinf_gain_2D)
%   86-103   settings
%   116-119  gam dpvar/lpi_ineq/objective -> fixed double GAM (see cx_Hinf_gain)
%   133-140  P = poslpivar(dd1,options1) [+ psatz term]     -> cx_hinf_lf
%   152-157  the 3 x 3 KYP operator, including the Tw terms
%   166-180  negativity: slack + lpi_eq 'symmetric'         -> cx_hinf_slack
% Tw ~= 0: -gam*Iw has an empty registry (R -> R) while Tw'*(P*Bw) carries
% the L_2 variable, and @cdopvar/plus demands identical registries, so -gam*Iw
% is restated with cx_on_registry first. Tw == 0: the Tw terms are zero
% operators and are dropped, which leaves the same LPI.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || ~isnumeric(gam)
    error('cx_Hinf_gain_coercive:gam','A fixed numeric gamma is required.')
end
PIE = initialize(PIE);
if PIE.dim==2,  prog = cx_Hinf_gain_2D(PIE,st,gam);     return,     end

Top = PIE.T;    Bwop = PIE.Bw;    Czop = PIE.Cz;
Tm  = opvar2copvar(Top);        Am  = opvar2copvar(PIE.A);
Bw  = opvar2copvar(Bwop);       Cz  = opvar2copvar(Czop);
Dzw = opvar2copvar(PIE.Dzw);
Iw  = opvar2copvar(mat2opvar(eye(size(Bwop,2)),Bwop.dim(:,2),PIE.vars,PIE.dom));
Iz  = opvar2copvar(mat2opvar(eye(size(Czop,1)),Czop.dim(:,1),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pm] = cx_hinf_lf(prog,Tm,PIE,st);

PB = Pm*Bw;     PA = Pm*Am;
K11 = -(gam*Iw);
K13 = PB'*Tm;               K31 = Tm'*PB;
if ~(PIE.Tw==0)                                 % boundary disturbance
    Twm = opvar2copvar(PIE.Tw);
    X   = Twm'*PB;
    K11 = cx_on_registry(K11,X) + X + PB'*Twm;
    K13 = K13 + Twm'*PA;    K31 = K31 + PA'*Twm;
end
Km = [K11,      Dzw',        K13;
      Dzw,      -(gam*Iz),   Cz;
      K31,      Cz',         PA'*Tm + Tm'*PA];
prog = cx_hinf_slack(prog,Km,PIE,st);
end
