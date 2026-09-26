function prog = cx_Hinf_gain_dual_coercive(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_HINF_GAIN_DUAL_COERCIVE(PIE,ST,GAM) builds, with copvar/cdopvar
% only, the LPI of 'PIETOOLS_Hinf_gain_dual_coercive' (dual KYP, coercive,
% storage T*P*T') at the FIXED gain GAM, unsolved. (The stock file's
% function line names itself PIETOOLS_Hinf_gain_dual; MATLAB runs it under
% its file name.)
%
% Mirrors executives/PIETOOLS_Hinf_gain_dual_coercive.m:
%   71-78    2-D PIE: forwarded (cx_Hinf_gain_dual_2D)
%   85-87    Tw ~= 0: error, as the executive
%   99-116   settings
%   132-135  gam dpvar/lpi_ineq/objective -> fixed double GAM (see cx_Hinf_gain)
%   148-155  P = poslpivar(dd1,options1) [+ psatz term]     -> cx_hinf_lf
%   169-174  the 3 x 3 KYP operator; products kept in the executive's
%            left-to-right order, (Cz*P)*T' etc.
%   183-197  negativity: slack + lpi_eq 'symmetric'         -> cx_hinf_slack
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || ~isnumeric(gam)
    error('cx_Hinf_gain_dual_coercive:gam','A fixed numeric gamma is required.')
end
PIE = initialize(PIE);
if PIE.dim==2,  prog = cx_Hinf_gain_dual_2D(PIE,st,gam);    return,     end
if ~(PIE.Tw==0)
    error('cx_Hinf_gain_dual_coercive:Tw',['Hinf-dual LPI cannot be solved for '...
          'systems with disturbances at the boundary (as the stock executive).'])
end

Top = PIE.T;    Bwop = PIE.Bw;    Czop = PIE.Cz;
Tm  = opvar2copvar(Top);        Am  = opvar2copvar(PIE.A);
Bw  = opvar2copvar(Bwop);       Cz  = opvar2copvar(Czop);
Dzw = opvar2copvar(PIE.Dzw);
Iw  = opvar2copvar(mat2opvar(eye(size(Bwop,2)),Bwop.dim(:,2),PIE.vars,PIE.dom));
Iz  = opvar2copvar(mat2opvar(eye(size(Czop,1)),Czop.dim(:,1),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pm] = cx_hinf_lf(prog,Tm,PIE,st);

Km = [-(gam*Iz),     Dzw,         Cz*Pm*Tm';
       Dzw',         -(gam*Iw),   Bw';
       Tm*Pm*Cz',    Bw,          Tm*Pm*Am' + Am*Pm*Tm'];
prog = cx_hinf_slack(prog,Km,PIE,st);
end
