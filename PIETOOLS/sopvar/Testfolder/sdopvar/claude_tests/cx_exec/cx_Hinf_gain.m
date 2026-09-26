function prog = cx_Hinf_gain(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_HINF_GAIN(PIE,ST,GAM) builds, with copvar/cdopvar only, the LPI
% of 'PIETOOLS_Hinf_gain' (primal KYP, NON-coercive Q-form) at the FIXED
% gain GAM, and returns it unsolved. ST is an lpisettings struct.
%
% Mirrors executives/PIETOOLS_Hinf_gain.m:
%   63-71    2-D PIE: forwarded (here to cx_Hinf_gain_2D)
%   78-86    Tw ~= 0: forwarded to the coercive executive
%   98-115   settings read as the executive reads them
%   128-131  gam dpvar + lpi_ineq(gam>=0) + objective: REPLACED by the fixed
%            double GAM. A dpvar cannot multiply a container (no dpvar
%            branch in @cdopvar/@copvar mtimes), and there is no container
%            lpi_ineq; cx_compare bisects GAM instead.
%   145-152  R = poslpivar(dd1,options1) [+ psatz term]     -> cx_hinf_lf
%   155-157  Q = lpivar(Top.dim, get_lpivar_degs(R,T)), T'Q - R = 0
%            -> cx_hinf_qdeg, lpivar_cdopvar (legacy [d1 d2 d3] mapping,
%               equal to lpivar in 1-D, test_lpivar_cdopvar), and
%               lpi_eq_cdopvar without 'symmetric'
%   170-175  the 3 x 3 KYP operator
%   184-198  negativity: slack + lpi_eq 'symmetric'         -> cx_hinf_slack
% eppos/eppos2 are read by the executive (104-105) but never used.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || ~isnumeric(gam)
    error('cx_Hinf_gain:gam','A fixed numeric gamma is required.')
end
PIE = initialize(PIE);
if PIE.dim==2,          prog = cx_Hinf_gain_2D(PIE,st,gam);         return,     end
if ~(PIE.Tw==0),        prog = cx_Hinf_gain_coercive(PIE,st,gam);   return,     end

Top = PIE.T;    Bwop = PIE.Bw;    Czop = PIE.Cz;
Tm  = opvar2copvar(Top);        Am  = opvar2copvar(PIE.A);
Bw  = opvar2copvar(Bwop);       Cz  = opvar2copvar(Czop);
Dzw = opvar2copvar(PIE.Dzw);    % zero D11 -> explicit zero block
Iw  = opvar2copvar(mat2opvar(eye(size(Bwop,2)),Bwop.dim(:,2),PIE.vars,PIE.dom));
Iz  = opvar2copvar(mat2opvar(eye(size(Czop,1)),Czop.dim(:,1),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Rm] = cx_hinf_lf(prog,Tm,PIE,st);
Qdeg = cx_hinf_qdeg(Rm);                        % = get_lpivar_degs(Rop,Top)
[sp,dm] = cx_space_list(Tm,'out');
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,Qdeg);
prog = lpi_eq_cdopvar(prog,Tm'*Qm - Rm);        % NOT symmetric (line 157)

Km = [-(gam*Iw),   Dzw',        Bw'*Qm;
       Dzw,        -(gam*Iz),   Cz;
       Qm'*Bw,     Cz',         Am'*Qm + Qm'*Am];
prog = cx_hinf_slack(prog,Km,PIE,st);
end
