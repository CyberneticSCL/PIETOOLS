function prog = cx_Hinf_gain_dual(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_HINF_GAIN_DUAL(PIE,ST,GAM) builds, with copvar/cdopvar only, the
% LPI of 'PIETOOLS_Hinf_gain_dual' (dual KYP, NON-coercive, R = T*Q) at the
% FIXED gain GAM, unsolved.
%
% Mirrors executives/PIETOOLS_Hinf_gain_dual.m:
%   72-79    2-D PIE: forwarded (cx_Hinf_gain_dual_2D)
%   86-88    Tw ~= 0: error, as the executive
%   100-117  settings
%   133-136  gam dpvar/lpi_ineq/objective -> fixed double GAM (see cx_Hinf_gain)
%   149-156  R = poslpivar(dd1,options1) [+ psatz term]     -> cx_hinf_lf
%   159-161  Q = lpivar(Top.dim, get_lpivar_degs(R,T)), T*Q - R = 0 (the
%            executive reuses the T'Q sizing rule for T*Q; mirrored as is)
%   175-180  the 3 x 3 dual KYP operator
%   189-203  negativity: slack + lpi_eq 'symmetric'         -> cx_hinf_slack
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || ~isnumeric(gam)
    error('cx_Hinf_gain_dual:gam','A fixed numeric gamma is required.')
end
PIE = initialize(PIE);
if PIE.dim==2,  prog = cx_Hinf_gain_dual_2D(PIE,st,gam);    return,     end
if ~(PIE.Tw==0)
    error('cx_Hinf_gain_dual:Tw',['Hinf-dual LPI cannot be solved for systems '...
          'with disturbances at the boundary (as the stock executive).'])
end

Top = PIE.T;    Bwop = PIE.Bw;    Czop = PIE.Cz;
Tm  = opvar2copvar(Top);        Am  = opvar2copvar(PIE.A);
Bw  = opvar2copvar(Bwop);       Cz  = opvar2copvar(Czop);
Dzw = opvar2copvar(PIE.Dzw);
Iw  = opvar2copvar(mat2opvar(eye(size(Bwop,2)),Bwop.dim(:,2),PIE.vars,PIE.dom));
Iz  = opvar2copvar(mat2opvar(eye(size(Czop,1)),Czop.dim(:,1),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Rm] = cx_hinf_lf(prog,Tm,PIE,st);
Qdeg = cx_hinf_qdeg(Rm);                        % = get_lpivar_degs(Rop,Top)
[sp,dm] = cx_space_list(Tm,'out');
[prog,Qm] = lpivar_cdopvar(prog,dm,sp,PIE.dom,Qdeg);
prog = lpi_eq_cdopvar(prog,Tm*Qm - Rm);         % NOT symmetric (line 161)

Km = [-(gam*Iz),   Dzw,         Cz*Qm;
       Dzw',       -(gam*Iw),   Bw';
       Qm'*Cz',    Bw,          Qm'*Am' + Am*Qm];
prog = cx_hinf_slack(prog,Km,PIE,st);
end
