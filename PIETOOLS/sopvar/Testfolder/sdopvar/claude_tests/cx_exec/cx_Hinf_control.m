function prog = cx_Hinf_control(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_HINF_CONTROL(PIE,ST,GAM) builds, with copvar/cdopvar only, the
% LPI of 'PIETOOLS_Hinf_control' (state-feedback synthesis, dual KYP in
% P and Z, u = Z*inv(P)*x) at the FIXED gain GAM, unsolved. The controller
% extraction (lpigetsol, getController) is not reproduced: there is no
% getsol for containers, and only the gamma threshold is compared.
%
% Mirrors executives/PIETOOLS_Hinf_control.m:
%   77-79    2-D PIE: error, as the executive
%   86-88    Tw ~= 0 or Tu ~= 0: error, as the executive
%   99-118   settings
%   130-133  gam dpvar/lpi_ineq/objective -> fixed double GAM (see cx_Hinf_gain)
%   147-154  P = poslpivar(dd1,options1) [+ psatz term]     -> cx_hinf_lf
%   157-158  P += blkdiag(eppos*I, eppos2*I), as an opvar2copvar identity
%   160      Z = lpivar(Buop.dim(:,[2,1]), ddZ): state -> R^nu, via
%            lpivar_cdopvar with lpivar's legacy [d1 d2 d3] degrees
%   172-177  the 3 x 3 KYP operator. Dzu (R^nu -> R^nz) has an empty
%            registry and Z carries the L_2 variable; @cdopvar/mtimes
%            demands identical registries, so Dzu is restated with
%            cx_on_registry before Dzu*Z (same operator, metadata only)
%   188-206  negativity: slack + lpi_eq 'symmetric'         -> cx_hinf_slack
% epneg is read (line 106) but never used.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || ~isnumeric(gam)
    error('cx_Hinf_control:gam','A fixed numeric gamma is required.')
end
PIE = initialize(PIE);
if PIE.dim==2
    error('cx_Hinf_control:dim','Optimal control of 2D PIEs is currently not supported.')
end
if ~(PIE.Tw==0) || ~(PIE.Tu==0)
    error('cx_Hinf_control:Tw',['Hinf-dual LPI cannot currently be solved for '...
          'systems with disturbances or inputs at the boundary.'])
end

Top = PIE.T;    Bwop = PIE.B1;    Buop = PIE.Bu;    Czop = PIE.C1;
Tm  = opvar2copvar(Top);        Am  = opvar2copvar(PIE.A);
Bw  = opvar2copvar(Bwop);       Bu  = opvar2copvar(Buop);
Cz  = opvar2copvar(Czop);       Dzw = opvar2copvar(PIE.D11);
Iw  = opvar2copvar(mat2opvar(eye(size(Bwop,2)),Bwop.dim(:,2),PIE.vars,PIE.dom));
Iz  = opvar2copvar(mat2opvar(eye(size(Czop,1)),Czop.dim(:,1),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pm] = cx_hinf_lf(prog,Tm,PIE,st);
% Strict positivity; eye(dim(1,:)) is the executive's own idiom (n x n).
Imat = blkdiag(st.eppos*eye(Top.dim(1,:)),st.eppos2*eye(Top.dim(2,:)));
Pm = Pm + opvar2copvar(mat2opvar(Imat,Top.dim(:,2),PIE.vars,PIE.dom));

% Z: state (Bu.dim(:,1)) -> R^nu (Bu.dim(:,2)); zero-dimension spaces dropped.
[dz,sz] = cx_hinf_dimsp(Buop.dim(:,[2,1]),PIE);
[prog,Zm] = lpivar_cdopvar(prog,dz,sz,PIE.dom,st.ddZ);
Dzu = cx_on_registry(opvar2copvar(PIE.Dzu),Zm);    % registry bridge only

CZ = Cz*Pm + Dzu*Zm;        AZ = Am*Pm + Bu*Zm;
Km = [-(gam*Iz),   Dzw,         CZ*Tm';
       Dzw',       -(gam*Iw),   Bw';
       Tm*CZ',     Bw,          AZ*(Tm') + Tm*AZ'];
prog = cx_hinf_slack(prog,Km,PIE,st);
end
