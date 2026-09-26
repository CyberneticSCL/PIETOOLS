function prog = cx_Hinf_estimator(PIE,st,gam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = CX_HINF_ESTIMATOR(PIE,ST,GAM) builds, with copvar/cdopvar only, the
% LPI of 'PIETOOLS_Hinf_estimator' (Luenberger estimator synthesis, primal
% KYP in P and Z, L = inv(P)*Z) at the FIXED gain GAM, unsolved. The
% observer extraction (lpigetsol, getObserver) is not reproduced: there is
% no getsol for containers, and only the gamma threshold is compared.
%
% Mirrors executives/PIETOOLS_Hinf_estimator.m:
%   71-81    2-D PIE: the executive forwards to PIETOOLS_Hinf_estimator_2D,
%            which has no transcription here (error)
%   97-116   settings
%   129-132  gam dpvar/lpi_ineq/objective -> fixed double GAM (see cx_Hinf_gain)
%   146-153  P = poslpivar(dd1,options1) [+ psatz term]     -> cx_hinf_lf
%   156-157  P += blkdiag(eppos*I, eppos2*I)
%   159      Z = lpivar(Cyop.dim(:,[2,1]), ddZ): R^ny -> state
%   172-177  the 3 x 3 KYP operator, Tw terms included when Tw ~= 0. Z*Dyw
%            and -gam*Iw + Tw'*(...) need matching registries; D21 and Iw
%            are restated with cx_on_registry (metadata only)
%   186-200  negativity: slack + lpi_eq 'symmetric'         -> cx_hinf_slack
% epneg is read (line 104) but never used.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || ~isnumeric(gam)
    error('cx_Hinf_estimator:gam','A fixed numeric gamma is required.')
end
PIE = initialize(PIE);
if PIE.dim==2
    error('cx_Hinf_estimator:dim','PIETOOLS_Hinf_estimator_2D has no container transcription.')
end

Top = PIE.T;    Bwop = PIE.B1;    Czop = PIE.C1;    Cyop = PIE.C2;
Tm  = opvar2copvar(Top);        Am  = opvar2copvar(PIE.A);
Bw  = opvar2copvar(Bwop);       Cz  = opvar2copvar(Czop);
Cy  = opvar2copvar(Cyop);       Dzw = opvar2copvar(PIE.D11);
Iw  = opvar2copvar(mat2opvar(eye(size(Bwop,2)),Bwop.dim(:,2),PIE.vars,PIE.dom));
Iz  = opvar2copvar(mat2opvar(eye(size(Czop,1)),Czop.dim(:,1),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pm] = cx_hinf_lf(prog,Tm,PIE,st);
Imat = blkdiag(st.eppos*eye(Top.dim(1,:)),st.eppos2*eye(Top.dim(2,:)));
Pm = Pm + opvar2copvar(mat2opvar(Imat,Top.dim(:,2),PIE.vars,PIE.dom));

% Z: R^ny (Cy.dim(:,1)) -> state (Cy.dim(:,2)).
[dz,sz] = cx_hinf_dimsp(Cyop.dim(:,[2,1]),PIE);
[prog,Zm] = lpivar_cdopvar(prog,dz,sz,PIE.dom,st.ddZ);
Dyw = cx_on_registry(opvar2copvar(PIE.D21),Zm);    % registry bridge only

PB = Pm*Bw + Zm*Dyw;        PA = Pm*Am + Zm*Cy;
K11 = -(gam*Iw);
K13 = -PB'*Tm;              K31 = -Tm'*PB;
if ~(PIE.Tw==0)                                     % boundary disturbance
    Twm = opvar2copvar(PIE.Tw);
    X   = Twm'*PB;
    K11 = cx_on_registry(K11,X) + X + PB'*Twm;
    K13 = K13 - Twm'*PA;    K31 = K31 - PA'*Twm;
end
Km = [K11,      -Dzw',       K13;
      -Dzw,     -(gam*Iz),   Cz;
      K31,      Cz',         PA'*Tm + Tm'*PA];
prog = cx_hinf_slack(prog,Km,PIE,st);
end
