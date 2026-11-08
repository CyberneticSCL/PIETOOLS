function PIEcl = closedLoopPIE(PIE, K, type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% closedLoopPIE.m     PIETOOLS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a PIE of the form
% T \dot x +Tw w + Tu u= A x + B1 w + B2 u
% z = C1 x + D11 w + D12 u
% y = C2 x + D21 w + D22 u
% where the PI operators are stored as fields of PIE object
% this function returns the closed loop pie_struct object that corresponds to the 
% PIE 
% if type='controller' (default)
% (T + Tu*K) \dot x +Tw w = (A + B2*K) x + B1 w 
% z = (C1 + D12*K) x + D11 w 
% y = (C2 + D22*K) x + D21 w 
% 
% if type='observer'
% [T 0][\dot x   ] +[Tw] w +[Tu] u = [A     0     ][x   ] +[B1   ] w +[B2   ] u 
% [0 T][\dot xhat]  [0 ]    [Tu ]     [-K*C2 A+K*C2][xhat] +[-K*D21]w+   [B2   ]u    % DB, 11/04/2026
% 
% z = [C1 0 ][x   ] +[D11] w +[D12] u
%     [0  C1][xhat]  [0  ]    [D12]    % DB, 11/04/2026
%
% if type='dynamic', the dynamic output feedback controller    % DB, 11/03/2026
% TK \dot xK = AK xK + BK y,    u = CK xK + DK y    % DB, 11/03/2026
% (as returned by PIETOOLS_Hinf_OFcontrol_sop), for D22 = 0, so that    % DB, 11/03/2026
% [T 0 ][\dot x ] = [A + B2*DK*C2   B2*CK][x ] + [B1 + B2*DK*D21] w    % DB, 11/03/2026
% [0 TK][\dot xK]   [BK*C2          AK   ][xK]   [BK*D21         ]    % DB, 11/03/2026
% z = [C1 + D12*DK*C2,  D12*CK][x; xK] + (D11 + D12*DK*D21) w    % DB, 11/03/2026
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE: a linear PIE, as a pie_struct object
% K : an opvar, opvar2d, or copvar object, mapping the PIE state to the
%     actuator input u (controller) or the observed output y to the PIE
%     state (observer); for type='dynamic', a struct with fields TK, AK,    % DB, 11/03/2026
%     BK, CK, DK, the operators (opvar, opvar2d or copvar objects) of the    % DB, 11/03/2026
%     controller PIE, or a struct with such a field G    % DB, 11/03/2026
% type: optional string input; 'controller' (default), 'observer' or    % DB, 11/03/2026
%     'dynamic'    % DB, 11/03/2026
% If Dyu=Tw=Tu=0, then given a state feedback controller K and an estimator gain L, DB, 11/03/2026
% the following call gives the closed-loop estimator-based control   DB, 11/03/2026
% Kc.TK = PIE.T;% DB, 11/03/2026
% Kc.AK = PIE.A + PIE.Bu*K + L*PIE.Cy;% DB, 11/03/2026
% Kc.BK = -L;% DB, 11/03/2026
% Kc.CK = K;% DB, 11/03/2026
% Kc.DK = 0*(K*L);        % DB, 11/03/2026
% PIEcl = closedLoopPIE(PIE,Kc,'dynamic'); % DB, 11/03/2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIEcl: an initialized pie_struct object, with the operators Tcl Twcl Acl
%        B1cl C1cl D11cl C2cl D21cl, as copvar objects. For the controller,
%        the actuator input u is eliminated, so that Tu, B2, D12 and D22
%        are zero. For the observer, the state is [x;xhat], and the
%        observed output is eliminated. For the dynamic controller, the    % DB, 11/03/2026
%        state is [x;xK], and both u and y are eliminated    % DB, 11/03/2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding SS - 5/30/2021
% added observer closed-loop construction, SS - 1/18/2024
% DJ, 12/07/2024: Make sure variables of output match those of input;
% DJ, 01/04/2025: Add support for 2D PIEs;
% DJ, 10/29/2026: Update to new pie_struct structure, PIEs with algebraic
%                   constraints are not supported;
% DB, 11/03/2026: Observer: the disturbance enters the xhat equation as -L*D21*w
%                   (xhat_t = A xhat + L(C2 xhat - y)), not +L*D21*w; the term
%                   L*D22*u (xhat row of Bu) was not checked;
% DB, 11/04/2026: Observer: the xhat row of Bu is B2 (not L*D22): the estimator
%                   xhat_t = A xhat + B2 u + L(C2 xhat + D22 u - y) has the same u
%                   path as the plant, so that u, including L*D22*u, cancels from
%                   the error xhat-x, as in the estimator executives;
%                   Likewise the xhat row of Tu is Tu, so that Tu*u_t cancels too;
% DB, 11/04/2026: Observer: the estimate zhat = C1 xhat + D12 u includes the term
%                   D12*u (xhat row of Dzu is D12), as in the user manual, Sec. 5.5.4;
% DB, 11/03/2026: Add type 'dynamic', the closed loop of the PIE and a
%                   dynamic output feedback controller with its own T operator
%                   (formerly closedLoopPIE_OF_sop);

if nargin==2
    type = 'controller';
elseif nargin==3
    % do nothing
end

if ~isa(PIE,'pie_struct')
    error('First input to closedLoopPIE() must be a pie_struct object');
end
PIE = initialize(PIE);
if size(PIE.K,1)>0
    error('Closed-loop PIEs are not supported for PIEs with an algebraic constraint (K, Kw, Ku).');
end

if ~any(strcmp(type,{'controller','observer','dynamic'}))    % DB, 11/03/2026
    error('Third input to closedLoopPIE() must be "controller", "observer" or "dynamic"');    % DB, 11/03/2026
end    % DB, 11/03/2026

% The gain, as a copvar
if strcmp(type,'dynamic')    % DB, 11/03/2026
    if ~isstruct(K)    % DB, 11/03/2026
        error('For type "dynamic", the second input to closedLoopPIE() should be a struct with fields TK, AK, BK, CK and DK');    % DB, 11/03/2026
    end    % DB, 11/03/2026
    if isfield(K,'G'),  K = K.G;    end    % DB, 11/03/2026
    for f = {'TK','AK','BK','CK','DK'}    % DB, 11/03/2026
        if ~isfield(K,f{1})    % DB, 11/03/2026
            error('For type "dynamic", the second input to closedLoopPIE() should have a field %s',f{1});    % DB, 11/03/2026
        end    % DB, 11/03/2026
        if isa(K.(f{1}),'opvar') || isa(K.(f{1}),'opvar2d')    % DB, 11/03/2026
            K.(f{1}) = op2copvar_sop(K.(f{1}));    % DB, 11/03/2026
        elseif ~isa(K.(f{1}),'copvar')    % DB, 11/03/2026
            error('The controller operators should be opvar, opvar2d or copvar class objects');    % DB, 11/03/2026
        end    % DB, 11/03/2026
    end    % DB, 11/03/2026
elseif isa(K,'opvar') || isa(K,'opvar2d')    % DB, 11/03/2026
    K = op2copvar_sop(K);
elseif ~isa(K,'copvar')
    error('Second input to closedLoopPIE() must be an opvar, opvar2d or copvar class object');
end

% Build controller closed-loop system by default
if strcmp(type,'controller') 
    if sum(K.dim_in)~=sum(PIE.space.x.dims)
        error('Input dimensions of controller gain should match size of the PIE state')
    elseif sum(K.dim_out)~=sum(PIE.space.u.dims)
        error('Output dimensions of controller gain should match size of the actuator input')
    end

    % (T + Tu*K) \dot x +Tw w = (A + Bu*K) x + Bw w, with the outputs
    % z = (Cz + Dzu*K) x + Dzw w, y = (Cy + Dyu*K) x + Dyw w, and no u
    PIEcl = PIE;
    PIEcl.T = PIE.T + PIE.Tu*K;     PIEcl.A = PIE.A + PIE.Bu*K;
    PIEcl.Cz = PIE.Cz + PIE.Dzu*K;  PIEcl.Cy = PIE.Cy + PIE.Dyu*K;
    PIEcl.Tu = [];  PIEcl.Bu = [];  PIEcl.Dzu = [];  PIEcl.Dyu = [];
    PIEcl.K = [];   PIEcl.Kw = [];  PIEcl.Ku = [];
    PIEcl.space.u = struct('dims',zeros(0,1),'depmat',false(0,numel(PIE.space.vars)));
elseif strcmp(type,'dynamic')    % DB, 11/03/2026
    % State [x;xK] of the plant and the controller, with u and y eliminated    % DB, 11/03/2026
    if ~(PIE.Dyu==0)    % DB, 11/03/2026
        error('The dynamic closed loop is implemented for D22 = 0');    % DB, 11/03/2026
    elseif ~(PIE.Tw==0) || ~(PIE.Tu==0)    % DB, 11/03/2026
        error('The dynamic closed loop is implemented for Tw = 0 and Tu = 0');    % DB, 11/03/2026
    end    % DB, 11/03/2026
    TK = K.TK;  AK = K.AK;  BK = K.BK;  CK = K.CK;  DK = K.DK;    % DB, 11/03/2026
    PIEcl = pie_struct();    % DB, 11/03/2026
    PIEcl.T  = blkdiag(PIE.T,TK);    % DB, 11/03/2026
    PIEcl.A  = [PIE.A+PIE.Bu*DK*PIE.Cy,     PIE.Bu*CK;    % DB, 11/03/2026
                BK*PIE.Cy,                  AK];    % DB, 11/03/2026
    PIEcl.Bw = [PIE.Bw+PIE.Bu*DK*PIE.Dyw;   BK*PIE.Dyw];    % DB, 11/03/2026
    PIEcl.Cz = [PIE.Cz+PIE.Dzu*DK*PIE.Cy,   PIE.Dzu*CK];    % DB, 11/03/2026
    PIEcl.Dzw = PIE.Dzw + PIE.Dzu*DK*PIE.Dyw;    % DB, 11/03/2026
else
    % Otherwise, build observer closed-loop system
    L = K;
    if sum(L.dim_in)~=sum(PIE.space.y.dims)
        error('Input dimensions of observer gain should match size of the observed output')
    elseif sum(L.dim_out)~=sum(PIE.space.x.dims)
        error('Output dimensions of observer gain should match size of the PIE state')
    end

    % State [x;xhat], with error dynamics for xhat, and the outputs z of x
    % and of xhat
    PIEcl = pie_struct();
    LCy = L*PIE.Cy;
    PIEcl.T = blkdiag(PIE.T,PIE.T);
    PIEcl.Tw = [PIE.Tw; 0*PIE.Tw];      PIEcl.Tu = [PIE.Tu; PIE.Tu];    % DB, 11/04/2026: xhat row of PIEcl.Tu should be Tu, not 0*Tu
    PIEcl.A = [PIE.A, 0*PIE.A; -LCy, PIE.A+LCy];
    PIEcl.Bw = [PIE.Bw; -L*PIE.Dyw];    PIEcl.Bu = [PIE.Bu; PIE.Bu];    % DB, 11/03/2026: sign of L*Dyw in PIEcl.Bw should be -L*Dyw; 
    % 11/04/2026: xhat row of PIEcl.Bu should be Bu, not Bu-L*Dyw
    PIEcl.Cz = blkdiag(PIE.Cz,PIE.Cz);
    PIEcl.Dzw = [PIE.Dzw; 0*PIE.Dzw];   PIEcl.Dzu = [PIE.Dzu; PIE.Dzu];    % DB, 11/04/2026
end

% The closed-loop PIE no longer corresponds to the original PDE, so the maps
% to the PDE are removed, and the PIE is initialized.
default_pie = pie_struct();
PIEcl.space_pde = default_pie.space_pde;
PIEcl.pie2pde = default_pie.pie2pde;
PIEcl.pde2pie = default_pie.pde2pie;
PIEcl = initialize(PIEcl);

end
