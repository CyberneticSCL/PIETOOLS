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
% [0 T][\dot xhat]  [0 ]    [0 ]     [-K*C2 A+K*C2][xhat] +[K*D21]    [K*D22]
% 
% z = [C1 0 ][x   ] +[D11] w +[D12] u
%     [0  C1][xhat]  [0  ]    [0  ]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE: a linear PIE, as a pie_struct object
% K : an opvar, opvar2d, or copvar object, mapping the PIE state to the
%     actuator input u (controller) or the observed output y to the PIE
%     state (observer)
% type: optional string input; 'controller' (default) or 'observer' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIEcl: an initialized pie_struct object, with the operators Tcl Twcl Acl
%        B1cl C1cl D11cl C2cl D21cl, as copvar objects. For the controller,
%        the actuator input u is eliminated, so that Tu, B2, D12 and D22
%        are zero. For the observer, the state is [x;xhat], and the
%        observed output is eliminated

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

% The gain, as a copvar
if isa(K,'opvar') || isa(K,'opvar2d')
    K = op2copvar_sop(K);
elseif ~isa(K,'copvar')
    error('Second input to closedLoopPIE() must be an opvar, opvar2d or copvar class object');
end

if ~strcmp(type,'controller') && ~strcmp(type,'observer')
    error('Third input to closedLoopPIE() must be "controller" or "observer"');
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
    PIEcl.Tw = [PIE.Tw; 0*PIE.Tw];      PIEcl.Tu = [PIE.Tu; 0*PIE.Tu];
    PIEcl.A = [PIE.A, 0*PIE.A; -LCy, PIE.A+LCy];
    PIEcl.Bw = [PIE.Bw; L*PIE.Dyw];     PIEcl.Bu = [PIE.Bu; L*PIE.Dyu];
    PIEcl.Cz = blkdiag(PIE.Cz,PIE.Cz);
    PIEcl.Dzw = [PIE.Dzw; 0*PIE.Dzw];   PIEcl.Dzu = [PIE.Dzu; 0*PIE.Dzu];
end

% The closed-loop PIE no longer corresponds to the original PDE, so the maps
% to the PDE are removed, and the PIE is initialized.
default_pie = pie_struct();
PIEcl.space_pde = default_pie.space_pde;
PIEcl.pie2pde = default_pie.pie2pde;
PIEcl.pde2pie = default_pie.pde2pie;
PIEcl = initialize(PIEcl);

end
