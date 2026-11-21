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
% [0 T][\dot xhat]  [0 ]    [Tu ]     [-K*C2 A+K*C2][xhat] +[-K*D21]w+   [B2   ]u
% 
% z = [C1 0 ][x   ] +[D11] w +[D12] u
%     [0  C1][xhat]  [0  ]    [D12]
%
% if type='dynamic', the dynamic output feedback controller
% TK \dot xK = AK xK + BK y,    u = CK xK + DK y
% for D22 = 0, Tw = 0 and Tu = 0, so that
% [T 0 ][\dot x ] = [A + B2*DK*C2   B2*CK][x ] + [B1 + B2*DK*D21] w
% [0 TK][\dot xK]   [BK*C2          AK   ][xK]   [BK*D21         ]
% z = [C1 + D12*DK*C2,  D12*CK][x; xK] + (D11 + D12*DK*D21) w
%
% PIEs with an algebraic constraint (K, Kw, Ku) are not supported.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE: a linear PIE, as a pie_struct object
% K : a copvar object, mapping the PIE state to the
%     actuator input u (controller) or the observed output y to the PIE
%     state (observer); for type='dynamic', a struct with fields TK, AK,
%     BK, CK, DK, the operators (copvar objects) of the
%     controller PIE, or a struct with such a field G
% type: optional string input; 'controller' (default), 'observer' or
%     'dynamic'
% If Dyu=Tw=Tu=0, then given a state feedback controller K and an estimator gain L,
% the following call gives the closed-loop estimator-based control
% Kc.TK = PIE.T;
% Kc.AK = PIE.A + PIE.Bu*K + L*PIE.Cy;
% Kc.BK = -L;
% Kc.CK = K;
% Kc.DK = 0*(K*L);
% PIEcl = closedLoopPIE(PIE,Kc,'dynamic');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIEcl: an initialized pie_struct object, with the operators Tcl Twcl Acl
%        B1cl C1cl D11cl C2cl D21cl, as copvar objects. For the controller,
%        the actuator input u is eliminated, so that Tu, B2, D12 and D22
%        are zero. For the observer, the state is [x;xhat], and the
%        observed output is eliminated. For the dynamic controller, the
%        state is [x;xK], and both u and y are eliminated
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - closedLoopPIE
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding SS - 5/30/2021

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

if ~any(strcmp(type,{'controller','observer','dynamic'}))
    error('Third input to closedLoopPIE() must be "controller", "observer" or "dynamic"');
end

% The gain, as a copvar
if strcmp(type,'dynamic')
    if ~isstruct(K)
        error('For type "dynamic", the second input to closedLoopPIE() should be a struct with fields TK, AK, BK, CK and DK');
    end
    if isfield(K,'G'),  K = K.G;    end
    for f = {'TK','AK','BK','CK','DK'}
        if ~isfield(K,f{1})
            error('For type "dynamic", the second input to closedLoopPIE() should have a field %s',f{1});
        end
        if ~isa(K.(f{1}),'copvar')
            error('The controller operators should be copvar class objects');
        end
    end
elseif ~isa(K,'copvar')
    error('Second input to closedLoopPIE() must be a copvar class object');
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
elseif strcmp(type,'dynamic')
    % State [x;xK] of the plant and the controller, with u and y eliminated
    if ~(PIE.Dyu==0)
        error('The dynamic closed loop is implemented for D22 = 0');
    elseif ~(PIE.Tw==0) || ~(PIE.Tu==0)
        error('The dynamic closed loop is implemented for Tw = 0 and Tu = 0');
    end
    TK = K.TK;  AK = K.AK;  BK = K.BK;  CK = K.CK;  DK = K.DK;
    PIEcl = pie_struct();
    PIEcl.T  = blkdiag(PIE.T,TK);
    PIEcl.A  = [PIE.A+PIE.Bu*DK*PIE.Cy,     PIE.Bu*CK;
                BK*PIE.Cy,                  AK];
    PIEcl.Bw = [PIE.Bw+PIE.Bu*DK*PIE.Dyw;   BK*PIE.Dyw];
    PIEcl.Cz = [PIE.Cz+PIE.Dzu*DK*PIE.Cy,   PIE.Dzu*CK];
    PIEcl.Dzw = PIE.Dzw + PIE.Dzu*DK*PIE.Dyw;
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
    PIEcl.Tw = [PIE.Tw; 0*PIE.Tw];      PIEcl.Tu = [PIE.Tu; PIE.Tu];
    PIEcl.A = [PIE.A, 0*PIE.A; -LCy, PIE.A+LCy];
    PIEcl.Bw = [PIE.Bw; -L*PIE.Dyw];    PIEcl.Bu = [PIE.Bu; PIE.Bu];
    PIEcl.Cz = blkdiag(PIE.Cz,PIE.Cz);
    PIEcl.Dzw = [PIE.Dzw; 0*PIE.Dzw];   PIEcl.Dzu = [PIE.Dzu; PIE.Dzu];
end

% The closed-loop PIE no longer corresponds to the original PDE, so the maps
% to the PDE are removed, and the PIE is initialized.
default_pie = pie_struct();
PIEcl.space_pde = default_pie.space_pde;
PIEcl.pie2pde = default_pie.pie2pde;
PIEcl.pde2pie = default_pie.pde2pie;
PIEcl = initialize(PIEcl);

end
