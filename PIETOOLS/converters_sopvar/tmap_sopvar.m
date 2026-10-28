function [Top, Pop_f2x, Pop_b2x, bc_state_tab, PDE] = tmap_sopvar(PDE)
% Compute 'copvar' operators Top.x, Top.u and Top.w such that the PDE state
% x can be expressed in terms of the fundamental state xf and the inputs as
%   x = Top.x * xf + Top.u * u + Top.w * w;
% where the boundary conditions of the PDE are satisfied by construction.
%
% INPUTS:
% - PDE:    An initialized 'pde_struct' with ordered state components, as
%           passed to the conversion routine after "reorder_comps".
%
% OUTPUTS:
% - Top:    Struct with 'copvar' fields x, u, w, mapping from xf, u and w to
%           the PDE state. The field u (w) is [] if there are no inputs u (w);
% - Pop_f2x, Pop_b2x, bc_state_tab:
%           The expansion of the PDE state in the fundamental and core
%           boundary state, see "ftc_expansion_sopvar";
% - PDE:    The PDE with BCs expanded in terms of the core boundary and
%           fundamental state;
%
% NOTES:
% - The core boundary state xb is expressed as a function of xf, u and w by
%   solving the BCs 0 = P_b2BC*xb + P_f2BC*xf + P_u2BC*u + P_w2BC*w. At
%   present this is only supported if the BCs and core boundary states are
%   finite-dimensional (e.g. 1D PDEs), in which case P_b2BC is a matrix. For
%   BCs that vary in space, a separate map along each variable is needed,
%   see the "Option 1" described in convert_PIETOOLS_PDE.
% - If the matrix P_b2BC is singular, as for generalized periodic boundary
%   conditions, the map from [v0;xf] to x is computed with additional moment
%   states v0, see "periodic_map_sopvar", and Top has the additional fields F,
%   K and periodic.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - tmap_sopvar
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
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding, 'copvar' version of the map from fundamental to PDE state
% in convert_PIETOOLS_PDE, 10/28/2026

vars = PDE.vars;
dom = PDE.dom;
nvars = size(vars,1);
x_tab = PDE.x_tab;

% % % STEP 1: x = Pop_b2x*xb + Pop_f2x*xf
[Pop_f2x, Pop_b2x, bc_state_tab] = ftc_expansion_sopvar(PDE);

% % % STEP 2: 0 = Pop_b2BC*xb + Pop_f2BC*xf + Pop_u2BC*u + Pop_w2BC*w
PDE = expand_PIETOOLS_PDE_BCs(PDE);
[Pop_b2BC, Pop_f2BC, Pop_u2BC, Pop_w2BC, layout] = bc_operators_sopvar(PDE,bc_state_tab);

% % % STEP 3: solve for xb, and set x = Top.x*xf + Top.u*u + Top.w*w.
no_u = isempty(Pop_u2BC) || is_zero_copvar(Pop_u2BC);
no_w = isempty(Pop_w2BC) || is_zero_copvar(Pop_w2BC);
Ginv = [];

if is_zero_copvar(Pop_f2BC) && no_u && no_w
    % The BCs impose no constraints that involve the core boundary state
    % (e.g. only conditions at the lower boundary), so that xb=0.
    Top.x = Pop_f2x;
else
    % The core boundary values are finite-dimensional if the BCs and core
    % boundary state are
    G = copvar2matrix(Pop_b2BC);
    if isempty(G)
        error(['The boundary conditions involve infinite-dimensional boundary values; ',...
                'conversion for such PDEs is not yet supported in the sopvar converter.'])
    end
    if size(G,1)<size(G,2)
        error('It appears insufficient boundary conditions have been specified; conversion to PIE is not supported.')
    elseif size(G,1)>size(G,2)
        error('It appears too many boundary conditions have been specified; either declare a higher order of differentiability of the PDE state, or reduce the number of boundary conditions.')
    end
    if abs(det(G))>eps
        % xb = -inv(G)*(P_f2BC*xf + P_u2BC*u + P_w2BC*w)
        entry = struct('r',0,'c',0,'rows',1:size(G,1),'cols',1:size(G,2), ...
                        'sub',ones(1,nvars),'C',polynomial(-pinv(G)));
        Ginv = assemble_copvar_blocks(entry,layout.b.masks,layout.b.dims,layout.BC.masks,layout.BC.dims,vars,dom);
        Top.x = Pop_b2x*(Ginv*Pop_f2BC) + Pop_f2x;
    elseif rank(G)<size(G,2)
        % The matrix is singular, as for generalized periodic boundary
        % conditions: introduce additional moment states.
        if any(x_tab(:,3)~=1)
            error('Currently supported generalized PBCs do not include coupled ODE states).')
        end
        for ii=1:numel(PDE.x)
            if isfield(PDE.x{ii},'tdiff') && PDE.x{ii}.tdiff~=1
                error('Currently supported generalized PBCs do not include higher-order temporal derivatives.')
            end
        end
        if ~no_u || ~no_w
            error('PIETOOLS:convert:SingularBCInputs', ...
                'Currently supported generalized PBCs do not include boundary inputs.')
        end
        [Tper,Fop,Kop,info] = periodic_map_sopvar(PDE,Pop_f2x,Pop_b2x,bc_state_tab,Pop_b2BC,Pop_f2BC);
        Top.x = Tper.x;
        Top.F = Fop;    Top.K = Kop;    Top.periodic = info;
    else
        error(['The PI operator defining the boundary conditions appears not to be invertible.',...
                ' Please make sure that your system is well-posed.'])
    end
end

Top.u = input_map(Pop_u2BC,Ginv,Pop_b2x,layout.u,layout.x,vars,dom);
Top.w = input_map(Pop_w2BC,Ginv,Pop_b2x,layout.w,layout.x,vars,dom);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function Pop = input_map(Pop_i2BC,Ginv,Pop_b2x,lay_i,lay_x,vars,dom)
% Map from an input to the PDE state: Pop_b2x*(Ginv*Pop_i2BC), where
% Ginv = -inv(P_b2BC).
if isempty(Pop_i2BC)
    Pop = [];
elseif isempty(Ginv) || is_zero_copvar(Pop_i2BC)
    empty_e = struct('r',{},'c',{},'rows',{},'cols',{},'sub',{},'C',{});
    Pop = assemble_copvar_blocks(empty_e,lay_x.masks,lay_x.dims,lay_i.masks,lay_i.dims,vars,dom);
else
    Pop = Pop_b2x*(Ginv*Pop_i2BC);
end
end

function tf = is_zero_copvar(P)
% True if all blocks of the container are structurally or numerically zero.
tf = true;
for k=1:numel(P.C)
    blk = P.C{k};
    if isempty(blk),    continue,   end
    for m=1:numel(blk.params)
        if any(any(blk.params{m}~=0)),  tf = false;     return,     end
    end
end
end

function G = copvar2matrix(P)
% Matrix representing a container of operators between finite-dimensional
% spaces; [] if any of the spaces is infinite-dimensional.
G = [];
if any(P.space_out(:)) || any(P.space_in(:)),   return,     end
rr = [0;cumsum(P.dim_out(:))];      cc = [0;cumsum(P.dim_in(:))];
G = zeros(rr(end),cc(end));
for i=1:numel(P.dim_out)
    for j=1:numel(P.dim_in)
        blk = P.C{i,j};
        if isempty(blk),    continue,   end
        G(rr(i)+1:rr(i+1),cc(j)+1:cc(j+1)) = full(blk.params{1});
    end
end
end
