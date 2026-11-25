function PIE = convert(PDE,out_type)
% PIE = convert(PDE,out_type)
% Convert a 'pde_struct' object PDE to an associated 'pie_struct' object
% PIE.
%
% INPUT:
% - PDE:        A pde_struct object defining a PDE in the terms format.
% - out_type:   char object specifying to what type of system the PDE is to
%               be converted. Only out_type='pie' is supported, as we only
%               convert PDEs to PIEs.
%               Alternatively, 'char' or cell of 'char' objects with the
%               flags 'silent' (suppress the display), 'nonlinear' (return
%               an 'npie_struct' also for a linear PDE) and/or 'Top'
%               (return only the operators T, Tw, Tu, K and the spaces);
%
% OUTPUT:
% - PIE:        A pie_struct object (or an npie_struct object, if the PDE
%               is nonlinear or the flag 'nonlinear' is specified) defining
%               the equivalent PIE representation of the input PDE, if this
%               representation exists.
%
% NOTES:
% - If no dynamics are specified for any of the states of the (linear) PDE,
%   it defines a boundary value problem (BVP), which is solved for the PDE
%   state. PIE is then a pie_struct object with only the
%   operators T, Tw, Tu and K, Kw, Ku specified, for which
%       x_pde = T*xf + Tw*wf + Tu*uf,       K*xf + Kw*wf + Ku*uf = 0,
%   for the fundamental state xf and inputs wf, uf. A PDE without dynamics
%   cannot have delays. Its static constraints 0=F(s) on the whole domain,
%   and other algebraic constraints, are rows of K.
% - Each state needs at least as many boundary conditions along each
%   direction as its order. A condition reducing to 0=0, such as
%   x(0)==x(0), is an error: states free along a direction are not
%   supported.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pde_struct/convert
%
% Copyright (C)2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% Initial coding DJ - 10/11/2022
% Call combines 1D-2D converter, DJ - 10/17/2024
% Display successful connversion message DB -12/22/24
% DJ, 06/01/2025: Move pre-processing of PDE to "convert_PIETOOLS_PDE";
% DJ, 11/02/2026: Accept the flags 'silent', 'nonlinear', 'Top' of
%                   "convert_PIETOOLS_PDE" as second input;
% DJ, 11/11/2026: Convert a PDE without dynamics (a BVP) to a PIE with
%                   only T and K operators, using "convert_PIETOOLS_BVP";
% DJ, 11/12/2026: Give an error for a PDE without dynamics that has static
%                   constraints;
% DJ (CC), 11/13/2026: Error for states sharing a moment, as in
%                   convert_PIETOOLS_PDE, and for delays in a PDE without
%                   dynamics; shorter error messages;
% MMP, 11/19/2026: Release comment pass (comment-only).
% DJ, 11/20/2026: Allow static constraints in a BVP (rows of K); set
%                   pie2pde.X, T, Tw, Tu; print a heading and summary;
% DJ, 11/25/2026: Order the components of xf of a BVP by their number of
%                   variables, as for a PDE with dynamics;

% We support only PDE to PIE conversion. The second argument is 'pie', or
% flags of the converter, with the 'pie' left out.
if nargin<2
    out_type = {};
end
flag = cellstr(out_type);
if ~all(ismember(flag,{'pie','silent','nonlinear','Top','return_Top'}))     % DJ, 11/02/2026
    error(['Only conversion from PDE to PIE is supported. The second argument must be ''pie'', ',...
        'or flags of convert_PIETOOLS_PDE (''silent'', ''nonlinear'', ''Top''), or be omitted.'])
end
flag = flag(~strcmp(flag,'pie'));

% Make sure the PDE is initialized
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end

% Check if any dynamics have been specified
is_BVP = ~isempty(PDE.x) && ~PDE.is_nonlinear && ~ismember('nonlinear',flag) &&...
            all(cellfun(@(xi) isfield(xi,'tdiff') && isnan(xi.tdiff),PDE.x(:)));

% Call the appropriate converter.
if is_BVP
    % Convert BVP to PIE
    PIE = BVP2pie_struct(PDE,ismember('silent',flag));
else
    % Convert PDE to PIE
    PIE = convert_PIETOOLS_PDE(PDE,flag);
end

% DB 12/22/24
if ~ismember('silent',flag)                                                 % DJ, 11/02/2026
    fprintf('\n --- Conversion to PIE was successful --- \n');
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function PIE = BVP2pie_struct(PDE,suppress_summary)
% PIE = BVP2PIE_STRUCT(PDE,SUPPRESS_SUMMARY) converts a PDE without
% dynamics, i.e. a boundary value problem, to a 'pie_struct' object with
% only the operators T, Tw, Tu and K, Kw, Ku specified, solving the BVP for
% the PDE state.
%
% INPUT
% - PDE:    initialized, linear 'pde_struct' object with nx states, with no
%           equation specified for any of them; its static constraints
%           (PDE.constraint_ID) become rows of K;
% - suppress_summary:   logical, true to not print the heading and how
%           the fundamental state and inputs relate to the PDE;
%
% OUTPUT
% - PIE:    initialized 'pie_struct' object, with T, Tw, Tu the nx x ncomp
%           block operators mapping the fundamental state xf (with ncomp
%           components) and inputs wf, uf to the PDE state, and K, Kw, Ku 
%           the constraint on them. Each
%           component of xf is derived from a single PDE state by 
%           PIE.pde2pie.x, and each component of wf and uf from a
%           PDE input by PIE.pde2pie.w and PIE.pde2pie.u;

if PDE.has_delay
    error('Delays are not supported in a PDE without dynamics.')
end
if ~suppress_summary
    fprintf('\n --- Solving the boundary value problem: parameterizing the PDE state by a fundamental state --- \n')
end
% Parameterize the states satisfying the BCs by the fundamental state xf
% and inputs wf, uf.
[Tops,Kops,Dmap,info] = convert_PIETOOLS_BVP(PDE);
% pde2pie describes each component of xf by a single PDE state, so split
% the states of coupled BCs that have been merged into one component.
if any(Dmap.x.idx==0)
    [Tops,Kops,Dmap,info] = convert_PIETOOLS_BVP(PDE,'blocks');
end
% Static constraints are rows of K, not outputs of the PDE in space_pde
PDE = static_constraints(PDE,'drop');
if any(Dmap.x.idx==0)                                                       % DJ (CC), 11/13/2026
    error("convert_PIETOOLS_PDE:blocks","Rank-deficient (e.g. periodic) boundary conditions that couple states in a shared moment are not supported.")
end
% List the components of xf by their number of variables, ascending (ODE
% states first), stable within each group, as for a PDE with dynamics. The
% columns of T and K, and the rows of the map to xf, are permuted.
[~,perm] = sort(sum(info.space.x.depmat,2)','ascend');                      % DJ, 11/25/2026
if ~isequal(perm,1:numel(perm))
    Tops.x = copvar(Tops.x.C(:,perm));
    if size(Kops.x,1)>0
        Kops.x = copvar(Kops.x.C(:,perm));
    end
    Dmap.x.idx = Dmap.x.idx(perm);
    Dmap.x.diff = Dmap.x.diff(perm,:,:);
    Dmap.x.C = copvar(Dmap.x.C.C(perm,:));
    info.space.x.dims = info.space.x.dims(perm);
    info.space.x.depmat = info.space.x.depmat(perm,:);
end
space_pie = info.space;
N = numel(space_pie.vars);     
nx = numel(PDE.x);      
ncomp = numel(space_pie.x.dims);

% % % Declare the operators and the function spaces of the PIE and the PDE
PIE = pie_struct();
PIE.T = Tops.x;     PIE.Tw = Tops.w;    PIE.Tu = Tops.u;
PIE.K = Kops.x;     PIE.Kw = Kops.w;    PIE.Ku = Kops.u;
PIE.space.vars = space_pie.vars;               
PIE.space.dom = space_pie.dom;
for ty = {'x','w','u'}
    PIE.space.(ty{1}) = space_pie.(ty{1});
end
space_pde = struct('vars',{space_pie.vars},'dom',space_pie.dom);
for ty = {'x','w','u','y','z'}
    tab = PDE.([ty{1},'_tab']);
    ID = cellfun(@(c) c.ID,PDE.(ty{1}));
    space_pde.(ty{1}) = struct('dims',tab(:,2),'depmat',logical(tab(:,3:2+N)),'ID',ID(:));
end
PIE.space_pde = space_pde;

% Spaces (structs with fields vars, dom, dim) of the PDE states, and the
% identity operators on them.
get_space = @(dims,dep,i) struct('vars',{space_pie.vars(dep(i,:))},'dom',space_pie.dom(dep(i,:),:),'dim',dims(i));
identity_op = @(s) mat2sopvar(eye(s.dim),struct('in',{s.vars},'out',{s.vars}),struct('in',s.dom,'out',s.dom));
x_sp = arrayfun(@(k) get_space(PDE.x_tab(:,2),logical(PDE.x_tab(:,3:2+N)),k),1:nx);

% % % Maps from the PDE to the PIE: component c of xf is C_c applied to
% % % the derivative D^{diff(c,:)} of the PDE state idx(c), and wf, uf
% % % consist of derivatives of the PDE inputs, as given by Dmap. The
% % % operator C of pde2pie is block diagonal,
% % % with the operators C_c, and the other blocks [] (zero, their spaces
% % % following from the diagonal, so that no ncomp^2 zero operators are
% % % constructed)
pde2pie = PIE.pde2pie;
idx_x = Dmap.x.idx(:);
pde2pie.x.idx = [idx_x,zeros(ncomp,2)];
pde2pie.x.diff = zeros(ncomp,N);
for c=1:ncomp
    pde2pie.x.diff(c,:) = reshape(Dmap.x.diff(c,idx_x(c),:),1,N);
end
pde2pie.x.tdiff = zeros(ncomp,1);       pde2pie.x.loc = zeros(ncomp,N);
pde2pie.x.delay = zeros(ncomp,1);
% The block-diagonal operator C, block c the C_c of the PDE state idx(c)
C_blocks = cell(ncomp);
for c=1:ncomp
    C_blocks{c,c} = Dmap.x.C.C{c,idx_x(c)};
end
pde2pie.x.C = copvar(C_blocks);
% The inputs: the PDE inputs, their derivatives, positions and operators
for ty = {'w','u'}
    ninp = numel(Dmap.(ty{1}).idx);
    pde2pie.(ty{1}).idx = Dmap.(ty{1}).idx(:);
    pde2pie.(ty{1}).diff = Dmap.(ty{1}).diff;
    pde2pie.(ty{1}).tdiff = zeros(ninp,1);
    pde2pie.(ty{1}).loc = Dmap.(ty{1}).loc;
    if ninp>0
        pde2pie.(ty{1}).C = Dmap.(ty{1}).C;
    end
end
PIE.pde2pie = pde2pie;

% % % Map from the PIE to the PDE state: x_pde = X_pde*(T*xf+Tw*wf+Tu*uf),
% % % with the rows of T the PDE states, so that X_pde is the identity,
% % % and pie2pde.T, Tw, Tu are T, Tw, Tu
pie2pde = PIE.pie2pde;
X_blocks = cell(nx);
for k=1:nx
    X_blocks{k,k} = identity_op(x_sp(k));
end
pie2pde.X = copvar(X_blocks);
pie2pde.T = Tops.x;     pie2pde.Tw = Tops.w;    pie2pde.Tu = Tops.u;
PIE.pie2pde = pie2pde;

PIE = initialize(PIE);

% % % Summary of the fundamental state and inputs in terms of the PDE; the
% % % auxiliary components (on fewer variables than their state, e.g. means
% % % for periodic BCs) are shown as F*x, F their block of pde2pie.x.C
if ~suppress_summary
    is_aux = arrayfun(@(c) ~isequal(PIE.space.x.depmat(c,:),space_pde.x.depmat(idx_x(c),:)),1:ncomp);
    for ty = {'x','w','u'}
        if isempty(PIE.space.(ty{1}).dims)
            continue
        elseif strcmp(ty{1},'x') && any(is_aux)
            aux_idcs = find(is_aux);
            print_convert_summary(PIE,'x',aux_idcs,arrayfun(@(c) C_blocks{c,c},aux_idcs,'UniformOutput',false));
            fprintf('\nThe PIE state must satisfy PIE.K*xf = 0.\n');
            fprintf('The states are defined by PIE.pde2pie.x.C.\n');
        else
            print_convert_summary(PIE,ty{1});
        end
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %