function PIE = convert_PIETOOLS_PDE(PDE,comp_order,flag)
% PIE = CONVERT_PIETOOLS_PDE(PDE,FLAG) converts a PDE in the terms format,
% possibly coupling ODE states with PDE states on several spatial
% variables, and possibly nonlinear, to an equivalent PIE
%   d/dt (T*xf + Tw*wf + Tu*uf) = A*xf + Bw*wf + Bu*uf;
%                             z = Cz*xf + Dzw*wf + Dzu*uf;
%                             y = Cy*xf + Dyw*wf + Dyu*uf;
%                          K*xf = 0,
% with fundamental state xf and fundamental inputs wf, uf. For a
% nonlinear PDE, A*xf, Cz*xf and Cy*xf are polynomials f, gz and gy.
%
% INPUT
% - PDE:    'pde_struct' object in the terms format. A PDE in a legacy
%           format (batch, or 'terms_legacy') is passed on to the legacy
%           converters;
% - flag:   (optional) char or cellstr, with 'silent' to suppress the
%           summary of how the PIE relates to the PDE (and the warning
%           for more BCs than the states require), 'Top' to stop
%           after computing T, Tw, Tu and K, and/or 'nonlinear' to return
%           an 'npie_struct' also for a linear PDE. For
%           compatibility, FLAG may be the third input, after an empty
%           (not supported) argument comp_order;
%
% OUTPUT
% - PIE:    'pie_struct' object with 'copvar' operators for a linear PDE,
%           'npie_struct' object with 'tensopmat' operators if the PDE is
%           nonlinear or the flag 'nonlinear' is set. The components of
%           xf are the blocks of the solution of the BCs (one for each
%           state and component of the fundamental state), ordered by
%           their number of variables, and those of wf and uf the
%           fundamental inputs of that solution. The fields space,
%           space_pde, pde2pie and pie2pde describe the function spaces of
%           the PIE and PDE, how xf, wf and uf are derived from the PDE,
%           and how the PDE components are obtained from the PIE. With
%           'Top', PIE only has the fields space, T, Tw, Tu and K;
%
% NOTES
% - Each PDE state has one main block of xf, its highest derivative in
%   each variable. The other (auxiliary) blocks xf_b are moments F*x_pde
%   (e.g. for BCs that are periodic in a variable): their rows of T are
%   F times the row of x_pde (the selector of xf_b if F involves no
%   derivative), their right-hand sides F times that of x_pde, and
%   K*xf=0 ties them to x_pde.
% - The inputs in the BCs are decomposed into fundamental inputs
%   (derivatives of the PDE inputs); the other inputs are their own
%   fundamental inputs. An input in a BC can be used elsewhere in the PDE
%   only if it is its own fundamental input (e.g. no spatial variables);
%   pie2pde.W, pie2pde.U only map the inputs that are not decomposed.
%   Inputs can only appear linearly.
% - Delays are modeled by transport states, and higher-order temporal
%   derivatives by additional states, described by pde2pie.x.delay and
%   pde2pie.x.tdiff. This is not supported for states with inputs in
%   their BCs.
% - Static constraints 0=E(s) (outputs listed in PDE.constraint_ID),
%   linear in the state, are added to K. ODE states and states of lower
%   dimension may enter the BCs of other states.
% - Conditions that are not a BC along a single direction (e.g.
%   0=x2(s)-x1(1), or 0=x1-2*x2+0.5*w for ODE states) are algebraic
%   constraints K*xf + Kw*wf + Ku*uf = 0; inputs in them are supported for
%   a linear PIE only.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - convert_PIETOOLS_PDE
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/02/2026: Initial coding (ND version);
% DJ, 11/06/2026: Display PIE states and inputs with a hat in the summary;
% DJ, 11/06/2026: Get the names of components and variables in the summary
%                   from 'varname2display';
% DJ, 11/08/2026: Drop the zero maps of the inputs by info.input_maps;
% DJ, 11/11/2026: Rename solve_PIETOOLS_BVP to convert_PIETOOLS_BVP;
% DJ (CC), 11/13/2026: Constraints and hotd for an 'npie_struct'; summary
%                   from space, space_pde, pde2pie; Iw, Iu from info; warn
%                   for surplus BCs; inline permute_cols; shorter errors;
% DJ, 11/15/2026: Bugfix in case row of f, gz or gy has no terms;
% DB, 11/18/2026: Move print_convert_summary to its own file;
% MMP, 11/19/2026: Release comment pass (comment-only).
% DJ, 11/20/2026: Take all constraints (K, Kw, Ku) from
%                   convert_PIETOOLS_BVP; new pie2pde fields;

% % % Pass PDEs in the legacy formats to the legacy converters, and
% % % convert other structures to pde_struct
if isfield(PDE,'n0') || isfield(PDE,'n1') || isfield(PDE,'n2')
    PIE = convert_PIETOOLS_PDE_batch(PDE);
    return
elseif isfield(PDE,'n') || isfield(PDE,'ODE') || isfield(PDE,'PDE')
    PIE = convert_PIETOOLS_PDE_terms_legacy(PDE);
    return
elseif isstruct(PDE) && isfield(PDE,'x')
    PDE = pde_struct(PDE);
elseif isa(PDE,'sys')
    PDE = PDE.params;
end

% % % Check the inputs
if nargin==1
    flag = {};
elseif nargin==2
    flag = comp_order;
elseif nargin==3 && ~isempty(comp_order)
    error("The order of the components in the PIE cannot be specified; the components of the PDE are not reordered.")
end
flag = cellstr(flag);
suppress_summary = ismember('silent',flag);
return_Top = ismember('Top',flag) || ismember('return_Top',flag);
% Flag 'nonlinear': return an 'npie_struct' also for a linear PDE
force_nonlinear = ismember('nonlinear',flag);
if ~isa(PDE,'pde_struct')
    error("The PDE must be specified as a 'pde_struct' object.")
end
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end
% DB 11/07/2026: outputs flagged as static constraints 0=E(s)
% are not outputs of the PIE; their operator is added to K.
if ~isempty(PDE.constraint_ID)                                              % DJ, 11/20/2026
    z_IDs = cellfun(@(zi) zi.ID,PDE.z);
    is_con = ismember(z_IDs(:),PDE.constraint_ID(:));
    if any(cellfun(@(zi) any(cellfun(@numel,zi.term)>1),PDE.z(is_con)))    % DJ (CC), 11/13/2026
        error("convert_PIETOOLS_PDE:nonlinearConstraint","Static constraints 0=F(s) must be linear in the state.")
    elseif any(cellfun(@(zi) any(cellfun(@(tm) any(isfield(tm,{'w','u'})),zi.term)),PDE.z(is_con)))
        error("Constraints 0=F(s) cannot involve the inputs w or u.")
    end
end
% Every state component needs a differential equation
x_noeq = find(cellfun(@(xi) isfield(xi,'tdiff') && isnan(xi.tdiff),PDE.x));
if ~isempty(x_noeq)
    error(['Conversion to PIE is not supported: no differential equation has been specified for state component(s) x{',strjoin(arrayfun(@num2str,x_noeq(:)','UniformOutput',false),'}, x{'),'}.',...
            ' A PDE without dynamics describes a boundary value problem, which cannot be converted to a PIE.'])
end
% Delays are modeled by transport states on additional variables, by
% expand_delays; the maps between the PDE and PIE refer to the original
% PDE0
PDE0 = PDE;     del_tab = [];
if PDE.has_delay
    if ~suppress_summary
        fprintf(['\n','Delayed states or inputs were encountered:\n --- Running "expand_delays" to remove the delays --- \n']);
    end
    [PDE,del_tab] = expand_delays(PDE);
end
% Higher-order temporal derivatives: the PIE is first derived with these,
% and expanded to first order afterwards (for a 'pie_struct' or an
% 'npie_struct')
tdiff_pde = ones(numel(PDE.x),1);
for i=1:numel(PDE.x)
    if isfield(PDE.x{i},'tdiff')
        tdiff_pde(i) = PDE.x{i}.tdiff;
    end
end
has_hotd = any(tdiff_pde>1);
is_nonlinear = PDE.is_nonlinear || force_nonlinear;
% Algebraic states are not supported in a nonlinear PIE, nor with higher-
% order temporal derivatives
if ~isempty(x_noeq) && ~return_Top && (has_hotd || is_nonlinear)
    error("convert_PIETOOLS_PDE:algebraicStates","States without dynamics are not yet supported in a nonlinear PIE, nor with higher-order temporal derivatives.")
end
N = size(PDE.vars,1);
nx = numel(PDE.x);
vars_pde = arrayfun(@(i) PDE.vars(i,1).varname{1},1:N,'UniformOutput',false);
if ~suppress_summary
    fprintf('\n --- Converting the PDE to an equivalent PIE --- \n')
end

% Function spaces are structs with fields vars, dom, dim. The anonymous
% functions give the space of a component, identity and zero operators, and
% the matrix of operators X(r,j) where hit(r,j) holds (zero elsewhere)
get_space = @(obj_tab,i) struct('vars',{vars_pde(logical(obj_tab(i,3:2+N)))},'dom',PDE.dom(logical(obj_tab(i,3:2+N)),:),'dim',obj_tab(i,2));
identity_op = @(s) mat2sopvar(eye(s.dim),struct('in',{s.vars},'out',{s.vars}),struct('in',s.dom,'out',s.dom));
zero_op = @(row_sp,col_sp) mat2sopvar(zeros(row_sp.dim,col_sp.dim),struct('in',{col_sp.vars},'out',{row_sp.vars}),struct('in',col_sp.dom,'out',row_sp.dom));
% The operators are 'copvar' objects (linear PDE) or 'tensopmat' objects
% (nonlinear PDE); get_blocks extracts the cell of 'sopvar' blocks of an
% operator.
if is_nonlinear
    cell2op = @tensopmat;
    get_blocks = @(P) cellfun(@(o) o.ops{1},P.ops,'UniformOutput',false);
else
    cell2op = @copvar;
    get_blocks = @(P) P.C;
end
select_cell = @(c,h) c{1+h};
op_matrix = @(row_sp,col_sp,hit,X) cell2op(arrayfun(@(r,j) select_cell({zero_op(row_sp(r),col_sp(j)),X(r,j)},hit(r,j)),...
            repmat((1:numel(row_sp))',1,numel(col_sp)),repmat(1:numel(col_sp),numel(row_sp),1),'UniformOutput',false));


% % % Derive the map from the fundamental state and inputs to the PDE
% % % state, with a block for each state
opts = {'all_inputs','blocks'};
% 'tensopmat' blocks for a nonlinear PDE, or if the flag 'nonlinear' is set
if is_nonlinear
    opts{end+1} = 'nonlinear';
end
% The derivatives of the states needed are those in the terms of the PDE
diffs = 'terms';
if return_Top
    diffs = [];
end
[Tops,Kops,Dmap,info] = convert_PIETOOLS_BVP(PDE,diffs,opts);               % DJ, 11/11/2026
% The static constraints are rows of Kops.x, not outputs of the PIE
PDE = static_constraints(PDE,'drop');                                       % DJ, 11/20/2026
PDE0 = static_constraints(PDE0,'drop');
if info.surplus && ~suppress_summary
    fprintf('\n Warning: There are more boundary conditions than the orders of the states require; the PDE may be overdetermined.\n')
end
blocks = info.blocks;       nblocks = numel(blocks);
spaces_pde = arrayfun(@(k) get_space(PDE.x_tab,k),1:nx);

% Each PDE state k has one main block (derived from this state alone, on
% all its variables); the other blocks are auxiliary states
% xf_b=F*x_pde, with x_idcs(b) the state of block b and F_cell{b} the
% block of F
if any(Dmap.x.idx==0)
    error("convert_PIETOOLS_PDE:blocks","Rank-deficient (e.g. periodic) boundary conditions that couple states in a shared moment are not supported.")
end
% List the blocks by the number of variables they depend on, ascending:
% the finite-dimensional blocks (ODE states, moments) first, so that e.g.
% an ODE state precedes the PDE states in the PIE,
% then the blocks depending on one variable, then on two, etc. The sort is
% stable, so that the order of the PDE states is kept within each group.
% The blocks are permuted in Tops, Kops, Dmap and info.space.
nvars_blk = arrayfun(@(b) numel(b.vars),blocks);
[~,perm] = sort(nvars_blk(:)','ascend');
if ~isequal(perm,1:nblocks)
    blocks = blocks(perm);
    Dmap.x.idx = Dmap.x.idx(perm);
    Dmap.x.diff = Dmap.x.diff(perm,:,:);
    for fn=fieldnames(info.space.x)'
        if size(info.space.x.(fn{1}),1)==nblocks
            info.space.x.(fn{1}) = info.space.x.(fn{1})(perm,:);
        end
    end
    % The block columns of the operators on xf: Tops.x is one operator
    % for 'Top', and a cell of operators for each state otherwise
    perm_cols = @(B) B(:,perm);
    permute_op = @(P) cell2op(perm_cols(get_blocks(P)));
    if return_Top
        Tops.x = permute_op(Tops.x);
    else
        for k=1:nx
            for l=find(~cellfun(@isempty,Tops.x{k}(:)))'
                Tops.x{k}{l} = permute_op(Tops.x{k}{l});
            end
        end
    end
    if size(Kops.x,1)>0
        Kops.x = permute_op(Kops.x);
    end
end
x_idcs = Dmap.x.idx(:)';
% Find the main block of each PDE state
main_idcs = zeros(nx,1);
for k=1:nx
    cand_idcs = find(x_idcs==k & arrayfun(@(b) isempty(setxor(b.vars,spaces_pde(k).vars)),blocks));
    if numel(cand_idcs)~=1
        error("convert_PIETOOLS_PDE:mainBlock","Could not identify a unique block of the fundamental state representing state x{%d}.",k)
    end
    main_idcs(k) = cand_idcs;
end
is_aux = true(nblocks,1);     is_aux(main_idcs) = false;
% The rows of T and A: the blocks of the states with dynamics
dyn_rows = find(~isnan(tdiff_pde(x_idcs(:))))';                             % DJ, 11/20/2026
F_blocks = get_blocks(Dmap.x.C);
F_blocks = F_blocks(perm,:);
F_cell =arrayfun(@(b) F_blocks{b,x_idcs(b)},1:nblocks,'UniformOutput',false);


% % % Collect the blocks, as 'sopvar' objects ([] for zero), of the maps to
% % % the state (T), and of the right-hand sides of the equations (x), the
% % % outputs z and y, from the columns (blocks of xf, wf, uf) to rows
col_spaces = struct('x',blocks,'w',info.inputs.w,'u',info.inputs.u);
op_blocks = struct('T',struct('x',[],'w',[],'u',[]),'x',struct('x',[],'w',[],'u',[]),...
           'z',struct('x',[],'w',[],'u',[]),'y',struct('x',[],'w',[],'u',[]));
poly_struct = struct('x',{cell(nx,1)},'z',{cell(numel(PDE.z),1)},'y',{cell(numel(PDE.y),1)});
% The map to x_pde (the derivative of order 0), as a row of blocks for each
% state if the derivatives were computed
for obj_type={'x','w','u'}
    % w, u: no maps if there are no inputs in the BCs
    if ~strcmp(obj_type{1},'x') && ~info.input_maps.(obj_type{1})           % DJ, 11/08/2026
        Tops.(obj_type{1}) = [];
    end
    if isempty(Tops.(obj_type{1}))
        continue
    end
    if return_Top
        op_blocks.T.(obj_type{1}) = get_blocks(Tops.(obj_type{1}));
    else
        T_rows = cellfun(@(c) get_blocks(c{1}),Tops.(obj_type{1}),'UniformOutput',false);
        op_blocks.T.(obj_type{1}) = vertcat(T_rows{:});
    end
end
% The maps to the derivatives, and from wf to the PDE inputs, for the terms
Tops_term = Tops;
Tops_term.Iw = info.Iw;     Tops_term.Iu = info.Iu;
% Sum the terms of each equation: the contributions of the state of a
% nonlinear PDE in poly_struct (polyopvar), the others as operators
for obj={'x','z','y'}
    if return_Top
        break
    end
    obj_tab = PDE.([obj{1},'_tab']);
    for obj_type={'x','w','u'}
        op_blocks.(obj{1}).(obj_type{1}) = cell(numel(PDE.(obj{1})),numel(col_spaces.(obj_type{1})));
    end
    % Convert each term of equation i to blocks acting on xf, wf and uf
    for i=1:numel(PDE.(obj{1}))
        vars_eq = vars_pde(logical(obj_tab(i,3:2+N)));
        for j=1:numel(PDE.(obj{1}){i}.term)
            [Pblk,Px] = convert_PDE_term(PDE.(obj{1}){i}.term{j},Tops_term,PDE,struct('nonlinear',is_nonlinear,'outvars',{vars_eq}));
            if ~isempty(Px)
                if isempty(poly_struct.(obj{1}){i})
                    poly_struct.(obj{1}){i} = Px;
                else
                    poly_struct.(obj{1}){i} = poly_struct.(obj{1}){i}+Px;
                end
            end
            for obj_type={'x','w','u'}
                if isempty(Pblk.(obj_type{1})) || (is_nonlinear && strcmp(obj_type{1},'x'))
                    continue
                end
                % Add the blocks of this term to those of the equation
                for k=find(~cellfun(@isempty,Pblk.(obj_type{1}).B))
                    if isempty(op_blocks.(obj{1}).(obj_type{1}){i,k})
                        op_blocks.(obj{1}).(obj_type{1}){i,k} = Pblk.(obj_type{1}).B{k};
                    else
                        op_blocks.(obj{1}).(obj_type{1}){i,k} = op_blocks.(obj{1}).(obj_type{1}){i,k}+Pblk.(obj_type{1}).B{k};
                    end
                end
            end
        end
    end
end


% % % Assemble the operators of the PIE, with a row for each block of the
% % % state, for the maps to and the right-hand side of the state equation:
% % % main block of state k: row k of the PDE; auxiliary block b: F_b*(row
% % % k of the PDE), for the state k that b is derived from, since
% % % xf_b=F_b*x_pde; for T also, unless it is the selector of xf_b:
% % % T*xf=[xf_b;x_pde].
if is_nonlinear
    PIE = npie_struct();
else
    PIE = pie_struct();
end
% The function spaces of the PIE: those of the fundamental state and
% inputs, and of the outputs of the PDE
space = info.space;
space.y = struct('dims',PDE.y_tab(:,2),'depmat',logical(PDE.y_tab(:,3:2+N)));
space.z = struct('dims',PDE.z_tab(:,2),'depmat',logical(PDE.z_tab(:,3:2+N)));
PIE.space = space;
row_spaces = struct('T',blocks,'x',blocks,'z',arrayfun(@(i) get_space(PDE.z_tab,i),1:numel(PDE.z)),...
             'y',arrayfun(@(i) get_space(PDE.y_tab,i),1:numel(PDE.y)));
op_names = struct('T',{{'T','Tw','Tu'}},'x',{{'A','Bw','Bu'}},'z',{{'Cz','Dzw','Dzu'}},'y',{{'Cy','Dyw','Dyu'}});
if is_nonlinear
    op_names.x{1} = '';
end   % f is set below
xwu_cell = {'x','w','u'};
% Top_pde.(x/w/u): maps from xf, wf, uf to all PDE states (pie2pde.T, Tw, Tu)
Top_pde = struct('x',[],'w',[],'u',[]);
% Declare the operators of the PIE, block by block
for obj={'T','x','z','y'}
    is_state_obj = any(strcmp(obj{1},{'T','x'}));
    for xwu_idx=1:3
        blocks_old = op_blocks.(obj{1}).(xwu_cell{xwu_idx});
        if isempty(op_names.(obj{1}){xwu_idx}) || isempty(blocks_old) || all(cellfun(@isempty,blocks_old(:)))
            continue
        end
        blocks_new = cell(numel(row_spaces.(obj{1})),size(blocks_old,2));
        % Block (ridx,cidx): the PDE row, or F times it for auxiliary ones
        for ridx=1:size(blocks_new,1)
            for cidx=1:size(blocks_new,2)
                block_op = [];
                if ~is_state_obj
                    block_op = blocks_old{ridx,cidx};
                elseif ~is_aux(ridx)
                    block_op = blocks_old{x_idcs(ridx),cidx};
                elseif strcmp(obj{1},'T') && xwu_idx==1 && ~any(blocks(ridx).Dorder)
                    if cidx==ridx
                        block_op = identity_op(blocks(ridx));
                    end   % selector
                elseif ~isempty(blocks_old{x_idcs(ridx),cidx})
                    block_op = F_cell{ridx}*blocks_old{x_idcs(ridx),cidx};
                end
                if isempty(block_op)
                    block_op = zero_op(row_spaces.(obj{1})(ridx),col_spaces.(xwu_cell{xwu_idx})(cidx));
                end
                blocks_new{ridx,cidx} = block_op;
            end
        end
        % T, A: only the rows of the states with dynamics; the main rows of
        % all PDE states of T, Tw, Tu give Top_pde
        rows_kept = 1:size(blocks_new,1);
        if is_state_obj
            rows_kept = dyn_rows;
            if strcmp(obj{1},'T')
                Top_pde.(xwu_cell{xwu_idx}) = cell2op(blocks_new(main_idcs(1:numel(PDE0.x)),:));
            end
        end
        try
            PIE.(op_names.(obj{1}){xwu_idx}) = cell2op(blocks_new(rows_kept,:));   % DJ, 11/20/2026
        catch ME
            error("convert_PIETOOLS_PDE:operator","Could not construct the operator %s: %s",op_names.(obj{1}){xwu_idx},ME.message)
        end
    end
end
% Constraints K*xf=0 on the fundamental state, for the blocks of xf
if size(Kops.x,1)>0
    PIE.K = Kops.x;
    % Constraints involving inputs: K*xf + Kw*wf + Ku*uf = 0
    if info.K_inputs
        PIE.Kw = Kops.w;    PIE.Ku = Kops.u;
    end
end
if return_Top
    return
end
% Polynomials of a nonlinear PIE, with the F-weighted rows for the
% auxiliary states of f
if is_nonlinear
    for obj={'x','z','y'}
        obj_tab = PDE.([obj{1},'_tab']);
        for i=find(cellfun(@isempty,poly_struct.(obj{1}))')
            poly_struct.(obj{1}){i} = constant2polyopvar(zeros(obj_tab(i,2),1));
            % Make sure the row maps to the correct space
            dep_i = logical(obj_tab(i,3:2+N));
            if any(dep_i)                                                   % DJ, 11/15/2026
                poly_struct.(obj{1}){i} = lift_output(poly_struct.(obj{1}){i},...
                    vars_pde(dep_i),PDE.dom(dep_i,:));
            end
        end
    end
    f_rows = poly_struct.x(x_idcs);
    for b=find(is_aux)'
        % F*f is not supported for f of degree>1: periodic BCs with a
        % nonlinear right-hand side error with a clear message
        try                                                                 % DJ (CC), 11/02/2026
            f_rows{b} = tensopmat(F_cell{b})*poly_struct.x{x_idcs(b)};
        catch
            error("convert_PIETOOLS_PDE:nonlinearAux","Rank-deficient (e.g. periodic) boundary conditions of a state with a nonlinear equation are not supported.")
        end
    end
    PIE.f = vertcat(f_rows{:});
    if ~isempty(poly_struct.z)
        PIE.gz = vertcat(poly_struct.z{:});
    end
    if ~isempty(poly_struct.y)
        PIE.gy = vertcat(poly_struct.y{:});
    end
end



% % % Function spaces of the PDE (before expanding delays), in its own
% % % variables; dep_pde is their dependence on the variables of the PIE
% (if the variables were merged and renamed, the k-th PDE variable is the
% k-th PIE variable)
nvars0 = size(PDE0.vars,1);     nx0 = numel(PDE0.x);        nw0 = numel(PDE0.w);
vars_pde0 = arrayfun(@(i) PDE0.vars(i,1).varname{1},1:nvars0,'UniformOutput',false);
space_pde = struct('vars',{vars_pde0},'dom',reshape(double(PDE0.dom),nvars0,2));
[~,vars_idcs0] = ismember(vars_pde0,vars_pde);
if any(vars_idcs0==0)
    vars_idcs0 = 1:nvars0;
end
dep_pde = struct();
for obj_type={'x','w','u','y','z'}
    obj_tab = PDE0.([obj_type{1},'_tab']);
    ID_arr = zeros(size(obj_tab,1),1);
    for i=1:numel(ID_arr)
        ID_arr(i) = PDE0.(obj_type{1}){i}.ID;
    end
    space_pde.(obj_type{1}) = struct('dims',obj_tab(:,2),'depmat',logical(obj_tab(:,3:2+nvars0)),'ID',ID_arr);
    dep_pde.(obj_type{1}) = false(size(obj_tab,1),N);
    dep_pde.(obj_type{1})(:,vars_idcs0) = space_pde.(obj_type{1}).depmat;
end
PIE.space_pde = space_pde;

% % % Maps from the PDE to the PIE, from the map Dmap of the BCs: 
% % %   xf_b = C_b*D^{diff}x_pde_idx(b), 
% % % with C_b the block F_b of Dmap.x.C, the identity for the main blocks
% % % and the integral for the auxiliary ones. The inputs are evaluated at
% % % the lower boundary of the variables they no longer depend on, C=I.
pde2pie_s = PIE.pde2pie;
% Components added for the delays are the history of a state or input of
% the original PDE (see del_tab): they depend on the delay variable, in
% which a derivative is a temporal derivative, and are evaluated at loc
idcs_xwu = zeros(nblocks,3);      loc = zeros(nblocks,N);      delay_idcs = zeros(nblocks,1);
diff_mat = zeros(nblocks,N);     tdiff = zeros(nblocks,1);
for b=1:nblocks
    diff_mat(b,:) = reshape(Dmap.x.diff(b,x_idcs(b),:),1,N);
end
% Source PDE component, and derivatives, of each block of xf
for b=1:nblocks
    x_idx = x_idcs(b);
    if x_idx<=nx0
        idcs_xwu(b,1) = x_idx;       continue
    end
    % A history block: the source (state, w or u) from del_tab
    del_idx = x_idx-nx0;              delay_src = del_tab(del_idx,1);
    if delay_src<=nx0
        idcs_xwu(b,1) = delay_src;
    elseif delay_src<=nx0+nw0
        idcs_xwu(b,2) = delay_src-nx0;
    else
        idcs_xwu(b,3) = delay_src-nx0-nw0;
    end
    % Delay variable: the first variable of the block that the source does
    % not depend on, where a derivative is a time derivative; positions and
    % derivatives in the original variables from del_tab
    col_idx = find(idcs_xwu(b,:));
    dep_idcs = find(space.x.depmat(b,:));
    rest_idcs = dep_idcs(~dep_pde.(xwu_cell{col_idx})(idcs_xwu(b,col_idx),dep_idcs));
    if isempty(rest_idcs)
        rest_idcs = dep_idcs;
    end
    delay_idcs(b) = rest_idcs(1);
    tdiff(b) = tdiff(b)+diff_mat(b,delay_idcs(b));
    diff_mat(b,:) = 0;
    loc(b,vars_idcs0) = del_tab(del_idx,2+(1:nvars0));
    diff_mat(b,vars_idcs0) = del_tab(del_idx,2+size(del_tab,2)/2-1+(1:nvars0));
end
% Map from the PDE state: sources, positions, delays, derivatives, and
% the operators C_b
pde2pie_s.x.idx = idcs_xwu;        pde2pie_s.x.loc = loc;        pde2pie_s.x.delay = delay_idcs;
pde2pie_s.x.diff = diff_mat;     pde2pie_s.x.tdiff = tdiff;
pde2pie_s.x.C = op_matrix(blocks,spaces_pde(x_idcs),eye(nblocks)==1,@(r,j) F_cell{r});
% Inputs: the derivatives, positions and maps C of the fundamental inputs
for obj_type={'w','u'}
    inputs = info.inputs.(obj_type{1});       ninp = numel(inputs);
    pde2pie_s.(obj_type{1}).idx = Dmap.(obj_type{1}).idx(:);
    pde2pie_s.(obj_type{1}).diff = Dmap.(obj_type{1}).diff;
    pde2pie_s.(obj_type{1}).tdiff = zeros(ninp,1);
    pde2pie_s.(obj_type{1}).loc = Dmap.(obj_type{1}).loc;
    if ninp>0
        pde2pie_s.(obj_type{1}).C = Dmap.(obj_type{1}).C;
    end
end
% Outputs: the PDE output of each PIE output
pde2pie_s.y.idx = PDE.y_tab(:,1);        pde2pie_s.z.idx = PDE.z_tab(:,1);
PIE.pde2pie = pde2pie_s;

% % % Maps from the PIE to the PDE: the PDE state k is the main block of
% % % the PIE state (the other blocks are tied to it by the dynamics), the
% % % outputs are the same, and the inputs not in the BCs are PIE inputs
pie2pde_s = PIE.pie2pde;
% X: from the rows of T to the PDE states (zero rows for algebraic states)
pie2pde_s.X = op_matrix(spaces_pde(1:nx0),blocks(dyn_rows),dyn_rows==main_idcs(1:nx0),@(k,b) identity_op(spaces_pde(k)));   % DJ, 11/20/2026
% T, Tw, Tu: from xf, wf, uf to all PDE states; the empty default is kept
% where zero
Tnames = struct('x','T','w','Tw','u','Tu');
for fn = {'x','w','u'}
    if ~isempty(Top_pde.(fn{1}))
        pie2pde_s.(Tnames.(fn{1})) = Top_pde.(fn{1});
    end
end
% Maps of the inputs and outputs: the same components, or direct inputs
for obj_type={'w','u','y','z'}
    obj_tab = PDE.([obj_type{1},'_tab']);      nr = size(obj_tab,1);
    row_spaces_pde = arrayfun(@(i) get_space(obj_tab,i),1:nr);
    is_input = any(strcmp(obj_type{1},{'w','u'}));
    if is_input
        col_spaces = info.inputs.(obj_type{1});
    else
        col_spaces = row_spaces_pde;
    end
    if nr==0 || isempty(col_spaces)
        continue
    end
    hit = eye(nr)==1;
    if is_input
        % The inputs that are their own fundamental input (not decomposed)
        hit = ([col_spaces.pde]==(1:nr)') & ~cellfun(@isempty,info.(['I',obj_type{1}])(:));
    end
    pie2pde_s.(upper(obj_type{1})) = op_matrix(row_spaces_pde,col_spaces,hit,@(r,j) identity_op(row_spaces_pde(r)));
end
PIE.pie2pde = pie2pde_s;

% The variables of the polynomials f, gz, gy of a nonlinear PIE are named
% as the PIE state PIE.x, by the default names x01_pie, ... assigned to
% the components on initialization
PIE = initialize(PIE);

% Expand higher-order temporal derivatives, introducing states for the
% temporal derivatives, so that the PIE is first-order in time. The order
% of the derivative in the equation of each component is that of its state.
tdiff_tab = [(1:nblocks)',zeros(nblocks,1)];
if has_hotd
    [PIE,tdiff_tab] = expand_tderivatives(PIE,tdiff_pde(x_idcs),true);
end
% Provide an overview of how the new state variables relate to the PDE
% variables.
if ~suppress_summary
    for obj_summary = {'x','u','w','y','z'}
        if isempty(PIE.space.(obj_summary{1}).dims)
            continue
        elseif strcmp(obj_summary{1},'x') && any(is_aux)
            % DB 11/07/2026: usual summary, where the auxiliary states are
            % displayed as the operator F applied to x, and a note on
            % K*xf=0 (tdiff_tab(:,1): the component each one derives from)
            aux_idcs = find(is_aux(tdiff_tab(:,1)))';
            print_convert_summary(PIE,'x',aux_idcs,F_cell(tdiff_tab(aux_idcs,1)));
            fprintf('\nThe PIE state must satisfy PIE.K*xf = 0.\n');
            fprintf('The states are defined by PIE.pde2pie.x.C.\n');
        else
            print_convert_summary(PIE,obj_summary{1});
        end
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
