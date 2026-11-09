function PIE = convert_PIETOOLS_PDE(PDE,comp_order,flag)
% PIE = CONVERT_PIETOOLS_PDE(PDE,FLAG) converts a PDE in the terms format,
% possibly coupling ODE states with PDE states on several spatial
% variables, and possibly nonlinear, to an equivalent PIE.
%
% INPUT
% - PDE:    'pde_struct' object defining a PDE in the terms format (see
%           also "@pde_struct/initialize"). A PDE in one of the legacy
%           formats (batch, or 'terms_legacy') is passed on to the legacy
%           converters of PIETOOLS, as in the PIETOOLS converter;
% - flag:   (optional) 'char' or cell of 'char' objects, containing
%           'silent' to suppress the display of an overview of how the
%           fundamental state and inputs of the PIE relate to those of the
%           PDE, and/or 'Top' (or 'return_Top') to stop the conversion
%           after the T operators (from fundamental to PDE state) have been
%           computed, and/or 'nonlinear' to return the PIE in the nonlinear
%           format, with polynomials f, gz, gy, also if the PDE is linear
%           (needed for the SOS programming of PIESOS). For compatibility,
%           FLAG may also be
%           specified as third input, following a (not supported) argument
%           comp_order;
%
% OUTPUT
% - PIE:    'pie_struct' object if the PDE is linear, and 'npie_struct'
%           object if it is nonlinear or the flag 'nonlinear' is
%           specified (see "@pie_struct/pie_struct"),
%           describing the PIE
%           d/dt (T*xf + Tw*wf + Tu*uf) = A*xf + Bw*wf + Bu*uf;
%                                        z = Cz*xf + Dzw*wf + Dzu*uf;
%                                        y = Cy*xf + Dyw*wf + Dyu*uf;
%                                K*xf = 0,
%           with the fundamental state xf, and the fundamental exogenous
%           and actuator inputs wf and uf. For a nonlinear PDE, A*xf,
%           Cz*xf, and Cy*xf are replaced by the polynomials f, gz,
%           and gy. The operators are 'copvar' objects (linear PDE) or
%           'tensopmat' objects (nonlinear PDE). The components of xf
%           are the blocks info.blocks of solve_PIETOOLS_BVP, one for each
%           state in each component of the fundamental state, and those of
%           wf and uf the fundamental inputs info.inputs.w and
%           info.inputs.u. The
%               fields space, space_pde, pie2pde, pde2pie
%           describe the function spaces of the PIE and PDE, and how the
%           PIE state and inputs are derived from those of the PDE: for
%           each component b of xf, the PDE state idx, the derivative
%           diff taken of it, and the operator C (see
%           "@pie_struct/pie_struct"). If flag contains 'Top', PIE only has
%           the fields space, T, Tw, Tu and K;
%
% NOTES
% - Each PDE state has one "main" block of xf, its highest derivative in
%   each variable. T has a row for each block: for a main block the row of
%   T defining x_pde, and for any other (auxiliary) block xf_b, a moment
%   of x_pde (xf_b=F*x_pde, or F*D^{a}x_pde if the moment is taken of a
%   derivative in the other variables, e.g. for BCs that are periodic in
%   one of several variables), the row of T defining F*x_pde: the selector
%   of xf_b if a=0, and F times the row of x_pde otherwise. The
%   right-hand side of that row is F times that of x_pde, and the
%   constraint K*xf=0 ties xf_b to x_pde.
% - The inputs in the boundary conditions are decomposed into fundamental
%   inputs wf and uf, as in solve_PIETOOLS_BVP: each is a derivative of an
%   element of an input of the PDE, evaluated at the lower boundary in some
%   variables (e.g. for a state with BCs along both s1 and s2, the corner
%   value of an input, and its derivatives along the boundaries). This is
%   stored in pde2pie (idx, diff, loc, and C selecting the elements), see
%   the description of Dmap in solve_PIETOOLS_BVP. This requires all
%   directions of the group of BCs to be admissible (e.g. no periodic
%   BCs), and each BC with an input to have a single element of a single
%   input as right-hand side. An input of fewer variables than its
%   boundary is constant in the others, which restricts the other inputs
%   (e.g. scalar inputs on all faces of a 2D state must be equal, so there
%   is one fundamental input). An input that is shared by the BCs of
%   states on different variables is one fundamental input (see
%   solve_PIETOOLS_BVP). An input that does not appear in the BCs is its
%   own fundamental input (see append_inputs). An input that appears both
%   in a BC and elsewhere in the PDE is not supported. The map pie2pde.w,
%   pie2pde.u from PIE inputs to PDE inputs is only specified for the
%   inputs not in the BCs.
% - The right-hand side of the PIE is derived by convert_PDE_term, so that
%   inputs can only appear linearly.
% - Delays are modeled by transport states (expand_delays), and
%   higher-order temporal derivatives (linear PDEs, any number of spatial
%   variables) by additional states (expand_tderivatives); pde2pie.x.
%   delay and pde2pie.x.tdiff describe these components. Expansion is not
%   supported for states with boundary inputs.
% - ODE states, and states of lower dimension, may appear in the BCs of
%   other states, see solve_PIETOOLS_BVP; the T operators then map the
%   fundamental state of the former to the latter.
%

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
% DJ, 11/02/2026: Initial coding (ND version);
% DJ, 11/06/2026: Display PIE states and inputs with a hat in the summary;
% DJ, 11/06/2026: Get the names of components and variables in the summary
%                   from 'varname2display';
% DJ, 11/08/2026: Drop the zero maps of the inputs by info.input_maps;;

% % % Pass PDEs in the legacy formats to the legacy converters, as the
% % % PIETOOLS converter does, and convert other structures to pde_struct
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
% DB 11/07/2026: outputs flagged as static constraints 0=E(s) (see
% initialize) are not outputs of the PIE; their operator is added to K: the PDE is
% converted without them, and with only them, which have the same
% fundamental state, and their rows Cz are the rows of K that are added
if ~isempty(PDE.constraint_ID) && ~return_Top
    z_IDs = cellfun(@(zi) zi.ID,PDE.z);
    is_con = ismember(z_IDs(:),PDE.constraint_ID(:));
    PDE_con = PDE;      PDE_con.constraint_ID = [];
    PDE_con.z = PDE.z(is_con);          PDE_con.z_tab = PDE.z_tab(is_con,:);
    PDE_con.z_tab(:,1) = (1:nnz(is_con))';
    PDE_main = PDE;     PDE_main.constraint_ID = [];
    PDE_main.z = PDE.z(~is_con);        PDE_main.z_tab = PDE.z_tab(~is_con,:);
    PDE_main.z_tab(:,1) = (1:nnz(~is_con))';
    PIE = convert_PIETOOLS_PDE(PDE_main,[flag(:)',{'constraint_main'}]);
    PIE_con = convert_PIETOOLS_PDE(PDE_con,[flag(:)',{'silent'}]);
    for ic=1:numel(PDE_con.z)
        for jc=1:numel(PDE_con.z{ic}.term)
            if any(isfield(PDE_con.z{ic}.term{jc},{'w','u'}))
                error("Constraints 0=F(s) cannot involve the inputs w or u.")
            end
        end
    end
    if isempty(PIE.K) || isempty(PIE.K.C)
        PIE.K = PIE_con.Cz;
    else
        PIE.K = [PIE.K; PIE_con.Cz];
    end
    return
end
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
% Higher-order temporal derivatives (linear PDEs): the PIE is first derived
% with these, and expanded to first order afterwards, see
% expand_tderivatives
tdiff_pde = ones(numel(PDE.x),1);
for i=1:numel(PDE.x)
    if isfield(PDE.x{i},'tdiff')
        tdiff_pde(i) = PDE.x{i}.tdiff;
    end
end
has_hotd = any(tdiff_pde>1);
is_nonlinear = PDE.is_nonlinear || force_nonlinear;
if has_hotd && is_nonlinear
    error("Higher-order temporal derivatives in nonlinear PDEs are not yet supported.")
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
% % % state, with a block for each state (see solve_PIETOOLS_BVP)
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
[Tops,Kops,Dmap,info] = solve_PIETOOLS_BVP(PDE,diffs,opts);
Gops = info.Gops;
blocks = info.blocks;       nblocks = numel(blocks);
spaces_pde = arrayfun(@(k) get_space(PDE.x_tab,k),1:nx);

% Each PDE state k has one main block (derived from this state alone, on
% all its variables); the other blocks are auxiliary states
% xf_b=F*x_pde, with x_idcs(b) the state of block b and F_cell{b} the
% block of F
if any(Dmap.x.idx==0)
    error("convert_PIETOOLS_PDE:blocks","Fundamental state components that cannot be split into one block per state are not supported.")
end
% List the blocks by the number of variables they depend on, ascending:
% the finite-dimensional blocks (ODE states, moments) first, as in the old
% converters, so that e.g. an ODE state precedes the PDE states in the PIE,
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
    Tops.x = permute_cols(Tops.x,perm,is_nonlinear);
    if size(Kops.x,1)>0
        Kops.x = permute_cols(Kops.x,perm,is_nonlinear);
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
if isfield(Gops,'Iw')
    Tops_term.Iw = Gops.Iw;     Tops_term.Iu = Gops.Iu;
end
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
            [Pterm,Pblk] = convert_PDE_term(PDE.(obj{1}){i}.term{j},Tops_term,PDE,struct('nonlinear',is_nonlinear,'outvars',{vars_eq}));
            if is_nonlinear && ~isempty(Pterm.x)
                if isempty(poly_struct.(obj{1}){i})
                    poly_struct.(obj{1}){i} = Pterm.x;
                else
                    poly_struct.(obj{1}){i} = poly_struct.(obj{1}){i}+Pterm.x;
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
        try
            PIE.(op_names.(obj{1}){xwu_idx}) = cell2op(blocks_new);
        catch ME
            error("convert_PIETOOLS_PDE:operator","Could not construct the operator %s: %s",op_names.(obj{1}){xwu_idx},ME.message)
        end
    end
end
% Constraints K*xf=0 on the fundamental state, for the blocks of xf
if size(Kops.x,1)>0
    PIE.K = Kops.x;
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
        end
    end
    f_rows = poly_struct.x(x_idcs);
    for b=find(is_aux)'
        % F*f is not supported for f of degree>1: periodic BCs with a
        % nonlinear right-hand side error with a clear message
%        f_rows{b} = tensopmat(F_cell{b})*poly_struct.x{x_idcs(b)};         % DJ (CC), 11/02/2026 (was)
        try                                                                 % DJ (CC), 11/02/2026
            f_rows{b} = tensopmat(F_cell{b})*poly_struct.x{x_idcs(b)};      % DJ (CC), 11/02/2026
        catch ME                                                            % DJ (CC), 11/02/2026
            error("convert_PIETOOLS_PDE:nonlinearAux",...
                "Moments of the state (periodic BCs) with a nonlinear right-hand side are not supported: %s",ME.message)    % DJ (CC), 11/02/2026
        end                                                                 % DJ (CC), 11/02/2026
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

% % % Maps from the PDE to the PIE, as Dmap of solve_PIETOOLS_BVP: 
% % %   xf_b = C_b*D^{diff}x_pde_idx(b), 
% % % with C_b the block F_b of Dmap.x.C, the identity for the main blocks 
% % % and the integral for the auxiliary ones. The inputs are evaluated at 
% % % the lower boundary of the variables they no longer depend on, and C=I.
pde2pie_s = PIE.pde2pie;
% Components added by expand_delays are the history of a state or input of
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
    del_idx = x_idx-nx0;              delay_src = del_tab(del_idx,1);
    if delay_src<=nx0
        idcs_xwu(b,1) = delay_src;
    elseif delay_src<=nx0+nw0
        idcs_xwu(b,2) = delay_src-nx0;
    else
        idcs_xwu(b,3) = delay_src-nx0-nw0;
    end
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
pde2pie_s.x.idx = idcs_xwu;        pde2pie_s.x.loc = loc;        pde2pie_s.x.delay = delay_idcs;
pde2pie_s.x.diff = diff_mat;     pde2pie_s.x.tdiff = tdiff;
pde2pie_s.x.C = op_matrix(blocks,spaces_pde(x_idcs),eye(nblocks)==1,@(r,j) F_cell{r});
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
pde2pie_s.y.idx = PDE.y_tab(:,1);        pde2pie_s.z.idx = PDE.z_tab(:,1);
PIE.pde2pie = pde2pie_s;

% % % Maps from the PIE to the PDE: the PDE state k is the main block of
% % % the PIE state (the other blocks are tied to it by the dynamics), the
% % % outputs are the same, and the inputs not in the BCs are PIE inputs
pie2pde_s = PIE.pie2pde;
pie2pde_s.x = op_matrix(spaces_pde(1:nx0),blocks,(1:nblocks)==main_idcs(1:nx0),@(k,b) identity_op(spaces_pde(k)));
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
        hit = ([col_spaces.pde]==(1:nr)') & ([col_spaces.group]==0);
    end
    pie2pde_s.(obj_type{1}) = op_matrix(row_spaces_pde,col_spaces,hit,@(r,j) identity_op(row_spaces_pde(r)));
end
PIE.pie2pde = pie2pde_s;

% The variables of the polynomials f, gz, gy of a nonlinear PIE are named
% as the PIE state PIE.x, by the default names x01_pie, ... that
% initialize assigns to the components (see convert_PDE_term)
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
        if strcmp(obj_summary{1},'x') && any(is_aux)
            % DB 11/07/2026: usual summary, where the auxiliary states are
            % displayed as the operator F applied to x, and a note on K*xf=0
            if size(PIE.x_tab,1)>=1
                print_convert_summary(PIE,'x',tdiff_tab(:,2),'',find(is_aux)',F_cell(is_aux));
            end
            fprintf('\nThe PIE state must satisfy PIE.K*xf = 0.\n');
            fprintf('The states are defined by PIE.pde2pie.x.C.\n');
        elseif size(PIE.([obj_summary{1},'_tab']),1)>=1
            print_convert_summary(PIE,obj_summary{1},tdiff_tab(:,2));
        end
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function print_convert_summary(PIE,obj,tdiff_list,aux_msg,aux_idcs,F_aux)
% PRINT_CONVERT_SUMMARY(PIE,obj,tdiff_list) prints information in the
% command window on how the state variables, inputs, and outputs in the
% computed PIE structure relate to those in the original PDE.
%
% INPUT
% - PIE:    'pie_struct' object defining a PIE;
% - obj:    char 'x', 'u', 'w', 'y', or 'z', indicating for which
%           object to display how they relate to the associated object in
%           the PDE;
% - tdiff_list: if obj='x', an nx x 1 array indicating for each of the nx
%           (vector-valued) fundamental state components what temporal
%           derivative of the corresponding PDE state component they
%           correspond to;
%
% OUTPUT
% Displays information in the command window on how the state components,
% inputs, and outputs in the PIE relate to those in the PDE.
% Set a name associated to each object.
if strcmp(obj,'x')
    object_name = 'fundamental state component';
elseif strcmp(obj,'y')
    object_name = 'observed output';
elseif strcmp(obj,'z')
    object_name = 'regulated output';
elseif strcmp(obj,'u')
    object_name = 'actuator input';
elseif strcmp(obj,'w')
    object_name = 'exogenous input';
end
obj_tab = PIE.([obj,'_tab']);
ncomps = size(obj_tab,1);
nvars = size(PIE.vars,1);
has_diffs = false;
if all(obj_tab(:,1)==(1:ncomps)') && (size(obj_tab,2)<=2+nvars || ~any(any(obj_tab(:,3+nvars:2+2*nvars))))
    % The components in the PIE correspond exactly to the components in the
    % PDE.
    return
elseif size(obj_tab,2)<=2+nvars || ~any(any(obj_tab(:,3+nvars:2+2*nvars)))
    % Otherwise, we list the new order of the components.
    fprintf(['\n','The ',object_name,'s have been reindexed as:\n']);
else
    % Indicate how the new components relate to those in the PDE
    has_diffs = true;
    fprintf(['\n','The following ',object_name,'s have been introduced:\n']);  
    if nargin>=4 && ~isempty(aux_msg)
        fprintf(aux_msg);       % DB 11/07/2026
    end
end

% Use UNICODE to add subscript indices to different components.
thin_space = char(8201);
partial = '\x2202';
sub_t = '\x209C';
% The PIE state and inputs are indicated as such in their names (a hat), unlike
% those of the PDE; the names are given by 'varname2display'
pie_tag = '';
if any(strcmp(obj,{'x','u','w'}))
    pie_tag = '_pie';
end
%sub_s = '\x209B';
%int = '\x222B';
%xdot = '\x1E8B';
%xdot = [partial,sub_t,' x'];
% Also use UNICODE for superscript orders of the temporal derivative
sup_num = cell(10,1);    % Superscripts
sup_num{1} = '\x2070';  % Superscript 0
sup_num{2} = '\xB9';
sup_num{3} = '\xB2';
sup_num{4} = '\xB3';
sup_num{5} = '\x2074';
sup_num{6} = '\x2075';
sup_num{7} = '\x2076';
sup_num{8} = '\x2077';
sup_num{9} = '\x2078';
sup_num{10} = '\x2079';

% Determine how many subscripts will be needed.
n_digits = 1+floor(log(ncomps)/log(10));


% % For the purpose of clearer display, we make an estimate of how long the
% % display will be. 
% First, determine which variables appear in the different state variables.
global_vars_obj = PIE.vars(any(PIE.([obj,'_tab'])(:,3:2+nvars),1),:);
lngth_varnames = 3;      % just three characters: (t)
if ~isempty(global_vars_obj) || has_diffs
    varname_cell = cell(nvars,1);
    if nvars==1
        varname_cell{1} = 's';
        lngth_varnames = lngth_varnames+2;    % add two characters: ,s
    else
        varname_cell{1} = varname2display('s1');
        lngth_varnames = lngth_varnames+2;    % add three characters: ,s1
        for kk=2:nvars
            varname_cell{kk} = varname2display(sprintf('s%d',kk));
            lngth_varnames = lngth_varnames+2;    % add three characters: ,sj
        end
    end
end
% Then, estimate the size of the string of characters denoting the
% components "PDE.x{ii}".
nvars_max = max(sum(PIE.([obj,'_tab'])(:,3:2+nvars),2));
lngth_varnames_mean = ceil(lngth_varnames*nvars_max/nvars);
LHS_length_max = 1+1 + n_digits + lngth_varnames_mean+3; % e.g. ' x13(t,s1,s2,s3)', 


% % For each of the components, display its size, and which variables it
% % depends on.
for ii=1:ncomps
    old_ID = obj_tab(ii,1);
    % Check if this component is derived from some higher-dimensional
    % component in the original PDE (e.g. a boundary value of some input)
    comp_idx = find(obj_tab(:,1)==old_ID,1,'last');

    % Establish the names of the variables on which the component depends.
    varnames_LHS_t = 't';
    LHS_length = 1;
    var_idcs = find(obj_tab(ii,3:2+PIE.dim));
    for kk=var_idcs
        varnames_LHS_t = [varnames_LHS_t,',',varname_cell{kk}];
        LHS_length = LHS_length+3;  % add 3 characters: ,s1
    end
    if ii~=comp_idx
        % In case the current component is derived from some other
        % component, it may not depend on the same variables
        var_idcs_RHS = find(obj_tab(comp_idx,3:2+PIE.dim));
        varnames_RHS_t = 't';
        for kk=var_idcs_RHS
            is_var = ismember(kk,var_idcs);
            if is_var
                % The new component still depends on this variable
                varnames_RHS_t = [varnames_RHS_t,',',varname_cell{kk}];
            else
                % The component is evaluated at the lower boundary w/r to
                % this spatial variable
                varnames_RHS_t = [varnames_RHS_t,',',num2str(PIE.dom(kk,1))];
            end
        end
    else
        varnames_RHS_t = varnames_LHS_t;
    end
    
    % Establish the index for the new component, displayed as a subscript.
    new_idx = ii;
    Lcomp_idx = '';
    if ncomps>1
        Lcomp_idx = num2str(new_idx);
        LHS_length = LHS_length + length(Lcomp_idx);
    end
    % Set the name of the component, including its dependence on spatial
    % variables.
    LHS_name = [' ',varname2display([obj,pie_tag,Lcomp_idx]),'(',varnames_LHS_t,')'];
    LHS_length = 1 + LHS_length + 3;

    % Establish the index for the old component
    Rcomp_idx = '';
    if ncomps>1
        Rcomp_idx = num2str(old_ID);
    end
    % Set the name of the component, including its dependence on spatial
    % variables.
    RHS_name = [thin_space,varname2display([obj,Rcomp_idx]),'(',varnames_RHS_t,')'];
    % For state components, also indicate what spatial/temporal derivative
    % of the PDE state they correspond to.
    if strcmp(obj,'x') && any(tdiff_list)
        tdiff = tdiff_list(ii);
        if tdiff==0
            diff_str = ['       ',thin_space];
        elseif tdiff==1
            diff_str = [' (',partial,'/',partial,'t)',thin_space];
        else
            t_sup = cell2mat(sup_num(str2num(num2str(tdiff)')+1)');
            diff_str = ['(',partial,'/',partial,'t)',t_sup,thin_space];
        end
    else
        diff_str = '';
    end
    if has_diffs
        Dvals = obj_tab(ii,2+nvars+1:2+2*nvars);
        for jj=1:nvars
            sdiff = Dvals(jj);
            if sdiff==0
                diff_str = ['        ',diff_str,thin_space];
            elseif sdiff==1
                diff_str = [diff_str,'(',partial,'/',partial,varname_cell{jj},')',thin_space];
            else
                s_sup = cell2mat(sup_num(str2num(num2str(sdiff)')+1)');
                diff_str = [diff_str,'(',partial,'/',partial,varname_cell{jj},')',s_sup,thin_space];
            end
        end
        RHS_name = [diff_str,RHS_name];
    end
    if nargin>=6 && ismember(ii,aux_idcs)
        % DB 11/07/2026: auxiliary state xf=F*x, displayed as the integral
        % form of F*x (not as a boundary value of x)
        F_str = 'F*x';
        try
            F_blk = F_aux{aux_idcs==ii};
            F_d = get_display(F_blk,varname2display(['x',Rcomp_idx]));
            if F_d.ok
                F_str = '';
                for kk=1:numel(F_d.terms)
                    if kk==1
                        sgn = repmat('-',1,F_d.signs(kk)<0);
                    else
                        sgn = [' ',char(43+2*(F_d.signs(kk)<0)),' '];
                    end
                    F_str = [F_str,sgn,F_d.terms{kk}];
                end
            end
        catch
        end
        RHS_name = [diff_str,thin_space,F_str];
    end
    % % % Finally, display:
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),' <--  ',RHS_name,'\n']);
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function P = permute_cols(P,perm,is_nonlinear)
% P = PERMUTE_COLS(P,PERM,IS_NONLINEAR) reorders the block columns of the
% operator P, or of each operator in a cell array of operators.
%
% INPUT
% - P:      'copvar' object (IS_NONLINEAR false), 'tensopmat' object
%           (IS_NONLINEAR true), [], or a cell array of such objects;
% - perm:   permutation of the block columns of P;
% - is_nonlinear:   true if the operators are 'tensopmat' objects;
%
% OUTPUT
% - P:      the same operators, with block column j of the output the
%           block column perm(j) of the input;
%
% DJ (CC), 11/02/2026: Initial coding;
if iscell(P)
    P = cellfun(@(Pk) permute_cols(Pk,perm,is_nonlinear),P,'UniformOutput',false);
elseif ~isempty(P)
    if is_nonlinear
        B = cellfun(@(o) o.ops{1},P.ops,'UniformOutput',false);
        P = tensopmat(B(:,perm));
    else
        P = copvar(P.C(:,perm));
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
