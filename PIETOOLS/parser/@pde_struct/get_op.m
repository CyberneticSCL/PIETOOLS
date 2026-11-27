function [Pop_x,Pop_w,Pop_u] = get_op(PDE,space_out)
% [POP_X,POP_W,POP_U] = GET_OP(PDE,SPACE_OUT) returns the PI operators
% defined by the free terms PDE, as 'opvar' objects acting on the states,
% exogenous inputs and controlled inputs,
%   PDE = Pop_x*[x_1;...;x_n] + Pop_w*[w_1;...;w_m] + Pop_u*[u_1;...;u_k],
% with x_j = PDE.x{j}, w_j = PDE.w{j}, u_j = PDE.u{j}.
%
% INPUT
% - PDE:        'pde_struct' object representing one or several free terms
%               (or a single PDE state or input variable), without
%               derivatives, evaluation at a position, delays or products
%               of states or inputs (nonlinear terms);
% - space_out:  (optional) the spatial variables on which the output
%               depends, as one of
%               - 1 x n 'polynomial' array of pvars, e.g. [s1,s2], or
%                 cellstr, e.g. {'s1','s2'}: the same variables for all
%                 free equations;
%               - M x 1 cell with a 'polynomial' array or cellstr per free
%                 equation, e.g. {s1; [s1,s3]}, for M free equations;
%               - 'pde_struct' object, a column of M PDE variables, whose
%                 spatial variables and domains give the output space of
%                 each free equation;
%               The domain of each variable is that of the states and
%               inputs, or of the PDE variables in space_out; a variable
%               of neither is a new variable on [0,1]. Defaults to the
%               variables on which each free equation depends.
%
% OUTPUT
% - Pop_x:      'opvar' object, one row per free equation of PDE, one
%               column per state component in PDE.x;
% - Pop_w:      'opvar' object for the components of PDE.w, [] if none;
% - Pop_u:      'opvar' object for the components of PDE.u, [] if none;
%
% NOTES
% - A free equation integrated over the full domain of s does not depend
%   on s, and maps into a space without s, unless the domain of
%   integration is given as [a,s,b] or s is in space_out.
% - In a term int(C*x,s,...), the kernel is C, with the variable of
%   integration as second argument.
% - A coefficient depending on a variable r of no state or input defines
%   a multiplier along a new variable r on [0,1]: for x in R, s*x gives
%   the operator R -> L2[0,1], x -> s*x. The time t and dummy variables
%   (names ending in '_dum') are not new variables.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_op
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
%
% Initial coding DJ - 11/22/2026

% % % Process the inputs
if ~is_pde_term(PDE)
    error("The input should correspond to one or several free terms in a PDE.")
end
M = numel(PDE.free);
dim_out = cellfun(@(eq) eq.size,PDE.free(:));

% % Registry of spatial variables and domains of the states and inputs.
objs = {'x','w','u'};
vnames_reg = cell(1,0);     dom_reg = zeros(0,2);
for kk=1:numel(objs)
    for jj=1:numel(PDE.(objs{kk}))
        [vnames_reg,dom_reg] = add_vars(vnames_reg,dom_reg,PDE.(objs{kk}){jj});
    end
end

% % Output space of each free equation: the variables on which it depends
% % by default, or the variables specified by space_out.
vnames_out = cell(M,1);
if nargin<2
    % Variables on which each free equation depends.
    for ii=1:M
        vars_ii = PDE.free{ii}.vars;
        vnames_out{ii} = arrayfun(@(k) vars_ii(k,1).varname{1},1:size(vars_ii,1),'UniformOutput',false);
    end
elseif isa(space_out,'pde_struct')
    % Variables and domains of a column of PDE variables.
    if numel(space_out.free)~=M
        error("The number of PDE variables in the output space should match the number of free equations.")
    end
    for obj = {'x','w','u','y','z'}
        for jj=1:numel(space_out.(obj{1}))
            [vnames_reg,dom_reg] = add_vars(vnames_reg,dom_reg,space_out.(obj{1}){jj});
        end
    end
    for ii=1:M
        if space_out.free{ii}.size~=dim_out(ii)
            error("The size of PDE variable %d of the output space should match that of free equation %d.",ii,ii)
        end
        vars_ii = space_out.free{ii}.vars;
        vnames_out{ii} = arrayfun(@(k) vars_ii(k,1).varname{1},1:size(vars_ii,1),'UniformOutput',false);
    end
else
    % A list of variables for all free equations, or a cell with a list
    % per free equation.
    if ischar(space_out) || iscellstr(space_out) || isa(space_out,'polynomial')
        space_cell = repmat({space_out},M,1);
    elseif iscell(space_out) && numel(space_out)==M
        space_cell = space_out(:);
    elseif iscell(space_out)
        error("The output space should give one list of variables per free equation (%d), not %d.",M,numel(space_out))
    else
        error("The output space should be specified as list of variables, cell with a list of variables per free equation, or 'pde_struct' object.")
    end
    for ii=1:M
        list_ii = space_cell{ii};
        if isa(list_ii,'polynomial')
            if ~isempty(list_ii) && ~ispvar(list_ii)
                error("The variables of the output space should be specified as pvar objects.")
            end
            list_ii = list_ii(:);
            vnames_out{ii} = arrayfun(@(k) list_ii(k).varname{1},1:size(list_ii,1),'UniformOutput',false);
        elseif ischar(list_ii) || iscellstr(list_ii) || isempty(list_ii)
            vnames_out{ii} = reshape(cellstr(list_ii),1,[]);
        else
            error("The variables of the output space should be specified as 'polynomial' array or cellstr.")
        end
    end
end
% A variable of the output space that is no variable of the states,
% inputs, or PDE variables in space_out is a new spatial variable on
% [0,1], along which the operator is a multiplier. The time t and dummy
% variables of integration cannot be spatial variables of the output.
for ii=1:M
    for vname = vnames_out{ii}(~ismember(vnames_out{ii},vnames_reg))
        if strcmp(vname{1},'t')
            error("Variable 't' is the time, and cannot be a spatial variable of the output space.")
        elseif endsWith(vname{1},'_dum')
            error("Variable '%s' appears outside of an integral over it; a dummy variable of integration cannot be a spatial variable of the output space.",vname{1})
        end
        vnames_reg = [vnames_reg, vname];                                   %#ok<AGROW>
        dom_reg = [dom_reg; 0,1];                                           %#ok<AGROW>
    end
end
[vnames_reg,idcs] = sort(vnames_reg);   dom_reg = dom_reg(idcs,:);


% % % Build the blocks, one term at a time.
blocks = struct('x',{cell(M,numel(PDE.x))},'w',{cell(M,numel(PDE.w))},'u',{cell(M,numel(PDE.u))});
for ii=1:M
    if ~isfield(PDE.free{ii},'term'),   continue;   end
    [~,loc_out] = ismember(vnames_out{ii},vnames_reg);
    dom_out = dom_reg(loc_out,:);
    for jj=1:numel(PDE.free{ii}.term)
        [obj,Ridx,blk] = term2sopvar(PDE.free{ii}.term{jj},PDE,vnames_out{ii},dom_out,dim_out(ii));
        if isempty(blocks.(obj){ii,Ridx})
            blocks.(obj){ii,Ridx} = blk;
        else
            blocks.(obj){ii,Ridx} = blocks.(obj){ii,Ridx}+blk;
        end
    end
end

% % % Assemble the containers.
space_out_mask = false(M,numel(vnames_reg));
for ii=1:M
    space_out_mask(ii,:) = ismember(vnames_reg,vnames_out{ii});
end
Pops = cell(1,3);
for kk=1:numel(objs)
    comps = PDE.(objs{kk});
    if isempty(comps)
        continue
    end
    N = numel(comps);
    space_in_mask = false(N,numel(vnames_reg));
    dim_in = zeros(N,1);
    for jj=1:N
        vars_jj = comps{jj}.vars;
        vnames_jj = arrayfun(@(k) vars_jj(k,1).varname{1},1:size(vars_jj,1),'UniformOutput',false);
        space_in_mask(jj,:) = ismember(vnames_reg,vnames_jj);
        dim_in(jj) = comps{jj}.size;
    end
    meta = struct('vars',{vnames_reg},'dom',dom_reg,'space_out',space_out_mask, ...
                  'space_in',space_in_mask,'dim_out',dim_out,'dim_in',dim_in);
    Pops{kk} = opvar(blocks.(objs{kk}),meta);
end
[Pop_x,Pop_w,Pop_u] = Pops{:};

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [obj,Ridx,blk] = term2sopvar(trm,PDE,vnames_out,dom_out,q)
% [OBJ,RIDX,BLK] = TERM2SOPVAR(TRM,PDE,VNAMES_OUT,DOM_OUT,Q) returns the
% 'sopvar' object BLK representing the term TRM as an operator on the
% component PDE.(OBJ){RIDX}, into the space of variables VNAMES_OUT.
%
% INPUT
% - trm:        'struct', a term of a free equation of a 'pde_struct';
% - PDE:        'pde_struct' object holding the states and inputs;
% - vnames_out: cellstr, the variables of the output space;
% - dom_out:    numel(vnames_out) x 2 array, their domains;
% - q:          integer, the number of rows of the free equation;
%
% OUTPUT
% - obj:        'char', 'x', 'w' or 'u';
% - Ridx:       integer, the index of the component in PDE.(obj);
% - blk:        q x p 'sopvar' object, with p the size of the component.

% Check that the term is a PI operator on a single state or input: a term
% with several factors (a product of states or inputs) is nonlinear.
if numel(trm)>1
    error("'get_op' is not supported for nonlinear terms.")
elseif is_LHS_term(trm)
    error("Outputs and temporal derivatives of states have no operator representation.")
end
obj = '';
for obj_k = {'x','w','u'}
    if isfield(trm,obj_k{1}) && ~isempty(trm.(obj_k{1}))
        obj = obj_k{1};
    end
end
Ridx = trm.(obj);
comp = PDE.(obj){Ridx};
vnames_in = arrayfun(@(k) comp.vars(k,1).varname{1},1:size(comp.vars,1),'UniformOutput',false);
dom_in = comp.dom;
N = numel(vnames_in);
if isfield(trm,'D') && any(trm.D)
    error("Derivatives of states are not PI operators on the states.")
elseif isfield(trm,'delay') && ~isempty(trm.delay) && any(double(trm.delay)~=0)
    error("Delayed states are not supported.")
end

% Integration domain and kernel variable along each variable of the
% component: multiplier [s,s], partial integral [a,s] or [s,b], or full
% integral [a,b]; the dummy variable becomes s_dum for a variable of the
% output space, and s itself for an input-only variable.
if ~isfield(trm,'C') || isempty(trm.C),   C = polynomial(eye(q,comp.size));
else,                                     C = polynomial(trm.C);   end
if all(size(C)==1) && comp.size>1
    C = C*eye(comp.size);
end
int_dom = cell(N,2);
for k=1:N
    vname = vnames_in{k};
    is_out = ismember(vname,vnames_out);
    if isfield(trm,'loc') && length(trm.loc)>=k
        loc_k = polynomial(trm.loc(k));
    else
        loc_k = pvar(vname);
    end
    if isdouble(loc_k)
        error("Evaluation at a position is not a PI operator on the states.")
    end
    has_int = isfield(trm,'I') && numel(trm.I)>=k && ~isempty(trm.I{k});
    if ~has_int
        if ~strcmp(loc_k.varname{1},vname)
            error("The term evaluates the state or input at '%s' in place of '%s', without an integral over '%s'.",loc_k.varname{1},vname,loc_k.varname{1})
        elseif ~is_out
            error("The term depends on variable '%s', which is not a variable of the output space.",vname)
        end
        int_dom(k,:) = {vname,vname};
        continue
    end
    lims = polynomial(trm.I{k});
    if isdouble(lims)
        int_dom(k,:) = {double(lims(1)),double(lims(2))};
    elseif ispvar(lims(2))
        int_dom(k,:) = {double(lims(1)),vname};
    else
        int_dom(k,:) = {vname,double(lims(2))};
    end
    if ~isdouble(lims) && ~is_out
        error("A partial integral in '%s' requires '%s' in the output space.",vname,vname)
    end
    % Rename the dummy variable of integration to the kernel convention.
    if ~strcmp(loc_k.varname{1},vname)
        if is_out
            C = subs(C,loc_k,pvar([vname,'_dum']));
        else
            C = subs(C,loc_k,pvar(vname));
        end
    end
end
if size(C,1)~=q || size(C,2)~=comp.size
    error("The coefficients of a term should be a %d x %d matrix.",q,comp.size)
end
% The coefficients may depend only on the variables of the output space,
% the dummy variables of integration of those, and the variables of
% integration of the component that are not in the output space.
vnames_int = vnames_in(~strcmp(int_dom(:,1)',int_dom(:,2)'));               % not multipliers
vnames_ok = [vnames_out, strcat(intersect(vnames_int,vnames_out),'_dum'), setdiff(vnames_int,vnames_out)];
C = polynomial(C);
vnames_C = setdiff(C.varname,vnames_ok);
if ~isempty(vnames_C)
    error("The coefficients of a term depend on '%s', which is neither a variable of the output space nor a variable of integration.",vnames_C{1})
end
vars_blk = struct('in',{vnames_in},'out',{vnames_out});
dom_blk = struct('in',dom_in,'out',dom_out);
blk = int2sopvar(C,vars_blk,dom_blk,int_dom);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [vnames_reg,dom_reg] = add_vars(vnames_reg,dom_reg,comp)
% [VNAMES_REG,DOM_REG] = ADD_VARS(VNAMES_REG,DOM_REG,COMP) adds the spatial
% variables of a PDE component COMP, with their domains, to the registry,
% checking that a variable has a single domain.
%
% INPUT
% - vnames_reg: 1 x nv cellstr, the registered variable names;
% - dom_reg:    nv x 2 array, their domains;
% - comp:       'struct', a component of a 'pde_struct', with fields vars
%               and dom;
%
% OUTPUT
% - vnames_reg, dom_reg: the updated registry.

if ~isfield(comp,'vars') || isempty(comp.vars)
    return
end
vnames = arrayfun(@(k) comp.vars(k,1).varname{1},1:size(comp.vars,1),'UniformOutput',false);
for k=1:numel(vnames)
    [is_reg,idx] = ismember(vnames{k},vnames_reg);
    if ~is_reg
        vnames_reg = [vnames_reg, vnames(k)];                               %#ok<AGROW>
        dom_reg = [dom_reg; comp.dom(k,:)];                                 %#ok<AGROW>
    elseif any(dom_reg(idx,:)~=comp.dom(k,:))
        error("Variable '%s' is defined on different domains.",vnames{k})
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
