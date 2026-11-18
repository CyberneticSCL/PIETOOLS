function [PIE,tdiff_tab] = expand_tderivatives(PIE,tdiff_list,suppress_summary)
% [PIE,TDIFF_TAB] = EXPAND_TDERIVATIVES(PIE,TDIFF_LIST,SUPPRESS_SUMMARY)
% takes a PIE with evolution equations of higher order in time, and adds
% components for the temporal derivatives of the state, so that the PIE is
% first order in time. This works for any number of spatial variables, and
% for PIEs with a constraint K*x=0, as it operates on the 'copvar'
% operators directly, with a block for each component of the state.
%
% INPUT
% - PIE:        'pie_struct' object with nb state components x_i, with
%               evolution equations of order p_i in time,
%                   (d/dt)^{p_i} (T_i*x + Tw_i*w + Tu_i*u)
%                       = A_i*x + Bw_i*w + Bu_i*u,
%               for T_i the i-th row of T, and similarly for the other
%               operators, and a constraint K*x=0;
% - tdiff_list: nb x 1 array of positive integers, the orders p_i of the
%               temporal derivatives in the equations;
% - suppress_summary:   (optional) logical, true to suppress the display
%               of how the state components in the output PIE relate to
%               those in the input PIE (false by default);
%
% OUTPUT
% - PIE:        'pie_struct' object, first order in time, describing the
%               same system with the state components of the input and
%               their temporal derivatives (see NOTES);
% - tdiff_tab:  nb_new x 2 array, specifying for each state component of
%               the output PIE what component of the input PIE it is the
%               temporal derivative of: component i is
%                   x_i = (d/dt)^{tdiff_tab(i,2)}*
%                         x_{old,tdiff_tab(i,1)};
%
% NOTES
% - A component x_j gets L_j temporal derivatives x_{j,l}=(d/dt)^l x_j,
%   l=1,...,L_j, for L_j the highest value of p_i-1 among the equations
%   with p_i>1 in which x_j appears (T_{i,j}~=0). The new components follow
%   the component they derive from, x_{j,0}=x_j.
% - The equation of x_j is moved to the last component x_{j,L_j}, with the
%   operators T_{i,c} of an equation of order p_i acting on x_{c,p_i-1}, so
%   that it is first order. For l<L_j, the new equations are
%   (d/dt) x_{j,l} = x_{j,l+1}.
%   As a result, T*x need no longer be the PDE state, and pie2pde.x is
%   adjusted, so that x_pde=X_pde*(T*x+Tw*w+Tu*u) still holds: with the old
%   operators X_old and T_old, the PDE state is x_pde=sum_i X_old(:,i)*
%   (T_old*x)_i, where (T_old*x)_i is the i-th row of the new
%   T*x+Tw*w+Tu*u if p_i=1, and sum_c T_old(i,c)*x_{c,0} otherwise.
% - The derivatives satisfy the constraint of x. For each constraint of K,
%   a constraint is added for each of the derivatives of the components in
%   it, acting on x_{c,l} for the same l. An error is thrown if a
%   constraint involves components with different numbers of derivatives.
% - The expansion is not supported for equations with p_i>1 that involve
%   inputs in T (boundary inputs), as these would require derivatives of
%   the inputs, or for constraints involving inputs, Kw and Ku.
% - The maps pde2pie of the new components are those of the component they
%   derive from, with the order of the temporal derivative increased.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - expand_tderivatives
%
% Copyright (C)2026  PIETOOLS Team
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
% Initial coding DJ, 06/18/2025;
% DJ, 10/29/2026: Update to new pie_struct structure (expansion of copvar
%                   operators, extension of the maps between PDE and PIE,
%                   and no support for PIEs with an algebraic constraint);
% DJ, 11/02/2026: Major overhaul to support ND PIEs and algebraic 
%                   constraints;

% % % Check the inputs
if nargin<2
    error("Insufficient input arguments.")
elseif nargin==2
    suppress_summary = false;
elseif nargin>3
    error("Too many input arguments.")
end
if ~isa(PIE,'pie_struct')
    error("The PIE must be specified as a 'pie_struct' object.")
end
if ~PIE.is_initialized
    PIE = initialize(PIE);
end
nb = numel(PIE.space.x.dims);
tdiff_list = tdiff_list(:);
if numel(tdiff_list)~=nb || any(tdiff_list<1) || any(tdiff_list~=round(tdiff_list))
    error("The orders of the temporal derivatives should be specified as an array of positive integers, one for each state component.")
end
if all(tdiff_list==1)
    tdiff_tab = [(1:nb)',zeros(nb,1)];
    if ~suppress_summary
        fprintf(['\n','No higher-order temporal derivatives were encountered; returning original PIE'])
    end
    return
end
% Blocks of an operator that are nonzero, and the identity on component j
has_op = @(B) ~isempty(B) && any(cellfun(@nnz,B.params(:)));
nz_blocks = @(P) cellfun(has_op,P.C);
space = PIE.space;
identity_op = @(j) mat2sopvar(eye(space.x.dims(j)),struct('in',{space.vars(space.x.depmat(j,:))},'out',{space.vars(space.x.depmat(j,:))}),struct('in',space.dom(space.x.depmat(j,:),:),'out',space.dom(space.x.depmat(j,:),:)));
% Metadata of the operator P with the rows (columns) of its spaces taken
% from the rows row_idcs (columns col_idcs) of P
meta_subset = @(P,row_idcs,col_idcs) struct('vars',{P.vars},'dom',P.dom,'space_out',P.space_out(row_idcs,:),'space_in',P.space_in(col_idcs,:),'dim_out',P.dim_out(row_idcs),'dim_in',P.dim_in(col_idcs));
is_exp = tdiff_list>1;
T_nz = nz_blocks(PIE.T);
A_nz = nz_blocks(PIE.A);
Tin_nz = [nz_blocks(PIE.Tw),nz_blocks(PIE.Tu)];
if any(any(Tin_nz(is_exp,:)))
    error("expand_tderivatives:boundaryInputs","Expansion of higher-order temporal derivatives is not supported for states with boundary inputs.")
end
if any(any([nz_blocks(PIE.Kw),nz_blocks(PIE.Ku)]))
    error("expand_tderivatives:constraintInputs","Expansion of higher-order temporal derivatives is not supported for constraints involving inputs.")
end


% % % Determine the components, and the new components for the temporal
% % % derivatives
% L(j): number of derivatives of component j; the new components are
% ordered per component j, with the derivatives of order 0,...,L(j),
% parent_idx and tdiff_ord the component and order of each, first(j) the
% new index of x_pie_j, and host(j) that of the last derivative, which
% gets the equation of x_pie_j
L = zeros(nb,1);
for i=find(is_exp)'
    cols = T_nz(i,:)';
    L(cols) = max(L(cols),tdiff_list(i)-1);
end
nnew = sum(L+1);
parent_idx = repelem((1:nb)',L+1,1);
tdiff_ord = cell2mat(arrayfun(@(l) (0:l)',L,'UniformOutput',false));
first = cumsum([1;L(1:end-1)+1]);
host = first+L;
tdiff_tab = [parent_idx,tdiff_ord];


% % % The operators T and A: the equation of x_pie_i moves to host(i),
% % % where T acts on x_pie_{c,p_i-1}, and the derivatives get the
% % % equations (d/dt) x_pie_{j,l} = x_pie_{j,l+1}
T_old = PIE.T.C;       A_old = PIE.A.C;
T_new = cell(nnew,nnew);       A_new = cell(nnew,nnew);
% Row i: move the equation of x_pie_i to host(i), add the chain equations
for i=1:nb
    tdiff_i = 0;
    if is_exp(i)
        tdiff_i = tdiff_list(i)-1;
    end
    for c=find(T_nz(i,:))
        T_new{host(i),first(c)+tdiff_i} = T_old{i,c};
    end
    for c=find(A_nz(i,:))
        A_new{host(i),first(c)} = A_old{i,c};
    end
    for l=0:L(i)-1
        T_new{first(i)+l,first(i)+l} = identity_op(i);
        A_new{first(i)+l,first(i)+l+1} = identity_op(i);
    end
end
PIE.T = copvar(T_new,meta_subset(PIE.T,parent_idx,parent_idx));
PIE.A = copvar(A_new,meta_subset(PIE.A,parent_idx,parent_idx));
% The rows of the operators on the inputs move to host(i) as well, and the
% columns of the operators on the state are those of x_pie_{c,0}
for fname = {'Tw','Tu','Bw','Bu'}
    P = PIE.(fname{1});
    blocks_new = cell(nnew,size(P.C,2));
    blocks_new(host,:) = P.C;
    PIE.(fname{1}) = copvar(blocks_new,meta_subset(P,parent_idx,1:size(P.C,2)));
end
for fname = {'Cz','Cy'}
    P = PIE.(fname{1});
    blocks_new = cell(size(P.C,1),nnew);
    blocks_new(:,first) = P.C;
    PIE.(fname{1}) = copvar(blocks_new,meta_subset(P,1:size(P.C,1),parent_idx));
end


% % % The constraint K*x_pie=0 holds for each of the derivatives, with the
% % % components of a constraint having the same number of derivatives
K = PIE.K;
nK = size(K.C,1);
if nK>0
    K_nz = nz_blocks(K);
    L_con = zeros(nK,1);
    for r=1:nK
        cols = find(K_nz(r,:));
        if isempty(cols)
            continue
        end
        if any(L(cols)~=L(cols(1)))
            error("expand_tderivatives:constraintOrders","A constraint on the state involves components with different numbers of temporal derivatives; this is not supported.")
        end
        L_con(r) = L(cols(1));
    end
    % The constraints (row of K, order of the derivative), by order
    K_rows = [repelem((1:nK)',L_con+1,1),cell2mat(arrayfun(@(l) (0:l)',L_con,'UniformOutput',false))];
    K_rows = sortrows(K_rows,[2,1]);
    K_new = cell(size(K_rows,1),nnew);
    for q=1:size(K_rows,1)
        for c=find(K_nz(K_rows(q,1),:))
            K_new{q,first(c)+K_rows(q,2)} = K.C{K_rows(q,1),c};
        end
    end
    PIE.K = copvar(K_new,meta_subset(K,K_rows(:,1),parent_idx));
else
    PIE.K = [];
end
PIE.Kw = [];        PIE.Ku = [];


% % % The function spaces, and the maps between the PDE and PIE: the new
% % % components are maps from the same PDE component, with a temporal
% % % derivative
space.x.dims = space.x.dims(parent_idx);
space.x.depmat = space.x.depmat(parent_idx,:);
pde2pie_s = PIE.pde2pie;
for fname = {'idx','diff','loc','delay'}
    if size(pde2pie_s.x.(fname{1}),1)==nb
        pde2pie_s.x.(fname{1}) = pde2pie_s.x.(fname{1})(parent_idx,:);
    end
end
if size(pde2pie_s.x.tdiff,1)==nb
    pde2pie_s.x.tdiff = pde2pie_s.x.tdiff(parent_idx) + tdiff_ord;
end
if isa(pde2pie_s.x.C,'copvar') && size(pde2pie_s.x.C.C,1)==nb
    blocks_new = cell(nnew,nnew);
    for k=1:nnew
        blocks_new{k,k} = pde2pie_s.x.C.C{parent_idx(k),parent_idx(k)};
    end
    pde2pie_s.x.C = copvar(blocks_new,meta_subset(pde2pie_s.x.C,parent_idx,parent_idx));
end
% The PDE state from the new T*x_pie+Tw*w+Tu*u, see NOTES
pie2pde_s = PIE.pie2pde;
X_old = pie2pde_s.x;
if isa(X_old,'copvar') && size(X_old.C,2)==nb
    X_nz = nz_blocks(X_old);
    X_new = cell(size(X_old.C,1),nnew);
    for i=1:nb
        if is_exp(i)
            cols = find(T_nz(i,:));       targets = first(cols)';     ops = T_old(i,cols);
        else
            targets = host(i);              ops = {[]};
        end
        for k=find(X_nz(:,i))'
            for m=1:numel(targets)
                P = X_old.C{k,i};
                if ~isempty(ops{m})
                    P = P*ops{m};
                end
                if isempty(X_new{k,targets(m)})
                    X_new{k,targets(m)} = P;
                else
                    X_new{k,targets(m)} = X_new{k,targets(m)} + P;
                end
            end
        end
    end
    pie2pde_s.x = copvar(X_new,meta_subset(X_old,1:size(X_old.C,1),parent_idx));
end
PIE.space = space;
PIE.pde2pie = pde2pie_s;
PIE.pie2pde = pie2pde_s;
PIE = initialize(PIE);
if ~suppress_summary
    print_tdiff_expansion_summary(PIE,tdiff_tab);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function print_tdiff_expansion_summary(PIE,tdiff_tab)
% PRINT_TDIFF_EXPANSION_SUMMARY(PIE,TDIFF_TAB) prints in the command window
% some information concerning how the state variables in the PIE structure
% 'PIE' relate to the state variables in the
% original PIE specified by the user in the call to "expand_tderivatives".
%
% INPUTS:
% - PIE:    A "pie_struct" class object defining a PIE.
% - tdiff_tab:  nx_new x 2 array of integers, indicating for each of the
%               j=1:nx_new state components what component i in the
%               original PIE it corresponds to, and what order temporal
%               derivative of that component:
%                   x_j = (d/dt)^{tdiff_tab(j,2)} x_{old,tdiff_tab(j,1)}.
%
% OUTPUTS:
% Displays information in the command window concerning how the state
% variables in the "PIE" structure relate to the fundamental state
% variables in the original PIE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncomps = size(tdiff_tab,1);
fprintf(['\n','Introduced the following (vector-valued) fundamental state components:\n'])

% Use UNICODE to add subscript indices to different components.
thin_space = char(8201);
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
partial = '\x2202';
sub_t = '\x209C';
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
nvars = numel(PIE.space.vars);
global_vars_obj = PIE.space.vars;
lngth_varnames = 3;      % just three characters: (t)
if ~isempty(global_vars_obj)
    varname_cell = cell(nvars,1);
    if nvars==1
        varname_cell{1} = 's';
        lngth_varnames = lngth_varnames+2;    % add two characters: ,s
    else
        varname_cell{1} = ['s',sub_num{2}];
        lngth_varnames = lngth_varnames+2;    % add three characters: ,s1
        for kk=2:nvars
            varname_cell{kk} = ['s',cell2mat(sub_num(str2num(num2str(kk)')+1)')];
            lngth_varnames = lngth_varnames+2;    % add three characters: ,sj
        end
    end
end
% Then, estimate the size of the string of characters denoting the components
% "PDE.x{ii}".
LHS_length_max = 1+1 + n_digits + lngth_varnames+2; % e.g. ' x13(t,s1,s2,s3)', 


% % For each of the state variables in the PIE structure, display
% % which of the state variables in the original PIE it corresponds to,
% % and what order temporal derivative.
for ii=1:ncomps
    Lstate_num = ii;   % State index associated to the new state.
    
    % Establish the names of the variables on which the component depends.
    varnames_ii_t = 't';
    var_idcs = find(PIE.space.x.depmat(ii,:));
    for kk=var_idcs
        varnames_ii_t = [varnames_ii_t,',',varname_cell{kk}];
    end
    
    % Establish the (subscript) index for the component.
    LHS_length = length(varnames_ii_t);
    if ncomps<=9
        % The state number consists of a single decimal.
        Lcomp_idx = sub_num{Lstate_num+1};
        LHS_length = LHS_length + 1;
    else
        % The state number consists of multiple decimals.
        Lcomp_idx = cell2mat(sub_num(str2num(num2str(Lstate_num)')+1)');
        LHS_length = LHS_length + length(num2str(Lstate_num));
    end
    % Set the name of the component, including its dependence on spatial
    % variables.
    LHS_name = [' x',Lcomp_idx,'(',varnames_ii_t,')'];
    LHS_length = 1 + LHS_length + 3;
        
    % For the added state components, indicate of which state component
    % they are the temporal derivative.
    Rstate_num = tdiff_tab(ii,1);
    Rcomp_idx = cell2mat(sub_num(str2num(num2str(Rstate_num)')+1)');

    % Also indicate the order of the derivative
    if tdiff_tab(ii,2)==0
        tdiff_str = ['   ',thin_space];
    elseif tdiff_tab(ii,2)==1
        tdiff_str = [' ',partial,sub_t,' '];
    else
        t_sup = cell2mat(sup_num(str2num(num2str(tdiff_tab(ii,2))')+1)');
        tdiff_str = [partial,sub_t,t_sup,' '];
    end
    RHS = [' <-- ',tdiff_str,'x',Rcomp_idx,'(',varnames_ii_t,')'];
    
    % % % Finally, display
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),RHS,'\n']);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %