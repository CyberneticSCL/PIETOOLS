function [PIE,tdiff_tab] = expand_tderivatives(PIE,tdiff_list,suppress_summary)
% [PIE,TDIFF_TAB] = EXPAND_TDERIVATIVES(PIE,TDIFF_LIST,SUPPRESS_SUMMARY)
% takes a nonlinear PIE with evolution equations of higher order in time,
% and adds components for the temporal derivatives of the state, so that the
% PIE is first order in time. This works for any number of spatial
% variables, and for PIEs with a constraint K*x=0, as it operates on the
% 'tensopmat' operators and 'polyopvar' polynomials with a block for each
% component of the state. See also `pie_struct/expand_tderivatives`.
%
% INPUT
% - PIE:        'npie_struct' object with nb state components x_i, with
%               evolution equations of order p_i in time,
%                   (d/dt)^{p_i} (T_i*x + Tw_i*w + Tu_i*u)
%                       = f_i(x) + Bw_i*w + Bu_i*u,
%               for T_i the i-th row of T, f_i the i-th row of the
%               polynomial f, and similarly for the other operators, and a
%               constraint K*x=0;
% - tdiff_list: nb x 1 array of positive integers, the orders p_i of the
%               temporal derivatives in the equations;
% - suppress_summary:   (optional) logical, true to suppress the display
%               of how the state components in the output PIE relate to
%               those in the input PIE (false by default);
%
% OUTPUT
% - PIE:        'npie_struct' object, first order in time, describing the
%               same system with the state components of the input and
%               their temporal derivatives (see NOTES);
% - tdiff_tab:  nb_new x 2 array, specifying for each state component of
%               the output PIE what component of the input PIE it is the
%               temporal derivative of: component i is
%                   x_i = (d/dt)^{tdiff_tab(i,2)}*
%                         x_{old,tdiff_tab(i,1)};
%
% NOTES
% - A component x_j gets L_j derivatives x_{j,l}=(d/dt)^l x_j, l=1,...,L_j,
%   for L_j the highest p_i-1 among the equations with p_i>1 in which x_j
%   appears (T_{i,j}~=0). The components x_{j,0}=x_j keep their position and
%   name, and the derivatives are added after them, ordered by j and l
%   (unlike for a 'pie_struct', as the names of components must be sorted).
% - The equation of x_j (the rows of T, f, Tw, Tu, Bw, Bu) moves to x_{j,L_j},
%   with T_{i,c} of an equation of order p_i acting on x_{c,p_i-1}. For
%   l<L_j, the equation (d/dt) x_{j,l} = x_{j,l+1} is added: the identity in
%   T, and the linear monomial x_{j,l+1} in f.
% - T*x need no longer be the PDE state, so pie2pde.x is adjusted for
%   x_pde=X_pde*(T*x+Tw*w+Tu*u) to hold: x_pde=sum_i X_old(:,i)*(T_old*x)_i,
%   with (T_old*x)_i the new row of T*x+Tw*w+Tu*u if p_i=1, and
%   sum_c T_old(i,c)*x_{c,0} otherwise.
% - The new components are named x<k>_pie for their index k (zero-padded as
%   in `initialize`), which does not show what they derive from: see
%   tdiff_tab, pde2pie.x.idx and pde2pie.x.tdiff. If the names would not be
%   unique and sorted (e.g. x10_pie < x2_pie), all components are renamed
%   x01_pie, ..., with a warning.
% - Each constraint of K is also imposed on the derivatives of its
%   components, acting on x_{c,l}. These components must have the same L_c.
% - Not supported (error): p_i>1 with inputs in T (boundary inputs, as these
%   need derivatives of the inputs), inputs in the constraint (Kw, Ku), and
%   a map pde2pie.x.C with blocks outside the diagonal.
% - The maps pde2pie of a derivative are those of its parent, with the order
%   of the temporal derivative increased.

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
% Initial coding DJ (CC), 11/09/2026, based on the pie_struct version;

% % % Check the inputs
if nargin<2
    error("Insufficient input arguments.")
elseif nargin==2
    suppress_summary = false;
elseif nargin>3
    error("Too many input arguments.")
end
if ~isa(PIE,'npie_struct')
    error("The PIE must be specified as an 'npie_struct' object.")
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
% Blocks of an operator that are nonzero (a zero operator may be a nonempty
% block)
nz_blocks = @(P) reshape(cellfun(@is_nz_blk,P.ops(:)),size(P.ops));
% The zero operator with the spaces of the rows "rows" (columns "cols") of P,
% in which the blocks are then set. It has a (zero) block in each row and
% column, so that the spaces of the rows and columns are not ambiguous
zero_op = @(P,rows,cols) mat2tensopmat(sparse(sum(P.dims.out(rows)),sum(P.dims.in(cols))),...
            struct('in',P.dims.in(cols),'out',P.dims.out(rows)),P.vars,P.dom,...
            struct('in',P.depmat.in(cols,:),'out',P.depmat.out(rows,:)));
% The identity on component j of the state, as a block of a 'tensopmat',
% with the variables of the function space sorted as for any operator
space = PIE.space;
[vsort,vidx] = sort(space.vars);
dom_sort = space.dom(vidx,:);
dep_sort = space.x.depmat(:,vidx);
identity_op = @(j) tensopvar(mat2sopvar(eye(space.x.dims(j)),...
            struct('in',{vsort(dep_sort(j,:))},'out',{vsort(dep_sort(j,:))}),...
            struct('in',dom_sort(dep_sort(j,:),:),'out',dom_sort(dep_sort(j,:),:))));
is_exp = tdiff_list>1;
T_nz = nz_blocks(PIE.T);
Tin_nz = [nz_blocks(PIE.Tw),nz_blocks(PIE.Tu)];
if any(any(Tin_nz(is_exp,:)))
    error("expand_tderivatives:boundaryInputs","Expansion of higher-order temporal derivatives is not supported for states with boundary inputs.")
end
if any(any([nz_blocks(PIE.Kw),nz_blocks(PIE.Ku)]))
    error("expand_tderivatives:constraintInputs","Expansion of higher-order temporal derivatives is not supported for constraints involving inputs.")
end
Cmap = PIE.pde2pie.x.C;
if size(Cmap.ops,1)==nb
    C_nz = nz_blocks(Cmap);
    C_nz(logical(eye(nb))) = false;
    if any(C_nz(:))
        error("expand_tderivatives:mapOffDiagonal","Expansion of higher-order temporal derivatives is not supported for a map pde2pie.x.C with blocks outside the diagonal.")
    end
end


% % % Determine the components, and the new components for the temporal
% % % derivatives
% L(j): number of derivatives of component j; the components x_pie_{j,0}=
% x_pie_j are the first nb components, followed by the derivatives, ordered
% per component j, parent_idx and tdiff_ord the component and order of each.
% ix(j,l) is the new index of x_pie_{j,l}, and host(j) that of the last
% derivative x_pie_{j,L(j)}, which gets the equation of x_pie_j
L = zeros(nb,1);
for i=find(is_exp)'
    cols = T_nz(i,:)';
    L(cols) = max(L(cols),tdiff_list(i)-1);
end
nnew = sum(L+1);
parent_idx = [(1:nb)'; repelem((1:nb)',L,1)];
tdiff_ord = [zeros(nb,1); cell2mat(arrayfun(@(l) (1:l)',L,'UniformOutput',false))];
dstart = nb+cumsum([0;L(1:end-1)]);
ix = @(c,l) c.*(l==0) + (dstart(c)+l).*(l>0);
host = ix((1:nb)',L);
tdiff_tab = [parent_idx,tdiff_ord];
% Names of the new components: those of the old components for x_pie_{c,0},
% and the default names assigned by initialize for the added components.
% The names must be unique and sorted in the order of the components, see
% initialize: if not (e.g. if the added names have more digits than the old
% ones), all components are renamed with the default names
names_old = space.x.varname;
names_new = [names_old; arrayfun(@(k) ['x',num2str(k,['%0',num2str(numel(num2str(nnew))),'d']),'_pie'],(nb+1:nnew)','UniformOutput',false)];
if ~issorted(names_new) || numel(unique(names_new))<nnew
    warning("expand_tderivatives:namesNotSorted","The names of the new components cannot be unique and sorted after those of the old components; all components are renamed with the default names.")
    names_new = arrayfun(@(k) ['x',num2str(k,['%0',num2str(numel(num2str(nnew))),'d']),'_pie'],(1:nnew)','UniformOutput',false);
end


% % % The operator T: the equation of x_pie_i moves to host(i), where T acts
% % % on x_pie_{c,p_i-1}, and the derivatives get the equations
% % % (d/dt) x_pie_{j,l} = x_pie_{j,l+1}
T_old = PIE.T.ops;
Tn = zero_op(PIE.T,parent_idx,parent_idx);
ops = Tn.ops;
for i=1:nb
    tdiff_i = 0;
    if is_exp(i)
        tdiff_i = tdiff_list(i)-1;
    end
    for c=find(T_nz(i,:))
        ops{host(i),ix(c,tdiff_i)} = T_old{i,c};
    end
    for l=0:L(i)-1
        ops{ix(i,l),ix(i,l)} = identity_op(i);
    end
end
Tn.ops = ops;       Tn.chkval = true;
PIE.T = Tn;
% The rows of the operators on the inputs move to host(i) as well. The
% operators Dzw, Dzu, Dyw, Dyu do not involve the state
for fname = {'Tw','Tu','Bw','Bu'}
    P = PIE.(fname{1});
    Pn = zero_op(P,parent_idx,1:size(P.ops,2));
    ops = Pn.ops;
    ops(host,:) = P.ops;
    Pn.ops = ops;       Pn.chkval = true;
    PIE.(fname{1}) = Pn;
end


% % % The polynomials: the state variables are renamed, with the
% % % components x_pie_c of the input becoming x_pie_{c,0}, and the new
% % % components added to the variables. The row of f for
% % % x_pie_i moves to host(i), and f gets the linear monomials x_pie_{j,l+1}
% % % in the rows of x_pie_{j,l}, l<L(j), the equations of the derivatives.
space.x.dims = space.x.dims(parent_idx);
space.x.depmat = space.x.depmat(parent_idx,:);
space.x.varname = names_new;
for fname = {'f','gz','gy'}
    Pol = PIE.(fname{1});
    % Index of the new component x_pie_{c,0} for each monomial factor
    [~,loc] = ismember(Pol.varname,names_old);
    [Zr,Zc,Zv] = find(Pol.Z);
    Z_new = sparse(Zr(:),Zc(:),loc(Zv(:)),size(Pol.Z,1),max(size(Pol.Z,2),1));
    Cn = Pol.C;
    if strcmp(fname{1},'f')
        ops_old = Pol.C.ops;
        nZ = size(ops_old,2);
        ch = setdiff((1:nnew)',host);       % rows with the equation of a derivative
        ops = cell(nnew,nZ+numel(ch));
        ops(host,1:nZ) = ops_old;
        for m=1:numel(ch)
            ops{ch(m),nZ+m} = identity_op(parent_idx(ch(m)));
        end
        Z_new = [Z_new; sparse(1:numel(ch),1,ix(parent_idx(ch),tdiff_ord(ch)+1),numel(ch),size(Z_new,2))];
        Cn = tensopmat(ops);
    end
    Pol.C = Cn;
    Pol.Z = Z_new;
    Pol.varname = names_new;
    Pol.dims.in = space.x.dims;
    Pol.vars.in = vsort;
    Pol.dom.in = dom_sort;
    Pol.depmat.in = space.x.depmat(:,vidx);
    PIE.(fname{1}) = initialize(Pol);
end


% % % The constraint K*x_pie=0 holds for each of the derivatives, with the
% % % components of a constraint having the same number of derivatives
K = PIE.K;
nK = size(K.ops,1);
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
    Kn = zero_op(K,K_rows(:,1),parent_idx);
    ops = Kn.ops;
    for q=1:size(K_rows,1)
        for c=find(K_nz(K_rows(q,1),:))
            ops{q,ix(c,K_rows(q,2))} = K.ops{K_rows(q,1),c};
        end
    end
    Kn.ops = ops;       Kn.chkval = true;
    PIE.K = Kn;
else
    PIE.K = [];
end
PIE.Kw = [];        PIE.Ku = [];


% % % The function spaces, and the maps between the PDE and PIE: the new
% % % components are maps from the same PDE component, with a temporal
% % % derivative
pde2pie_s = PIE.pde2pie;
for fname = {'idx','diff','loc','delay'}
    if size(pde2pie_s.x.(fname{1}),1)==nb
        pde2pie_s.x.(fname{1}) = pde2pie_s.x.(fname{1})(parent_idx,:);
    end
end
if size(pde2pie_s.x.tdiff,1)==nb
    pde2pie_s.x.tdiff = pde2pie_s.x.tdiff(parent_idx) + tdiff_ord;
end
if size(Cmap.ops,1)==nb
    Cn = zero_op(Cmap,parent_idx,parent_idx);
    ops = Cn.ops;
    for k=1:nnew
        if ~isempty(Cmap.ops{parent_idx(k),parent_idx(k)})
            ops{k,k} = Cmap.ops{parent_idx(k),parent_idx(k)};
        end
    end
    Cn.ops = ops;       Cn.chkval = true;
    pde2pie_s.x.C = Cn;
end
% The PDE state from the new T*x_pie+Tw*w+Tu*u, see NOTES
pie2pde_s = PIE.pie2pde;
X_old = pie2pde_s.x;
if size(X_old.ops,2)==nb
    nX = size(X_old.ops,1);
    X_nz = nz_blocks(X_old);
    X_acc = cell(nX,nnew);
    for i=1:nb
        if is_exp(i)
            cols = find(T_nz(i,:));       targets = cols;     ops_T = T_old(i,cols);
        else
            targets = host(i);              ops_T = {[]};
        end
        for k=find(X_nz(:,i))'
            for m=1:numel(targets)
                P = X_old.ops{k,i};
                if ~isempty(ops_T{m})
                    P = P*ops_T{m};
                end
                if isempty(X_acc{k,targets(m)})
                    X_acc{k,targets(m)} = P;
                else
                    X_acc{k,targets(m)} = X_acc{k,targets(m)} + P;
                end
            end
        end
    end
    % Blocks that are not set stay empty (zero), unlike for the other
    % operators: a product of X_pde with a polynomial, e.g. f_pde = X_pde*f,
    % composes all nonempty blocks, and the composition of a zero block with
    % the monomials of f of degree>1 is not supported
    Xn = tensopmat();
    Xn.ops = X_acc;
    Xn.dims = struct('in',X_old.dims.in(parent_idx),'out',X_old.dims.out);
    Xn.vars = X_old.vars;       Xn.dom = X_old.dom;
    Xn.depmat = struct('in',X_old.depmat.in(parent_idx,:),'out',X_old.depmat.out);
    Xn.chkval = true;
    pie2pde_s.x = Xn;
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
function is_nz = is_nz_blk(B)
% IS_NZ = IS_NZ_BLK(B) checks whether a block of a 'tensopmat' is nonzero.
%
% INPUT
% - B:      block of a 'tensopmat': empty, a 'tensopvar' or 'intop' object,
%           or a matrix of type 'double', 'polynomial' or 'dpvar';
%
% OUTPUT
% - is_nz:  logical, false if the block is empty, or zero up to a tolerance
%           of 1e-12, true otherwise;
if isempty(B)
    is_nz = false;
elseif isa(B,'tensopvar') || isa(B,'intop')
    [~,is_zero] = clean_op(B,1e-12);
    is_nz = ~is_zero;
else
    [~,is_zero] = clean_dpvar(B,1e-12);
    is_nz = ~is_zero;
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
