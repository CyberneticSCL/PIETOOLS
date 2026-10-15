function Zop = monopmat(dims,vars,dom,depmat,degs,opts)
% ZOP = MONOPMAT(DIMS,VARS,DOM,DEPMAT,DEGS,OPTS) constructs a basis of
% tensor-PI operators in the specified variables,
%       ZOP = blkdiag(Zop{1},...,Zop{q})
% where each Zop{i} is a tensor-PI operator,
%   Zop{i}: L2^{ni}[Omega_in{i}] --> L2^{mi}[Omega_out{i}]
% for
%   Omega_in = [a1,b1]^{d_in(i,1)} x ... x [aN,bN]^{d_in(i,N)}
%   Omega_out = [a1,b1]^{d_out(i,1)} x ... x [aM,bM]^{d_out(i,M)}
%
% INPUTS
% - dims:   q x 1 array specifying the number of state variables (ni) on
%           which each of the blocks in Zop acts. May also be a struct
%           with field 'in' containing this array;
% - vars:   struct with fields
%       'in':   1 x N cellstr specifying the names of the input spatial
%               variables to the operator;
%       'out':  1 x M cellstr specifying the names of the output spatial
%               variables to the operator;
% - dom:    struct with fields
%       'in':   N x 2 array specifying the intervals on which the input
%               variables are defined;
%       'out':  M x 2 array specifying the intervals on which the output
%               variables are defined;
% - depmat: struct with fields
%       'in':   q x N integer array specifying the degree to which each of
%               the states depends on each of the input variables (d_in);
%       'out':  q x M logical array specifying whether or not each of the
%               outputs depends on each of the output variables (d_out).
%               A state may map to output variables on which its input
%               does not depend: its output is constant in these
%               variables, and only the monomials in degs.out depend on
%               them;
% - degs:   (optional) struct with fields
%       'in':   1 x q cell, specifying the degrees of the input variables
%               in each operator block. Each element degs.in{i} can be a
%               di x N cell, where di depends on depmat.in(i,:),
%               specifying the degrees of the monomial basis in each of
%               the factors of each of the input variables in block i.
%               Element degs.in{i} can also be a di x N array, specifying
%               the maximal monomial degree in each of the factors of each
%               variable. If degs.in{i} is scalar, the same degree is used
%               for all variables;
%       'out':  1 x q cell, specifying the degrees of the output variables
%               in each operator block. Each element degs.out{i} can be a
%               1 x M cell, specifying the degrees of the monomial basis
%               in each of the output variables in block i, or a 1 x M
%               array specifying the maximal monomial degree in each of
%               the variables.
%       If DEGS is not a struct, but a numeric array or a cell, the same
%       DEGS is used for 'in' and 'out'. Defaults to 0;
% - opts:   (optional) struct with fields
%       'exclude':  1 x q cell, specifying for each of the operator blocks
%               whether any terms need to be excluded. Element
%               opts.exclude{i} can be a 3^N3(i) logical array, where
%               N3(i) is the number of 3-PI variables (appearing both in
%               the input and output) of the operator block i, with
%               element k along direction j set to true to exclude the
%               multiplier term (k=1), lower integral term (k=2), or upper
%               integral term (k=3), along spatial variable j in the list
%               of 3-PI variables (intersect(vars.in,vars.out)). A single
%               logical array may also be specified, to be used for all
%               blocks. By default, no terms are excluded;
%       'sep':  1 x q cell, specifying for each of the operator blocks
%               whether to enforce separability along each of the
%               variables. Element opts.sep{i} is a 1 x N3(i) logical
%               array, specifying for each of the 3-PI variables whether
%               the operator is separable along that variable, meaning
%               that the kernel in the upper and lower integrals is the
%               same. A single logical array (or scalar) may also be
%               specified, to be used for all blocks. By default, no
%               separability is enforced;
%       'otype':  integer in {1,2} specifying whether to take the proper
%               tensor product (otype = 2) or the pointwise product
%               (otype = 1) of the 3-PI operators in constructing the
%               tensor-PI operator basis. If otype==1, then d_out =
%               depmat.out. If otype==2, then d_out = d*depmat.out, where
%               d = max(depmat.in(i,:)). Defaults to 1;
%
% OUTPUTS
% - Zop:    m x n 'tensopmat' object representing a basis of operators of
%           specified degrees on the declared function space, where n =
%           sum(dims) and where m depends on the selected monomial
%           degrees, if every block has a single factor, i.e.
%           max(depmat.in(i,:))==1 for all i;
%
% NOTES
% - A block i with dj>1 factors (dj = max(depmat.in(i,:))) is a tensor
%   product of dj operators acting on L2^{dims(i)}, and acts on
%   L2^{dims(i)^dj}. Its column dimension is dims(i)^dj, instead of
%   dims(i), and its row dimension is a product over the dj factors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - monopmat
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
% DJ, 10/15/2026: Initial coding;


%% Process the inputs
% Check that the state sizes are properly specified
if isstruct(dims) && isfield(dims,'in')
    % For monomial basis operator, only input dimensions are relevant
    dims = dims.in;
end
if ~isnumeric(dims) || any(round(abs(dims))~=dims)
    error("Input dimensions of the operator should be specified as q x 1 array of nonnegative integers.")
end
q = numel(dims);

% Check that the variables are properly specified
if ~isstruct(vars) || ~isfield(vars,'in') || ~isfield(vars,'out')
    error("Variables to use in the operator should be specified as a struct with fields 'in' and 'out'.")
end
vars_in = vars.in;          vars_out = vars.out;
if isa(vars_in,'polynomial')
    vars_in = pvar2varname(vars_in);
elseif ~iscellstr(vars_in)
    error("Input variable names should be specified as 1 x N cellstr array.")
end
if isa(vars_out,'polynomial')
    vars_out = pvar2varname(vars_out);
elseif ~iscellstr(vars_out)
    error("Output variable names should be specified as 1 x M cellstr array.")
end
vars_in = reshape(vars_in,1,[]);    vars_out = reshape(vars_out,1,[]);
if numel(unique(vars_in))~=numel(vars_in) || numel(unique(vars_out))~=numel(vars_out)
    error("Variable names should be distinct.")
end
N = numel(vars_in);         M = numel(vars_out);
% Keep track of how many variables map from input to output
[vars_shared,N3_idcs_in,N3_idcs_out] = intersect(vars_in,vars_out,'stable');
N3 = numel(vars_shared);

% Check that the domains are properly specified
if ~isstruct(dom) || ~isfield(dom,'in') || ~isfield(dom,'out')
    error("Domains of the variables should be specified as a struct with fields 'in' and 'out'.")
end
dom_in = dom.in;            dom_out = dom.out;
if ~all(size(dom_in)==[N,2]) || ~all(size(dom_out)==[M,2])
    error("Number of intervals should match number of variables.")
end

% Check that the dependency array is appropriate
if ~isstruct(depmat) || ~isfield(depmat,'in') || ~isfield(depmat,'out')
    error("Dependency array should be specified as a struct with fields 'in' and 'out'.")
end
depmat_in = depmat.in;      depmat_out = depmat.out;
if (~isnumeric(depmat_in) && ~islogical(depmat_in)) || ~all(size(depmat_in)==[q,N])
    error("Input dependency matrix should be specified as q x N array of integers, for q inputs and N variables.")
elseif (~isnumeric(depmat_out) && ~islogical(depmat_out)) || ~all(size(depmat_out)==[q,M])
    error("Output dependency matrix should be specified as q x M array of integers, for q outputs and M variables.")
elseif any(depmat_out(:)>1)
    error("Output dependency matrix must be boolean; distributed output degree will be determined by input degree and type of tensor product used.")
end
% Create logical array indicating on what vars each input/output depends
is_dep_in = logical(depmat_in);     is_dep_out = logical(depmat_out);
% Keep track of number of vars on which each input/output depends
N_arr = sum(is_dep_in,2);           M_arr = sum(is_dep_out,2);
% Also keep track of for which vars a 3-PI operator is used
is_dep3 = depmat_in(:,N3_idcs_in) & depmat_out(:,N3_idcs_out);
N3_arr = sum(is_dep3,2);

% Check that the monomial degrees are appropriate
if nargin<=4
    degs = 0;
end
if isnumeric(degs) || isa(degs,'cell')
    % Assume the same maximal degree is used for input and output
    degval = degs;
    degs = struct();
    degs.in = degval;       degs.out = degval;
elseif ~isa(degs,'struct') || ~isfield(degs,'in') || ~isfield(degs,'out')
    error("Options 'degs' must be specified as a struct with fields 'in' and 'out'.")
end



%% Process the options
% Extract exclude options
% opts is the 6th input: use the defaults if it is omitted
if nargin<=5  
    opts = struct();
end
if ~isfield(opts,'exclude')
    % Include all types of operators by default
    exclude_arr = false([3*ones(1,N3),1]);
else
    exclude_arr = opts.exclude;
end
if N3>0
    if islogical(exclude_arr) || (isnumeric(exclude_arr) && all(ismember(exclude_arr(:),[0;1])))
        if numel(exclude_arr)==3
            % Assume same exclude along all dimensions
            exclude_arr = exclude_arr(:);
            excludeL = exclude_arr;
            new_shape = ones(1,N3+1);
            for i=2:N3
                new_shape(i) = 3;
                exclude_arr = exclude_arr & reshape(excludeL,new_shape);
                new_shape(i) = 1;
            end
        elseif numel(exclude_arr)~=3^N3
            error("Option 'exclude' must be specified as 3^N = 3 x ... x 3 logical array for N spatial variables.")
        end
        exclude_arr = reshape(exclude_arr,[3*ones(1,N3),1]);
        % Assume same exclude array for each state
        exclude_cell = cell(1,q);
        indcs_full = repmat({1:3},1,N3);
        for j=1:q
            indcs_tmp = indcs_full;
            indcs_tmp(~is_dep3(j,:)) = {1};
            indcs_tmp = [indcs_tmp,{1}];
            % Column of size 3^Ni, one dimension per 3-PI variable of the state;
            % squeeze would give a 1 x 3 row for Ni=1 at position 2, which the sep
            % update below cannot index
            exclude_cell{j} = reshape(exclude_arr(indcs_tmp{:}),[3*ones(1,N3_arr(j)),1,1]);
        end
    elseif isa(exclude_arr,'cell')
        if numel(exclude_arr)~=q
            error("Option 'exclude' should be specified as a 1 x q cell for q coupled states.")
        end
        exclude_cell = exclude_arr;
        for j=1:q
            Ni = N3_arr(j);
            exclude_arr = exclude_cell{j};
            if ~islogical(exclude_arr) || (numel(exclude_arr)~=3 && numel(exclude_arr)~=3^Ni)
                error("Each element of the option 'exclude' must be specified as a 3 x ... x 3 cell array, corresponding to the spatial dimensionality of the state.")
            end
            if numel(exclude_arr)==3
                % Assume same exclude along all dimensions
                exclude_arr = exclude_arr(:);
                excludeL = exclude_arr;
                new_shape = ones(1,Ni+1);
                for i=2:Ni 
                    new_shape(i) = 3;
                    exclude_arr = exclude_arr & reshape(excludeL,new_shape);
                    new_shape(i) = 1;
                end
            end
            if Ni==0
                exclude_cell{j} = false;
            else
                exclude_cell{j} = reshape(exclude_arr,[3*ones(1,Ni),1]);
            end
        end
    else
        error("Option 'exclude' must be specified as cell or logical array.")
    end
else
    exclude_cell = repmat({false},1,q);
end

% Extract separability options
if ~isfield(opts,'sep')
    % Default to no separability
    is_sep = false(1,N3);
else
    is_sep = opts.sep;
    if isnumeric(is_sep)
        is_sep = logical(is_sep);
    end
end
if islogical(is_sep)
    % Separability is the same for all state vectors
    if isscalar(is_sep)
        is_sep = repmat(is_sep,1,N3);
    elseif numel(is_sep)~=N3
        error("Option 'sep' must be specified as a 1 x N logical arrray, for N spatial variables.")
    end
    is_sep = reshape(is_sep,1,N3);
    sep_cell = cell(1,q);
    for j=1:q
        sep_cell{j} = is_sep(is_dep3(j,:));
    end
elseif isa(is_sep,'cell')
    % A different separability option is specified for each state
    if numel(is_sep)~=q
        error("Option 'sep' must be specified as 1 x q cell for q coupled states.")
    end
    sep_cell = is_sep;
    for j=1:q
        if ~islogical(sep_cell{j}) || numel(sep_cell{j})~=N3_arr(j)
            error("Each field of option 'sep' must be specified as a 1 x Ni logical array, for Ni the spatial dimensionality of the state.")
        end
    end
end
% Update the exclude array to account for separability
for j=1:q
    exclude_arr = exclude_cell{j};
    indcs_full = [repmat({1:3},1,N3_arr(j)),{1}];
    for i=find(sep_cell{j})
        % If upper integral is excluded, lower integral must be too
        indcs_tmp2 = indcs_full;    indcs_tmp2{i} = 2;   
        indcs_tmp3 = indcs_full;    indcs_tmp3{i} = 3;
        exclude_arr(indcs_tmp2{:}) = exclude_arr(indcs_tmp2{:}) | exclude_arr(indcs_tmp3{:});
        % Exclude upper integral if we enforce separability
        exclude_arr(indcs_tmp3{:}) = true;
    end
    exclude_cell{j} = exclude_arr;
end

% Check that the type of product is appropriate
if ~isfield(opts,'otype')
    % By default, use pointwise product
    otype = [1,2];
else
    otype = opts.otype;
end
if isscalar(otype) && ismember(otype,[1;2])
    otype = [otype,2];
elseif numel(otype)~=2 || any(~ismember(otype,[1;2]))
    error("Type of product used in declaring the operator must be specified as either 1 or 2.")
end
if otype(2)==1
    error("Pointwise product of the inputs is not supported.")
end


%% Declare the actual operator

% Extract the full list of monomials to use in defining the operator
degs.in = process_degs(degs.in,N_arr,depmat_in,is_dep_in,otype(2));
degs.out = process_degs(degs.out,M_arr,depmat_out,is_dep_out,otype(1));

% Declare a tensor-PI operator acting on each input
Zop = tensopmat();
for j=1:q
    
    % Set the input and output domains of each of the PI operators in the
    % tensor-PI operator j
    % --> they act on the same state, and have the same input/output space
    vars_tmp = struct();
    vars_tmp.in = vars_in(is_dep_in(j,:));           
    vars_tmp.out = vars_out(is_dep_out(j,:));
    dom_tmp = struct();
    dom_tmp.in = dom_in(is_dep_in(j,:),:);
    dom_tmp.out = dom_out(is_dep_out(j,:),:);
    var_out_idcs = find(is_dep_out(j,:));
    var_in_idcs = find(is_dep_in(j,:));
    % Positions, in the input and output variables of this state, of the
    % 3-PI variables (in input order), of the output-only and of the
    % input-only variables, from the dependence of this state,
    % NOT the global lists
    [~,N3_idcs_j_in,N3_idcs_j_out] = intersect(vars_tmp.in,vars_tmp.out,'stable');
    N1_idcs_j = setdiff(1:numel(vars_tmp.out),N3_idcs_j_out);
    N2_idcs_j = setdiff(1:numel(vars_tmp.in),N3_idcs_j_in);


    % Extract monomials in terms of which each of the operators is defined
    ZRj = degs.in{j};               ZLj = degs.out{j};
    dj = size(ZRj,1);
    if size(ZLj,1)==1 && dj>1
        % Use primary monomial only for first factor
        ZLj = [ZLj;num2cell(zeros(dj-1,M_arr(j)))];
    end
    % Check the number of monomials included in the operators
    nZL_arr = zeros(size(ZLj));
    nZR_arr = zeros(size(ZRj));
    for k=1:numel(ZLj)
        nZL_arr(k) = numel(ZLj{k});
    end
    for k=1:numel(ZRj)
        nZR_arr(k) = numel(ZRj{k});
        if ~ismember(0,ZRj{k})
            error("Monomial basis in independent variables must include constant monomial.")
        end
        % Make sure constant monomial is in first index
        ZRj{k} = sort([ZRj{k}]);
    end
    nZL_list = prod(nZL_arr,2);     % total number of monomials in all vars
    nZR_list = prod(nZR_arr,2);
    
    % Initialize empty parameters for the operators
    nj = dims(j);
    params = cell([3*ones(1,N3_arr(j)),1]);

    % Keep track of which parameters we exclude, which parameters are
    % multipliers, and which parameters are equal due to separability
    include_j = ~exclude_cell{j};
    is_sep_j = sep_cell{j};
    is_mult_params = false(numel(params),N3_arr(j));
    params_rtn = params(include_j);     % excluded parameters will be set to 0
    rtn2full_idcs = reshape(1:3^N3_arr(j),size(params));
    for i=1:N3_arr(j)
        if is_sep_j(i)
            % Set parameter of upper integral equal to that of lower
            % integral
            idx_old = [repmat({1:3},1,N3_arr(j)),{1}];
            idx_old{i} = 3;
            idx_new = [repmat({1:3},1,N3_arr(j)),{1}];
            idx_new{i} = 2;
            rtn2full_idcs(idx_old{:}) = rtn2full_idcs(idx_new{:});
        end
        % Indicate which parameters are multiplier along variable i
        mult_idcs = (1:3^(i-1))' + (0:3^(N3_arr(j)-i)-1)*3^i;
        is_mult_params(mult_idcs,i) = true;
    end
    is_mult_params = is_mult_params(include_j(:),:);    

    % Get the total row dimension of the operator and the parameters
    nr_op_list = zeros(dj,1);                       % number of rows in each operator
    nr_term_arr = zeros(dj,numel(params_rtn));      % number of rows in each term (distinct integral) in each operator
    for k=1:dj
        nr_term_arr(k,:) = prod(nZL_arr(k,N1_idcs_j),2)*prod(nZL_arr(k,N3_idcs_j_out).*(nZR_arr(k,N3_idcs_j_in).^(1-is_mult_params)),2)*prod(nZR_arr(k,N2_idcs_j),2);
        nr_op_list(k) = nj*sum(nr_term_arr(k,:));
    end
    nr_op_tot = prod(nr_op_list);
    if otype(1)==1 && dj>1
        nr_sopvar_rows = repmat(nr_op_tot,dj,1);
    else
        nr_sopvar_rows = nr_op_list;
    end
    nr_param_list = nr_sopvar_rows.*nZL_list;

    % Set the operator in each factor
    for k=1:dj
        % Determine the number of monomials in each parameter
        nZL_k = nZL_list(k);       
        nZR_k = nZR_list(k);

        % Determine the row dimension of the preceding, current and
        % following operators
        sz_prev = prod([1;nr_op_list(1:k-1)]);
        sz_post = prod([1;nr_op_list(k+1:end)]);
        sz_now = nr_op_list(k);

        % Set indices for all right-monomials, excluding potentially
        % introduced constant monomial
        if any(nZR_arr(k,:))
            zr_idcs_full = permute(reshape(1:nZR_list(k),[nZR_arr(k,end:-1:1),1]),[N_arr(j):-1:1,N_arr(j)+1]);
            rtn_idcs_ZR = [cellfun(@(a)1:a,num2cell(nZR_arr(k,:)),'UniformOutput',false),{1}];
        else
            zr_idcs_full = 1;
            rtn_idcs_ZR = {1,1};
        end

        % Declare all the parameters to include        
        r_strt = 0;
        for ll=1:numel(params_rtn)
            % Set row and column indices to include all monomials that
            % appear in the parameters, accounting for exclusion of
            % right-monomials in multiplier terms
            is_mult_ll = is_mult_params(ll,:);
            mult_mask = false(1,N_arr(j));
            mult_mask(N3_idcs_j_in(is_mult_ll)) = true;
            rtn_idcs_ZR_ll = rtn_idcs_ZR;
            rtn_idcs_ZR_ll([mult_mask,false]) = {1};     % Keep only first monomial
            zr_idcs_ll = zr_idcs_full(rtn_idcs_ZR_ll{:});
            zl_idcs_ll = (1:nZL_k)';        % Include all left monomials
            % Set row and column indices of sparse matrix C, defining a
            % vector of monomials as
            %   Zfull = kron(eye(nZ),ZL)' * C * ZR
            c_idcs_ll = repmat(zr_idcs_ll(:),nZL_k,1);
            r_idcs_ll = (0:numel(c_idcs_ll)-1)'*nZL_k + repelem(zl_idcs_ll,numel(zr_idcs_ll),1);
            % Account for number of state variables
            %   kron(n,Zfull) = kron(eye(n*nZ),ZL)' * C * kron(eye(n),ZR)
            c_idcs0 = repmat(c_idcs_ll,nj,1) + repelem((0:nj-1)'*nZR_k,numel(c_idcs_ll),1);
            r_idcs0 = repmat(r_idcs_ll,nj,1) + repelem((0:nj-1)'*nr_term_arr(k,ll)*nZL_k,numel(r_idcs_ll),1);
            
            % Account for elementwise product of operators
            if otype(1)==1 && dj>1
                % Keep track of row index and left monomial associated with
                % each coefficient
                rr_idcs0 = ceil(r_idcs0(:)./nZL_k);       
                zL_idcs0 = r_idcs0(:) - (rr_idcs0-1)*nZL_k;
                % Repeat the column indices to match the total size
                c_idcs = repmat(c_idcs0,sz_post*sz_prev,1);
                % Perform repetition along the row direction
                rr_idcs = rr_idcs0 + r_strt/nZL_k;
                rr_idcs = (rr_idcs-1)*sz_post + (1:sz_post);
                zL_idcs = repmat(zL_idcs0(:),1,sz_post);
                r_idcs = (rr_idcs(:)-1)*nZL_k + zL_idcs(:);
                r_idcs = repmat(r_idcs,sz_prev,1) + repelem((0:sz_prev-1)'*sz_now*sz_post,numel(r_idcs),1);
            else
                % Account for the vertical concatenation of the parameters
                r_idcs = r_idcs0 + r_strt;
                c_idcs = c_idcs0;
            end
            r_strt = r_strt + r_idcs0(end);
            % Set the actual parameter
            params_rtn{ll} = sparse(r_idcs,c_idcs,1,nr_param_list(k),nj*nZR_k);
        end
        
        % Declare the jth operator in the factor
        params(include_j) = params_rtn;
        params(:) = params(rtn2full_idcs);      % account for separability
        Zop_k = sopvar(params,vars_tmp,ZLj(k,:),ZRj(k,:),dom_tmp,[nr_sopvar_rows(k),nj]);

        % Take tensor product of the operators to match distributed
        % monomial degree of state
        if k==1
            % Cell route, so that the variables are stored sorted
            Zop_j = tensopvar({Zop_k});
            Zop_j.otype = otype;
        else
            Zop_j = otimes(Zop_j,Zop_k,otype);
        end
    end
    % Take block-diagonal cocatenation to obtain the full operator
    Zop = blkdiag(Zop,tensopmat(Zop_j));
end


end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function degs = process_degs(degs,N_arr,depmat,is_deps,otype)
% Check the specified monomial degrees and augment to account for the
% number of coupled state variables, number of factors in the tensor-PI 
% operators, and the number of spatial variables in terms of which each
% operator is defined. Return the full list of monomial degrees in each of
% the variables (rather than just the maximal monomial degrees);

[q,N] = size(depmat);

% If degrees are specified for only one state, repeat for all states
if isnumeric(degs)
    % The same maximal monomial degree is used for all states
    degs = repmat({degs},1,q);
elseif isa(degs,'cell') && isscalar(degs)
    % The same full monomial list is used for all states
    degs = repmat(degs,1,q);
elseif ~isa(degs,'cell') || numel(degs)~=q
    error("Degrees must be specified as integer array or 1 x q cell for q inputs.")
end

% Check the monomial degrees for each of the state variables
for j=1:q
    degval = degs{j};
    if ~any(depmat(j,:))
        degs{j} = cell(1,0);
        continue
    else
        nfctrs = max(depmat(j,:))^(otype-1);
    end
    if isnumeric(degval)
        % The maximal monomial degrees are specified
        if isscalar(degval)
            % Same cumulative degree for each of the N variables
            degval = degval*ones(1,N_arr(j));
        elseif size(degval,2)==N
            % Retain only cumulative degree for relevant variables
            degval = degval(:,is_deps(j,:));
        elseif size(degval,2)~=N_arr(j)
            error("Column dimension of maximal degree array must match total number of spatial variables.")
        end
        % Make sure monomial degrees are specified for each factor
        if size(degval,1)==1
            % Set monomial degree adding up to cumulative degree
            deg_max = degval;
            degval = zeros(nfctrs,N_arr(j));
            for k=1:nfctrs
                degval(k,:) = ceil(deg_max./(nfctrs-(k-1)));
                deg_max = deg_max - degval(k,:);
            end
        elseif size(degval,1)==max(depmat(:).^(otype-1))
            % Discard rows exceeding the degree of this monomial
            degval = degval(1:nfctrs,:);
        elseif size(degval,1)~=nfctrs
            error("Row dimension of maximal degree array must match the distributed monomial degree of the factors.")
        end
        % Convert to true monomial degrees
        degs{j} = cellfun(@(a)(0:a)',num2cell(degval),'UniformOutput',false);

    else
        % The actual monomial degrees are specified
        if size(degval,2)==N
            % Retain only relevant monomials
            degval = degval(:,is_deps(j,:));
        elseif size(degval,2)~=N_arr(j)
            error("Number of monomial vectors should match number of variables in terms of which the operator is defined.")
        end
        % Make sure monomials are specified for each factor
        if size(degval,1)==1
            % Same monomials for each of the factors
            degval = repmat(degval,nfctrs,1);
        elseif size(degval,1)==max(depmat(:).^(otype-1)) 
            % Discard monomials exceeding degree of the factor
            degval = degval(1:nfctrs,:);
        elseif size(degval,1)~=nfctrs
            error("Number of monomial vectors should match number of factors in the tensor-PI operator.")
        end
        degs{j} = degval;
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %