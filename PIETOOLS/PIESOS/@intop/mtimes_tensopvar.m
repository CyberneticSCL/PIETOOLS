function Cop = mtimes_tensopvar(Aop,Bop)
% COP = MTIMES_TENSOPVAR(AOP,BOP) returns the functional COP representing
% the composition of the functional AOP with the tensor-PI operator BOP.
%
% INPUTS
% - Aop:    p x m 'intop' object representing a PI functional;
% - Bop:    m x n 'tensopvar' object representing a tensor-PI operator.
%           The factors of Bop may be 'tensopvar' objects themselves, e.g.
%           the pointwise product of operators as a factor in a Kronecker
%           product, see get_params;
%
% OUTPUTS
% - Cop:    p x n 'intop' object representing the composition of Aop and
%           Bop;
%
% NOTES
% - Elementwise products of the input states of Bop (Bop.otype(2)=1) are
%   not supported.
% - The coefficients of Aop are taken directly from the sparse matrices
%   Aop.params, in the monomials Aop.Z.in. These form a tensor product over
%   the variables, but not over the dummy variables of each variable, so
%   they are first embedded in the tensor product over all dummy variables.
% - The composition is computed on the coefficients, represented as structs
%   with fields i,j,v,m,n (the nonzero values and their row and column
%   indices). The coefficients of Cop are declared as sparse matrices in
%   the monomials of the dummy variables of Cop, with the decision
%   variables Aop.Zd.
% - Monomials that are not used by any of the integrals of Cop are removed
%   from the monomials Cop.Z.in.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes_tensopvar
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
% DJ, 10/21/2026: Compare empty variable lists regardless of shape;
% DJ, 10/21/2026: Name variable groups S1, S2 as in 'sopvar';

% Check that the functional is properly specified.
if ~isa(Aop,'intop')
    error("Functional must be specified as an object of type 'intop'.")
elseif ~Aop.chkval || ~isvalid(Aop)
    error("Operation is not supported: The functional is not properly declared.")
end

% Check that the tensor-PI operator is properly specified.
if ~isa(Bop,'tensopvar')
    error("Tensor-PI operator must be specified as an object of type 'tensopvar'.")
end
if ~Bop.chkval
    Bop = initialize(Bop);
end
if Bop.otype(2)==1
    error("Elementwise product of the input states is not supported.")
end

% Make sure the input domain of the functional matches the output of the
% operator
if size(Aop,2)~=size(Bop,1)
    error("Number of columns of functional must match number of rows of tensor-PI operator.")
% Compare the variables as columns, so that empty cells of different shape (0x0, 1x0) are equal
elseif ~isequal(Aop.vars.in(:),Bop.vars.out(:)) || ~isequal(Aop.dom.in,Bop.dom.out) || ~isequal(Aop.depmat.in,Bop.deps.out)
    error("Input function space of the functional must match output space of the operator.")
end

% Composition with an identity operator leaves the functional unchanged
% (Aop maps from Bop's output space, so for identity Bop, Cop = Aop, which
% also maps from Bop's input space; compared explicitly rather than assumed)
if isequal(Bop.vars.in,Aop.vars.in) && isequal(Bop.dom.in,Aop.dom.in) && isequal(Bop.deps.in,Aop.depmat.in) && is_id_op(Bop)
    Cop = Aop;
    return
end

% Decompose the tensor-PI operator
vars1 = Bop.vars.out;       M = numel(vars1);
vars2 = Bop.vars.in;        N = numel(vars2);
dom = Bop.dom;
deps1 = Bop.deps.out;         % to what degree does the output depend on each var
deps2 = Bop.deps.in;
% Determine which variables are passed through Bop
% Names as in 'sopvar': S1 = vars.in\vars.out (integrated out),
% S2 = vars.out\vars.in (introduced), S3 = common (passed through)
[vars_S3,idcs_S3_out,idcs_S3_in] = intersect(vars1,vars2,'stable');
idcs_S2 = setdiff(1:M,idcs_S3_out);
idcs_S1 = setdiff(1:N,idcs_S3_in);
N1 = numel(idcs_S1);
N2 = numel(idcs_S2);
N3 = numel(idcs_S3_out);
dom_S3 = dom.out(idcs_S3_out,:);

% Clean zero terms/integrals from both operands before the decomposition
% below
Aop = clean_op(Aop);
[Bop,is_zeroB] = clean_op(Bop);
if is_zeroB
    % Bop is zero up to tolerance -- return the zero functional directly.
    intmat = cell(1,N1+N3);
    for i=1:N1+N3
        intmat{i} = repmat([0,deps2(i)+1],1,deps2(i));
    end
    nd = numel(Aop.Zd);
    params = {sparse(size(Aop,1)*(nd+1),size(Bop,2))};
    Zin = cell(1,N1+N3);
    for i=1:N1+N3
        Zin{i} = zeros(1,deps2(i));
    end
    vars = Bop.vars;        vars.out = cell(1,0);
    dom = Bop.dom;          dom.out = zeros(0,2);
    Cop = intop(params,intmat,vars,dom,Zin,Aop.Zd,true);
    return
end

% Set the names of the variables as they appear in the parameters defining
% Bop, and in the monomials of the integrated variables
varnames_B = struct();
varnames_B.out = cell(1,M);
for i=1:M
    outvars_Bi = cell(1,deps1(i));
    for k=1:deps1(i)
        outvars_Bi{k} = ['t',num2str(i),'_',num2str(k)];
    end
    varnames_B.out{i} = outvars_Bi;
end
varnames_B.in = cell(1,N);
for i=1:N
    invars_Bi = cell(1,deps2(i));
    for k=1:deps2(i)
        invars_Bi{k} = ['t',num2str(i),'_dum',num2str(k)];
    end
    varnames_B.in{i} = invars_Bi;
end
varnames_B_full = struct('in',{hcat_cellstr(varnames_B.in)},'out',{hcat_cellstr(varnames_B.out)});
varnames_Z =struct('in',{varnames_B.out([idcs_S2,idcs_S3_out(:)'])},'out',{cell(1,0)});
varnames_Z_full = struct('in',{hcat_cellstr(varnames_B.out)},'out',{cell(1,0)});

% Express parameters defining A in terms of full list of variables;
% order as [vars_S2, vars_S3]
intmat_A = Aop.intmat;
Aparams = permute(Aop.params,[idcs_S2,idcs_S3_out(:)',M+1,M+2]);
sz_Aparams = size(Aparams,1:N1+N2+N3);
full_idcsA = cellfun(@(a) 1:a,num2cell(sz_Aparams),'UniformOutput',false);

% Get intervals of integration and associated parameters defining Bop,
% in terms of the full list of variables, [vars_S2, vars_S3, vars_S1]
% Start with vars introduced by Bop
% --> just whatever integral Aop takes
intmat_S2 = intmat_A(idcs_S2);
% Next, vars removed by Bop
% --> must be integral over full domain
intmat_S1 = cell(1,N1);
for j=1:N1
    d2j = deps2(idcs_S1(j));
    intmat_S1{j} = repmat([0,d2j+1],1,d2j);
end
% Finally, vars passed through Bop
% --> need to check for each factor in Bop what integrals are taken
dtot_arr = deps1(idcs_S3_out)+deps2(idcs_S3_in);
Bvars = Bop.vars;
Bvars.S3 = vars_S3;
[Bparams,intmatB_S3] = get_params(Bop,Bvars,varnames_B,dtot_arr);

intmatA_S3 = intmat_A(idcs_S3_out);
for i=1:N3
    intmat_i = intmatA_S3{i};
    intmatA_S3{i}(intmat_i==deps1(idcs_S3_out(i))+1) = dtot_arr(i)+1;
end

% Express the parameters in terms of the same basis of monomials
%   Aparams{k} = Acoeff{i} * (Im o ZAR)
% The sparse parameters of Aop are in the same monomials Aop.Z.in already. These
% form a tensor product over the variables, but not over the dummy variables of
% each variable: embed in the tensor product over all dummy variables ZAR, by
% relabeling the columns
[Zgrid_A,Zdegs_A] = grid_basis(Aop.Z.in);
zmapA = Zmap(Aop.Z.in,Zgrid_A);
nZgrid_A = prod(cellfun(@(a) size(a,1),Zgrid_A));
Acoeffs = cell(size(Aparams));
for k=1:numel(Aparams)
    Acoeffs{k} = sparse2struct(Aparams{k},zmapA,nZgrid_A);
end
ZAR = quadPoly(1);
if ~isempty(Zdegs_A)
    sz_coeffs = prod(cellfun(@(a) numel(a),Zdegs_A));
    ZAR = quadPoly(speye(sz_coeffs),cell(1,0),Zdegs_A,[sz_coeffs,1],cell(1,0),varnames_Z_full.in);
end
%   Bparams{k} = (Im o ZBL) * Bcoeffs{k} * (In o ZBR)
[Bcoeffs,ZBL,ZBR] = quad_expand_quadPoly(Bparams,varnames_B_full);
Z_inner = ZAR*ZBL;
[C_inner,~,Z_inner] = quad_expand_quadPoly(Z_inner,varnames_Z_full);
for k=1:numel(Acoeffs)
    Acoeffs{k} = mtimes_kronR(Acoeffs{k},C_inner{1});
end

% Take the integrals specified in Aop of the S2 variables introduced by Bop
for i=1:N2
    % Determine the variables and domain associated with var i
    varnames_i = varnames_B.out{idcs_S2(i)};
    dom_i = dom.out(idcs_S2(i),:);
    domvals_i = [{dom_i(1)},varnames_i,{dom_i(2)}];
    % Perform the integrals specified by intmat_S2{i}
    intmat_i = intmat_S2{i};
    params_i = cell(size(intmat_i,1),1);
    for j=1:size(intmat_i,1)
        params_i{j} = int_params(Z_inner,intmat_i(j,:),domvals_i);
    end
    % Express coefficient maps in terms of a unified basis
    varnames_Z.in(1) = [];
    varnames_Z_full.in = hcat_cellstr(varnames_Z.in);
    [Ci,~,Z_inner] = quad_expand_quadPoly(params_i,varnames_Z_full);
    % Add parameters of which we have just evaluated the integrals
    sz_Aparams(i) = 1;
    Acoeffs_new = repmat({0},[sz_Aparams,1]);
    for j=1:size(intmat_i,1)
        full_idcsA{i} = j;
        Acoeffs_old = Acoeffs(full_idcsA{:});
        for k=1:numel(Acoeffs_old)
            Acoeffs_new{k} = plus_spstruct(Acoeffs_new{k},mtimes_kronR(Acoeffs_old{k},Ci{j}));
        end
    end
    full_idcsA{i} = 1;
    Acoeffs = Acoeffs_new;
end
% Remove dimensions along variables S2
sz_Aparams = sz_Aparams(N2+1:end);
full_idcsA = full_idcsA(N2+1:end);
Acoeffs = reshape(Acoeffs,[sz_Aparams,ones(1,2-numel(sz_Aparams))]);

% Take all combinations of integrals in Aop and Bop for S3 vars
param_idcsA = cell(1,N3);       param_idcsB = cell(1,N3);
intmat_S3 = cell(1,N3);
int_orderA = Aop.int_order(idcs_S3_out);
var_orderA = cell(1,N3);
for i=1:N3
    intmatA_i = intmatA_S3{i};      ntrms_A = size(intmatA_S3{i},1);
    intmatB_i = intmatB_S3{i};      ntrms_B = size(intmatB_S3{i},1);
    intmat_S3{i} = [repelem(intmatA_i,ntrms_B,1), repmat(intmatB_i,ntrms_A,1)];
    % Keep track of which parameters in A and B correspond to each integral
    param_idcsA{i} = repelem((1:ntrms_A)',ntrms_B,1);
    param_idcsB{i} = repmat((1:ntrms_B)',ntrms_A,1);
    % Also keep track of the order of integration
    int_orderA_i = int_orderA{i};
    [~,var_orderA_i] = sort(int_orderA_i,2);    
    var_orderA{i} = repelem(var_orderA_i,ntrms_B,1);
end

% Take the appropriate integrals for the S3 vars
for i=1:N3
    % Extract the intervals of integration for variable i
    intmat_i = intmat_S3{i};
    ntrms_i = size(intmat_S3{i},1);
    di_A = size(intmatA_S3{i},2)/2;
    var_orderA_i = var_orderA{i};
    [var_orderA_i,~,var_orderA_idx] = unique(var_orderA_i,'rows');
    n_ordsA = size(var_orderA_i,1);
    Zcell = repmat({Z_inner},n_ordsA,1);
    Ccell_trms = num2cell(ones(ntrms_i,1));
    param_idcsA_i = reshape(param_idcsA{i},[],1);
    param_idcsB_i = reshape(param_idcsB{i},[],1);
    % Get the dummy variables associated with variable i
    varnames_i = varnames_B.out{idcs_S3_out(i)};
    varnames_i_in = varnames_B.in{idcs_S3_in(i)};
    varnames_Z_full.in = [varnames_i_in, varnames_Z_full.in];
    domvals_i = [dom_S3(i,1),varnames_i,varnames_i_in,dom_S3(i,2)];
    % Loop over all variables, performing the necessary integrals of
    % Z over the variables in A
    for k=1:di_A
        ntrms_k = size(intmat_i,1);                                       
        intmat_i_tmp = intmat_i;                                        
        param_idcsA_tmp = param_idcsA_i;                                 
        param_idcsB_tmp = param_idcsB_i;                                
        var_orderA_idx_tmp = var_orderA_idx;                            
        Ccell_tmp = Ccell_trms;                                          
        intmat_blocks = cell(ntrms_k,1);                                  
        param_idcsA_blocks = cell(ntrms_k,1);                             
        param_idcsB_blocks = cell(ntrms_k,1);                            
        var_orderA_idx_blocks = cell(ntrms_k,1);                         
        Ccell_blocks = cell(ntrms_k,1);
        for j=1:ntrms_k       % could be parallelized, if we split terms into a cell
            % Take the (remaining) variable in A that must be integrated
            % first
            trm = j;                                                     
            varsA = var_orderA_i(var_orderA_idx_tmp(trm),1:k);           
            varA = varsA(end);
            int_dom_varA = intmat_i_tmp(trm,[2*varA-1,2*varA]);         
            row_out = intmat_i_tmp(trm,:);                              

            % Determine which variables in B are integrated over a domain
            % defined by varA
            is_varlimA = row_out==varA;    
            is_varlimA(2*varsA(1:end-1)-1) = false;
            is_varlimA(2*varsA(1:end-1)) = false;
            if int_dom_varA(1)==int_dom_varA(2)
                % Aop acts as multiplier along varA: varA = int_dom_varA(1)
                % --> replace varA with int_dom_varA(1)
                row_out(is_varlimA) = int_dom_varA(1);
                row_out([2*varA-1,2*varA]) = int_dom_varA;  
            elseif ~any(is_varlimA)
                % This variable does not appear in the limits of
                % integration of any other variable
                % --> we can just take whatever integral has been
                % specified for varA
            elseif any(is_varlimA(1:2:end) & is_varlimA(2:2:end))
                % The variable appears in a multiplier term
                % --> reorder
                varB = find(is_varlimA(1:2:end) & is_varlimA(2:2:end),1,'last');
                % Integrate varB over domain of varA
                row_out([2*varB-1,2*varB]) = int_dom_varA;
                % Evaluate others multiplier vars at varB
                is_varlimA([2*varB-1,2*varB]) = false;
                row_out(is_varlimA) = varB;
                % Evaluate varA==varB
                row_out([2*varA-1,2*varA]) = [varB,varB];
            else
                % The variable appears as limit of integration of at least
                % one other variable
                % --> use the function "switch_ints" to get all possible
                %     integrals that are consistent with these limits
                intmat_tmp = switch_ints(row_out,varsA);
                ntrms_ll = size(intmat_tmp,1);                            
                intmat_blocks{j} = intmat_tmp;
                param_idcsA_blocks{j} = param_idcsA_tmp(trm)*ones(ntrms_ll,1);
                param_idcsB_blocks{j} = param_idcsB_tmp(trm)*ones(ntrms_ll,1);
                var_orderA_idx_blocks{j} = repmat(var_orderA_idx_tmp(trm),ntrms_ll,1);
                Ccell_blocks{j} = repmat(Ccell_tmp(trm),ntrms_ll,1); 
                continue
            end
            intmat_blocks{j} = row_out;                                   
            param_idcsA_blocks{j} = param_idcsA_tmp(trm);               
            param_idcsB_blocks{j} = param_idcsB_tmp(trm);              
            var_orderA_idx_blocks{j} = var_orderA_idx_tmp(trm);  
            Ccell_blocks{j} = Ccell_tmp(trm);
        end
        intmat_i = vertcat(intmat_blocks{:});                            
        param_idcsA_i = vertcat(param_idcsA_blocks{:});                   
        param_idcsB_i = vertcat(param_idcsB_blocks{:});                  
        var_orderA_idx = vertcat(var_orderA_idx_blocks{:});
        Ccell_trms = vertcat(Ccell_blocks{:});

        % Determine which variables from Aop need to be integrated,
        % and over which domain  
        intdoms_A = intmat_i(:,1:2*di_A);
        var_listA = var_orderA_i(:,k);
        for idx = 1:n_ordsA     % could be parallelized
            % Determine which integrals/coefficients correspond to the
            % variable order var_orderA_i(idx,:)
            varA = var_listA(idx);
            varnameA = varnames_i(varA);
            rtn_idcs = var_orderA_idx==idx;
            Ccell_idx = Ccell_trms(rtn_idcs);
            % Consider all unique integral w/r to varA that must be taken
            intdoms_Ak = intdoms_A(rtn_idcs,[2*varA-1,2*varA]);            
            [intdoms_Ak,~,idcs2] = unique(intdoms_Ak,'rows');
            % Compute the actual integral
            Zcell_tmp = cell(size(intdoms_Ak,1),1);
            for ll=1:size(intdoms_Ak,1)
                Lval = domvals_i(intdoms_Ak(ll,1)+1);
                Uval = domvals_i(intdoms_Ak(ll,2)+1);
                if ~isequal(Lval,Uval)
                    Zcell_tmp{ll} = int(Zcell{idx},varnameA,Lval,Uval);
                elseif ~isequal(varnameA,Lval)
                    Zcell_tmp{ll} = subs(Zcell{idx},varnameA,Lval);
                end
            end
    
            % Express integrated monomials in terms of the same monomial basis
            % could probably be improved, only consider vars i
            % We should really also split of monomials in variables not i
            varnames_Z_i = varnames_Z_full;
            varnames_Z_i.in = setdiff(varnames_Z_i.in,varnames_i(var_orderA_i(idx,1:k)));
            [Ccell_tmp,~,Z_tmp] = quad_expand_quadPoly(Zcell_tmp,varnames_Z_i);
            Ccell_tmp = Ccell_tmp(idcs2);       % account for repeated integrals
            for ll=1:numel(Ccell_tmp)
                Ccell_idx{ll} = mtimes_spstruct(Ccell_idx{ll},Ccell_tmp{ll});
            end
            Ccell_trms(rtn_idcs) = Ccell_idx;
            Zcell{idx} = Z_tmp;
        end
    end

    % Set the intervals of integration (in variable i) of the output Cop
    intmat_C = intmat_i(:,2*di_A+1:end);
    intmat_C(intmat_C>0) = intmat_C(intmat_C>0) - di_A;
    intmat_C = canonical_intmat(intmat_C);
    % Retain only unique integrals
    [intmat_aug,~,idcs2] = unique([param_idcsB_i(:),intmat_C],'rows');
    intmat_C = intmat_aug(:,2:end);
    param_idcsB{i} = intmat_aug(:,1);
    intmat_S3{i} = intmat_C;

    % Express integrated monomials in terms of the same monomial basis
    % could probably be improved, only consider vars i
    varnames_Z.in{i} = varnames_i_in;
    varnames_Z_full.in = hcat_cellstr(varnames_Z.in);
    [Ccell_ords,~,Z_inner] = quad_expand_quadPoly(Zcell,varnames_Z_full);

    % Multiply the coefficients defining A with the appropriate factors
    ntrms_old = numel(Ccell_trms);    
    ntrms_new = size(intmat_C,1);
    sz_Aparams(i) = ntrms_new;
    full_idcsA{i} = 1:ntrms_new;
    tmp_idcsA = full_idcsA;
    Acoeffs_new = num2cell(zeros([sz_Aparams,1]));
    for j=1:ntrms_old             % can be parallelized
        % Get coefficients for A associated with term j
        idxA = param_idcsA_i(j);
        tmp_idcsA{i} = idxA;
        Aparams_old = Acoeffs(tmp_idcsA{:});
        % Multiply with the appropriate factor
        tmp_idcsA{i} = idcs2(j);
        Aparams_new = Acoeffs_new(tmp_idcsA{:});
        Ck = mtimes_spstruct(Ccell_trms{j},Ccell_ords{var_orderA_idx(j)});
        for k=1:numel(Aparams_new)
            Aparams_new{k} = plus_spstruct(Aparams_new{k},mtimes_kronR(Aparams_old{k},Ck));
        end
        Acoeffs_new(tmp_idcsA{:}) = Aparams_new;
    end
    Acoeffs = Acoeffs_new;
    clear Acoeffs_new;
end

% Take the product of the monomials with the coefficients defining B
% --> compute Bnew = (B o I_nZ), so that Bnew*(Im o Z) = (In o Z)*B
nZ = size(Z_inner,1);
% (B o I_nZ) is not formed: A*(B o I_nZ) is computed from B itself, by
% mtimes_kronL below
% Multiply with right-monomials from B parameters
ZR = IkronZ(size(ZBR,1),Z_inner)*ZBR;

% Terms with the same intervals of integration are merged, as unique_intmat does
% below, but here on the coefficient structs, before each term is converted to a
% sparse matrix
sz_C = [sz_Aparams,ones(1,2-numel(sz_Aparams))];
sz_new = sz_C;      map_S3 = cell(1,N3);
for i=1:N3
    [P_i,intmat_u] = uniquerows_integerTable(intmat_S3{i});
    if size(intmat_u,1)<size(intmat_S3{i},1)
        map_S3{i} = full(P_i*(1:size(P_i,2))');      % new index of each old term
        intmat_S3{i} = intmat_u;
        sz_new(i) = size(intmat_u,1);
    else
        map_S3{i} = (1:size(intmat_S3{i},1))';
    end
end
sz_new_c = [sz_new(1:N3),prod(sz_new(N3+1:end))];
Csum = cell(sz_new);        Csum(:) = {0};

% Take the product of the coefficients in A and B
for k=1:numel(Acoeffs)      % can be parallelized
    idcs_k = cell(1,N3+1);
    [idcs_k{:}] = ind2sub(sz_C,k);
    idcs_B = cellfun(@(a,b) a(b),param_idcsB,idcs_k(1:end-1),'UniformOutput',false);
    AB = mtimes_kronL(Acoeffs{k},Bcoeffs{idcs_B{:}},nZ);
    % Add the product to the term it merges into
    idcs_new = idcs_k;
    for i=1:N3
        idcs_new{i} = map_S3{i}(idcs_k{i});
    end
    if N3==0
        t = 1;
    else
        t = sub2ind(sz_new_c,idcs_new{:});
    end
    Csum{t} = plus_spstruct(Csum{t},AB);
end
% Declare the coefficient matrices of the integrals as sparse matrices, with
% the monomials of the dummy variables of Cop given by ZR
[ridcs,cidcs,vals] = find(ZR.C);
B_C = struct('i',ridcs(:),'j',cidcs(:),'v',vals(:),'m',size(ZR.C,1),'n',size(ZR.C,2));
Cparams = cell(sz_new);
for t=1:numel(Csum)
    C_C = mtimes_kronR(Csum{t},B_C);
    Cparams{t} = sparse(C_C.i,C_C.j,C_C.v,C_C.m,C_C.n);
end
clear Csum

% Reorder the variables
vars_in = [Bvars.out(idcs_S3_out),Bvars.in(idcs_S1)];
[vars_in,new_order] = sort(vars_in);
doms_in = [dom.out(idcs_S3_out,:); dom.in(idcs_S1,:)];
doms_in = doms_in(new_order,:);
varsC = struct('in',{vars_in},'out',{cell(1,0)});
domC = struct('in',doms_in,'out',zeros(0,2));
Cparams = permute(Cparams,[new_order,numel(new_order)+1:ndims(Cparams)]);
intmatC = [intmat_S3,intmat_S1];
intmatC = intmatC(new_order);

% Declare the functional
% Monomials of Cop along each variable, ordered as the variables of Cop.
% ZR.Zt specifies the degrees of the monomials in each dummy variable,
% for the variables in the order of Bop.vars.in
cs = [0,cumsum(deps2)];
ZinB = cell(1,N);
for i=1:N
    ZinB{i} = cartesian_grid(ZR.Zt(cs(i)+1:cs(i+1)));
end
[~,idxB] = ismember(vars_in,Bvars.in);
Zin_C = ZinB(idxB);
if ~isequal(idxB,1:N)
    % Variables of Cop are ordered differently than those of Bop: permute the
    % directions of the monomials
    zperm = Zpermmap(ZinB,idxB);
    nZ_C = prod(cellfun(@(a) size(a,1),ZinB));
    for t=1:numel(Cparams)
        Cparams{t} = remap_cols(Cparams{t},zperm,nZ_C);
    end
end

% Declare the functional
Cop = intop(Cparams,intmatC,varsC,domC,Zin_C,Aop.Zd,true);
Cop = unique_intmat(Cop);
Cop = prune_basis(Cop);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [params,intmat] = get_params(Rop,vars,varnames,dtot,degs_out,degs_in,prune)
% Extract the parameters and intervals of integration defining the
% tensor-PI operator Rop
%
% INPUTS
% - Rop:        'tensopvar' object representing a tensor-PI operator;
% - vars:       struct with fields
%     in:       1 x N cellstr, specifying the names of the input variables
%               to the tensor-PI operator;
%    out:       1 x M cellstr, specifying the names of the output variables
%               of the tensor-PI operator;
%     S3:       1 x N3 cellstr, specifying the names of the variables for
%               which the operator acts as 3-PI operator;
% - varnames:   struct with fields
%         in:   1 x N cell, with each element a cellstr specifying the
%               names of the dummy variables associated with each spatial
%               input variable to the operator;
%        out:   1 x M cell, with each element a cellstr specifying the
%               names of the dummy variables associated with each spatial
%               output variable to the operator;
% - dtot:       1 x N3 array specifying the total degree of each of the
%               3-PI variables as it appears in the operator
% - degs_out,degs_in:   (optional) 1 x M and 1 x N arrays of integers, with
%               for each output and input variable the index of the dummy
%               variable associated with the first factor in Rop. Default
%               to ones, but differ for a nested operator (see below);
% - prune:      (optional) boolean, set to false to retain the parameters
%               for all of the integrals defining the factors, also if they
%               are zero in each term. Defaults to true;
%
% NOTES
% - A factor of Rop may be a 'tensopvar' object itself, e.g. the pointwise
%   product (otype=[1,2]) of operators applied to different copies of the
%   state, as in the nonlinear terms of a PIE, which is a factor in an
%   operator taking a Kronecker product (otype=[2,2]). Its parameters are
%   then extracted by recursion, continuing the numbering of the dummy
%   variables for the factors preceding it: the factors of a pointwise
%   product share the output dummy variable, but have a different input
%   dummy variable, which is what the products of the parameters use;
%
% OUTPUTS
% - params:     struct (quadPolyCell-style; see otimes_params.m) with
%               fields 'C' (n_ints{1} x ... x n_ints{N3} cell of sparse
%               coefficient matrices), 'Zs'/'Zt' (shared monomial bases),
%               'dim', 'ns'/'nt';
% - intmat:     1 x N3 cell indicating the intervals of integration along
%               each of the 3-PI type variables in the operator. 
%               Element intmat{i} is a n_ints{i} x 2*dtot(i) array of 
%               integers, specifying the upper and lower bounds of the
%               interval of integration for each variable
%

% Extract the individual factors defining the operator
ops = Rop.ops;
otypeR = Rop.otype;
if otypeR(2)==1
    error("Pointwise product of operators along the input dimension is not supported.")
end
[ntrms,nfctrs] = size(ops);
N3 = numel(vars.S3);

% Declare initial intmats for each variable
intmat = repmat({zeros(1,0)},1,numel(vars.S3));

% Declare starting degrees for each variable
if nargin<6 || isempty(degs_in)
    degs_in = ones(1,numel(vars.in));
end
if nargin<5 || isempty(degs_out)
    degs_out = ones(1,numel(vars.out));
end
if nargin<7 || isempty(prune)
    prune = true;
end

% Determine, for each factor column j, which 3-PI "slots" (per S3 direction 
% k) are structurally zero in EVERY term
is_zero_shared = cell(1,nfctrs);                                          
for j=1:nfctrs                                                           
    is_zero_j = [];                                                      
    is_consistent = true;                                                  
    for i=1:ntrms                                                        
        if isa(ops{i,j},'tensopvar') || isa(ops{i,j},'intop')      
            % Pruning below only applies to the plain-operator branch;  
            % unsafe/unmodeled to aggregate sparsity for this column.  
            is_consistent = false;                             
            break                                                    
        end                                                    
        vars_ij_S3_pp = ops{i,j}.vars_S3;                       
        [~,~,idcs_S3_pp] = intersect(vars_ij_S3_pp,vars.S3,'stable');   
        params_tmp_pp = ops{i,j}.params;                             
        sz_params_pp = size(params_tmp_pp,1:numel(idcs_S3_pp));       
        sz_full_pp = ones(1,max(N3,1));                             
        sz_full_pp(idcs_S3_pp) = sz_params_pp;                   
        params_tmp_pp = reshape(params_tmp_pp,[sz_full_pp,1]); 
        is_zero_ij_pp = cellfun(@(a) nnz(a)==0,params_tmp_pp);        
        if isempty(is_zero_j)                                   
            is_zero_j = is_zero_ij_pp;                            
        elseif ~isequal(size(is_zero_j),size(is_zero_ij_pp))  
            % Different terms disagree on this factor's own structure --  
            % should not happen for a well-formed tensopvar
            error("Same factor in different terms of 'tensopvar' should map between the same function spaces.")                                            
        else                                                        
            is_zero_j = is_zero_j & is_zero_ij_pp;   
        end                                               
    end                    
    if is_consistent             
        is_zero_shared{j} = is_zero_j;
    else                                
        is_zero_shared{j} = [];   % signals "prune nothing" below    
    end                                                  
end             

if ~prune
    % Keep all integrals (pruning of a nested operator could differ between terms)
    is_zero_shared = cell(1,nfctrs);
end

% Determine the integrals and parameters defining the operator
for i=1:ntrms
    degs_in_j = degs_in;
    degs_out_j = degs_out;
    for j=1:nfctrs
        if isa(ops{i,j},'tensopvar')
            % % The operator is a nested 'tensopvar' object
            opN = ops{i,j};
            if ~opN.chkval
                opN = initialize(opN);
            end
            if opN.otype(2)==1 || (otypeR(1)==1 && opN.otype(1)==2 && size(opN.ops,2)>1)
                error("Pointwise product of the input dimensions, or of a tensor product of operators, in a nested tensor-PI operator is not supported.")
            end
            % Extract the parameters and integrals by recursion, starting at the
            % dummy variables for the factors before. A single term in the
            % operator can retain the pruning of zero integrals of the nested
            % operator, as there are no other terms to combine it with
            [params_ij,intmat_ij] = get_params(opN,vars,varnames,dtot,degs_out_j,degs_in_j,prune && ntrms==1);
            if i==1
                % Combine the integrals with those of the preceding factors
                for k=1:N3
                    intmat{k} = [repelem(intmat{k},size(intmat_ij{k},1),1), repmat(intmat_ij{k},size(intmat{k},1),1)];
                end
            end
            % Move to the next degree of the variables, past those of the nested operator
            degs_N = opN.deps;
            [~,idcs_N,idcs_out] = intersect(opN.vars.out,vars.out,'stable');
            if otypeR(1)==2
                degs_out_j(idcs_out) = degs_out_j(idcs_out)+degs_N.out(idcs_N);
            end
            [~,idcs_N,idcs_in] = intersect(opN.vars.in,vars.in,'stable');
            if otypeR(2)==2
                degs_in_j(idcs_in) = degs_in_j(idcs_in)+degs_N.in(idcs_N);
            end
        elseif isa(ops{i,j},'intop')
            error("This case is currently not supported.")
        else
            % % The operator is specified as 'sopvar' object
            % Extract the variables defining the operator
            opdim = size(ops{i,j});
            vars_ij = ops{i,j}.vars;
            vars_ij_S3 = ops{i,j}.vars_S3;
            [~,~,idcs_out] = intersect(vars_ij.out,vars.out,'stable');
            [~,~,idcs_in] = intersect(vars_ij.in,vars.in,'stable');
            [~,~,idcs_S3] = intersect(vars_ij_S3,vars.S3,'stable');

            % Set the output variable names
            varnames_out_ij = cell(1,numel(vars_ij.out));
            ZL_ij = ops{i,j}.ZL;
            for k=1:numel(vars_ij.out)
                idx = idcs_out(k);
                dk_out = degs_out_j(idx);
                varnames_out_ij{k} = varnames.out{idx}{dk_out};
            end

            % Set the input variables and monomials
            varnames_in_ij = cell(1,numel(vars_ij.in));
            ZR_ij = ops{i,j}.ZR;
            for k=1:numel(vars_ij.in)
                idx = idcs_in(k);
                dk_in = degs_in_j(idx);
                varnames_in_ij{k} = varnames.in{idx}{dk_in};
            end

            % Reshape the parameters into a cell of expected dimensions
            params_tmp = ops{i,j}.params;
            sz_params = size(params_tmp,1:numel(idcs_S3));
            sz_full = ones(1,max(N3,1));
            sz_full(idcs_S3) = sz_params;
            
            % Store parameters as a struct
            params_tmp = reshape(params_tmp,[sz_full,1]);
            params_ij = struct();
            params_ij.C = params_tmp;           params_ij.dim = opdim;
            params_ij.Zs = ZL_ij;               params_ij.Zt = ZR_ij;
            params_ij.ns = varnames_out_ij;     params_ij.nt = varnames_in_ij;
            params_ij = set_vars_quadPolyCell(params_ij,[varnames.out{:}],[varnames.in{:}]);

            % Set the associated intmats
            % Prune terms which are defined by all-zero parameters
            is_zero_shared_j = is_zero_shared{j};
            if isempty(is_zero_shared_j)
                % Structure differs across terms (or this factor is itself
                % a nested tensopvar/intop) -- unsafe to prune, keep every
                % slot.
                is_zero_shared_j = false(size(params_tmp)); 
            end
            for k=1:N3
                % Check which parameters are nonzero along direction k in 
                % EVERY term
                rtn_idcs = ~all(is_zero_shared_j,setdiff(1:N3,k));
                sz_param = size(params_ij.C,1:N3);
                params_ij.C = params_ij.C(rtn_idcs & true([sz_param,1]));
                sz_param(k) = nnz(rtn_idcs);
                params_ij.C = reshape(params_ij.C,[sz_param,1]);
                if sz_param(k)==0
                    error("get_params:inconsistent_zero_slot", ...  
                        ['Internal error: factor column %d is structurally ' ...
                         'zero along every slot of direction %d across all terms -- ' ...
                         'this should already have been removed by clean_op.'],j,k) 
                end
                if i==1
                    % Structural intmat is shared across terms -- 
                    % build it once
                    if ismember(k,idcs_S3)
                        % Operator acts as 3-PI operator along this variable
                        dk_out = degs_out_j(ismember(vars.out,vars.S3(k)));
                        intmat_k = [dk_out,dk_out;       % multiplier
                                    0,dk_out;            % lower integral
                                    dk_out,dtot(k)+1]; % upper integral
                        intmat_k = intmat_k(rtn_idcs,:);
                        intmat{k} = [repelem(intmat{k},size(intmat_k,1),1), repmat(intmat_k,size(intmat{k},1),1)];
                    elseif ismember(vars.S3(k),vars_ij.in)
                        % Operator acts as integral along this variable
                        intmat_k = [0,dtot(k)+1];
                        intmat{k} = [repelem(intmat{k},size(intmat_k,1),1), repmat(intmat_k,size(intmat{k},1),1)];
                    end
                end
            end

            % Move to the next degree of the variables
            if otypeR(1)==2
                degs_out_j(idcs_out) = degs_out_j(idcs_out)+any(idcs_out);
            end
            if otypeR(2)==2
                degs_in_j(idcs_in) = degs_in_j(idcs_in)+any(idcs_in);
            end
        end
        % Take tensor product with previous parameters
        if j==1
            params_i = params_ij;
        else
            params_i = otimes_params(params_i,params_ij,otypeR);
        end
    end
    % Combine parameter associated with same integral
    if i==1
        params = params_i;
    elseif isstruct(params)

        if ~isequal(params.Zs,params_i.Zs) || ~isequal(params.Zt,params_i.Zt)
            % Bases differ (in content and/or size) -- unify before summing.
            varnames_ik = struct('in',{params.nt},'out',{params.ns});
            Zs_u = params.Zs;   Zt_u = params.Zt;
            for k=1:numel(params.C)
                qA = quadPoly(params.C{k},  params.Zs,  params.Zt,   ...
                    params.dim,  params.ns,   params.nt);
                qB = quadPoly(params_i.C{k},params_i.Zs,params_i.Zt, ...
                    params_i.dim,params_i.ns, params_i.nt);
                [Csum,Zu] = quad_expand({qA,qB},varnames_ik);
                params.C{k} = Csum{1} + Csum{2};
                Zs_u = Zu.out;   Zt_u = Zu.in;
            end
            params.Zs = Zs_u;   params.Zt = Zt_u;
        else                                                               
            % The bases already agree: no basis union is needed.
            for k=1:numel(params.C)
                params.C{k} = params.C{k} + params_i.C{k};
            end
        end 
    else
        error("get_params:inconsistent_params_type", ...
            "Internal error: 'params' should always be struct-typed at this point.")
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function P = int_params(P,intmat,pvars)
% Integrate P over the domain specified by intmat.
%
% INPUTS
% - P:      'quadPoly' object to integrate
% - intmat: 1 x 2*d array of integers, specifying intervals of integration
%           in `intop` format;
% - pvars:  1 x d+2 cellstr array, specifying the names of the variables
%           defining the limits of integration, where pvars(2:end-1) are
%           the actual variable names, and pvars(1), pvars(end) define the
%           lower and upper boundaries of the domain on which the variable
%           exist;
%
% OUTPUTS
% - P:      'quadPoly' object representing the desired integral of the
%           input P. Integration is performed with respect to each of the
%           variables listed in pvars(2:end-1), where variable pvars(i+1)
%           is integrated over domain pvars(intmat(2*i-1)) up to 
%           pvars(intmat(2*i)). If the upper and lower limits match, the
%           variable is evaluated at pvars(intmat(2*i)) instead (as per
%           'intop' format);

% Split lower and upper boundaries of the domains
int_order = get_int_order(intmat);
Llist = intmat(:,1:2:end);
Ulist = intmat(:,2:2:end);
d = numel(Llist);

[~,var_order] = sort(int_order);
for j=1:d
    % Determine the interval of integration for variable si_k
    idx = var_order(j);
    var = pvars(idx+1);
    L = pvars(Llist(idx)+1);
    U = pvars(Ulist(idx)+1);
    % Perform the integration
    if isequal(Llist(idx),Ulist(idx)) && ~ismember(Llist(idx),[0;d+1])
        % Multiplier operator
        P = subs(P,var,U);
    else
        % Integral operator
        P = int(P,var,L,U);
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [coeffs,ZL,ZR] = quad_expand_quadPoly(params,varnames)
% Express all parameters in the cell in terms of the same basis of
% monomials Z, so that params{k} = ZL*coeffs{k}*ZR;
%
% INPUTS
% - params:     cell of 'quadPoly' objects, or a struct representing a
%               'quadPolyCell' object, with the parameters;
% - varnames:   struct with fields 'in' and 'out', specifying the names of
%               the input and output variables of the parameters;
%
% OUTPUTS
% - coeffs:     cell of structs with fields i,j,v,m,n, representing the
%               sparse coefficient matrices by their nonzero values and the
%               row and column indices in which they appear;
% - ZL,ZR:      'quadPoly' objects representing the monomial bases in the
%               output variables, as a row, and in the input variables, as
%               a column;
%

% Get the coefficients and monomial degrees defining the operator
if isa(params,'cell')
    [coeffs,Zdegs] = quad_expand(params,varnames);
else
    % Parameters are already expressed in terms of a unified basis
    coeffs = params.C;
    ns = params.ns;     nt = params.nt;
    if ~isequal(reshape(ns,1,[]),reshape(varnames.out,1,[])) ...
            || ~isequal(reshape(nt,1,[]),reshape(varnames.in,1,[]))
        error("Input and output variables of parameters should match specified variables.")
    end
    Zs = params.Zs;     Zt = params.Zt;
    Zdegs = struct('in',{Zt},'out',{Zs});
end
if ~isa(coeffs,'cell')
    coeffs = {coeffs};
end

% Store the coefficients only by the nonzero values
for k=1:numel(coeffs)
    [ridcs,cidcs,vals] = find(coeffs{k});
    [m,n] = size(coeffs{k});
    coeffs{k} = struct();
    coeffs{k}.i = ridcs;
    coeffs{k}.j = cidcs;
    coeffs{k}.v = vals;
    coeffs{k}.m = m;
    coeffs{k}.n = n;
end

% Declare 'quadPoly' objects representing the left- and right-bases
ZR = quadPoly(1);
if ~isempty(Zdegs.in)
    n = 1;
    nZ_arr = cellfun(@(a) numel(a),Zdegs.in);
    nZ = prod(nZ_arr);
    sz_coeffs = n*nZ;
    ZR = quadPoly(speye(sz_coeffs),cell(1,0),Zdegs.in,[sz_coeffs,n],cell(1,0),varnames.in);
end
ZL = quadPoly(1);
if ~isempty(Zdegs.out)
    m = 1;
    nZ_arr = cellfun(@(a) numel(a),Zdegs.out);
    nZ = prod(nZ_arr);
    sz_coeffs = m*nZ;
    ridcs = repmat(1:m,nZ,1);
    cidcs = (1:nZ)'+(0:nZ-1)'*nZ +(0:m-1)*nZ^2;
    Cmat = sparse(ridcs(:),cidcs(:),1);
    ZL = quadPoly(Cmat,cell(1,0),Zdegs.out,[m,sz_coeffs],cell(1,0),varnames.out);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function IZ = IkronZ(n,Z)
% Compute the Kronecker product In o Z for a 'quadPoly' object Z,
% representing a RIGHT monomial basis, so that Z.C is the identity

nZ = prod(cellfun(@(a) numel(a),Z.Zt));
IZ = quadPoly(speye(n*nZ),Z.Zs,Z.Zt,n*size(Z),Z.ns,Z.nt);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function B = mtimes_kronR(A,C)
% Compute the product A*(Im o C) for A a coefficient matrix or 'dpvar'
% object, C a sparse coefficient matrix, and o the Kronecker product
% NOTE: A and C should be specified as structs, representing sparse
% matrices (or dpvar objects) by the nonzero values and the row/column
% indices in which they appear

% Get the dimension of the coefficient matrices and identity matrix
n = C.m;    p = C.n;
r = A.m;    mn = A.n;
if r==1 && mn==1
    % Separately deal with scalar multiplication
    B = A;
    B.i = repmat(C.i(:),numel(A.i),1);
    B.j = repmat(C.j(:),numel(A.j),1);
    B.v = reshape(C.v(:).*A.v(:)',[],1);
    B.m = n;
    B.n = p;
    return
end
m = mn/n;
if round(m)~=m
    error("Dimensions don't match.")
end

% At most one nonzero in each row of C (e.g. the map from the monomials of the
% parameters of B to the unified basis): the product only relabels the columns
% of A, so the join below is skipped
if ~isempty(C.i) && max(accumarray(C.i(:),1,[n,1]))<=1
    colC = zeros(n,1);      valC = zeros(n,1);
    colC(C.i) = C.j;        valC(C.i) = C.v;
    cidcs = A.j(:);         blk_idcs = ceil(cidcs/n);
    zidcs = cidcs - (blk_idcs-1)*n;
    is_rtn = colC(zidcs)>0;         % rows of C without nonzeros give 0
    ridcs = A.i(:);         vals = A.v(:);
    B = A;
    B.i = ridcs(is_rtn);
    B.j = (blk_idcs(is_rtn)-1)*p + colC(zidcs(is_rtn));
    B.v = vals(is_rtn).*valC(zidcs(is_rtn));
    B.m = r;
    B.n = p*m;
    return
end

% Split A into blocks of columns of size n
ridcs = A.i(:);
cidcs = A.j(:);
blk_idcs = ceil(cidcs/n);
% Vertically concatenate the blocks
ridcs = ridcs + (blk_idcs-1)*r;
cidcs = cidcs - (blk_idcs-1)*n;
A.i = ridcs;
A.j = cidcs;
A.m = r*m;
A.n = n;
% Multiply each block with C
AC = mtimes_spstruct(A,C);
ridcs = AC.i;
cidcs = AC.j;
vals = AC.v;
% Re-split into blocks of rows
blk_idcs = ceil(ridcs/r);
% Horizontally concatenate the blocks
ridcs = ridcs - (blk_idcs-1)*r;
cidcs = cidcs + (blk_idcs-1)*p;

% Return the sparse coefficient matrix/dpvar object as a struct
B = A;
B.i = ridcs(:);
B.j = cidcs(:);
B.v = vals(:);
B.m = r;
B.n = p*m;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function C = mtimes_kronL(A,B,n)
% Compute the product A*(B o I_n) of sparse coefficient matrices, specified
% as structs, with o the Kronecker product, without forming (B o I_n),
% which has nnz(B)*n nonzeros: column (b-1)*n+z of A (block b, monomial z)
% is mapped to the columns (j-1)*n+z for the nonzeros (b,j) of B, so the
% join is over the rows of B only

% Multiplication with 0 returns 0
if isempty(A.v) || isempty(B.v)
    C = A;
    C.i = zeros(0,1);
    C.j = zeros(0,1);
    C.v = zeros(0,1);
    C.m = A.m;
    C.n = n*B.n;
    return
end
if A.n~=n*B.m
    error("Dimensions don't match.")
end

% Column vectors throughout: the indices may be stored as rows, and combining
% rows with columns would broadcast
ridcsA = A.i(:);    cidcsA = A.j(:);    valsA = A.v(:);
ridcsB = B.i(:);    cidcsB = B.j(:);    valsB = B.v(:);
blk_idcs = ceil(cidcsA/n);
zidcs = cidcsA - (blk_idcs-1)*n;
[rtn_idcs1,rtn_idcs2] = find_common_idcs(blk_idcs,ridcsB,B.m);

% Declare the product as a struct
C = A;
C.i = ridcsA(rtn_idcs1);
C.j = (cidcsB(rtn_idcs2)-1)*n + zidcs(rtn_idcs1);
C.v = valsA(rtn_idcs1).*valsB(rtn_idcs2);
C.m = A.m;
C.n = n*B.n;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function C = mtimes_spstruct(A,B)
% Take the matrix product of two sparse matrices, specified as structs with
% fields i, j, v, m, n;

% The first factor can be 1, the initial value of a product
if isequal(A,1)
    C = B;
    return
end
% Multiplication with 0 returns 0
if isempty(A.v) || isempty(B.v)
    C = A;
    C.i = zeros(0,1);
    C.j = zeros(0,1);
    C.v = zeros(0,1);
    C.m = A.m;
    C.n = B.n;
    return
end

% Extract the nonzero values of each matrix and associated row and column 
% indices
ridcsA = A.i(:);    cidcsA = A.j(:);    valsA = A.v(:);
ridcsB = B.i(:);    cidcsB = B.j(:);    valsB = B.v(:);

% Find all values for which column indices of A match row indices of B
[rtn_idcs1, rtn_idcs2] = find_common_idcs(cidcsA,ridcsB,A.n);
ridcs = ridcsA(rtn_idcs1);
cidcs = cidcsB(rtn_idcs2);
vals = valsA(rtn_idcs1).*valsB(rtn_idcs2);

% Declare the product as a struct
C = A;
C.i = ridcs(:);
C.j = cidcs(:);
C.v = vals(:);
C.m = A.m;
C.n = B.n;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function C = plus_spstruct(A,B)
% Take the sum of two sparse matrices, specified as structs with
% fields i, j, v, m, n;

% The first term can be 0, the initial value of a sum
if isequal(A,0)
    C = B;
    return
end

% Extract the nonzero values of each matrix and associated row and column
% indices
ridcsA = A.i(:);    cidcsA = A.j(:);    valsA = A.v(:);
ridcsB = B.i(:);    cidcsB = B.j(:);    valsB = B.v(:);

% Sum coefficients in the same row and column
% If A.n is at most 20 times the number of nonzeros, accumulate in a sparse
% matrix; otherwise, sort the linear indices (m*n < 2^53)
ridcs = [ridcsA;ridcsB];    cidcs = [cidcsA;cidcsB];    vals = [valsA;valsB];
% Sum repeated positions only once the stored values outnumber those after the
% last sum by a factor 2; until then, the entries are only concatenated
nref = numel(valsA);
if isfield(A,'nref'),   nref = min(A.nref,nref);    end
if numel(vals)<=2*max(nref,numel(valsB))
    C = A;      C.i = ridcs;    C.j = cidcs;    C.v = vals;     C.nref = nref;
    return
end
if A.n<=20*numel(vals)
    [ridcs,cidcs,vals] = find(sparse(ridcs,cidcs,vals,A.m,A.n));
else
    [lin,~,idcs_lin] = unique((cidcs-1)*A.m + ridcs);
    vals = accumarray(idcs_lin,vals);
    ridcs = mod(lin-1,A.m) + 1;
    cidcs = (lin-ridcs)/A.m + 1;
    rtn_idcs = vals~=0;
    ridcs = ridcs(rtn_idcs);        cidcs = cidcs(rtn_idcs);        vals = vals(rtn_idcs);
end

% Declare the sum as a struct
C = A;
C.i = ridcs(:);
C.j = cidcs(:);
C.v = vals(:);
C.nref = numel(C.v);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function intmat = switch_ints(intmat,varsA)
% INTMAT = SWITCH_INTS(INTMAT,VARA) takes a set of intervals of
% integration for some variables s1,...,sd in [a,b], and declares a new set
% of intervals that add up to the same integral, but now with none of the
% limits of integration involving varsA(end)
%
% INPUTS
% - intmat: 1 x 2*d array of integers ranging from 0 to d+1, with elements
%           (i,i+1) for in 1:2:2*d-1 specifying the lower and upper limits
%           of the interval of integration of variable i. Here
%           intmat(2*i-1,2*i) = [l,u] implies si is integrated over
%           [sl,su], where we set s0=a and sd = b;
% - varsA:  1 x m array for m<=d specifying which variables should not
%           appear in the limit of integration. If m>1, only the last
%           variable is actively removed from the limits of integration,
%           and the variables varA(1:end-1) are assumed to have been
%           removed in an earlier call to this function, and are therefore
%           ignored here;
%
% OUTPUTS
% - intmat: nprms x 2*d array of integers, specifying intervals of
%           integration, so that the sum of the integrals (with the same
%           kernel) returns the same integral as specified by the input
%           intmat, but now with none of the limits of integration
%           involving varA(end);
%           


% Determine which variables are bigger and smaller than varA in the 
% specified integral
varA = varsA(end);
d = size(intmat,2)/2;
Llims = intmat(1:2:end);    Ulims = intmat(2:2:end);
is_geqA = Llims==varA;      is_leqA = Ulims==varA;
is_geqA(varA) = true;       is_leqA(varA) = true;
is_geqA(varsA(1:end-1)) = false;
is_leqA(varsA(1:end-1)) = false;
is_geq_i = is_geqA;         is_leq_i = is_leqA;
for i=1:d
    if ~any(is_geq_i)
        break
    end
    % Check which variables greater than varA have variable upper limit
    Ulims_i = Ulims(is_geq_i);
    is_varU_i = Ulims_i<=d;
    % Variable upper limit must be greater than varA as well
    is_geq_i(Ulims_i(is_varU_i)) = true;
    is_geq_i = is_geq_i & ~is_geqA;
    is_geqA = is_geqA | is_geq_i;    
end
for i=1:d
    if ~any(is_leq_i)
        break
    end
    % Check which variables smaller than varA have variable upper limit
    Llims_i = Llims(is_leq_i);
    is_varL_i = Llims_i>=1;
    % Variable lower limit must be smaller than varA as well
    is_leq_i(Llims_i(is_varL_i)) = true;
    is_leq_i = is_leq_i & ~is_leqA;
    is_leqA = is_leqA | is_leq_i;    
end
is_geqA(varA) = false;      is_leqA(varA) = false;

% Determine which variables must be involved in the expansion of the
% integrals
is_varA = false(1,d);
is_varA(varA) = true;
expnd_vars = (is_varA | is_leqA | is_geqA) & Llims~=Ulims;
expnd_vars_i = expnd_vars;
for i=1:d
    % Check which upper and lower limits are variable
    if ~any(expnd_vars_i) || all(expnd_vars)
        break
    end
    Llims_expnd = Llims(expnd_vars_i);
    Ulims_expnd = Ulims(expnd_vars_i);
    is_varL_i = Llims_expnd>=1;
    is_varU_i = Ulims_expnd<=d;
    % Add variable upper and lower limits to the list to expand
    expnd_vars_i(union(Llims_expnd(is_varL_i),Ulims_expnd(is_varU_i))) = true;
    expnd_vars_i = expnd_vars_i & ~expnd_vars;
    expnd_vars = expnd_vars | expnd_vars_i;    
end
var_nums = find(expnd_vars);
n = numel(var_nums);

% We consider the variables specified by expnd_vars, and integrate them all
% over the full domain, [a,b], represented by [0,di+1].
% We split this full domain into n! subdomains, considering all possible
% orderings of the n variables.
% We keep track of the ordering by assigning a magnitude to each variable,
% 'ord', so that ord(i)>=ord(j) implies si>=sj. For example, if n=3,
%   ord(i,:) = [3,1,2] implies s2 <= s3 <= s1 in term i
ords = enumerate_valid_orders(Llims,Ulims,expnd_vars,var_nums,n);

% Determine the lower and upper domains bounding each variable
nprms = size(ords,1);
[~,var_idcs] = sort(ords,2);                                    % var_idcs(i,:) = [k1,...,kn] implies s_k1 <= ... <= s_kn in term i
Ulims_tmp = [var_nums(var_idcs(:,2:end)),(d+1)*ones(nprms,1)];  % integrate s_ki up to s_k{i+1}
% Account for the order of the variables
cidcs = ords;
ridcs = repmat((1:nprms)',1,n);
lidcs = sub2ind(size(ords),ridcs,cidcs);
Ulims = Ulims_tmp(lidcs);
% Initialize a new intmat, integrating up to the obtained upper limits
intmat(:,2*var_nums-1) = 0;
intmat = repmat(intmat,nprms,1);
intmat(:,2*var_nums) = Ulims;
% If any variable is integrated up to upper limit varA, 
% switch the order of integration
%   int_{0}^{sC} dsA int_{0}^{sA} dsB
%   = int_{0}^{sC} dsB int_{sB}^{sC} dsA
[trm_nums,var_idx_leqA] = find(Ulims==varA);
var_num_leqA = var_nums(var_idx_leqA);
intmat(trm_nums,2*varA-1) = var_num_leqA;       % set new lower domain of sA
Ulims_leqA = intmat(trm_nums,2*varA);
lidcs = sub2ind(size(intmat),trm_nums(:),2*var_num_leqA(:));
intmat(lidcs) = Ulims_leqA;                     % set new upper domain of sB

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function intmat = canonical_intmat(intmat)
% Express the integrals as int_{a}^{si}
%
% INPUTS
% - intmat: n x 2*d array of integers specifying intervals of integration
%           per 'intop' format'. At least one of the integrals must be over
%           the full spatial interval [a,b] of the variable 
%           (intmat(i,2*j-1)=0 and intmat(i,2*j)=d+1 for some j=j(i) in 
%           each row i);
%
% OUTPUTS
% - intmat: n x 2*d array of integers, specifying the same full domain of
%           integration but now expressed solely in terms of integrals with
%           lower boundary corresponding to the lower boundary a of the
%           interval;
%
% NOTES
% - Only works for the particular setting considered in `mtimes_tensopvar`,
% may not work for more general cases;
%
% EXAMPLE
% For intmat = [0,4,0,3,1,4], we have the integral
%   int_{a}^{b} int_{s1}^{b} int_{a}^{s3} ds2 ds3 ds1
% We can express this equivalently as
%   int_{a}^{b} int_{a}^{s3} int_{a}^{s3} ds2 ds1 ds3
% defined by intmat = [0,3,0,3,0,4]

[ntrms,d] = size(intmat);
d = d/2;
for i=1:ntrms
    int_order = get_int_order(intmat(i,:));
    [~,var_order] = sort(int_order);
    for j=1:d
        var = var_order(d-j+1);
        Llim = intmat(i,2*var-1);
        Ulim = intmat(i,2*var);
        if Llim==Ulim
            % Impose the delta function: var==Ulim
            intmat(i,intmat(i,:)==var) = Ulim;
            continue
        elseif Llim==0
            % We are already integrating from a up to some variable
            continue
        end
        % Get the domain of integration of variable Llim
        Llim_L = intmat(i,2*Llim-1);
        Ulim_L = intmat(i,2*Llim);
        if Ulim_L==Ulim
            % If the upper limits match, we switch the order of integration
            % int_{a}^{y} int_{x}^{y} ... dz dx
            % = int_{a}^{y} int_{a}^{z} ... dx dz
            intmat(i,2*var-1) = Llim_L;
            intmat(i,2*Llim) = var;
        end
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function out = hcat_cellstr(c)
% OUT = HCAT_CELLSTR(C) horizontally concatenates a 1 x K cell of cellstr
% objects, C{1},...,C{K}, always returning a cellstr -- unlike the plain
% '[c{:}]', which degrades to a 0x0 double (not a 1x0 cell) whenever C
% itself has zero elements, since '[ ]' with no expanded arguments
% defaults to double.
%
% INPUTS
% - c:      1 x K cell, with each element a 1 x ni cellstr object;
%
% OUTPUTS
% - out:    1 x N cellstr, where N = sum_{i=1}^{K} ni, corresponding to the
%           horizontal concatenation of all the cellstr objects in c;

if isempty(c)
    out = cell(1,0);
else
    out = [c{:}];
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function ords = enumerate_valid_orders(Llims,Ulims,expnd_vars,var_nums,n)
% ORDS = ENUMERATE_VALID_ORDERS(LLIMS,ULIMS,EXPND_VARS,VAR_NUMS,N) returns
% all possible orderings of the variables VARS(EXPND_VARS) that are 
% consistent with the domains of integration specified by LLIMS and ULIMS,
% in the function `switch_ints`.
%
% INPUTS
% - Llims,Ulims:  1 x d arrays, switch_ints' own (unmodified) lower/upper
%                 limit arrays;
% - expnd_vars:   1 x d logical, switch_ints' own (unmodified) closure;
% - var_nums:     1 x n array, find(expnd_vars)
% - n:            numel(var_nums);
%
% OUTPUTS
% - ords:         nprms x n array as described above; nprms==0 is possible
%                 (an unsatisfiable/contradictory constraint set), matching
%                 the original, which would likewise filter every row away.
%
% DETAILS
% We consider the variables specified by expnd_vars, and integrate them all
% over the full domain, [a,b], represented by [0,di+1].
% We split this full domain into n! subdomains, considering all possible
% orderings of the n variables.
% We keep track of the ordering by assigning a magnitude to each variable,
% 'ord', so that ord(i)>=ord(j) implies si>=sj. For example, if n=3,
%   ord(i,:) = [3,1,2] implies s2 <= s3 <= s1 in term i

% For each expand variable (by its position 1..n in var_nums), find the
% position of its Llims/Ulims reference WITHIN var_nums, when that
% reference is itself an expand variable (0 otherwise, meaning "no
% constraint" -- the reference is a fixed/non-expand variable or a domain
% endpoint, exactly the case the original's ords_lwr/ords_upr default-fill
% of 0/(d+1) always satisfies).
d = numel(expnd_vars);
pred = zeros(1,n);      succ = zeros(1,n);
for c = 1:n
    lv = Llims(var_nums(c));
    if lv>=1 && lv<=d && expnd_vars(lv)
        pred(c) = find(var_nums==lv,1);
    end
    uv = Ulims(var_nums(c));
    if uv>=1 && uv<=d && expnd_vars(uv)
        succ(c) = find(var_nums==uv,1);
    end
end

% Build a "must come before" adjacency and each node's indegree, for
% Kahn-style backtracking: a node is "ready" to be placed next once all its
% prerequisites are already placed. pred(c)=p and succ(p)=c can assert the
% IDENTICAL constraint "p before c" (discovered from two different sides:
% c's own Llims, and separately p's own Ulims) -- collect each node's
% prerequisite set and deduplicate with unique() before computing indegree,
% or such a pair double-counts c's indegree and c can never become "ready"
% (confirmed by direct trace: d=2 chain Llims=[0,1] Ulims=[2,3] gives
% pred(2)=1 AND succ(1)=2, the same "1 before 2" constraint from both
% sides).
parents = cell(1,n);
for c = 1:n
    if pred(c)~=0
        parents{c}(end+1) = pred(c);
    end
    if succ(c)~=0
        s = succ(c);
        parents{s}(end+1) = c;
    end
end
children = cell(1,n);
indeg0 = zeros(1,n);
for c = 1:n
    % unique([]) returns a 0x1 (column), not 1x0 -- 'for p = (0x1 empty)'
    % iterates ONCE with p empty (only a 1x0 row iterates zero times), so
    % force row shape or an empty prerequisite set spuriously adds an edge
    % to an empty index (confirmed: this exact issue caused every closure
    % with an unconstrained variable to find zero valid orderings, i.e.
    % every real case -- every variable has at least one free bound).
    parents{c} = reshape(unique(parents{c}),1,[]);
    indeg0(c) = numel(parents{c});
    for p = parents{c}
        children{p}(end+1) = c;
    end
end

rows = cell(0,1);
seq = zeros(1,n);
placed = false(1,n);
indeg = indeg0;
backtrack(1);

if isempty(rows)
    ords = zeros(0,n);
else
    ords = vertcat(rows{:});
end

    function backtrack(pos)
        % Nested function: shares seq/placed/indeg/rows with the parent.
        if pos>n
            row = zeros(1,n);
            row(seq) = 1:n;
            rows{end+1,1} = row;
            return
        end
        avail = find(~placed & indeg==0);
        for v = avail
            seq(pos) = v;
            placed(v) = true;
            ch = children{v};
            indeg(ch) = indeg(ch)-1;
            backtrack(pos+1);
            indeg(ch) = indeg(ch)+1;
            placed(v) = false;
        end
    end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [Zgrid,Zdegs] = grid_basis(Zin)
% [ZGRID,ZDEGS] = GRID_BASIS(ZIN) returns for the cell ZIN of arrays of
% monomial degrees (one array for each variable, one column for each of its
% dummy variables) the smallest tensor product over the dummy variables
% that contains all monomials: ZDEGS{k} is the column vector of the degrees
% appearing in dummy variable k (in order, for all variables), and ZGRID{i}
% the array of all combinations of the degrees ZDEGS of the dummy variables
% of variable i, the first dummy variable changing slowest.

N = numel(Zin);
Zgrid = cell(1,N);
Zdegs = cell(1,0);
for i=1:N
    di = size(Zin{i},2);
    Dk = cell(1,di);
    for k=1:di
        Dk{k} = unique(Zin{i}(:,k));
    end
    Zgrid{i} = cartesian_grid(Dk);
    Zdegs = [Zdegs,Dk];
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Zgrid = cartesian_grid(Dk)
% ZGRID = CARTESIAN_GRID(DK) returns the array of all combinations of the
% degrees in the cell DK of column vectors, one column for each of the
% cells, the first changing slowest.

Zgrid = zeros(1,0);
for k=1:numel(Dk)
    Dk_k = Dk{k}(:);
    Zgrid = [repelem(Zgrid,numel(Dk_k),1), repmat(Dk_k,size(Zgrid,1),1)];
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function A = sparse2struct(P,zmap,nZgrid)
% A = SPARSE2STRUCT(P,ZMAP,NZGRID) represents the sparse coefficient matrix
% P, with columns (j-1)*nZ+z, by the struct of its nonzeros with fields
% i,j,v,m,n, with the monomial z moved to the monomial ZMAP(z) in a basis
% of NZGRID monomials.

nZ = numel(zmap);
[i,j,v] = find(P);
jb = ceil(j/nZ);
z = j - (jb-1)*nZ;
A = struct('i',i(:),'j',(jb(:)-1)*nZgrid + zmap(z(:)),'v',v(:),'m',size(P,1),'n',(size(P,2)/nZ)*nZgrid);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
