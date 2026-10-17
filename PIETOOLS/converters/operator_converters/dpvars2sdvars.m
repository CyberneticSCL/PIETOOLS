function [A,B,Z,dvars] = dpvars2sdvars(params,vars,Z,chk_vars)
% [COEFFS,Z,DVARS] = quad_expand(PARAMS,VARS,Z,CHK_VARS)
% takes a cell of parameters and expresses them all in terms of a common
% set of monomials and decision variables in the specified variables in
% sdopvar format,
%   PARAMS{k} = (Im o ZL{1}(s1) o ... o ZL{M}(sM))^T
%                   unvec(A{k} + B{k}'*DVARS')
%                       (In o ZR{1}(t1) o ... o ZR{N}(tN))
% 
% INPUTS
% - params:     cell of elements of type 'double', 'polynomial', 'dpvar',
%               or 'quadPoly';
% - varnames:   struct with fields
%         in:   1 x N cellstr, specifying the names of the input variables
%               appearing in the parameters;
%        out:   1 x M cellstr, specfiying the names of the output variables
%               appearing in the parameters;
% - Z:          (optional) struct with fields
%         in:   1 x N cell, with each element a di x 1 array specifying the
%               degrees of the monomials in varnames.in{i} in terms of
%               which to express the parameters;
%        out:   1 x M cell, with each element a pi x 1 array specifying the
%               degrees of the monomials in varnames.out{i} in terms of
%               which to express the parameters;
%               If not specified, monomial degrees will be determined based
%               on what monomials are encountered in the parameters.
% - chk_vars:   (optional) logical values indicating whether (true) or not
%               (false) to update the parameters to express them in terms
%               of the same variables. If all parameters are already
%               expressed in terms of the same variables, set to false to
%               reduce computational complexity. Defaults to true;
%
% OUTPUTS
% - A:          cell of the same dimensions as params, expressing for each
%               of the parameters the coefficient matrix defining the
%               constant (no decision variables) term in the quadratic 
%               format;
% - B:          cell of the same dimensions as params, expressing for each
%               of the parameters the coefficient matrix defining the
%               decision variable term in the quadratic format;
% - Z:          struct with fields
%         in:   1 x N cell, with each element a di x 1 array specifying the
%               degrees of the monomials in varnames.in{i} in terms of
%               which the parameters are expressed;
%        out:   1 x M cell, with each element a pi x 1 array specifying the
%               degrees of the monomials in varnames.out{i} in terms of
%               which the parameters are expressed;
% - dvars:      K x 1 cellstr specifying the names of the decision
%               variables defining the parameters, so that d = dvars;
%
% NOTES
% - In the quadratic form, the parameters are expressed as
%   params{k} = (Im o s1.^Z.out{1} o ... o sM.^Z.out{M})^T * Ccell{k} *
%                   (In o t1.^Z.in{1} o ... o tN.^Z.in{N})
% - The elements of Ccell will be sparse matrices, unless the input params
% involve decision variables, in which case the elements will be 'dpvar'
% objects.
%   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - quad_expand
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
% DJ, 10/05/2026: Initial coding;

% Check that the parameters are appropriate
use_cell = true;
if isa(params,'double') || isa(params,'polynomial') || isa(params,'dpvar') || isa(params,'quadPoly')
    params = {params};
    use_cell = false;
elseif ~isa(params,'cell')
    error("Parmeters should be specified as 'cell' of 'dpvar' or 'quadPoly' objects.")
end

% Check that the variables are properly specified
if ~isa(vars,'struct') || ~isfield(vars,'in') || ~isfield(vars,'out')
    error("Variable names must be specified as struct with fields 'in' and 'out'.")
end

% Extract the variable names
varnames1 = vars.out;       M = numel(varnames1);
varnames2 = vars.in;        N = numel(varnames2);
if isempty(varnames1)
    varnames1 = cell(1,0);
end
if isempty(varnames2)
    varnames2 = cell(1,0);
end
if ~iscellstr(varnames1) || ~iscellstr(varnames2)
    error("Input and output variable names should be specified as cellstr objects.")
end
varnames_full = [varnames1,varnames2];

% Express all parameters in terms of the same set of variables, and
% build the full left and right monomial bases
if nargin>=4 && ~chk_vars
    % Assume all parameters are specified in terms of the same variables,
    % and monomials are specified
    ZL = Z.out;
    ZR = Z.in;
elseif nargin>=3
    % Use monomial basis specified by user
    ZL = Z.out;
    ZR = Z.in;
    % Make sure all parameters are expressed in terms of the same variables
    for k=1:numel(params)
        param_k = params{k};
        if isa(param_k,'double')
            % Include only degree-0 monomial for constant terms
            degmat_full = zeros(1,M+N);
            params{k} = polynomial(reshape(param_k,1,[]),degmat_full,varnames_full,size(param_k));
        elseif isa(param_k,'quadPoly')
            % Express the parameter in terms of full set of variables
            param_k = set_vars(param_k,varnames1,varnames2);
            params{k} = param_k;
        elseif isa(param_k,'polynomial') || isa(param_k,'dpvar')
            % Express the parameter in terms of full set of variables
            [~,idcs1_new,idcs1_old] = intersect(varnames1,param_k.varname','stable');
            [~,idcs2_new,idcs2_old] = intersect(varnames2,param_k.varname','stable');
            degmat_full = zeros(size(param_k.degmat,1),M+N);
            degmat_full(:,idcs1_new) = param_k.degmat(:,idcs1_old);
            degmat_full(:,M+idcs2_new) = param_k.degmat(:,idcs2_old);
            if isa(param_k,'polynomial')
                params{k} = polynomial(param_k.C,degmat_full,varnames_full,size(param_k),false);
            else
                params{k} = dpvar(param_k.C,degmat_full,varnames_full,param_k.dvarname,size(param_k));
            end
        else
            error("Parameters must be specified as objects of type 'polynomial', 'dpvar', or 'quadPoly'.")
        end
    end
else
    % Declare a basis of monomials including all monomials appearing in the
    % parameters
    ZL = repmat({0},1,numel(varnames1));
    ZR = repmat({0},1,numel(varnames2));
    for k=1:numel(params)
        param_k = params{k};
        if isa(param_k,'double')
            % Include only degree-0 monomial for constant terms
            degmat_full = zeros(1,M+N);
            params{k} = polynomial(reshape(param_k,1,[]),degmat_full,varnames_full,size(param_k));
        elseif isa(param_k,'quadPoly')
            % Express the parameter in terms of full set of variables
            [param_k,s_idcs,t_idcs] = set_vars(param_k,varnames1,varnames2);
            params{k} = param_k;
            % Extract monomial degrees in common variables
            for i=1:numel(s_idcs{1})
                degs_i = param_k.Zs{s_idcs{1}(i)};
                ZL{s_idcs{1}(i)} = unique([ZL{s_idcs{1}(i)}; degs_i]);
            end
            for i=1:numel(t_idcs{1})
                degs_i = param_k.Zt{t_idcs{1}(i)};
                ZR{t_idcs{1}(i)} = unique([ZR{t_idcs{1}(i)}; degs_i]);
            end
        elseif isa(param_k,'polynomial') || isa(param_k,'dpvar')
            % Express the parameter in terms of full set of variables
            [~,idcs1_new,idcs1_old] = intersect(varnames1,param_k.varname','stable');
            [~,idcs2_new,idcs2_old] = intersect(varnames2,param_k.varname','stable');
            degmat_full = zeros(size(param_k.degmat,1),M+N);
            degmat_full(:,idcs1_new) = param_k.degmat(:,idcs1_old);
            degmat_full(:,M+idcs2_new) = param_k.degmat(:,idcs2_old);
            if isa(param_k,'polynomial')
                params{k} = polynomial(param_k.C,degmat_full,varnames_full,size(param_k),false);
            else
                params{k} = dpvar(param_k.C,degmat_full,varnames_full,param_k.dvarname,size(param_k));
            end
            % Extract monomial degrees in common variables
            for i=1:numel(idcs1_new)
                degs_i = unique(param_k.degmat(:,idcs1_old(i)));
                ZL{idcs1_new(i)} = unique([ZL{idcs1_new(i)}; degs_i]);
            end
            for i=1:numel(idcs2_new)
                degs_i = unique(param_k.degmat(:,idcs2_old(i)));
                ZR{idcs2_new(i)} = unique([ZR{idcs2_new(i)}; degs_i]);
            end
        else
            error("Parameters must be specified as objects of type 'polynomial', 'dpvar', or 'quadPoly'.")
        end    
    end
end

% Determine how each of the monomial bases relates to the full bases 0:d
nZL_arr = ones(1,M);   ZL_maps = cell(1,M);
for i=1:M
    nZL_arr(i) = numel(ZL{i});
    is_mon = (0:max(ZL{i}))'==ZL{i}';
    ZL_maps{i} = is_mon*(1:nZL_arr(i))'; 
end
nZR_arr = ones(1,N);   ZR_maps = cell(1,N);
for i=1:N
    nZR_arr(i) = numel(ZR{i});
    is_mon = (0:max(ZR{i}))'==ZR{i}';
    ZR_maps{i} = is_mon*(1:nZR_arr(i))'; 
end
nnZL_arr = cumprod([nZL_arr,1],'reverse');
nnZR_arr = cumprod([nZR_arr,1],'reverse');

% Determine the full set of decision variables across ALL parameters, so
% that every B{k} below is expressed against the same 'dvars'.
dvars = cell(0,1);
for k=1:numel(params)
    if isa(params{k},'dpvar')
        dvars = [dvars; setdiff(params{k}.dvarname(:),dvars,'stable')];
    end
end
ndvars = numel(dvars);

% Finally, express the operators in terms of the unified monomial basis
A = cell(size(params));
B = cell(size(params));
for k=1:numel(params)       % can be parallelized
    param_k = params{k};
    [m,n] = size(param_k);
    if isa(param_k,'double')
        % Convert double to polynomial
        param_k = polynomial(reshape(param_k,1,[]),zeros(1,M+N),varnames_full,[m,n]);
    end
    if isa(param_k,'polynomial')
        % Determine the row index and column index associated with each
        % coefficient defining the polynomial
        [rridcs,ccidcs,vals] = find(param_k.C);
        cidcs = ceil(ccidcs(:)/m);
        ridcs = ccidcs(:) - (cidcs-1)*m;

        % Determine the left-monomial index associated with each
        % coefficient
        degmat = param_k.degmat;
        ZLidcs = ones(numel(ridcs),1);
        for i=1:M
            % Check which monomial in variable i is considered
            Zi_vals = degmat(rridcs,i);
            ZLidcs_i = ZL_maps{i}(Zi_vals+1);
            % Account for Kronecker product with other monomials
            ZLidcs = ZLidcs + (ZLidcs_i-1)*nnZL_arr(i+1);
        end

        % Determine the right-monomial index associated with each
        % coefficient
        ZRidcs = ones(numel(cidcs),1);
        for i=1:N
            % Check which monomial in variable i is considered
            Zi_vals = degmat(rridcs,i+M);
            ZRidcs_i = ZR_maps{i}(Zi_vals+1);
            % Account for Kronecker product with other monomials
            ZRidcs = ZRidcs + (ZRidcs_i-1)*nnZR_arr(i+1);
        end

        % Declare a coefficient matrix acting on the appropriate monomials
        rridcs = (ridcs-1)*nnZL_arr(1) + ZLidcs;
        ccidcs = (cidcs-1)*nnZR_arr(1) + ZRidcs;
        lidcsA = (ccidcs-1)*nnZL_arr(1)*m + rridcs;
        A{k} = sparse(lidcsA,1,vals,m*nnZL_arr(1)*n*nnZR_arr(1),1);
        % We store A as m*n row vector
        
        % Also declare a 0 coefficient matrix acting on the decision 
        % variable terms
        B{k} = sparse([],[],[],ndvars,m*nnZL_arr(1)*n*nnZR_arr(1));
        % We store B rather than Bt, which is a nd x m*n matrix; sized     
        % against the GLOBAL 'dvars'

    elseif isa(param_k,'dpvar')
        % Decompose the operator. ndvars_k/didcs index param_k's OWN
        % dvarname list (that is the block stride of param_k.C); glob_didcs
        % maps that local index into the unified 'dvars' computed above.     
        dvars_k = param_k.dvarname;     ndvars_k = numel(dvars_k);
        [~,glob_didcs] = ismember(dvars_k,dvars);
        degmat = param_k.degmat;        nZ = size(param_k.degmat,1);
        % Determine the row index and column index associated with each
        % coefficient defining the polynomial
        [rridcs,ccidcs,vals] = find(param_k.C);
        ridcs = ceil(rridcs(:)/(ndvars_k+1));          % row number      
        cidcs = ceil(ccidcs(:)/nZ);                    % column number
        didcs = rridcs(:) - (ridcs-1)*(ndvars_k+1);    % decision variable index  
        Zidcs = ccidcs(:) - (cidcs-1)*nZ;              % monomial index

        % Determine the left-monomial index associated with each
        % coefficient
        ZLidcs = ones(numel(rridcs),1);
        for i=1:M
            % Check which monomial in variable i is considered
            Zi_vals = degmat(Zidcs,i);
            ZLidcs_i = ZL_maps{i}(Zi_vals+1);
            % Account for Kronecker product with other monomials
            ZLidcs = ZLidcs + (ZLidcs_i-1)*nnZL_arr(i+1);
        end

        % Determine the right-monomial index associated with each
        % coefficient
        ZRidcs = ones(numel(ccidcs),1);
        for i=1:N
            % Check which monomial in variable i is considered
            Zi_vals = degmat(Zidcs,i+M);
            ZRidcs_i = ZR_maps{i}(Zi_vals+1);
            % Account for Kronecker product with other monomials
            ZRidcs = ZRidcs + (ZRidcs_i-1)*nnZR_arr(i+1);
        end

        % Declare a coefficient matrix associated with the fixed term
        isA = didcs==1;         isB = ~isA;
        rridcsA = (ridcs(isA)-1)*nnZL_arr(1) + ZLidcs(isA);
        ccidcsA = (cidcs(isA)-1)*nnZR_arr(1) + ZRidcs(isA);
        lidcsA = (ccidcsA-1)*nnZL_arr(1)*m + rridcsA;
        valsA = vals(isA);
        A{k} = sparse(lidcsA,1,valsA,m*nnZL_arr(1)*n*nnZR_arr(1),1);
        % We store A as m*n row vector

        % Declare a coefficient matrix acting on the decision variable terms
        rridcsB = (ridcs(isB)-1)*nnZL_arr(1) + ZLidcs(isB);
        ccidcsB = (cidcs(isB)-1)*nnZR_arr(1) + ZRidcs(isB);
        lidcsB = (ccidcsB-1)*nnZL_arr(1)*m + rridcsB;
        valsB = vals(isB);
        glob_didcsB = glob_didcs(didcs(isB)-1);
        B{k} = sparse(glob_didcsB,lidcsB,valsB,ndvars,m*nnZL_arr(1)*n*nnZR_arr(1));
        % We store B rather than Bt, which is a nd x m*n matrix, 
        % against the GLOBAL 'dvars' 

    elseif isa(param_k,'quadPoly')
        % Decompose the operator
        Zs = param_k.Zs;        nZs_arr = cellfun(@(a) numel(a),Zs);
        Zt = param_k.Zt;        nZt_arr = cellfun(@(a) numel(a),Zt);
        nnZs = cumprod([nZs_arr,1],'reverse');
        nnZt = cumprod([nZt_arr,1],'reverse');

        % Determine the row index and column index associated with each
        % coefficient defining the polynomial
        [rridcs,ccidcs,vals] = find(param_k.C);
        ridcs = ceil(rridcs(:)/nnZs(1));               % row number
        cidcs = ceil(ccidcs(:)/nnZt(1));               % column number

        % Determine the left-monomial index associated with each
        % coefficient
        ZLidcs = ones(numel(rridcs),1);
        Zsidcs_i = ridcs;
        for i=1:M
            % Remove contribution from previous monomials
            rridcs = rridcs(:) - (Zsidcs_i-1)*nnZs(i);
            % Remove contribution from subsequent monomials
            Zsidcs_i = ceil(rridcs/nnZs(i+1));
            % Account for new monomial basis
            Zi_vals = Zs{i}(Zsidcs_i);
            ZLidcs_i = ZL_maps{i}(Zi_vals+1);
            % Account for Kronecker product with other monomials
            ZLidcs = ZLidcs + (ZLidcs_i-1)*nnZL_arr(i+1);
        end

        % Determine the right-monomial index associated with each
        % coefficient
        ZRidcs = ones(numel(ccidcs),1);
        Ztidcs_i = cidcs;
        for i=1:N
            % Remove contribution from previous monomials
            ccidcs = ccidcs(:) - (Ztidcs_i-1)*nnZt(i);
            % Remove contribution from subsequent monomials
            Ztidcs_i = ceil(ccidcs/nnZt(i+1));
            % Check which monomial in variable i is considered
            Zi_vals = Zt{i}(Ztidcs_i);
            ZRidcs_i = ZR_maps{i}(Zi_vals+1);
            % Account for Kronecker product with other monomials
            ZRidcs = ZRidcs + (ZRidcs_i-1)*nnZR_arr(i+1);
        end

        % Declare a coefficient matrix acting on the appropriate monomials
        rridcs = (ridcs-1)*nnZL_arr(1) + ZLidcs;
        ccidcs = (cidcs-1)*nnZR_arr(1) + ZRidcs;
        lidcsA = (ccidcs-1)*nnZL_arr(1)*m + rridcs;
        A{k} = sparse(lidcsA,1,vals,m*nnZL_arr(1)*n*nnZR_arr(1),1);
        % We store A as m*n row vector
        
        % Also declare a 0 coefficient matrix acting on the decision 
        % variable terms
        B{k} = sparse([],[],[],ndvars,m*nnZL_arr(1)*n*nnZR_arr(1));

    else
        error("Parameters must be specified as objects of type 'polynomial', 'dpvar', or 'quadPoly'.")
    end
end
% Declare the outputs
Z = struct('in',{ZR},'out',{ZL});
if ~use_cell
    A = A{1};
    B = B{1};
end

end