function [intmat,params] = split_ints(intrvals,params,opts)
% [INTMAT,PARAMS] = SPLIT_INTS(INTRVALS,PARAMS,OPTS) takes a set of
% intervals of integration, and determines an associated set of indicator
% functions to represent these intervals, splitting terms where necessary
% and updating the kernels accordingly.
%
% INPUTS
% - intrvals:   1 x N cell, specifying the intervals over which the N
%               variables must be integrated. Element i is a ntrmsi x
%               2*di array of integers, with columns (2*j-1,2*j)
%               indicating the interval over which variable si_j must be
%               integrated. In particular,
%               intrvals{i}(trm,2*j-1:2*j) = [p,q] indicates that variable
%               si_j is integrated over the domain [si_p,si_q], where
%               si_0 = ai and si_{di+1}=bi for [ai,bi] the interval on
%               which si exists;
% - params:     ntrms1 x ... x ntrmsN cell, specifying for each of the
%               integrals to take the kernel to use in this integral;
% - opts:       (optional) 'full' to not recombine the intervals obtained
%               from splitting a term, even if they involve the same
%               kernel. By default, such intervals are recombined;
%
% OUTPUTS
% - intmat:     1 x N cell, representing the indicator functions defining
%               the desired integrals;
% - params:     ntrms_new1 x ... x ntrms_newN cell, specifying the updated
%               kernels of the integrals, so that the output represents
%               the same integral as the input;

% DJ (CC), 10/15/2026: Initial coding;

intmat = intrvals;
int_order = intmat2int_order(intmat);
N = numel(intmat);

% full_idcs{k}: index range along direction k when slicing params in another
% direction's pass; updated at the end of each outer iteration, since the
% split grows direction i
sz_params = [size(params,1:N),1];
full_idcs = [cellfun(@(a)1:a,num2cell(sz_params),'UniformOutput',false)];

if nargin>=3 && strcmp(opts,'full')
    use_full_expand = true;
else
    use_full_expand = false;
end

for i=1:N
    % Get the intervals of integration along variable i
    intmat_i = intmat{i};
    [ntrms,di] = size(intmat_i);
    di = di/2;
    if di<=1
        continue
    end
    repdim = ones(1,N+1);
    % Loop over all terms and variables, in reverse order to which they
    % would be integrated
    [~,var_order] = sort(int_order{i},2);
    tmp_idcs = full_idcs;
    trm = 1;
    for cntr=1:ntrms
        % Get the parameters associated with the current term
        intmat_j = intmat_i(trm,:);
        tmp_idcs{i} = trm;
        params_j = params(tmp_idcs{:});
        param_fctrs = 1;
        % Get the order in which the variables are integrated in this term
        var_order_j = var_order(cntr,:);
        % Keep track of the ordering of the variables within the domain
        % (a <= si ... <= sj <= b)
        var_idcs = [0,di+1];
        % Get all possible orders of all variables satisfying the specified
        % integrals
        for j=1:di
            % Get the upper and lower domains of integration of var j
            var = var_order_j(end-j+1);            
            Llim = intmat_i(trm,2*var-1);
            Ulim = intmat_i(trm,2*var);
            if Llim==Ulim
                % Enforce the delta function, var k == Ulim
                intmat_j(intmat_j==var) = Ulim;
                intmat_i(trm,intmat_i(trm,:)==var) = Ulim;
                continue
            end
            var_idcs_j = zeros(0,size(var_idcs,2)+1);
            for l=1:size(var_idcs,1)
                % Determine the placement of Llim and Ulim in the order
                idxL = find(var_idcs(l,:)==Llim);
                idxU = find(var_idcs(l,:)==Ulim);
                ntrms_j = abs(idxU-idxL);
                if idxL>idxU
                    % Switch the limits of the integral
                    param_fctrs(l) = -param_fctrs(l);
                    idxL_new = idxU;
                    idxU = idxL;
                    idxL = idxL_new;
                end
                param_fctrs = [param_fctrs(1:l-1); param_fctrs(l)*ones(ntrms_j,1); param_fctrs(l+1:end)];
                % Consider every possible placement of var j between 
                % Llim and Ulim
                int_lims = [var_idcs(l,idxL:idxU-1)', var_idcs(l,idxL+1:idxU)'];
                intmat_l = repmat(intmat_j(l,:),ntrms_j,1);
                intmat_l(:,2*var-[1,0]) = int_lims;
                intmat_j = [intmat_j(1:l-1,:); intmat_l; intmat_j(l+1:end,:)];
                % Place var j in the associated positions between 
                % Llim and Ulim
                Llims = repmat([var_idcs(l,idxL:idxU),0],ntrms_j,1);
                Ulims = repmat([0,var_idcs(l,idxL:idxU)],ntrms_j,1);
                var_idcs_l = tril(Llims) +[zeros(ntrms_j,1),diag(var*ones(1,ntrms_j)),zeros(ntrms_j,1)] + triu(Ulims,2);
                var_idcs_l = [repmat(var_idcs(l,1:idxL-1),ntrms_j,1), var_idcs_l, repmat(var_idcs(l,idxU+1:end),ntrms_j,1)];
                var_idcs_j = [var_idcs_j; var_idcs_l];
            end
            var_idcs = var_idcs_j;
        end
        if ~use_full_expand
            % Recombine intervals if they involve the same parameter
            % (the same factor param_fctr)
            for j=1:di
                % Consider the variables in the order in which they must be
                % integrated
                var = var_order_j(j);
                if intmat_j(1,2*var)==intmat_j(1,2*var-1) || size(intmat_j,1)==1
                    % The variable is substituted for another
                    % --> ignore
                    continue
                end
                % Consider all terms involving the same parameters and the same
                % integrals in all vars but var j            
                intmat_var = intmat_j(:,2*var-[1,0]);
                intmat_novar = intmat_j(:,setdiff(1:2*di,2*var-[1,0]));
                [~,~,idcs2] = unique([param_fctrs,intmat_novar],'rows');
                % Find all integrals in varj for which the upper limit of one
                % matches the lower limit of another
                [Lidcs,Uidcs] = find_common_idcs(intmat_var(:,1),intmat_var(:,2));
                combine_idcs = idcs2(Lidcs)==idcs2(Uidcs);
                Lidcs = Lidcs(combine_idcs);
                Uidcs = Uidcs(combine_idcs);
                % Keep all integrals that cannot be combined
                no_combine_idcs = setdiff(1:size(intmat_j,1),union(Lidcs,Uidcs));
                intmat_nocombine = intmat_j(no_combine_idcs,:);
                param_fctr_nocombine = param_fctrs(no_combine_idcs);
                % Cobmine integrals in var j that can be combined
                intmat_combine = intmat_j(Lidcs,:);
                intmat_combine(:,2*var-1) = intmat_j(Uidcs,2*var-1);    % replace lower limit of interval
                param_fctr_combine = param_fctrs(Lidcs);
                % Update the intervals
                intmat_j = [intmat_nocombine; intmat_combine];
                param_fctrs = [param_fctr_nocombine; param_fctr_combine];
            end
        end
        % Add the intervals to the full list for variable i
        intmat_i = [intmat_i(1:trm-1,:); intmat_j; intmat_i(trm+1:end,:)];
        % Multiply the parameters with the appropriate factor
        ntrms_j = size(intmat_j,1);
        if ntrms_j>1
            repdim(i) = ntrms_j;
            params_j = repmat(params_j,repdim);
            if any(param_fctrs<0)
                tmp_idcs{i} = find(param_fctrs<0);
                params_j(tmp_idcs{:}) = cellfun(@(a) -a, params_j(tmp_idcs{:}),'UniformOutput',false);
            end
            % Add to the full list of parameters
            tmp_idcs1 = tmp_idcs;
            tmp_idcs2 = tmp_idcs;
            tmp_idcs1{i} = 1:trm-1;
            tmp_idcs2{i} = trm+1:sz_params(i);
            params = cat(i,params(tmp_idcs1{:}),params_j,params(tmp_idcs2{:}));
            sz_params(i) = sz_params(i) + ntrms_j-1;
        elseif param_fctrs<0
            % Multiply parameter with a factor of -1
            params(tmp_idcs{:}) = cellfun(@(a) -a, params_j,'UniformOutput',false);
        end
        % Move to the next integral specified by the user
        trm = trm+ntrms_j;
    end
    % Set the indicator functions representing the desired integral
    intmat{i} = intmat_i;
    full_idcs{i} = 1:size(intmat_i,1);
end

end