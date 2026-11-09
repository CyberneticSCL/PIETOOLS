function intmat = remove_badints(intmat)
% INTMAT = REMOVE_BADINTS(INTMAT) removes from the intmat cell of an
% 'intop' structure, specifying domains of integration, all terms
% corresponding to "bad" integrals, e.g. an integral in which s1 is
% integrated over [a,s2] and s2 over [a,s1].
%
% INPUTS
% - intmat: 1 x N cell of qi x 2*di arrays specifying a set of intervals
%           over which to integrate the di degrees of variable i in each of
%           the qi different terms, as per the 'intop' format;
%
% OUTPUTS
% - intmat: 1 x N cell specifying the same domains of integration, but
%           with inconsistent integrals removed;
%
% NOTES
% - A term is valid iff it admits an integration order (see is_consistent()
%   in isvalid.m and get_int_order.m).

% DJ (CC), 10/15/2026: Initial coding;

for i=1:numel(intmat)
    % Extract the domains of integration for the ith variable
    intmat_i = intmat{i};

    % Loop over all specified integrals
    [nterms,di] = size(intmat_i);
    di = di/2;
    rtn_trms = 1:nterms;
    for j=1:nterms
        % First check if any of the variables is integrated over the full
        % domain
        is_a = intmat_i(j,1:2:end)==0;
        is_b = intmat_i(j,2:2:end)==di+1;
        if ~any(is_a & is_b)
            rtn_trms = setdiff(rtn_trms,j);
            continue
        end
        % Get rid of variables which are evaluated at another variable
        intmat_j = intmat_i(j,:);
        mult_vars = find(intmat_j(1:2:end)==intmat_j(2:2:end));
        if any(mult_vars==intmat_j(2*mult_vars))
            % Evaluating vark==vark makes no sense
            rtn_trms = setdiff(rtn_trms,j);
            continue
        end
        % Reject mutual delta cycles (vari fixed at varj, varj fixed at vari):
        % the elimination loop below only checks for self-reference up front
        is_cyclic = false;
        for v0 = reshape(mult_vars,1,[])
            visited = v0;
            w = intmat_j(2*v0);
            while ismember(w,mult_vars)
                if ismember(w,visited)
                    is_cyclic = true;
                    break
                end
                visited = [visited,w];
                w = intmat_j(2*w);
            end
            if is_cyclic
                break
            end
        end
        if is_cyclic
            rtn_trms = setdiff(rtn_trms,j);
            continue
        end
        for k=1:numel(mult_vars)
            % Enforce vark = varU
            var = mult_vars(k);
            intmat_j(intmat_j==var) = intmat_j(2*var);
            mult_vars(mult_vars==var) = intmat_j(2*var);
            % Remove vark
            intmat_j(intmat_j>var) = intmat_j(intmat_j>var)-1;
            mult_vars(k+1:end) = mult_vars(k+1:end)-1;
            intmat_j(2*var-[1,0]) = [];            
        end
        d = numel(intmat_j)/2;
        if d==0
            rtn_trms = setdiff(rtn_trms,j);
            continue
        end

        % Define a 1 x d cell specifying for each of the d variables what other
        % variables are used to define the associated domain of integration
        var_list0 = mat2cell(intmat_j,1,2*ones(1,d));
        % Keep track of which variables we know are integrated over all space
        safe_vars = [0,d+1];
        % Keep track of remaining variables to check
        rem_vars = 1:d;
        % Check if the domain of integration of any variable involves the
        % variable itself
        var_list = var_list0;
        unsafe_idcs = cellfun(@(a,b) ismember(a,b),num2cell(rem_vars),var_list);
        if any(unsafe_idcs)
            % One of the domains of integration is inappropriate
            rtn_trms = setdiff(rtn_trms,j);
            continue
        end
        for k=1:d        
            % Remove variables which are integrated over the full domain
            var_list = cellfun(@(a) setdiff(a,safe_vars),var_list,'UniformOutput',false);
            % Update the list of safe variables
            safe_idcs = cellfun(@(a) isempty(a),var_list);
            if all(safe_idcs)
                % All variables are good
                break
            end
            safe_vars = [safe_vars,find(safe_idcs)];
            % For remaining variables, determine on what interval they are defined
            rem_vars = rem_vars(~safe_idcs);
            var_list = var_list(~safe_idcs);
            var_list = cellfun(@(a) cell2mat(var_list0(a)),var_list,'UniformOutput',false);
            % Check if any of the variables maps to itself
            unsafe_idcs = cellfun(@(a,b) ismember(a,b),num2cell(rem_vars),var_list);
            if any(unsafe_idcs)
                % The domain of integration of one variable involves the
                % variable itself
                rtn_trms = setdiff(rtn_trms,j);
                break
            end
        end
        if ~all(safe_idcs)
            rtn_trms = setdiff(rtn_trms,j);
            continue
        end
    end
    intmat{i} = intmat_i(rtn_trms,:);
end

end