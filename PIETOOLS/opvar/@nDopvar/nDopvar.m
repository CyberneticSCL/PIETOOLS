classdef (InferiorClasses={?polynomial,?dpvar})nDopvar
    properties
        vars_in = {};
        vars_out = {};
        dims = [1,1];
        params = dictionary([],[]);
        % alternatively, use two arrays, keysArray, paramsCell to store
        % params

        % for faster key operations, we can store key as uint16
    end

    methods
        function P = nDopvar(varsin,varsout,dims,params)
            P.vars_in = varsin(:).';
            P.vars_out = varsout(:).';
            P.dims = dims;
            P.params = params;
        end

        function vars = vars(P)
            vars = union(P.vars_in,P.vars_out);
        end
    end

    methods(Static)
        function idx = keys2index(key,n)
            key = key(:)-1;
            pow4 = 4.^(0:n-1);
            idx  = mod(floor(key ./ pow4), 4);    % (numel(key)-by-n)
        end

        function key = index2keys(idx,n)
            pow4 = 4.^(0:n-1);
            key = idx*pow4(:)+1;
            key = key.';
        end
    end
end