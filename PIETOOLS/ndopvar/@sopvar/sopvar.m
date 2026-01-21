classdef (InferiorClasses={?polynomial,?dpvar})sopvar
    % Represents PI maps from L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
    % note [Si,...,Sn] may not have all indices from i to n, likewise
    % for j to m.
    properties
        vars_in = {};   % cell string  {'Si',...,'Sn'}
        vars_out = {};  % cell strings {'Sj',...,'Sm'}
        dims = [1,1];  % vector-space dimension [q,p]
        params = {quadPoly(0,[],[],[1,1],{},{})};
        % each value in the params cell is a quadPoly
        % for maps to and from \R spaces, add an extra dimension, s0?
    end

    methods
        function P = sopvar(varsin,varsout,dims,params)
            P.vars_in = varsin(:).';
            P.vars_out = varsout(:).';
            P.dims = dims;

            if nargin==3
                vars = union(P.vars_out, P.vars_in);
                n= numel(vars);
                celldim = repmat(3,1,n);
                [~,idx] = ismember(vars, intersect(P.vars_in,P.vars_out));
                celldim(~idx) = 1;
                if n == 0
                    P.params = repmat({zeros(dims(1),dims(2))},1,1);
                else
                    P.params = repmat({zeros(dims(1),dims(2))},celldim);
                end
            else
                P.params = params;
            end
        end
    end

    methods(Static)
        out = randOpvar(varin,varout,dim,deg,den);
    end
end