classdef (InferiorClasses={?polynomial,?dpvar})nDopvar
    properties
        vars_in = [];
        vars_out = [];
        dims = [1,1];
        params = dictionary([],[]);
    end

    methods
        function P = nDopvar(varargin)
            P.vars_in = varargin{1};
            P.vars_out = varargin{2};
            P.dims = varargin{3};
            P.params = varargin{4};
        end

        function vars = vars(P)
            vars.vars_in = P.vars_in;
            vars.vars_out = P.vars_out;
        end
    end
end