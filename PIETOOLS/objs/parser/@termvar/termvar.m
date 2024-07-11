classdef (InferiorClasses={?signals}) termvar < opvar & signals
    properties
        % this does not require any properties; it is just an interface
        % between state objects and opvar operators
    end
    methods
        function obj = termvar(operator,state)
            obj@signals(state);
            obj@opvar(operator);
        end
    end
end