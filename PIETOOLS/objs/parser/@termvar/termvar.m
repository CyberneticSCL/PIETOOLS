classdef (InferiorClasses={?signals}) termvar
    properties
        operator; 
        state;
    end
    methods
        function obj = termvar(operator,state)
            obj.state = state;
            obj.operator = operator; 
        end
        function val = length(obj)
            val = size(obj.operator,1);
        end
    end
end