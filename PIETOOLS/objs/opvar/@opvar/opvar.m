classdef (InferiorClasses={?polynomial,?dpvar}) opvar
    properties
        dom;
        vars;
        N;
    end
    properties (Dependent)
        dim;
    end
    methods
    end
end