classdef sys
    properties
        equation;
        type char {mustBeMember(type,{'pde','dde','ddf','pie'})} = 'pde'
    end
    properties (Dependent)
        states;
        params;
    end
    methods
        obj = addequation(obj,eqn);
        obj = parse(obj);    
        obj = convert(obj,convertto);
        obj = removeequation(obj,eqnNumber);
    end
end
