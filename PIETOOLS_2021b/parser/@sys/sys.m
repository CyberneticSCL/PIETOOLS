classdef sys
    properties
        equation {validateEquation(equation)} = {};
        type char {mustBeMember(type,{'pde','dde','ddf','pie'})} = 'pde';
        params;
    end
    properties (Dependent)
        states;
    end
    methods
        function prop = get.states(obj)
            prop = getStatesFromEquations(obj);
        end
        function prop = get.params(obj)
            prop = getParams(obj);
        end
        
        obj = addequation(obj,eqn);
        obj = getParams(obj);    
        obj = convert(obj,convertto);
        obj = removeequation(obj,eqnNumber);
    end
end
function validateEquation(prop)
if ~iscell(prop)
    error('Equations must be stored in cell column array');
end
eqntype = cellfun(@(x) ~isa(x,'terms'), prop,'un',0);
eqntype = cell2mat(eqntype);
if ~isempty(eqntype)&&any(eqntype(:))
    error("Equation entries should be terms type object.");
end
end