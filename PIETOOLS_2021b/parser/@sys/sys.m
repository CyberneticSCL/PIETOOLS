classdef sys
    properties
        equation {validateEquation(equation)} = {};
        type char {mustBeMember(type,{'pde','dde','ddf','pie'})} = 'pde';
        params {mustBeA(params,'pdeparams')} = pdeparams();
        ControlledInputs;
        ObservedOutputs;
    end
    properties (Dependent)
        states;
    end
    methods
        function prop = get.states(obj)
            prop = getStatesFromEquations(obj);
        end
        function prop = get.params(obj)
            if isempty(obj.params)
                prop = getParams(obj);
            else
                prop = obj.params;
            end
        end
        function obj = set.params(obj,params)
            obj.params = params;
        end
        function out = get.ObservedOutputs(obj)
            if isempty(obj.ObservedOutputs)
                statelist = getStatesFromEquations(obj);
                if ~isempty(statelist)
                    out = zeros(length(statelist.veclength),1);
                else
                    out = [];
                end
            else
                out = obj.ObservedOutputs;
            end
        end
        function obj = set.ObservedOutputs(obj,OO)
            obj.ObservedOutputs = OO;
        end
        function out = get.ControlledInputs(obj)
            if isempty(obj.ControlledInputs)
                statelist = getStatesFromEquations(obj);
                if ~isempty(statelist)
                    out = zeros(length(statelist.veclength),1);
                else
                    out = [];
                end
            else
               out = obj.ControlledInputs;
            end
        end
        function obj = set.ControlledInputs(obj,CI)
            obj.ControlledInputs = CI;
        end
        
        prop = getStatesFromEquations(obj);
        obj = addequation(obj,eqn);
        obj = getParams(obj);    
        obj = convert(obj,convertto);
        obj = removeequation(obj,eqnNumber);
        obj = setControl(obj,input);
        obj = removeControl(obj,input);
        obj = setObserve(obj,output);
        obj = removeObserve(obj,output);
        
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