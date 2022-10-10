classdef sys
    properties
        equation {validateEquation(equation)} = [];
        type char {mustBeMember(type,{'pde','dde','ddf','pie'})} = 'pde';
        params {mustBeA(params,{'pde_struct','pie_struct'})} = pde_struct();
        ControlledInputs;
        ObservedOutputs;
    end
    properties (Dependent)
        states;
    end
    methods
        function obj = sys(type)
            if nargin==0
                type = 'pde';
            end
            obj.equation = [];
            obj.params = pde_struct();
            obj.ControlledInputs = [];
            obj.ObservedOutputs = [];
            obj.type = type;
            fprintf('Initialized sys() object of type "%s"\n',obj.type);
        end
        function prop = get.states(obj)
            prop = getStatesFromEquations(obj);
        end
        function prop = get.params(obj)
            if isempty(obj.params.x)
                obj = getParams(obj);
                prop = obj.params;
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
% if ~iscell(prop)
%     error('Equations must be stored in cell column array');
% end
eqntype = (~isa(prop,'terms'))&& (~isempty(prop));
% eqntype = cell2mat(eqntype);
if ~isempty(eqntype)&&any(eqntype(:))
    error("Equation entries should be terms type object.");
end
end