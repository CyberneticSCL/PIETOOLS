classdef (InferiorClasses={?signals,?termvar}) sys
    properties
        equation {validateEquation(equation)} = [];
        type char {mustBeMember(type,{'pde','dde','ddf','pie'})} = 'pde';
        params {mustBeA(params,{'pde_struct','pie_struct'})} = pde_struct();
        ControlledInputs = [];
        ObservedOutputs = [];
    end
    properties (Dependent)
        states;
        dom;
    end

    methods
        function obj = sys(varargin)
            if nargin==0
                type = 'pde';
            else
                type = varargin{1};
            end
            obj.equation = [];
            obj.params = pde_struct();
            obj.ControlledInputs = [];
            obj.ObservedOutputs = [];
            obj.type = type;
            tmpMsg = ['Initialized sys() object of type ' obj.type];
            disp(tmpMsg);


            if nargin==2
                if strcmp(obj.type,'pde')&&isa(varargin{2},'termvar')
                    self.addequation(varargin{2});
                elseif isa(varargin{2},'pde_struct')||isa(varargin{2},'pie_struct')
                    obj.params = varargin{2};
                else
                    tmpMsg = 'Unknown second input to the sys() type object.';
                    error(tmpMsg);
                end
            end
        end
        function prop = get.states(obj)
            prop = getStatesFromEquations(obj);
        end
        function prop = get.params(obj)
            if isempty(obj.params.dom)
                prop = getParams(obj);
            else
                prop = obj.params;
            end
        end
        function obj = set.params(obj,params)
            obj.params = params;
        end
        function obj = set.dom(obj,dom)
            obj.dom = dom;
        end
        function out = get.ObservedOutputs(obj)
            if isempty(obj.ObservedOutputs)
                statelist = obj.states;
                if ~isempty(statelist)
                    out = zeros(length(statelist.len),1);
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
                statelist = obj.states;
                if ~isempty(statelist)
                    out = zeros(length(statelist.len),1);
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
        disp(obj);
    end
end
function validateEquation(prop)
eqntype = (~isa(prop,'termvar'))&& (~isempty(prop));
if ~isempty(eqntype)&&any(eqntype(:))
    tmpMsg ="Equation entries should be terms type object.";
    error(tmpMsg);
end
end