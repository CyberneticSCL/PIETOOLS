classdef (InferiorClasses={?signals,?termvar}) sys
    properties
        equation {validateEquation(equation)} = [];
        params {mustBeA(params,{'pde_struct','tds_struct','pie_struct'})} = pde_struct(); 
    end
    properties
        CInames; % statename list of controlled inputs
        OOnames; % statename list of observed outputs
    end
    properties (Dependent)
        type;
        states;
        dom;
        ControlledInputs;  % this property is redundant, may need to be removed in later releases
        ObservedOutputs;  % this property is redundant, may need to be removed in later releases
    end
    properties (Hidden, SetAccess=private)
        tampered=0;
    end

    methods
        function obj = sys(varargin)
            if nargin==0
                type = 'pde';
            else
                type = varargin{1};
            end
            obj.equation = [];
            switch type
                case 'pde'
                    obj.params = pde_struct();
                case 'pie'
                    obj.params = pie_struct();
                case 'tds'
                    obj.params = tds_struct();
                otherwise
                    error('Unknown system type.');
            end
            tmpMsg = ['Initialized sys() object of type ' obj.type];
            disp(tmpMsg);


            if nargin==2
                if strcmp(obj.type,'pde')&&isa(varargin{2},'termvar')
                    obj = addequation(obj,varargin{2});
                elseif isa(varargin{2},'pde_struct')||isa(varargin{2},'pie_struct')
                    obj.tampered = 1;
                    warning('Parameters of "sys" object is being directly modified. Only "convert()" method is supported in this scenario. Data such as equations, inputs, outputs, etc., cannot be modified.');
                    obj.params = varargin{2};
                else
                    tmpMsg = 'Unknown second input to the sys() type object.';
                    error(tmpMsg);
                end
            end
        end
        function prop = get.type(obj)
            if isempty(obj.params)||isa(obj.params,'pde_struct')
                prop = 'pde';
            elseif isa(obj.params,'pie_struct')
                prop = 'pie';
            else
                prop = 'tds';
            end
        end
        function prop = get.states(obj)
            if obj.tampered
                error('This property/method is unavailable: "sys()" object parameters were manually defined/edited.');
            end
            prop = getStatesFromEquations(obj);
        end
        function prop = get.dom(obj)
            if obj.tampered
                error('This property/method is unavailable: "sys()" object parameters were manually defined/edited.');
            end
            prop = getdomain(obj.params);
        end
        function out = get.ObservedOutputs(obj)
            if obj.tampered
                error('This property/method is unavailable: "sys()" object parameters were manually defined/edited.');
            end
            statelist = obj.states;
            if ~isempty(statelist)
                out = zeros(length(statelist.len),1);
            else
                out = [];
            end
            out(ismember(out,obj.OOnames)) = 1;
        end
        function out = get.ControlledInputs(obj)
            if obj.tampered
                error('This property/method is unavailable: "sys()" object parameters were manually defined/edited.');
            end
            statelist = obj.states;
            if ~isempty(statelist)
                out = zeros(length(statelist.len),1);
            else
                out = [];
            end
            out(ismember(out,obj.CInames)) = 1;
        end
        function prop = get.params(obj)
            if ~isempty(obj.equation)
                prop = getParams(obj);
            else
                prop = obj.params;
            end
        end

        function obj = set.params(obj,params)
            obj.tampered = 1;
            warning('Parameters of "sys" object is being directly modified. Only "convert()" method is supported in this scenario. Data such as equations, inputs, outputs, etc., cannot be modified.');
            obj.params = params;
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