classdef (InferiorClasses={?polynomial,?dpvar}) signals
    properties
        type; % N x 1 vector of strings
        len; % N x 1 vector; length of signal
        var; % N x 1 cell; independent variables in the signal
        dom; % N x 1 cell; cell elements themselves are arrays
        maxdiff; % N x 1 cell; cell elements themselves are arrays
        diffOrder; % N x 1 cell; cell elements themselves are arrays
    end
    properties (Hidden, Dependent)
        isInfinite; % N x 1 vector; 1 if the signal is infinitedimensional
        dim; % N x 1 vector; spatial dimension of signal
    end
    properties (Hidden, SetAccess=private)
        statename; % N x 1 vector; global id for easy identification of the signal
    end

    methods (Access=private)
        function obj = setID(obj,val)
            obj.statename = val;
        end
    end

    methods
        function obj = signals(varargin)
            % no input constructor
            if nargin == 0
                obj = odevar();
                return
            end

            % if input is a signals object do nothing
            if nargin==1 && isa(varargin{1},'signals')
                obj = varargin{1};
                return;
            end

            % set default values
            if nargin>0
                type = varargin{1};
                len = [1];
                var = [pvar('t')];
                maxdiff = [inf];
                diffOrder = [0];
                dom = [];
            end
            if nargin>1
                len = varargin{2};
                if strcmp(type,"pde")
                    var = [pvar('t'),pvar('s1')];
                    dom = [0,1];
                else
                    var = [pvar('t')];
                    dom = [];
                end
                maxdiff = inf(length(var),1);
                diffOrder = zeros(length(var),1);
            end
            if nargin>2
                var = varargin{3};
                dom = repmat({[0,1]},length(var)-1,1);
                maxdiff = inf(length(var),1);
                diffOrder = zeros(length(var),1);
            else
                if strcmp(type,"pde")
                    var = [pvar('t'),pvar('s1')];
                end
            end
            if nargin>3
                dom = varargin{4};
            else
                if strcmp(type,"pde")
                    dom = [0,1];
                end
            end
            if nargin>4
                maxdiff = varargin{5};
            else
                if strcmp(type,"pde")
                    maxdiff = [inf,inf];
                end
            end
            if nargin>5
                diffOrder = varargin{6};
            else
                if strcmp(type,"pde")
                    diffOrder = [0,0];
                end
            end
            obj.type = type;
            obj.var = var;
            obj.len = len;
            obj.dom = dom;
            obj.maxdiff = maxdiff;
            obj.diffOrder = diffOrder;
            obj.statename = stateNameGenerator();
        end

        function val = get.isInfinite(obj)
            val = 0;
            if strcmp(obj.type,"pde") || (length(obj.var)>=2)
                val = 1;
            end
        end

        function val = get.dim(obj)
            val = length(obj.var);
            val = val-1;
        end

        function obj = setInfinite(obj,N)
            if any(strcmp(obj.type,"ode"))
                tmpMsg = 'Input to the method is "ode" type. Cannot be set as infinite';
                error(tmpMsg);
            end
            if nargin==1
                N = 1;
            end
            warning("Setting an object as ND resets the properties to default. Check metadata related to differentiability and domain.");
            obj.var = repmat(pvar('t'),1,N+1);
            obj.maxdiff = inf(1,N+1);
            obj.diffOrder = zeros(1,N+1);
            obj.dom = repmat([0,1],N,1);
            for i=2:N+1
                obj.var(i) = pvar(['s',num2str(i-1)]);
            end
        end

        function obj = setFinite(obj)
            if any(strcmp(obj.type,"pde"))
                tmpMsg = 'Input to the method is "pde" type. Cannot be set as finite';
                error(tmpMsg);
            end
            warning("Setting an object as 0D resets the properties to default. Check metadata related to differentiability and domain.");
            obj.var = [pvar('t')];
            obj.maxdiff = [inf];
            obj.diffOrder = [0];
            obj.dom = [];
        end
    end
end