classdef (InferiorClasses={?polynomial,?dpvar}) signals
    properties
        isInfinite; % N x 1 vector; 1 if the signal is infinitedimensional
        type; % N x 1 vector of strings
        len; % N x 1 vector; length of signal
        var; % N x 1 cell; independent variables in the signal
        dom; % N x 1 cell; cell elements themselves are arrays
        maxdiff; % N x 1 cell; cell elements themselves are arrays
        diffOrder; % N x 1 cell; cell elements themselves are arrays
        statename; % N x 1 vector; global id for easy identification of the signal
    end
    methods
        function obj = signals(varargin)
            if narargin==1 && isa(varargin{1},'signals')
                obj = varargin{1};
                return;
            end
            isInfinite = 0;
            type = ['ode'];
            len = [1];
            var = {[pvar('t')]};
            maxdiff = {[inf]};
            diffOrder = {[0]};
            dom = {[]};
            if nargin>0
                isInfinite = varargin{1};
            end
            if nargin>1
                type = varargin{2};
            end
            if nargin>2
                len = varargin{3};
            end
            if nargin>3
                var = varargin{4};
            else
                if isInfinite
                    var = {[pvar('t'),pvar('s1')]};
                end
            end
            if nargin>4
                dom = varargin{5};
            else
                if isInfinite
                    dom = {[0,1]};
                end
            end
            if nargin>5
                maxdiff = varargin{6};
            else
                if isInfinite
                    maxdiff = {[inf,inf]};
                end
            end
            if nargin>6
                diffOrder = varargin{7};
            else
                if isInfinite
                    diffOrder = {[0,0]};
                end
            end
            obj.isInfinite = isInfinite;
            obj.type = type;
            obj.var = var;
            obj.len = len;
            obj.dom = dom;
            obj.maxdiff = maxdiff;
            obj.diffOrder = diffOrder;
            obj.statename = stateNameGenerator();
        end

        function obj = setID(obj,val)
            obj.statename = val;
        end
        function obj = setInfinite(obj)
            obj.isInfinite = 1;
            obj.var = {[pvar('t'),pvar('s1')]};
            obj.maxdiff = {[inf,inf]};
            obj.diffOrder = {[0,0]};
            obj.dom = {[0,1]};
        end
        function obj = setFinite(obj)
            obj.isInfinite = 0;
            obj.var = {[pvar('t')]};
            obj.maxdiff = {[inf]};
            obj.diffOrder = {[0]};
            obj.dom = {[]};
        end


        obj = combine(A,B);
        obj = eq(A,B);
        obj = isequal(A,B);
        obj = ismember(A,B);
        obj = plus(A,B);
        obj = subsref(A,B);
        % almost trivial methods below
        obj = horzcat(A,B);
        obj = uplus(A,B);
        obj = uminus(A,B);
        obj = minus(A,B);
        obj = ne(A,B);

        % mathematical operations
        obj = diff(A,var,order);
        obj = int(A,var,lim);
    end
end