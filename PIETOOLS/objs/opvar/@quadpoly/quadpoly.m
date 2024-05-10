classdef (InferiorClasses={?polynomial})quadpoly
    % data structure storing a polynomial in a quadratic form given by
    % f(s,theta,d) = (I x Z_d(d))^T (I x Z_s(s))^T C (I X Z_th(theta))
    % where d is decision variable vector, s is left variable vector, theta
    % is right dummy variable vector corresponding to s. Z_d has max
    % degree one.
    properties (Access=public)
        lvarname = {}; % variable names of left monomial vector
        rvarname = {}; % variable names of right monomial vector
        dvarname = {}; % decision variable names (if a decision polynomial)
        ldegmat = sparse([]);
        rdegmat = sparse([]);
        C = sparse([]);
        matdim = [0,0];
        isdecVar = 0;  % using a single object for decision variables
    end

    methods
        function obj = quadpoly(varargin)
            if nargin == 0
                % Return quadpoly object with default values
                return
            elseif iscellstr(varargin) % if input is a cell array of strings create decision quadpoly with those names
                if nargout==0
                    for i=1:nargin
                        Cf = sparse([0;1]);
                        dmat = zeros(1,0);
                        vname = {};
                        dvname = {varargin{i}};
                        matdim = [1,1];
                        isdVar = 1;
                        assignin('caller', varargin{i}, quadpoly(Cf,dmat,dmat,vname,vname,dvname,matdim,isdVar));
                        clear obj;
                    end
                elseif nargin==1 && nargout==1
                    Cf = sparse([0;1]);
                    dmat = zeros(1,0);
                    vname = {};
                    dvname = {varargin{1}};
                    matdim = [1,1];
                    obj = quadpoly(Cf,dmat,dmat,vname,vname,dvname,matdim,isdVar);
                else
                    error('For cell char inputs, number of outputs should be at most 1');
                end
            elseif nargin==1 % single input is either polynomial, quadpoly or double
                if isa(varargin{1},'quadpoly')
                    % If input is already a quadpoly
                    obj = varargin{1};
                elseif isa(varargin{1},'polynomial')
                    % If input is a poly, convert to dpvar
                    obj = poly2quadpoly(varargin{1});
                elseif isa(varargin{1},'double') && isreal(varargin{1}) && ismatrix(varargin{1})
                    % If input is a real double, then convert to a poly
                    szc = size(varargin{1});
                    obj.C = varargin{1};
                    obj.ldegmat = sparse(1,0);
                    obj.rdegmat = sparse(1,0);
                    obj.lvarname  = {};
                    obj.rvarname  = {};
                    obj.dvarname = {};
                    obj.matdim = szc;
                else
                    error(['For single input, argument must be a quadpoly, '...
                        'real double, or a polynomial']);
                end
            elseif nargin<6
                errstr1 = 'Invalid number of inputs for the "quadpoly" command.';
                error([errstr1]);
            elseif nargin<=8 % more than one input means fields of the dpvar object
                if nargin==6   % adding an option to skip the matdim for scalars
                    varargin{7} = [1 1];
                    varargin{8} = 0;
                elseif nargin==7
                    varargin{8} = 0;
                end
                % convert row cells to columns cells and ensure inputs are
                % row/column vectors
                if size(varargin{4},1)==1
                    varargin{4} = varargin{4}';
                end
                if size(varargin{5},1)==1
                    varargin{5} = varargin{5}';
                end
                if size(varargin{6},1)==1 % convert row cells to columns cells
                    varargin{6} = varargin{6}';
                end

                obj.lvarname = varargin{4};
                obj.rvarname = varargin{5};
                obj.dvarname = varargin{6};
                obj.C = varargin{1};
                obj.ldegmat = varargin{2};
                obj.rdegmat = varargin{3};
                obj.matdim = varargin{7}; % why (:)' -DJ????
                obj.isdecVar = varargin{8}; % this is optional which specifies if error check needs to be skipped
            else
                error('Too many inputs for "quadpoly" class object');
            end
        end
    end
end