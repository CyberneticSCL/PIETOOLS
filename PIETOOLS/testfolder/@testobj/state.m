classdef (InferiorClasses={?polynomial,?dpvar})state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab class used to define symbolic objects to manipulate and define
% ode-pde equations in pietools. These objects will be used as ode
% state, pde state, output or input variables and can be vector-valued. To
% define a state object, use the syntax: x = state(type,length) where
% type: 'ode','pde','out','in' specifies the type of the variable
% length: specifies length of the variable (default=1)
% Once, x, is defined as a state variable, it can be used in various
% operations such as +,-,*,concatenation, diff(),int(),subs() to define
% dynamical systems
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications

    properties (SetAccess=protected)
        type = {'finite'; 'ode'}; % first cell element tells signal dimension 'finite' or 'infinite', second element tells signal type 'ode', 'pde', 'in', 'out'
                                  % third cell element may optionally be present in case of inputs and outputs qualifying additional properties such
                                  % as 'observed' or 'regulated' for outputs; 'disturbance' or 'control' for inputs
        Length = 1; % length of signal vector
        Vars = [pvar('t')]; % row vector of pvars corresponding to independent variables; First is always time t, then space is stored sequentially s1, s2, s3, etc.
        Diff = [0]; % row vector of non-negative integers showing the order of differentiation; i-th column corresponds to i-th variable in Vars
        Subs = [pvar('t')]; % row vector of pvars showing the substituted value; i-th column corresponds to i-th variable in Vars
        Int = [];   % Matrix of size (length(Vars)-1)x3; First value of row i states if there is an integral on i-th spatial variable, second value is kernel of lower integral, third value is kernel of upper integral
        Dom = []; % Matrix of size (length(Vars)-1)x2; domain of spatial variables, i-th row corresponds to domain of the spatial variable si in Vars.
    end
    properties (Hidden)
        MaxDiff = [inf]; % max allowable order of differentiation of each variable in Vars
    end
    properties (Hidden, SetAccess=protected)
        statename;
    end
    methods (Access = {?terms, ?sys, ?state})
        [out, varargout] = combine(varargin)
    end
    
    methods
        function obj = state(varargin) %constructor
            if nargout==0
                for i=1:nargin
                    obj = state();
                    obj.statename = stateNameGenerator();
                    assignin('caller', varargin{i}, obj);
                end
            else
                if nargin>=1
                    if strcmp(varargin{1},'in')
                        obj.type = {'finite'; 'in'; 'disturbance'};
                    elseif strcmp(varargin{1},'out')
                        obj.type = {'finite'; 'out'; 'regulated'};
                    elseif strcmp(varargin{1},'pde')
                        obj.type = {'infinite'; 'pde'};
                        obj.Vars = [pvar('t'),pvar('s1')];
                        obj.Diff = [0,0];
                        obj.Subs = obj.Vars;
                        obj.Int = [0,0,0];
                    end
                end
                if nargin>=2
                    
                elseif nargin==4 % internal use only, dont use this for constructing state vectors
                    obj.type = varargin{1};
                    obj.veclength = varargin{2};
                    obj.var = varargin{3};
                    obj.statename = varargin{4};
                    obj.diff_order = zeros(1,length(varargin{3}));
                elseif nargin>3
                    error('State class definition only takes 3 inputs');
                end
                if nargin>=3
                    obj.Length = varargin{3};
                end
                obj.statename = stateNameGenerator();
            end
        end
        
        function obj = set.dom(obj,dom)
            obj.dom = dom;
        end
        function obj = set.maxdiff(obj,maxdiff)
            obj.maxdiff = maxdiff;
        end

        % other class methods
        obj = subs(obj,var,var_val);
        obj = diff(obj,var,order);
        logval = eq(objA,objB);
        logval = isequal(objA,objB);
        obj = horzcat(varargin);
        obj = int(obj,var,limits);
        [logval,idx] = ismember(objA,objB)
        obj = minus(objA,objB);
        obj = mtimes(obj,K);
        logval = ne(objA,objB);
        obj = plus(objA,objB);
        obj = uplus(obj);
        obj = uminus(obj);
    end
end