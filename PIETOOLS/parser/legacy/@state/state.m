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
%
% Initial coding SS - 02/10/2022
% MMP, 11/19/2026: Release comment pass (comment-only).

    properties
        type {mustBeMember(type,{'ode','pde','in','out'})} = 'ode';
        veclength {mustBeInteger,mustBePositive}=1;
        var {mustBeVector,mustBeA(var,["polynomial","double"])} = [pvar('t')];
        diff_order {mustBeInteger,mustBeVector,mustBeNonnegative}= [0];
    end
    properties (Hidden)
        maxdiff = "undefined";
        dom = [0,1];
    end
    properties (Hidden, SetAccess=protected)
        statename;
    end
    methods (Access = {?terms, ?sys, ?state})
        objterms = state2terms(obj,operator,var,val);
        [out, varargout] = combine(varargin)
    end
    
    methods
        function obj = state(varargin) %constructor
            % OBJ = STATE(TYPE,VECLENGTH,VAR) constructs a 'state' object;
            % without an output argument, STATE('x','y',...) assigns a new
            % 'ode' state of length 1 to each named variable in the caller
            % workspace.
            %
            % INPUT
            % - type:       (optional) char, 'ode' (default), 'pde', 'in'
            %               or 'out';
            % - veclength:  (optional) positive integer, the length
            %               (default 1);
            % - var:        (optional) row vector of 'polynomial'
            %               variables of the state (default t; t and s for
            %               'pde');
            % - statename:  (optional, internal) the ID of the state
            %               (default: a new ID).
            %
            % OUTPUT
            % - obj:        'state' object, all derivative orders 0.

            % Command syntax: one new state per named variable
            if nargout==0
                for i=1:nargin
                    obj = state();
                    obj.statename = stateNameGenerator();
                    assignin('caller', varargin{i}, obj);
                end
            else
                % Type, length, variables and ID from the inputs
                if nargin==1
                    obj.type = varargin{1};
                    if strcmp(varargin{1},'pde')
                        obj.var = [pvar('t'),pvar('s')];
                        obj.diff_order = [0,0];
                    end
                    obj.statename = stateNameGenerator();
                elseif nargin==2
                    % Type and length; a 'pde' state in t and s
                    obj.type = varargin{1};
                    obj.veclength = varargin{2};
                    if strcmp(obj.type,'pde')
                        obj.var = [pvar('t'),pvar('s')];
                        obj.diff_order = [0,0];
                    end
                    obj.statename = stateNameGenerator();
                elseif nargin==3
                    % Type, length and a row of variables
                    if size(varargin{3},1)~=1
                        error('var must be a row vector');
                    end
                    obj.type = varargin{1};
                    obj.veclength = varargin{2};
                    obj.var = varargin{3};
                    obj.statename = stateNameGenerator();
                    obj.diff_order = zeros(1,length(varargin{3}));
                elseif nargin==4 % internal use only, dont use this for constructing state vectors
                    obj.type = varargin{1};
                    obj.veclength = varargin{2};
                    obj.var = varargin{3};
                    obj.statename = varargin{4};
                    obj.diff_order = zeros(1,length(varargin{3}));
                elseif nargin>3
                    error('State class definition only takes 3 inputs');
                end
            end
        end

        % OBJ = SET.DOM(OBJ,DOM) and OBJ = SET.MAXDIFF(OBJ,MAXDIFF) store
        % the domain ('double' 1x2) and the maximal derivative orders of
        % the 'state' object OBJ, and return the object.
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