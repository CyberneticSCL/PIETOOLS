classdef (InferiorClasses={?polynomial,?dpvar}) opvar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This defines the class of operator variables P: R^p x L2^q to R^m x L2^n
% Date: 7/20/19
% Version: 1.0
% 
% CLASS properties
% P.dim: an 2x2 array with entries [m,p]
%                                  [n,q]
% P.P: a mxp matrix
% P.Q1: a mxq matrix valued polynomial in s
% P.Q2: a nxp matrix valued polynomial in s
% P.R.R0: a nxq matrix valued poynomial in s
% P.R.R1: a nxq matrix valued poynomial in s, theta
% P.R.R2: a nxq matrix valued poynomial in s, theta
% P.I: domain    
% P.var1, P.var2: polynomial variables s, theta (internal property, recommend not to modify)
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - opvar
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 7_26_2019
% Changed var1 and var2 to be pvars instead of chars. SS - 9/6/2019
% DJ - 09/28/23: Add conversion cell to opvar.



    properties
        P = [];
        Q1 = polynomial([]);
        Q2 = polynomial([]);
        R = struct('R0',polynomial([]),'R1',polynomial([]),'R2',polynomial([]));
        I(1,2) double = [0,1];
        var1 polynomial = pvar('s');
        var2 polynomial = pvar('theta');
        dim(2,2) double = [0 0; 0 0];
    end
    properties (Hidden)
        internal_call = 0;
    end
    properties (Dependent)
        dimdependent;
    end
    methods
        function [P] = opvar(varargin) %constructor
            if iscellstr(varargin)
            for i=1:nargin
                if ischar(varargin{i})
                    if nargout==0
                        assignin('caller', varargin{i}, opvar());
                    end
                else
                    error("Input must be strings");
                end
            end
            else
                P = opvar;
                A = varargin{1};
                if ~isa(A,"cell") || (~all(size(A)==[2,2]) && ~all(size(A)==[1,3]))
                    error('First argument to opvar must be 2x2 cell specifying parameters defining a 4-PI operator')
                end
                if all(size(A)==[1,3])
                    A = {[],[]; [],A};
                end
                for j=1:3
                    if ~isa(A{j},'double') && ~isa(A{j},'polynomial')
                        error('Parameters of opvar should be specified as arrays of type ''double'' or ''polynomial''')
                    end
                end
                % Set values of R->R, L2->R and R->L2 parameters
                if ~isempty(A{1,1})
                    P.P = A{1};
                end
                if ~isempty(A{2,1})
                    if P.dim(1,2)~=0 && size(A{2,1},2)~=P.dim(1,2)
                        error('Input dimensions of specified parameters do not match')
                    end
                    P.Q2 = polynomial(A{2,1});
                end
                if ~isempty(A{1,2})
                    if P.dim(1,1)~=0 && size(A{1,2},1)~=P.dim(1,1)
                        error('Output dimensions of specified parameters do not match')
                    end
                    P.Q1 = polynomial(A{1,2});
                end
                % Set values of 3-PI parameters
                if isa(A{2,2},'double') || isa(A{2,2},'polynomial')
                    if ~isempty(A{2,2})
                        if P.dim(2,2)~=0 && size(A{2,2},2)~=P.dim(2,2)
                            error('Input dimensions of specified parameters do not match')
                        elseif P.dim(2,1)~=0 && size(A{2,2},1)~=P.dim(2,1)
                            error('Output dimensions of specified parameters do not match')
                        end
                        P.R.R0 = polynomial(A{2,2});
                    end
                elseif isa(A{2,2},'cell') && numel(A{2,2})<=3
                    fset = {'R0','R1','R2'};
                    A22 = A{2,2};
                    for j=1:numel(A22)
                        if isempty(A22{j})
                            continue
                        end
                        if ~isa(A22{j},'double') && ~isa(A22{j},'polynomial')
                            error('Parameters of opvar should be specified as arrays of type ''double'' or ''polynomial''')
                        elseif P.dim(2,2)~=0 && size(A22{j},2)~=P.dim(2,2)
                            error('Input dimensions of specified parameters do not match')
                        elseif P.dim(2,1)~=0 && size(A22{j},1)~=P.dim(2,1)
                            error('Output dimensions of specified parameters do not match')
                        end
                        P.R.(fset{j})=A22{j};
                    end
                else
                    error('Parameters defining 3-PI subcomponent should be specified as 1x3 cell')
                end
                if nargin==1
                    return
                end
                % Set the domain of the operator
                if ~isa(varargin{2},'double') || ~all(size(varargin{2})==[1,2])
                    error('The domain of the PI operator should be specified as 1x2 array')
                elseif varargin{2}(1)>=varargin{2}(2)
                    error('The specified domain is invalid: make sure a<b in the interval [a,b]')
                else
                    P.I = varargin{2};
                end
                if nargin==2
                    return
                end
                % Set the first spatial variable.
                var1 = varargin{3};
                if ischar(varargin{3})
                    var1 = {var1};
                end
                if ~isa(var1,'polynomial') && ~iscellstr(var1)
                    error('Variables of the operator should be specified as ''polynomial'' or ''cellstr''')
                else
                    var1 = polynomial(var1);
                end
                if length(var1)==1
                    P.var1 = var1;
                elseif length(var1)==2 && nargin==3
                    P.var1 = var1(1);   P.var2 = var1(2);
                else
                    error('At most two spatial variables can be specified')
                end
                if nargin==3
                    return
                end
                % Set the second spatial variable.
                var2 = varargin{4};
                if ischar(varargin{4})
                    var2 = {var2};
                end
                if ~isa(var2,'polynomial') && ~iscellstr(var2)
                    error('Variables of the operator should be specified as ''polynomial'' or ''cellstr''')
                else
                    var2 = polynomial(var2);
                end
                if length(var2)>1
                    error('At most two spatial variables can be specified')
                end
                P.var2 = var2(1);
                if nargin>4
                    error('opvar accepts at most 4 arguments')
                end
            end
        end
        function [obj] = set.P(obj,P) 
            obj.P = P;
            if ~obj.internal_call
                obj = opvar_autofill(obj);
            end
        end
        function [obj] = set.Q1(obj,Q1)
            obj.Q1 = Q1;
            if ~obj.internal_call
                obj = opvar_autofill(obj);
            end
        end
        function [obj] = set.Q2(obj,Q2)
            obj.Q2 = Q2;
            if ~obj.internal_call
                obj = opvar_autofill(obj);
            end
        end
        function [obj] = set.R(obj,R)
            obj.R = R;
            if ~obj.internal_call
                obj = opvar_autofill(obj);
            end
        end
        function [d] = get.dimdependent(obj)
            % for consistent dimensions following vectors should have
            % all values in each of them to be equal or zero.
            M = [size(obj.P,1); size(obj.Q1,1)];
            P = [size(obj.P,2); size(obj.Q2,2)];
            N = [size(obj.Q2,1); size(obj.R.R0,1); size(obj.R.R1,1); size(obj.R.R2,1)];
            Q = [size(obj.Q1,2); size(obj.R.R0,2); size(obj.R.R1,2); size(obj.R.R2,2)];
            
            M = M(M~=0);
            if isempty(M)
                m=0;
            elseif all(M/max(M)==1)
                m = max(M);
            else
                m = nan;
            end
            
            P = P(P~=0); %#ok<*PROP>
            if isempty(P)
                p=0;
            elseif all(P/max(P)==1)
                p = max(P);
            else 
                p = nan;
            end
            
            Q = Q(Q~=0);
            if isempty(Q) 
                q=0;
            elseif all(Q/max(Q)==1)
                q = max(Q);
            else
                q = nan;
            end
            
            N = N(N~=0);
            if isempty(N) 
                n=0;
            elseif all(N/max(N)==1)
                n = max(N);
            else 
                n = nan;
            end
            
            d = [m p; n q];
        end
        function [val] = get.dim(obj)
            val = obj.dimdependent;
        end
        function [obj] = set.dim(obj,val)
            obj.dim = val;
            obj.internal_call=1;
            if isempty(obj.P)
                obj.P = zeros(val(1,:));
            end
            if isempty(obj.Q1)
                obj.Q1 = zeros(val(1,1),val(2,2));
            end
            if isempty(obj.Q2)
                obj.Q2 = zeros(val(2,1),val(1,2));
            end
            if isempty(obj.R.R0)
                obj.R.R0 = zeros(val(2,1),val(2,2));
            end
            if isempty(obj.R.R1)
                obj.R.R1 = zeros(val(2,1),val(2,2));
            end
            if isempty(obj.R.R2)
                obj.R.R2 = zeros(val(2,1),val(2,2));
            end
            obj.internal_call=0;
        end
        function obj = opvar_autofill(obj)
            obj.dim = obj.dim;
        end
    end
end