classdef (InferiorClasses={?polynomial,?dpvar,?opvar}) dopvar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This defines the class of operator variables P: R^p x L2^q to R^m x L2^n
% Date: 7/3/21
% Version: 1.0
% 
% CLASS properties
% P.dim: an 2x2 array with entries [m,p]
%                                  [n,q]
% P.P: a mxp matrix or dpvar
% P.Q1: a mxq matrix valued polynomial in s or dpvar
% P.Q2: a nxp matrix valued polynomial in s or dpvar
% P.R.R0: a nxq matrix valued poynomial in s or dpvar
% P.R.R1: a nxq matrix valued poynomial in s, theta or dpvar
% P.R.R2: a nxq matrix valued poynomial in s, theta or dpvar
% P.I: domain    
% P.var1, P.var2: polynomial variables s, theta (internal property, recommend not to modify)
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - dopvar
%
% Copyright (C)2021  M. Peet, S. Shivakumar
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




    properties
        P = dpvar();
        Q1 = dpvar();
        Q2 = dpvar();
        R = struct('R0',dpvar(),'R1',dpvar(),'R2',dpvar()); 
        I(1,2) double = [0,1];
        var1 polynomial = pvar('s');
        var2 polynomial = pvar('theta');
        dim(2,2) double = [0 0; 0 0];
    end
    properties (Dependent)
        dimdependent;
    end
    methods
        function obj = dopvar(varargin) %constructor            
            for i=1:nargin
                if ischar(varargin{i})
                    if nargout==0
                        assignin('caller', varargin{i}, dopvar());
                    end
                else
                    error("Input must be strings");
                end
            end
        end
        function [obj] = set.P(obj,P) 
            obj.P = P;
        end
        function [obj] = set.Q1(obj,Q1)
            obj.Q1 = Q1;
        end
        function [obj] = set.Q2(obj,Q2)
            obj.Q2 = Q2;
        end
        function [obj] = set.R(obj,R)
            obj.R = R;
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
            if isempty(obj.P)
                obj.P = dpvar(zeros(val(1,:)));
            end
            if isempty(obj.Q1)
                obj.Q1 = dpvar(zeros(val(1,1),val(2,2)));
            end
            if isempty(obj.Q2)
                obj.Q2 = dpvar(zeros(val(2,1),val(1,2)));
            end
            if isempty(obj.R.R0)
                obj.R.R0 = dpvar(zeros(val(2,1),val(2,2)));
            end
            if isempty(obj.R.R1)
                obj.R.R1 = dpvar(zeros(val(2,1),val(2,2)));
            end
            if isempty(obj.R.R2)
                obj.R.R2 = dpvar(zeros(val(2,1),val(2,2)));
            end
        end
    end
end