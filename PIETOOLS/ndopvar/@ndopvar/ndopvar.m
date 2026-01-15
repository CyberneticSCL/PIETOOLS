classdef (InferiorClasses={?polynomial,?dpvar,?nopvar}) ndopvar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This defines the class of PI operator variables,
%
%   Pop(r): L2^n[[a1,b1]x...x[aN,bN]] --> L2^m[[a1,b1]x...x[aN,bN]],
%
% taking the form
%
% (Pop(r) x)(s) = sum_{j in {0,1,2}^N} int_{[a,b]} I_j(s,t)*R(s,t;r)*x(t) dt 
%
% for s=(s1,...,sN) in [a,b]=[a1,b1]x...x[aN,bN], where
%
%   I_j(s,t) = I_j1(s1,t1)*...*I_jN(sN,tN)
%
% for I_k(si,thi) indicator functions such that
%
%                                       { f(si),                    k=0;
%   int_{ai}^{bi} I_k(si,ti)*f(ti) dti ={ int_{ai}^{si} f(ti) dti,  k=1;
%                                       { int_{si}^{bi} f(ti) dti,  k=2;
%
% and parameterized by polynomial variables Rj(s,t;r), defined by decision
% variables r. We represent the parameters in a quadratic form as
%
%   Rj(s,t;r) = (Im o Zd(s))^T (Ik o [1;d])^T Cj (In o Zd(t)),
%
% where o denotes the Kronecker product, and where
%
%   Zd(s) = Zd1(s1) o ... o ZdN(sN),
%
% for Zdi(si) the vector of monomials of degree at most di in variable si.
%
%
% CLASS properties
% - P.C:    3^N cell, with element {i1,i2,...,iN} a sparse matrix of
%           dimension m*prod(deg+1)*(q+1) x n*prod(deg+1), representing
%           the coefficient matrix Cj for (j1,...,jN)=(i1+1,...,iN+1);
% - P.deg:  Nx1 array specifying the maximal monomial degrees d1,...,dN;
% - P.dim:  1x2 array specifying the dimensions [m,n] of the operator;
% - P.dom:  Nx2 array with each row dom(i,:) = [ai,bi] representing the
%           spatial interval along the ith direction on which the operator
%           is defined;
% - P.vars: kx2 'polynomial' array, specifying the primary variables s
%           (first column) and dummy variables t (second column) defining
%           the parameters R(s,t), as pvar objects.
% - P.dvarnames: qx1 cell of strings, specifying the names of the decision
%                variables, r, parameterizing the operator;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ndopvar
%
% Copyright (C) 2026 PIETOOLS Team
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
% DJ, 01/06/2206: Initial coding

properties
    C = {};
    deg = zeros(0,1);
    dom = zeros(0,2);
    dim = [0,0];
    dvarname = {};
    vars = polynomial(zeros(0,2));
end

methods
    function [P] = ndopvar(varargin) %constructor
        if nargout==0
            % Declare operators as
            %   ndopvar P1 P2 P3
            for i=1:nargin
                if ischar(varargin{i})
                    if nargout==0
                        assignin('caller', varargin{i}, ndopvar());
                    end
                else
                    error("Input must be strings");
                end
            end
        end
    end
    function [obj] = set.C(obj,C) 
        obj.C = C;
    end
    function [obj] = set.deg(obj,deg)
        obj.deg = deg;
    end
    function [obj] = set.dom(obj,dom)
        obj.dom = dom;
    end
    function [obj] = set.vars(obj,vars)
        obj.vars = vars;
    end
    function [obj] = set.dvarname(obj,dvarname)
        obj.dvarname = dvarname;
    end
    function [dim] = get.dim(obj)
        % % Determine the dimensions of the operator, m x n, from the
        % % dimensions of the coefficient matrices,
        % %     m*prod(deg+1)*(q+1) x n*prod(deg+1).

        % Get the number of monomials and decision variables
        degs = obj.deg;
        if isempty(degs)
            degs = 0;
        end
        nZ = prod(obj.deg+1);
        N = numel(obj.deg);
        q = numel(obj.dvarname);

        % Check the dimensions of the coefficient matrices
        m_min = inf;    n_min = inf;
        m_max = 0;      n_max = 0;
        sz_C = [size(obj.C),1];
        % if sz_C(2)==1
        %     sz_C = sz_C(1);
        % end
        for ii=1:numel(obj.C)
            [m,n] = size(obj.C{ii});
            m_min = min(m_min,m/(nZ*(q+1)));
            m_max = max(m_max,m/(nZ*(q+1)));
            % Get the number of monomials and decision variables
            idcs = cell(1,N);
            [idcs{:}] = ind2sub(sz_C,ii);
            idcs = cell2mat(idcs);
            nZ_t = prod(degs(logical(idcs-1))+1);
            n_min = min(n_min,n/(nZ_t));
            n_max = max(n_max,n/(nZ_t));
        end
        
        % Set the dimensions
        dim = [nan,nan];
        if m_min==m_max && round(m_min)==m_min
            dim(1) = m_min;
        end
        if n_min==n_max && round(n_min)==n_min
            dim(2) = n_min;
        end
    end
end

end