classdef (InferiorClasses={?polynomial,?dpvar,?nopvar}) mopvar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This defines the class of mixed-domain PI operators, each element of
% which is a nopvar object
%
%   Pop(r): [   L2^d_1^n_1   ] -->  [ L2^d_1^n_1    ],
                    ...         
%               L2^d_N^n_N   ]      [ L2^d_M^n_M    ]
%
% CLASS properties
% - P.C:    MxN matrix of nopvar objects;
% - P.in.vars:  Nx1 array specifying the number of spatial variables (n_i) in L2^d_i^n_i; -- Note: n_i can be 0
% - P.in.dim:   Nx1 array specifying the number of components (d_i) in L2^d_i^n_i;
% - P.in.dom:   Nx1 array specifying the domains of the spatial variables ;
% - P.dom:  Nx2 array with each row dom(i,:) = [ai,bi] representing the
%           spatial interval along the ith direction on which the operator
%           is defined;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - nopvar
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
% DJ, 01/15/2026: Initial coding

properties
    C = {};
    deg = zeros(0,1);
    dom = zeros(0,2);
    dim = [0,0];
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
    function [dim] = get.dim(obj)
        % % Determine the dimensions of the operator, m x n, from the
        % % dimensions of the coefficient matrices,
        % %     m*prod(deg+1) x n*prod(deg+1).

        % Get the number of monomials and decision variables
        degs = obj.deg;
        if isempty(degs)
            degs = 0;
        end
        nZ = prod(obj.deg+1);
        N = numel(obj.deg);

        % Check the dimensions of the coefficient matrices
        m_min = inf;    n_min = inf;
        m_max = 0;      n_max = 0;
        sz_C = [size(obj.C),1];
        % if sz_C(2)==1
        %     sz_C = sz_C(1);
        % end
        for ii=1:numel(obj.C)
            if isempty(obj.C{ii})
                continue
            end
            [m,n] = size(obj.C{ii});
            m_min = min(m_min,m/nZ);
            m_max = max(m_max,m/nZ);
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