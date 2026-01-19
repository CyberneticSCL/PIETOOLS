classdef (InferiorClasses={?polynomial,?dpvar,?nopvar,?ndopvar})tensopvar
% This function defines the class of 'tensopvar' objects, representing
% coefficient operators acting on distributed monomials to defined 
% distributed polynomials.
%
% CLASS properties
% - C.ops:  1 x nZ cell (for nZ the number of distributed monomials) 
%           specfiying the operators acting on the different distributed
%           monomials. Each element is itself a m x d cell of 'nopvar' 
%           objects, where m is the number of terms, and d the number of 
%           factors in the monomial. In particular, if 
%               Z_{i}(x)=x1^d1 o ... o xp^dp
%           where o denotes the tensor product and d = d1+...+dp,
%           and C.ops{ii} = P, we have
%               (C(i)*Z_{i}(x))(s) = 
%                       (P{1}x1)(s)*(P{2}x1)(s)*...*(P{d}xp)(s)
%           Note that Z_{i}(x) will always be scalar-valued, but each of
%           the operators P{j} may be matrix-valued, in which case the
%           products such as (P{1}x1)(s)*(P{2}x1)(s) will be elementwise;
% - C.dim:  1x2 array of integers, specifying the dimensions of the
%           matrix-valued distributed polynomial defined by C;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - tensopvar
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
% DJ, 01/16/2026: Initial coding

properties
    ops = {};
    dim = [0,0];
end

methods

    function [C] = tensopvar(varargin) %constructor
        if nargout==0
            % Declare state variables as
            %   polyopvar x1 x2 x3
            for i=1:nargin
                if ischar(varargin{i})
                    if nargout==0
                        assignin('caller', varargin{i}, ndopvar());
                    end
                else
                    error("Input must be strings");
                end
            end
        elseif nargout==1
            % Declare state variable as
            %   C = tensopvar(Pop);
            if nargin==0
                return
            end
            if nargin==1
                if isa(varargin{1},'nopvar') || isa(varargin{1},'ndopvar')
                    C.ops{1} = varargin{1};
                elseif isa(varargin{1},'double') && all(size(varargin{1})==[1,2])
                    C.ops = cell(varargin{1});
                else
                    error("Input must be 'nopvar' or 'ndopvar' object.")
                end
            elseif nargin==2
                if isa(varargin{1},'double') && isa(varargin{2},'double')
                    C.ops = cell(varargin{1},varargin{2});
                end
            else
                error("Too many input arguments")
            end
        end
    end

    function [dim] = get.dim(obj)
        % % Determine the dimensions of the matrix-valued operator
        % % Determine the dimensions of the matrix-valued operator, m x n, 
        % % from the individual operators

        % Check the dimensions of the individual operators
        nr = size(obj.ops,1);
        m_min = inf*ones(nr,1);     n_min = inf*ones(nr,1);
        m_max = zeros(nr,1);        n_max = zeros(nr,1);
        for ii = 1:size(obj.ops,1)
        for jj=1:numel(obj.ops(ii,:))
            if isempty(obj.ops{ii,jj})
                continue
            end
            if isa(obj.ops{ii,jj},'nopvar')
                m = obj.ops{ii,jj}.dim(1);
                n = obj.ops{ii,jj}.dim(2);
                m_min(ii) = min(m_min(ii),m);   n_min(ii) = min(n_min(ii),n);
                m_max(ii) = max(m_max(ii),m);   n_max(ii) = max(n_max(ii),n);
            elseif isa(obj.ops{ii,jj},'cell')
                for kk=1:numel(obj.ops{ii,jj})
                    m = obj.ops{ii,jj}{kk}.dim(1);
                    n = obj.ops{ii,jj}{kk}.dim(2);
                    m_min(ii) = min(m_min(ii),m);   n_min(ii) = min(n_min(ii),n);
                    m_max(ii) = max(m_max(ii),m);   n_max(ii) = max(n_max(ii),n);
                end
            end
        end
        end
        
        % Set the dimensions
        dim = [nan,nan];
        if all(m_min==m_max) && all(round(m_min)==m_min)
            dim(1) = sum(m_min);
        end
        if all(n_min==n_max) && all(round(n_min)==n_min)
            if all(n_min==n_min(1)*ones(size(n_min)))
                dim(2) = n_min(1);
            end
        end
    end

end

end