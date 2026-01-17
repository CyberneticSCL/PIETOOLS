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
        dim = [1,1];
    end

end

end