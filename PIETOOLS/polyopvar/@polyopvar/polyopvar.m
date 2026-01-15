classdef (InferiorClasses={?polynomial,?dpvar,?nopvar,?ndopvar})polyopvar
% This function defines the class of 'polyopvar' objects, representing
% polynomial functions on distributed states.
%
%
% CLASS properties
% - P.varname:  p x 1 cellstr, specifying the names of the p distributed
%               state variables on which the polynomial is defined,
%                   (x1,x2,...,xp);
% - P.degmat:   nZ x p array of integers, representing the basis of
%               distributed monomials in terms of which the polynomial is
%               defined. If P.degmat(i,:) = [d1,...,dp], then the ith
%               distributed monomial is defined by
%                   [Z(x)]_{i} = x1^d1 o ... o xp^dp,
%               where o denotes the tensor product, and where
%                   xi^di = xi o xi o ... o xi } <-- di products
% - P.C:        1 x nZ cell specifying the coefficients acting on the
%               different distributed monomials, so that
%                   P = C{1}*[Z(x)]_{1} + ... + C{nZ}*[Z(x)]_{nZ}
%               If [Z(x)]_{i} = xk for some k in 1,...,p, then C{i} can
%               just be an 'nopvar' or 'ndopvar' object, representing an
%               operator acting on the state xk. Otherwise, if
%                   [Z(x)]_{i} = x1^d1 o ... o xp^dp,
%               then C{i} is itself a m x d cell for d=d1+...+dp, with each
%               element an 'nopvar' or 'ndopvar' object
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - polyopvar
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
        C = cell(0,0);
        degmat = zeros(0,0);
        varname = cell(0,1);
        pvarname = cell(0,1);
        varmat = zeros(0,0);
        dom = zeros(0,2);
        dim = [0,0];
    end

    methods

        function [P] = polyopvar(varargin) %constructor
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
                %   x = polyopvar('x',s1,[0,1]);
                if nargin==0
                    return
                end
                varname = varargin{1};
                % Check which spatial variables x depends on
                if nargin<=1
                    pvarname = {};
                else
                    if isa(varargin{2},'cellstr')
                        pvarname = varargin{2};
                    elseif isa(varargin{2},'polynomial')
                        if ~ispvar(varargin{2})
                            error("Spatial variables must be specified as array of 'pvar' objects.")
                        end
                        pvarname = cell(numel(varargin{2}));
                        for jj=1:numel(pvarname)
                            pvarname(jj) = varargin{2}(jj).varname;
                        end
                    else
                        error("Spatial variables must be specified as array of 'pvar' objects.")
                    end
                end
                N = numel(pvarname);
                % Check the domain on which the spatial variables exist
                if nargin<=2
                    dom = [zeros(N,1),ones(N,1)];
                else
                    if all(size(varargin{3})==[1,2])
                        dom = repmat(varargin{3},[N,1]);
                    elseif all(size(varargin{3}==[N,2]))
                        dom = varargin{3};
                    else
                        error("The spatial domain should be specified as px2 array for 'p' variables.")
                    end
                end
                % Declare an identity operator mapping x to x
                C = nopvar();
                C.dom = dom;
                C.vars = polynomial(pvarname);
                C.deg = zeros(N,1);
                C.C = cell([3*ones(1,N),1]);
                C.C{1} = 1;
                % Declare the distributed polynomial variable
                P.C = {C};
                P.degmat = 1;
                P.varname = varname;
                P.pvarname = pvarname;
                P.varmat = true(1,N);
                P.dom = dom;
            end
        end

        function [dim] = get.dim(obj)
            % % Determine the dimensions of the matrix-valued distributed
            % polynomial
            dim = [1,1];
        end

    end
end