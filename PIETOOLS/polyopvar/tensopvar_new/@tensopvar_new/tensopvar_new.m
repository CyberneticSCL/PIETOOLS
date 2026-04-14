classdef (InferiorClasses={?polynomial,?dpvar,?nopvar,?ndopvar})tensopvar_new
% This function defines the class of 'tensopvar' objects, representing
% linear combinations of tensor products of PI operators,
%
%   (C*(x1 o ... o xd))(s) = sum_{i=1}^{m} (Ri1*x1)(s1) kron ... kron (Rid*xd)(sd)
%
% where the variables s1 through sd may overlap, and where o denotes the
% tensor product, and kron the Kronecker product
%
% CLASS properties
% - C.ops:  m x d cell with each element {i,j} a 'nopvar' or 'dopvar'
%           object representing the operator Rij, in the ith term and jth
%           factor. Note that this operator can map between different
%           function spaces, and can also be specified as a tensopvar
%           object itself;
% - C.var1: M x p 'pvar' array specifying the spatial variables s on which
%           the output of the operator depends. Here, variables in the same
%           row are defined on the same interval, as occurs when taking the
%           tensor product of a state variable with itself,
%               v^{o n}(s1,...,sn)=v(s1)*...*v(sn)
%           Variables in different rows may exist on different intervals,
%           as may be the case for multivariate PDEs;
% - C.var2: N x q 'pvar' array specifying the dummy variables used in the
%           definition of the operator; 
% - C.dom1: M x 2 array of type 'double', specifying for each of the
%           spatial variables the interval on which it is defined;
% - C.dom2: N x 2 array of type 'double', specifying for each of the
%           dummy variables the interval on which it is defined;
% - C.dim:  1 x 2 array specifying the matrix dimension of the operator;
% - C.use_kron: 1 x 2 binary array specifying whether the products of the
%               different are computed pointwise (false) or using the
%               Kronecker product (true) along the row and column
%               dimensions. If use_kron(1)==false, then the row dimensions
%               of all the operators Rij must match. If use_kron(2)==false,
%               then the column dimensions of all the operators Rij must
%               match;

% DEPENDENT properties
% - C.varnames: d x 2 cell of which each row j specifies the primary (first
%               column and dummy (second column) variables in which the
%               operators Rij for i=1:m are defined. Each element is a
%               cellstr object of variable names.
% - C.doms:     d x 2 cell of which each row j specifies the domains on
%               which the primary (first column) and dummy (second column)
%               variables appearing in Rij are defined;


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
% DJ, 04/14/2026: Initial coding

properties
    ops = {};
    var1 = polynomial(zeros(0,0));
    var2 = polynomial(zeros(0,0));
    dom1 = zeros(0,0);
    dom2 = zeros(0,0);
    dim = zeros(1,2);
    use_kron = true(1,2);
end
properties (Dependent,Hidden)
    varnames
    doms
    dims
end

methods

    function [C] = tensopvar_new(varargin) %constructor
        if nargout==0
            % Declare coefficient operator as
            %   tensopvar C1 C2 C3
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
            % Declare coefficient operator as
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


    function varnames = get.varnames(obj)
        % Determine the spatial and dummy variables on which each factor
        % in the tensor-PI operator depends
        d = size(obj.ops,2);
        varnames = cell(d,2);
        for j=1:d
            op = obj.ops{1,j};
            varnames{j,1} = pvar2varname(op.var1);
            varnames{j,2} = pvar2varname(op.var2);
        end
    end

    function doms = get.doms(obj)
        % Determine the spatial domains of the variables on which each
        % factor in the tensor-PI operator depends
        d = size(obj.ops,2);
        doms = cell(d,2);
        for j=1:d
            op = obj.ops{1,j};
            % Input and output domain of 'nopvar' objects is the same
            doms{j,1} = op.dom;
            doms{j,2} = op.dom;
        end
    end

    function dims = get.dims(obj)
        % Determine the matrix dimensions of the operators in each factor
        % of the tensor-PI operator
        d = size(obj.ops,2);
        dims = zeros(d,2);
        for j=1:d
            op = obj.ops{1,j};
            dims(j,:) = op.dim;
        end
    end
end

end