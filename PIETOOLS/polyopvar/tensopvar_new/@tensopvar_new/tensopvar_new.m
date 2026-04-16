classdef (InferiorClasses={?polynomial,?dpvar,?nopvar,?ndopvar})tensopvar_new
% This function defines the class of 'tensopvar' objects, representing
% linear combinations of tensor products of PI operators,
%
%   C = sum_{i=1}^{m} Ri1 \otimes ... \otimes Rid
%
% where \otimes denotes the tensor product
%
% CLASS properties
% - C.ops:  m x d cell with each element {i,j} a 'nopvar' or 'dopvar'
%           object representing the operator Rij, in the ith term and jth
%           factor. Note that this operator can map between different
%           function spaces, and can also be specified as a tensopvar
%           object itself;
% - C.dim:  1 x 2 array specifying the matrix dimensions of the operator;
% - C.vars: N x 2 'pvar' array specifying the spatial variables (first 
%           column) and dummy variables (second column) in terms of which
%           the operator is defined;
% - C.dom:  N x 2 'double' array specifying the intervals on which each of
%           the spatial variables are defined;
% - C.depmat1:  N x d 'double' array with element (i,j) specifying whether
%               the operators in factor j map to spatial variable i. If
%               factor j is a tensopvar object, element (i,j) may be bigger
%               than 1, to indicate the operator maps to a higher-degree
%               monomial in the spatial variable;
% - C.depmat2:  N x d 'double' array with element (i,j) specifying whether
%               the operators in factor j maps from spatial variable i. If
%               factor j is a tensopvar object, element (i,j) may be bigger
%               than 1, to indicate the operator maps from a higher-degree
%               monomial in the spatial variable;
% - C.type: 1 x 2 logical array specifying what type of tensor product is
%           taken, set to true to indicate that a proper tensor product is
%           taken along the row and column dimensions, or false to indicate
%           that a pointwise product is taken. For example, if Rop is
%           defined by Pop1,Pop2:L2^n[0,1] --> L2^m[0,1], then different
%           values of "type" imply
%               - (1,1):    Rop: L2^{n^2}[[0,1]^2] --> L2^{m^2}[[0,1]^2]
%               - (0,1):    Rop: L2^{n^2}[[0,1]^2] --> L2^{m}[0,1]
%               - (1,0):    Rop: L2^{n}[0,1] --> L2^{m^2}[[0,1]^2]
%               - (0,0):    Rop: L2^{n}[0,1] --> L2^{m}[0,1]
%
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
% DJ, 04/15/2026: Initial coding

properties
    ops = {};
    vars = polynomial(zeros(0,2));
    dom = zeros(0,2);
    depmat1 = zeros(0,0);
    depmat2 = zeros(0,0);
    type = [true,true];
end
properties (Dependent)
    dim
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
                    C = ndopvar2tensopvar_new(varargin{1});
                else
                    error("Input must be 'nopvar' or 'ndopvar' object.")
                end
            else
                error("Too many input arguments")
            end
        end
    end

    function dim = get.dim(obj)
        % Determine the matrix dimensions of the operator
        if isempty(obj.ops)
            dim = [0,0];
            return
        end
        % Establish the dimensions based on those of each factor
        dim = [nan,nan];
        dim_arr = obj.dims;
        if obj.type(1) 
            % Kronecker product
            dim(1) = prod(dim_arr(:,1));
        elseif max(dim_arr(:,1))==min(dim_arr(:,1))
            % Elementwise product
            dim(1) = dim_arr(1,1);
        end
        if obj.type(2) 
            % Kronecker product
            dim(2) = prod(dim_arr(:,2));
        elseif max(dim_arr(:,2))==min(dim_arr(:,2))
            % Elementwise product
            dim(2) = dim_arr(1,2);
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