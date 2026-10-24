classdef (InferiorClasses={?polynomial,?dpvar,?nopvar,?ndopvar,?sopvar})tensopvar
% Class of 'tensopvar' objects, representing linear combinations of tensor
% products of PI operators,
%
%   C = sum_{i=1}^{K} Ri1 \otimes ... \otimes Rid
%
% where \otimes denotes the tensor product. This operator will map
%   
%   C : L2[[a_in1,b_in1]^{d_in1} x ... x L2[[a_inN,b_inN]^{d_inN}
%      --> L2[[a_out1,b_out1]^{d_out1} x ... x L2[[a_outM,b_outM]^{d_outM}
%
% CLASS properties
% - C.ops:  K x d cell with each element {i,j} a 'sopvar' (or 'intop')
%           object representing the operator Rij, in the ith term and jth
%           factor. This operator can map between different function
%           spaces, and can also be a 'tensopvar' object itself (a nested
%           operator);
% - C.otype: 1 x 2 array of integers 1 or 2 specifying what type of tensor
%           product is taken along the row (otype(1)) and column
%           (otype(2)) dimensions: 2 indicates that a proper tensor
%           product is taken, and 1 that a pointwise product is taken.
%           Defaults to [2,2]. For example, if Rop is defined by
%           Pop1,Pop2:L2^n[0,1] --> L2^m[0,1], then different values of
%           otype imply
%               - (2,2):    Rop: L2^{n^2}[[0,1]^2] --> L2^{m^2}[[0,1]^2]
%               - (1,2):    Rop: L2^{n^2}[[0,1]^2] --> L2^{m}[0,1]
%               - (2,1):    Rop: L2^{n}[0,1] --> L2^{m^2}[[0,1]^2]
%               - (1,1):    Rop: L2^{n}[0,1] --> L2^{m}[0,1]
%
% PROTECTED properties
% These properties are determined from the operators in C.ops, and should
% not be set by the user
%
% - C.vars:     struct with fields
%       in:     1 x N cellstr specifying the variable names of the
%               function which the operator maps from;
%      out:     1 x M cellstr specifying the variable names of the
%               function which the operator maps to;
% - C.dom:      struct with fields
%       in:     N x 2 'double' array specifying the intervals on which
%               each of the input variables is defined;
%      out:     M x 2 'double' array specifying the intervals on which
%               each of the output variables is defined;
% - C.dims:     d x 2 array specifying the matrix dimensions of each factor
%               in the tensor-PI operator, from which C.matdim follows.
%               C.dims.in and C.dims.out return C.dims(:,2) and
%               C.dims(:,1), the input and output dimensions of the
%               factors;
% - C.depmat:   struct with fields
%         in:   d x N 'double' array with element (i,j) specifying whether
%               the operator in factor i maps FROM spatial variable j. If
%               factor i is a 'tensopvar' object, element (i,j) may be
%               bigger than 1, to indicate the operator maps from a
%               higher-degree monomial in the spatial variable;
%        out:   d x M 'double' array with element (i,j) specifying whether
%               the operator in factor i maps TO spatial variable j. If
%               factor i is a 'tensopvar' object, element (i,j) may be
%               bigger than 1, to indicate the operator maps to a
%               higher-degree monomial in the spatial variable;
%
% DEPENDENT properties
% These properties are computed from the protected properties
%
% - C.matdim:   1 x 2 array specifying the matrix dimensions of the
%               tensor-PI operator, computed from C.dims: along each
%               dimension, the product of the dimensions of the factors
%               for a Kronecker product (otype=2), or their common value
%               for a pointwise product (otype=1);
% - C.deps:     struct with fields
%       in:     1 x N array of integers, specifying for each of the
%               input variables sk to what degree d_ink the input depends
%               on this variable;
%      out:     1 x M array of integers, specifying for each of the
%               output variables sk to what degree d_outk the output
%               depends on this variable;
%               these are computed from C.depmat, along each dimension as
%               the sum over the factors for a Kronecker product
%               (otype=2), or the maximum over the factors for a pointwise
%               product (otype=1);
%
% NOTES
% - vars.in and vars.out always store the spatial variables in the same
%   (sorted) order; this order corresponds positionally to depmat.in/
%   depmat.out, dom.in/dom.out, and dims.

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
% DJ, 10/15/2026: Initial coding;

properties
    ops = {};
    otype = [2,2];
end
properties (Access=protected)
    % These properties are determined from those of the factors in 'ops'
    vars = struct('in',{cell(1,0)},'out',{cell(1,0)});
    dom = struct('in',zeros(0,2),'out',zeros(0,2));
    dims = zeros(0,2);
    depmat = struct('in',zeros(0,0),'out',zeros(0,0));
end
properties (Dependent)
    matdim
    deps
end
properties (Hidden)
    chkval = false; % indicates whether we have checked whether the object is properly declared
end

methods

function [C] = tensopvar(varargin) %constructor
    % C = TENSOPVAR(POP) returns a 'tensopvar' object C representing the
    % tensor-PI operator defined by the 'sopvar' object, or the cell of
    % 'sopvar' objects, POP.
    %
    % C = TENSOPVAR() returns an empty 'tensopvar' object.
    %
    % TENSOPVAR C1 C2 C3 declares the empty 'tensopvar' objects C1, C2, C3
    % in the calling workspace.
    %
    % INPUTS
    % - varargin:   (optional) for C = TENSOPVAR(POP), a single input
    %               POP: 'sopvar' object, or K x d cell of operators
    %               representing the K terms and d factors of the operator
    %               (as in C.ops). For TENSOPVAR C1 C2 C3, the names of the
    %               objects to declare, specified as strings;
    %
    % OUTPUTS
    % - C:      'tensopvar' object representing the tensor-PI operator
    %           defined by POP;
    %
    if nargout==0
        % Declare coefficient operator as
        %   tensopvar C1 C2 C3
        for i=1:nargin
            if ischar(varargin{i})
                if nargout==0
                    assignin('caller', varargin{i}, tensopvar());           % DJ, 10/14/2026
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
            if isa(varargin{1},'sopvar') || isa(varargin{1},'cell')
                C = sopvar2tensopvar(varargin{1});
            else
                error("Input must be a 'sopvar' object, or a cell of 'sopvar' objects.")
            end
        else
            error("Too many input arguments")
        end
    end
end



function matdim = get.matdim(obj)
    % Determine the matrix dimensions of the operator
    if isempty(obj.ops)
        matdim = [0,0];
        return
    end
    % Establish the dimensions based on those of each factor
    matdim = [nan,nan];
    dim_arr = obj.dims;
    if obj.otype(1)==2
        % Kronecker product
        matdim(1) = prod(dim_arr(:,1));
    elseif max(dim_arr(:,1))==min(dim_arr(:,1))
        % Elementwise product
        matdim(1) = dim_arr(1,1);
    end
    if obj.otype(2)==2
        % Kronecker product
        matdim(2) = prod(dim_arr(:,2));
    elseif max(dim_arr(:,2))==min(dim_arr(:,2))
        % Elementwise product
        matdim(2) = dim_arr(1,2);
    end
end



function deps = get.deps(obj)
    % Determine the degrees of the monomial in each of the variables
    % in the input and outputs of the tensor-PI operator
    deps = struct();
    if obj.otype(1)==1
        % Pointwise product along primary variable
        deps.out = max(obj.depmat.out,[],1);
    else
        % Tensor product along primary variable
        deps.out = sum(obj.depmat.out,1);
    end
    if obj.otype(2)==1
        % Pointwise product along dummy variable
        deps.in = max(obj.depmat.in,[],1);
    else
        % Tensor product along dummy variable
        deps.in = sum(obj.depmat.in,1);
    end
end

end

end