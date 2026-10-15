classdef (InferiorClasses={?tensopvar,?ndopvar,?nopvar,?sopvar,?dpvar,?polynomial})tensopmat
% Class of 'tensopmat' objects, representing a matrix structure wherein
% each element is a tensor-PI operator
%   Cop = [C11, ..., C1d]
%         [ : ,  . ,  : ]
%         [Cp1, ..., Cpd]
%
% CLASS properties
% - C.ops:  p x d cell with each element {i,j} a 'tensopvar' object Cij
%           representing a tensor-PI operator between different function
%           spaces, or an 'intop' object representing a functional, or a
%           matrix M of type 'double', 'polynomial' or 'dpvar' acting on R
%           (see NOTES). An empty element represents the zero operator;
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
% - C.dims:     struct with fields
%       in:     d x 1 array specifying the dimensions of the input to each
%               column-block (matdim(2) = sum(dims.in));
%      out:     p x 1 array specifying the dimensions of the output of each
%               row-block (matdim(1) = sum(dims.out));
% - C.depmat:   struct with fields
%         in:   d x N 'double' array with element (j,k) specifying to what
%               degree the operators in column-block j map from spatial
%               variable k;
%        out:   p x M 'double' array with element (i,k) specifying to what
%               degree the operators in row-block i map to spatial
%               variable k. If a block is a 'tensopvar' object, element
%               (i,k) may be bigger than 1, to indicate that the operator
%               maps to a higher-degree monomial in the spatial variable;
%
% DEPENDENT properties
% - C.matdim:   1 x 2 array specifying the matrix dimensions of the
%               tensor-PI operator (computed from C.dims);
%
% NOTES
% - vars.in and vars.out always store the spatial variables in the same
%   (sorted) order; this order corresponds positionally to depmat.in/
%   depmat.out and dom.in/dom.out.
% - A block C.ops{i,j} = M of type 'double', 'polynomial' or 'dpvar' is an
%   m x n matrix, representing the operator from R^n to L2^m[S] that maps x
%   to the function M*x, where S are the spatial variables appearing in M
%   (the decision variables of a 'dpvar' are not spatial variables). Such a
%   block can only act on R, so column j has depmat.in(j,:)=0, and no other
%   block in the column can act on a function space. A row consisting only
%   of such blocks maps to the variables appearing in them; if the row has
%   other blocks, these determine the output space, and the variables in M
%   must be among them (a constant is a constant function).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - tensopmat
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
    ops = cell(0,0);
    dims = struct('in',zeros(0,1),'out',zeros(0,1));
    vars = struct('in',{cell(1,0)},'out',{cell(1,0)});
    dom = struct('in',zeros(0,2),'out',zeros(0,2));
    depmat = struct('in',zeros(0,0),'out',zeros(0,0));
end
properties (Dependent)
    matdim
end
properties (Hidden)
    chkval = false; % indicates whether we have checked whether the object is properly declared
end

methods
    function obj = tensopmat(Cop,vars,dom) 
        % OBJ = TENSOPMAT(COP,VARS,DOM) returns a 'tensopmat' object OBJ
        % representing a block matrix of tensor-PI operators.
        %
        % OBJ = TENSOPMAT(COP) converts the object COP, a single operator
        % or a cell of operators, to a 'tensopmat' object.
        %
        % OBJ = TENSOPMAT(COP,VARS,DOM) additionally specifies the
        % variables and their domains as structs with fields 'in' and
        % 'out', see `mat2tensopvar`, for blocks of type 'double',
        % 'polynomial' or 'dpvar'.
        %
        % OBJ = TENSOPMAT() returns an empty 'tensopmat' object.
        %
        % INPUTS
        % - Cop:    (optional) 'tensopvar', 'intop', 'sopvar', 'opvar' or
        %           'opvar2d' object, matrix of type 'double',
        %           'polynomial' or 'dpvar' acting on R, or p x d cell of
        %           such objects, representing the blocks of OBJ (empty
        %           elements of a cell represent zero operators);
        % - vars:   (optional) struct with fields 'in' (1 x N cellstr) and
        %           'out' (1 x M cellstr) specifying the names of the
        %           input and output variables. Only with dom, and only
        %           for Cop a cell or a matrix;
        % - dom:    (optional) struct with fields 'in' (N x 2 array) and
        %           'out' (M x 2 array) specifying the intervals on which
        %           the input and output variables are defined. Only with
        %           vars;
        %
        % OUTPUTS
        % - obj:    'tensopmat' object representing the block matrix
        %           defined by Cop;
        %
        % NOTES
        % - The domain of each variable appearing in a block of type
        %   'double', 'polynomial' or 'dpvar' is taken from the other
        %   blocks, or from vars and dom.
        if nargin==0
            return
        end
        if nargin==2 || nargin>3  
            error("Variables and domains must be specified together, as inputs 'vars' and 'dom'.") 
        elseif nargin==3 && ~(isa(Cop,'cell') || is_native_op(Cop))         
            error("Variables and domains can only be specified for a cell of blocks, or a matrix of type 'double', 'polynomial' or 'dpvar'.") 
        end 
        if is_native_op(Cop) 
            Cop = {Cop};    % single block acting on R 
        end
        if isa(Cop,'cell')
            % Declare a 'tensopmat' defined by a cell of 'tensopvar' 
            % objects
            obj = tensopmat();
            for j=1:numel(Cop)
                if isempty(Cop{j})
                    continue
                elseif isa(Cop{j},'opvar') || isa(Cop{j},'opvar2d')
                    Cop{j} = sopvar2tensopvar(opvar2sopvar(Cop{j}));
                elseif isa(Cop{j},'sopvar') || isa(Cop{j},'cell')
                    Cop{j} = sopvar2tensopvar(Cop{j});
                elseif ~isa(Cop{j},'tensopvar') && ~isa(Cop{j},'intop') && ~is_native_op(Cop{j})
                    error("Elements of cell must be of type 'tensopvar', 'intop', 'double', 'polynomial' or 'dpvar' for conversion to a 'tensopmat'.")
                end
            end
            obj.ops = Cop;
            if nargin==3 
                % Variables and domains of the blocks of type 
                % double/polynomial/dpvar, taken as given by initialize if
                % not defined by the other blocks  
                if ~isstruct(vars) || ~isfield(vars,'in') || ~isfield(vars,'out') || ~isstruct(dom) || ~isfield(dom,'in') || ~isfield(dom,'out') 
                    error("Variables and domains must be specified as structs with fields 'in' and 'out'.") 
                elseif ~iscellstr(vars.in) || ~iscellstr(vars.out) || size(dom.in,1)~=numel(vars.in) || size(dom.out,1)~=numel(vars.out)
                    error("Variables must be specified as 'cellstr' arrays, with a domain (row of an N x 2 array) for each.") 
                end                                                         
                obj.vars = struct('in',{reshape(vars.in,1,[])},'out',{reshape(vars.out,1,[])});
                obj.dom = struct('in',reshape(dom.in,[],2),'out',reshape(dom.out,[],2)); 
            end                                                         
            obj = initialize(obj);
        elseif isa(Cop,'opvar') || isa(Cop,'opvar2d')
            obj = tensopmat(opvar2sopvar(Cop));
        elseif isa(Cop,'sopvar')
            obj = tensopmat(sopvar2tensopvar(Cop));
        elseif isa(Cop,'tensopvar') || isa(Cop,'intop')
            obj.ops = {Cop};
            obj.dims.in = size(Cop,2);  obj.dims.out = size(Cop,1);
            obj.vars = Cop.vars;
            obj.dom = Cop.dom;
            obj.depmat.in = Cop.deps.in;
            obj.depmat.out = Cop.deps.out;
            obj.chkval = true;
        else
            error("Input must be object of type 'tensopvar'.")
        end
    end

    function matdim = get.matdim(obj)
        % Determine the matrix dimensions of the operator
        matdim = [sum(obj.dims.out),sum(obj.dims.in)];
    end
end

end
