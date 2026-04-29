classdef (InferiorClasses={?tensopvar})tensopmat
% This function defines the class of 'tensopmat' objects, representing
% a matrix structure wherein each element is a tensor-PI operator
%   Cop = [C11, ..., C1d]
%         [ : ,  . ,  : ]
%         [Cp1, ..., Cpd]
%
% CLASS properties
% - C.ops:  p x d cell with each element {i,j} a 'tensopvar' object Cij
%           representing a tensor-PI operator between different function
%           spaces.

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
% DJ, 04/15/2026: Initial coding

    properties
        ops = cell(0,0);
        vars = polynomial(zeros(0,2));
        dom = zeros(0,2);
        depmat1 = zeros(0,0);
        depmat2 = zeros(0,0);
    end
    properties(Dependent)
        dim
    end


    methods
        function obj = tensopmat(Cop)
            %TENSOPMAT Construct an instance of this class
            % Converts a single 'tensopvar' or 'intop' object to a
            % tensopmat object
            if nargin==0
                return
            end
            if isa(Cop,'cell')
                % Declare a 'tensopmat' defined by a cell of 'tensopvar' 
                % objects
                obj = tensopmat();
                for j=1:numel(Cop)
                    if isempty(Cop{j})
                        continue
                    elseif isa(Cop{j},'nopvar') || isa(Cop{j},'ndopvar') || isa(Cop{j},'cell')
                        Cop{j} = ndopvar2tensopvar(Cop{j});
                    elseif ~isa(Cop{j},'tensopvar')
                        error("Elements of cell must be of type 'tensopvar' for conversion to a 'tensopmat'.")
                    end
                end
                obj.ops = Cop;
                % Determine the spatial variables and domain
                [vars,dom] = get_vars(obj);
                obj.vars = vars;
                obj.dom = dom;
                % Determine the dependency arrays
                [depmat1,depmat2] = get_deps(obj);
                obj.depmat1 = depmat1;
                obj.depmat2 = depmat2;
            elseif isa(Cop,'nopvar')
                obj = tensopmat(ndopvar2tensopvar(Cop));
            elseif isa(Cop,'tensopvar')
                obj.ops = {Cop};
                obj.vars = Cop.vars;
                obj.dom = Cop.dom;
                obj.depmat1 = Cop.dep(1,:);
                obj.depmat2 = Cop.dep(2,:);
            elseif isa(Cop,'intop')
                obj.ops = {Cop};
                obj.vars = polynomial([Cop.pvarname(1,1),[Cop.pvarname(1,1),'_dum']]);
                obj.dom = Cop.dom;
                obj.depmat1 = 0;
                obj.depmat2 = size(Cop.pvarname,2);
            else
                error("Input must be object of type 'tensopvar'.")
            end
        end

        function dim = get.dim(obj)
            % Get matrix dimensions of each tensor-PI operator in the
            % structure
            dim1 = zeros(size(obj.ops));
            dim2 = zeros(size(obj.ops));
            for j=1:numel(obj.ops)
                [dim1(j),dim2(j)] = size(obj.ops{j});
            end
            dim1 = max(dim1,[],2);
            dim2 = max(dim2,[],1)';
            dim = {dim1,dim2};
        end
    end
end