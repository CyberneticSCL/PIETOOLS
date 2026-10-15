classdef (InferiorClasses={?polynomial,?dpvar,?nopvar,?ndopvar,?sopvar,?tensopvar,?tensopmat})polyopvar
% This function defines the class of 'polyopvar' objects, representing
% polynomial functions on distributed states.
%
% CLASS properties
% - P.C:        m x n 'tensopmat' object, representing the
%               tensor-PI operators acting on the distributed monomial
%               basis, so that
%                   P(x) = C*Z(x);
% - P.Z:        nZ x d sparse array of integers, representing the
%               distributed monomial basis, where d is the maximal number
%               of factors of the monomials, and nZ the number of
%               monomials. Then
%                   Z(i,:) = [k1,...,kd]
%               implies that Z(i,:) represents the distributed monomial
%                   x_{k1} o ... o x_{kd}
%               where o denotes the tensor product.
%               Elements of Z(i,:) may be zero. The degree of the monomial
%               is the number of nonzeros in Z(i,:). If the tensor-PI
%               operators acting on monomial i have no more factors than
%               this degree, zeros are just padding, and initialize moves
%               them to the right. If they have more factors, then factor
%               j acts on state x_{Z(i,j)} if Z(i,j)>0, and is a
%               multiplier from R (a known function, e.g. m(s) = 1-s^2)
%               to a function space, i.e. an operator with input
%               dimension 1 and no input variables, if Z(i,j)=0 (also for
%               j>d). For example, Z(i,:)=[1,0,2] represents
%               x_1 o m o x_2, on which a tensor-PI operator with three
%               factors with pointwise output, the identity, multiplier
%               m, and the identity, acts to produce x_1(s)*m(s)*x_2(s).
%               num_factors returns the number of factors acting on each
%               monomial. A monomial of degree zero has the multiplier(s)
%               from R as its only factor(s). combine_factors combines
%               consecutive multipliers from R, e.g. Z(i,:)=[1,0,0,2] to
%               [1,0,2], by multiplying them.
% - P.varname:  p x 1 cellstr, specifying the names of the p distributed
%               state variables on which the polynomial is defined,
%                   (x_{1},x_{2},...,x_{p});
% - P.matdim:   1 x 2 array of integers, specifying the row and column
%               dimensions of the matrix-valued polynomial;
% - P.dims:     struct with fields
%       in:     p x 1 array specifying the vector sizes of each of the
%               state variables;
%      out:     nr x 1 array specifying for each row of P.C the row
%               dimension of this row in the polynomial;
% - P.vars:     struct with fields
%       in:     1 x N cellstr specifying the variable names in terms of
%               which the state variables are defined;
%      out:     1 x M cellstr specifying the spatial variable names in
%               terms of which the polynomial is defined;
% - P.dom:      struct with fields
%       in:     N x 2 'double' array specifying the intervals on which
%               each of the input variables is defined;
%      out:     M x 2 'double' array specifying the intervals on which
%               each of the output variables is defined;
% - P.depmat:   struct with fields:
%         in:   p x N 'logical' array with element (i,j) indicating
%               whether or not state variable x_i depends on spatial
%               variable s_j
%        out:   nr x N 'double' array with element (i,j) indicating to
%               what degree the ith row of the output depends on variable
%               s_j.
%
% NOTES
% - vars.in and vars.out always store the spatial variables in the same
%   (sorted) order; this order corresponds positionally to depmat.in/
%   depmat.out and dom.in/dom.out.
% - varname always stores the state variables in the same (sorted) order;
%   this order corresponds positionally to dims.in/depmat.in, and to the
%   values (not rows) of Z.
% - The monomials in Z are NOT necessarily ordered in any particular way.
%   Functions like clean_poly may order the monomials.
% - chkval is set to false when a field is assigned outside of the
%   methods of the class, e.g. f.C = Cnew, see subsasgn. Use initialize
%   to check the fields and set chkval to true.

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
% DJ, 10/15/2026: Initial coding;

properties
    C = tensopmat();
    Z = zeros(1,0);
    varname = cell(0,1);
    vars = struct('in',{cell(1,0)},'out',{cell(1,0)});
    dom = struct('in',zeros(1,0),'out',zeros(1,0));
    dims = struct('in',zeros(1,0),'out',zeros(1,0));
    depmat = struct('in',false(1,0),'out',zeros(1,0));
end
properties(Hidden)
    chkval = false;     % indicates whether we have checked the polynomial is properly declared
end
properties(Dependent)
    matdim;
    degmat;
end
properties(Dependent,Hidden)
    vartab;
end
    
methods

function [P] = polyopvar(varargin) %constructor
    if nargout==0
        % Declare state variables as
        %   polyopvar x1(s1,s2) x2(s1) x3(s3)
        for i=1:nargin
            if ischar(varargin{i})
                obj = varargin{i};
                % Split "x1(s1,s2)" into "x1" and "(s1,s2)"
                is_bracket = ismember(obj,'(');
                if any(is_bracket)
                    if sum(is_bracket)>1 || ~strcmp(obj(end),')')
                        error("Distributed polynomial variables must be specified as e.g. 'x1(s1,s2)'")
                    end
                    % Determine the name of the state variable
                    xname = obj(1:find(is_bracket,1)-1);
                    % Determine the names of the indepdnent
                    % variables
                    pnames = obj(find(is_bracket,1)+1:end-1);
                    pnames = strsplit(pnames,{',',' '});
                    % Declare the independent variables as pvar
                    % objects in the workspace
                    j = 1;
                    while j<=numel(pnames)
                        pname = pnames{j};
                        if isempty(pname) || strcmp(pname,'t')
                            % assume t is temporal variable
                            pnames(j) = [];
                            continue
                        elseif ~isnan(str2double(pname))
                            error("Independent variable names must include a letter.")
                        else
                            assignin('caller',pname,pvar(pname));
                            j = j+1;
                        end
                    end
                    % Declare the state variable as polyopvar
                    % object in the workspace
                    assignin('caller', xname, polyopvar(xname,pnames'));
                else
                    assignin('caller', obj, polyopvar(obj));
                end
            else
                error("Input must be strings");
            end
        end
    elseif nargout==1
        % Declare state variable as
        %   x = polyopvar(varname,dims,vars,dom,depmat);
        % with inputs
        %   - varname:  p x 1 cellstr, names of the state variables;
        %   - dims:     p x 1 array of positive integers, vector size of each
        %               state variable;
        %   - vars:     cellstr or array of 'pvar', the N spatial variables;
        %   - dom:      N x 2 array, interval of each spatial variable
        %               (1 x 2 is repeated for all variables);
        %   - depmat:   p x N logical, element (i,j) true if x_i depends on
        %               spatial variable j;
        % If nargin<5, assume the last inputs are omitted, unless
        % varargin{2} is not numeric (e.g. cellstr or 'polynomial'), in
        % which case dims is omitted instead. For omitted fields, set
        % default values
        %   - dims = ones(p,1)  (scalar-valued states), p = numel(varname);
        %   - depmat = true(p,N);
        %   - dom = repmat([0,1],N,1);
        %   - vars = cell(1,0);
        if nargin==0
            return
        end
        varname = varargin{1};
        if ischar(varname) && size(varname,1)==1
            varname = {varname};
        elseif ~iscellstr(varname)
            error("Variable names must be specified as string of char objects.")
        end
        varname = varname(:);       % store the names as a column
        p = numel(varname);
        % If the 2nd input is not numeric, it specifies vars rather than
        % dims, so insert empty dims to keep positions fixed below.
        if nargin>=2 && ~isnumeric(varargin{2})
            varargin = [varargin(1),{[]},varargin(2:end)]; 
        end
        nargs = numel(varargin);
        % Check the vector sizes of the state variables
        if nargs<2 || isempty(varargin{2})
            dims = ones(p,1);
        elseif numel(varargin{2})~=p || any(mod(varargin{2}(:),1)~=0 | varargin{2}(:)<1) 
            error("State dimensions must be specified as array of %d positive integers.",p) 
        else                     
            dims = reshape(varargin{2},[],1);
        end                      
        % Check which spatial variables x depends on
        if nargs<=2
            vars = struct('in',{cell(1,0)},'out',{cell(1,0)});
            M = 0;
        else
            if iscellstr(varargin{3})
                pvarname = reshape(varargin{3},1,[]);   % vars.in is 1 x N
                M = numel(pvarname);
                vars = struct('in',{pvarname},'out',{pvarname});
            elseif isa(varargin{3},'polynomial')
                if ~ispvar(varargin{3})
                    error("Spatial variables must be specified as array of 'pvar' objects.")
                end
                pvarname = reshape(pvar2varname(varargin{3}),1,[]);
                M = numel(pvarname);
                vars = struct('in',{pvarname},'out',{pvarname});
            elseif isa(varargin{3},'struct') 
                vars = varargin{3};                                        
                M = numel(vars.in);     % only vars.in is used, not vars.out
            else
                error("Spatial variables must be specified as array of 'pvar' objects.")
            end
        end
        % Check the domain on which the spatial variables exist
        if nargs<=3
            dom = [zeros(M,1),ones(M,1)];
            dom = struct('in',dom,'out',dom);
        elseif isa(varargin{4},'struct')
            dom = varargin{4};
        else
            if all(size(varargin{4})==[1,2])
                dom = repmat(varargin{4},[M,1]);
                dom = struct('in',dom,'out',dom);
            elseif all(size(varargin{4})==[M,2])
                dom = varargin{4};
                dom = struct('in',dom,'out',dom);
            else
                error("The spatial domain should be specified as px2 array for 'p' variables.")
            end
        end
        % Check on which spatial variables each state variable depends
        if nargs<=4 || isempty(varargin{5})
            depmat = true(p,M);
        elseif ~(isnumeric(varargin{5}) || islogical(varargin{5})) || ~isequal(size(varargin{5}),[p,M]) 
            error("Dependency array must be specified as a %d x %d logical array.",p,M) 
        else                   
            depmat = logical(varargin{5});
        end
        % % Declare an identity operator mapping x_j to x_j, for each state,
        % % restricted to the spatial variables x_j depends on, and of size
        % % dims(j)
        Ccell = cell(numel(varname),numel(varname));
        for j=1:numel(varname)
            vars_j = struct('in',{vars.in(depmat(j,:))},'out',{vars.in(depmat(j,:))}); 
            dom_j = struct('in',dom.in(depmat(j,:),:),'out',dom.in(depmat(j,:),:));
            Cj = mat2sopvar(eye(dims(j)),vars_j,dom_j);
            Ccell{j,j} = tensopvar(Cj);
        end
        Cop = tensopmat(Ccell);
        % Declare the distributed polynomial variable
        P.C = Cop;
        % Canonicalize order of state variables (varname), so declaring
        % the same variables in a different order produces the same
        % varname/dims.in/depmat.in.
        [varname_sorted,var_sort_idx] = sort(varname);
        [~,var_inv_idx] = sort(var_sort_idx);
        P.Z = sparse(var_inv_idx(:));
        P.varname = varname_sorted;
        P.dims.in = dims(var_sort_idx);     P.dims.out = Cop.dims.out;
        % Canonicalize order of input vars so states over the same        
        % variables declared in different orders agree on vars.in/dom.in.
        [vars_in,sort_idx] = sort(vars.in);
        dom_in = dom.in(sort_idx,:);
        P.vars.in = vars_in;        P.vars.out = Cop.vars.out;
        P.dom.in = dom_in;          P.dom.out = Cop.dom.out;
        P.depmat.in = depmat(var_sort_idx,sort_idx);
        P.depmat.out = Cop.depmat.out;
        P.chkval = true;
    end
end

function matdim = get.matdim(obj)
    % % Determine the dimensions of the matrix-valued distributed
    % % polynomial
    [nr,nc] = size(obj.C);
    degmat1 = obj.degmat;
    nZ = size(degmat1,1);
    szZ = prod(repmat(obj.dims.in(:)',nZ,1).^degmat1,2);
    matdim = [nr,nc/sum(szZ)];
end


function vartab = get.vartab(obj)
    % % Return a polyopvar object specifying the variables that
    % % appear in the distributed polynomial
    % Pass dims and depmat, so states keep their size and spatial dependency
    vartab = polyopvar(obj.varname,obj.dims.in,obj.vars.in,obj.dom.in,obj.depmat.in);
end


function degmat = get.degmat(obj)
% Get the degree of each variable in each of the monomials
nZ = size(obj.Z,1);
nvars = numel(obj.varname);
degmat = spalloc(nZ,nvars,nnz(obj.Z));
for i=1:nvars
    % Establish how many times variable i appears in each monomial
    is_var_i = obj.Z==i;
    degmat(:,i) = sum(is_var_i,2);
end

end

end
end