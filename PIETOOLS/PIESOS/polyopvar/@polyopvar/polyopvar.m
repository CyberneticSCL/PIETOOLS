classdef (InferiorClasses={?polynomial,?dpvar,?sopvar,?tensopvar,?tensopmat})polyopvar
% Class of 'polyopvar' objects, representing polynomials in distributed
% states x = (x_1,...,x_p),
%   P(x) = C*Z(x),
% for Z(x) a vector of distributed monomials in x, and C a matrix of
% tensor-PI operators acting on these monomials.
%
% CLASS properties
% - P.C:        'tensopmat' object with nr row-blocks and a column-block
%               for each monomial (row of P.Z), the tensor-PI operators
%               acting on the monomials;
% - P.Z:        nZ x d sparse array of integers, the nZ distributed
%               monomials: Z(i,:) = [k1,...,kd] represents
%                   x_{k1} o ... o x_{kd},
%               for o the tensor product; its degree is the number of
%               nonzeros in Z(i,:). If the operators acting on monomial i
%               have no more factors than this degree, zeros are padding,
%               which initialize moves to the right. If they have more,
%               factor j acts on state x_{Z(i,j)} if Z(i,j)>0, and is a
%               multiplier from R (an operator with input dimension 1 and
%               no input variables, e.g. a known function m(s) = 1-s^2)
%               if Z(i,j)=0, also for j>d. E.g. Z(i,:)=[1,0,2], acted on
%               by the factors identity, m and identity with pointwise
%               output, gives x_1(s)*m(s)*x_2(s). A monomial of degree
%               zero has only multipliers from R as factors;
% - P.varname:  p x 1 'cellstr', the names of the p state variables
%               x_1,...,x_p;
% - P.matdim:   (dependent) 1 x 2 array, the row and column dimensions of
%               the polynomial: [m,1], or [1,k] for an adjoint with a
%               weight with k columns (the number of rows of P.C for a
%               scalar weight);
% - P.dims:     'struct' with fields
%       in:     p x 1 array, the vector size of each state;
%      out:     nr x 1 array, the row dimension of each row-block of P.C;
% - P.vars:     'struct' with fields
%       in:     1 x N 'cellstr', the spatial variables of the states;
%      out:     1 x M 'cellstr', the output variables of the polynomial;
% - P.dom:      'struct' with fields
%       in:     N x 2 'double' array, the domain of each variable in
%               vars.in;
%      out:     M x 2 'double' array, the domain of each variable in
%               vars.out;
% - P.depmat:   'struct' with fields
%       in:     p x N 'logical' array, element (i,j) true if state x_i
%               depends on variable vars.in{j};
%      out:     nr x M 'double' array, element (i,j) the degree of
%               row-block i in variable vars.out{j}, i.e. the number of
%               copies of the variable;
% - P.degmat:   (dependent) nZ x p sparse array, the degree of each
%               monomial in each state;
% - P.vartab:   (dependent, hidden) 'polyopvar' object [x_1;...;x_p] of
%               the states;
% - P.is_adjoint: 'logical', true if P represents the L2-adjoint U'*W of
%               the polynomial U defined by C and Z, with weight W, so
%               that U'*V = innerprod(U,W,V) for a polynomial V. Defaults
%               to false;
% - P.weight:   n x k 'double', 'polynomial' or 'dpvar' array, the weight
%               W of an adjoint P = U'*W for U with n rows; a scalar w
%               represents w*eye(n). Defaults to 1, the only value if
%               P.is_adjoint is false (a matrix times a polynomial is
%               multiplied into C);
%
% NOTES
% - Only scalar or column-vector-valued polynomials, and their adjoints
%   (rows), are supported.
% - An adjoint is a functional, not a polynomial in the states. Only
%   U'*V, U'*W, sums of adjoints with the same weight, the negative of an
%   adjoint, and the horizontal concatenation of adjoints are supported;
%   other operations throw an error.
% - vars.in and vars.out are sorted, and dom and depmat follow their
%   order. varname is sorted, and dims.in, depmat.in and the values of Z
%   follow its order. The rows of Z are not ordered.
% - The hidden property chkval is set to false when a field is assigned
%   outside of the methods of the class, e.g. f.C = Cnew. initialize
%   checks the fields and sets chkval to true.

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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/24/2026: Initial coding;

properties
    C = tensopmat();
    Z = zeros(1,0);
    varname = cell(0,1);
    vars = struct('in',{cell(1,0)},'out',{cell(1,0)});
    dom = struct('in',zeros(1,0),'out',zeros(1,0));
    dims = struct('in',zeros(1,0),'out',zeros(1,0));
    depmat = struct('in',false(1,0),'out',zeros(1,0));
    is_adjoint = false;
    weight = 1;
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
% P = POLYOPVAR(VARNAME,DIMS,VARS,DOM,DEPMAT) returns the vector
% P = [x_1;...;x_p] of the distributed states named VARNAME, each acted
% on by the identity operator. DIMS may be omitted:
%   P = POLYOPVAR(VARNAME,VARS,DOM,DEPMAT).
% P = POLYOPVAR() returns an empty 'polyopvar' object.
% Without output, POLYOPVAR('x1(s1,s2)','x2(s1)',...), or
%   polyopvar x1(s1,s2) x2(s1)
% declares each state, and its spatial variables as 'pvar' objects with
% domain [0,1], in the workspace of the caller. A variable t is taken to
% be time, and ignored.
%
% INPUTS
% - varname:    p x 1 'cellstr' (or char), the names of the states;
% - dims:       (optional) p x 1 array of positive integers, the vector
%               size of each state. Defaults to ones(p,1);
% - vars:       (optional) 1 x N 'cellstr' or 'pvar' array, or 'struct'
%               with field 'in', the spatial variables. Defaults to none;
% - dom:        (optional) N x 2 array (or 1 x 2 for all variables), or
%               'struct' with field 'in', the domains of the variables.
%               Defaults to [0,1];
% - depmat:     (optional) p x N 'logical' array, true if x_i depends on
%               variable j. Defaults to true(p,N);
%
% OUTPUTS
% - P:          sum(dims) x 1 'polyopvar' object [x_1;...;x_p].
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
        if numel(unique(varname))<p
            error("The names of the state variables must be distinct.")
        end
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
        if any(dom.in(:,1)>=dom.in(:,2))
            error("The lower limit of each spatial domain must be smaller than the upper limit.")
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
        % Declare the distributed polynomial variable, with the output space
        % of Cop, see set.C
        P.C = Cop;
        % Canonicalize order of state variables (varname), so declaring
        % the same variables in a different order produces the same
        % varname/dims.in/depmat.in.
        [varname_sorted,var_sort_idx] = sort(varname);
        [~,var_inv_idx] = sort(var_sort_idx);
        P.Z = sparse(var_inv_idx(:));
        P.varname = varname_sorted;
        P.dims.in = dims(var_sort_idx);
        % Canonicalize order of input vars so states over the same        
        % variables declared in different orders agree on vars.in/dom.in.
        [vars_in,sort_idx] = sort(vars.in);
        dom_in = dom.in(sort_idx,:);
        P.vars.in = vars_in;
        P.dom.in = dom_in;
        P.depmat.in = depmat(var_sort_idx,sort_idx);
        P.chkval = true;
    end
end

function obj = set.C(obj,c)
% OBJ = SET.C(OBJ,C) sets obj.C = C, and the output space of OBJ to that
% of C. Called on every assignment to obj.C, also within the methods of
% the class.
%
% INPUTS
% - obj:    'polyopvar' object;
% - c:      'tensopmat' object, the new value of obj.C;
%
% OUTPUTS
% - obj:    'polyopvar' object with obj.C = c. If c.chkval is true,
%           obj.dims.out, obj.vars.out, obj.dom.out and obj.depmat.out
%           are set to those of c; otherwise, nothing else is changed.
%           The input space (of the states) is not changed;
    obj.C = c;
    if isa(c,'tensopmat') && c.chkval
        obj.dims.out = c.dims.out;
        obj.vars.out = c.vars.out;
        obj.dom.out = c.dom.out;
        obj.depmat.out = c.depmat.out;
    end
end

function matdim = get.matdim(obj)
    % MATDIM = GET.MATDIM(OBJ) returns the matrix dimensions of the
    % distributed polynomial OBJ.
    %
    % INPUTS
    % - obj:      'polyopvar' object;
    %
    % OUTPUTS
    % - matdim:   1 x 2 array, the row and column dimensions of obj;
    %             [1,k] for an adjoint U'*W, with k the number of columns
    %             of W (of rows of U for a scalar W).
    [nr,nc] = size(obj.C);
    degmat1 = obj.degmat;
    nZ = size(degmat1,1);
    szZ = prod(repmat(obj.dims.in(:)',nZ,1).^degmat1,2);
    matdim = [nr,nc/sum(szZ)];
    if obj.is_adjoint
        % The row U'*W, with W of size n x k (or scalar, W*eye(n)) and U of
        % size n x 1, has size 1 x k
        matdim = [matdim(2),matdim(1)];
        if ~all(size(obj.weight)==1)
            matdim(2) = size(obj.weight,2);
        end
    end
end


function vartab = get.vartab(obj)
    % VARTAB = GET.VARTAB(OBJ) returns the states of OBJ as a 'polyopvar'
    % object, with their sizes, spatial variables and domains.
    %
    % INPUTS
    % - obj:      'polyopvar' object;
    %
    % OUTPUTS
    % - vartab:   'polyopvar' object [x_1;...;x_p] of the states of obj.
    vartab = polyopvar(obj.varname,obj.dims.in,obj.vars.in,obj.dom.in,obj.depmat.in);
end


function degmat = get.degmat(obj)
% DEGMAT = GET.DEGMAT(OBJ) returns the degree of each monomial of OBJ in
% each state.
%
% INPUTS
% - obj:      'polyopvar' object;
%
% OUTPUTS
% - degmat:   nZ x p sparse array, the degree of each monomial in each
%             state.
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