function [logval,errmsg] = isvalid(Cop)
% [LOGVAL,ERRMSG] = ISVALID(COP) takes a 'tensopmat' object COP
% representing a block matrix of tensor-PI operators, and checks whether
% all the fields are properly specified.
%
% INPUTS
% - Cop:    'tensopmat' object with field Cop.ops specified as a p x d cell
%           array of 'tensopvar' or 'intop' objects, representing a
%           concatenation of tensor-PI operators, or of matrices of type
%           'double', 'polynomial' or 'dpvar' acting on R;
%
% OUTPUTS
% - logval: boolean value specifying whether (true) or not (false) the
%           object is properly specified;
% - errmsg: 'string' used to declare an error message indicating what is
%           wrong with the object if logval = false;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - isvalid
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

% Initialize
logval = true;  % innocent until proven guilty
errmsg = "Tensopmat object is not properly specified: ";
Cop.chkval = false;

% Check that all fields are of appropriate type
if ~isa(Cop.ops,'cell')
    errmsg = errmsg + "Field 'ops' must be specified as a cell.";
    logval = false;
elseif ~isstruct(Cop.dims) || ~isfield(Cop.dims,'in') || ~isfield(Cop.dims,'out')
    errmsg = errmsg + "Field 'dims' must be specified as a struct with fields 'in' and 'out'.";
    logval = false;
elseif ~isstruct(Cop.vars) || ~isfield(Cop.vars,'in') || ~isfield(Cop.vars,'out')
    errmsg = errmsg + "Field 'vars' must be specified as a struct with fields 'in' and 'out'.";
    logval = false;
elseif ~isstruct(Cop.dom) || ~isfield(Cop.dom,'in') || ~isfield(Cop.dom,'out')
    errmsg = errmsg + "Field 'dom' must be specified as a struct with fields 'in' and 'out'.";
    logval = false;
elseif ~isstruct(Cop.depmat) || ~isfield(Cop.depmat,'in') || ~isfield(Cop.depmat,'out')
    errmsg = errmsg + "Field 'depmat' must be specified as a struct with fields 'in' and 'out'.";
    logval = false;
end
if ~logval
    return
end

% Check that the dimensions are properly specified
dims = Cop.dims;
q = numel(dims.in);     p = numel(dims.out);
if any(size(dims.in)~=[q,1]) || any(size(dims.out)~=[p,1])
    errmsg = errmsg + "Input and output dimensions must be specified as column vectors.";
    logval = false;
    return
elseif any(round(abs(dims.in))~=dims.in) || any(round(abs(dims.out))~=dims.out)
    errmsg = errmsg + "Input and output dimensions must be specified as real nonnegative integers.";
    logval = false;
    return
end

% Check that the variables are properly specified
vars = Cop.vars;
if ~iscellstr(vars.in) || ~iscellstr(vars.out)
    errmsg = errmsg + "Input and output variables must be specified as cellstr objects.";
    logval = false;
    return
end
N = numel(vars.in);     M = numel(vars.out);
if ~all(size(vars.in)==[1,N]) || ~all(size(vars.out)==[1,M])
    errmsg = errmsg + "Input and output variables must be specified as row arrays.";
    logval = false;
elseif ~isequal(vars.in,unique(vars.in)) || ~isequal(vars.out,unique(vars.out))
    errmsg = errmsg + "Input and output variable names must be unique and sorted.";
    logval = false;
end
if ~logval
    return
end

% Check that the variable domains are appropriate
dom = Cop.dom;
if ~all(size(dom.in)==[N,2])
    errmsg = errmsg + "Field 'dom.in' must be specified as a N x 2 array of type 'double', for N input variables.";
    logval = false;
elseif ~all(size(dom.out)==[M,2])
    errmsg = errmsg + "Field 'dom.out' must be specified as a M x 2 array of type 'double', for M input variables.";
    logval = false;
elseif any(dom.in(:,1)>=dom.in(:,2)) || any(dom.out(:,1)>=dom.out(:,2))
    errmsg = errmsg + "Lower boundary of domains cannot exceed upper boundary.";
    logval = false;
end

% Check that the dependency arrays are appropriate
depmat = Cop.depmat;
if ~all(size(depmat.in)==[q,N])
    errmsg = errmsg + "Field 'depmat.in' must be specified as a q x N array of type 'double', for N input variables and q coupled inputs.";
    logval = false;
elseif ~all(size(depmat.out)==[p,M])
    errmsg = errmsg + "Field 'depmat.out' must be specified as a p x M array of type 'double', for M input variables and p coupled outputs.";
    logval = false;
elseif any(round(abs(depmat.in(:)))~=depmat.in(:)) || any(round(abs(depmat.out(:)))~=depmat.out(:))
    errmsg = errmsg + "Dependency arrays must be specified as real nonnegative integers.";
    logval = false;
end


% Check that the operators are properly specified
for i=1:p
    for j=1:q
        % Check that the operator is of appropriate type
        op_ij = Cop.ops{i,j};
        if isempty(op_ij)
            continue
        elseif ~isa(op_ij,'tensopvar') && ~isa(op_ij,'intop') && is_native_op(op_ij)
            % Matrix acting on R: no input variables, output variables among those of the row
            if dims.in(j)~=size(op_ij,2) || dims.out(i)~=size(op_ij,1)
                errmsg = errmsg + "Row and column dimensions of blocks in same row-block and column-block must match.";
                logval = false;
                return
            elseif any(depmat.in(j,:))
                errmsg = errmsg + "Blocks of type 'double', 'polynomial' or 'dpvar' act on R, so cannot share a column with an operator acting on a function space.";
                logval = false;
                return
            elseif ~all(ismember(native_vars(op_ij),vars.out(logical(depmat.out(i,:)))))
                errmsg = errmsg + "Variables in blocks of type 'double', 'polynomial' or 'dpvar' must be output variables of the row they are in.";
                logval = false;
                return
            end
            continue
        elseif ~isa(op_ij,'tensopvar') && ~isa(op_ij,'intop')
            errmsg = errmsg + "Blocks in 'tensopmat' structure must be of type 'tensopvar', 'intop', 'double', 'polynomial' or 'dpvar'.";
            logval = false;
            return
        elseif isa(op_ij,'intop') && ~op_ij.chkval
            [logval1,errmsg1] = isvalid(op_ij);
            if ~logval1
                errmsg_ij = ['Operator P.ops{',num2str(i),',',num2str(j),'} is not properly declared:'];
                errmsg = errmsg + errmsg_ij + errmsg1;
                logval = false;
                return
            end
        end

        % Check that the dimensions match
        if dims.in(j)~=size(op_ij,2)
            errmsg = errmsg + "Column dimensions of operators in same column-block must match.";
            logval = false;                                                
            return                                                         
        end
        if dims.out(i)~=size(op_ij,1)
            errmsg = errmsg + "Row dimensions of operators in same row-block must match.";  
            logval = false;                                                
            return 
        end

        % Check that the input variables match
        vars_ij = op_ij.vars;      deps_ij = op_ij.deps;
        [~,idcs1_in,idcs2_in] = unique([vars.in,vars_ij.in],'stable');
        if numel(idcs1_in)>N
            errmsg = errmsg + "Operators in the same column should be defined in terms of the same input variables.";
            logval = false;
            return
        end
        idcs2_in = idcs2_in(N+1:end);
        if ~isequal(depmat.in(j,idcs2_in),deps_ij.in)
            errmsg = errmsg + "Operators in the same column should depend to the same degree on the same input variables";
            logval = false;
            return
        end
        % Check that the output variables match
        [~,idcs1_out,idcs2_out] = unique([vars.out,vars_ij.out],'stable');
        if numel(idcs1_out)>M
            errmsg = errmsg + "Operators in the same row should be defined in terms of the same output variables.";
            logval = false;
            return
        end
        idcs2_out = idcs2_out(M+1:end);
        if ~isequal(depmat.out(i,idcs2_out),deps_ij.out)
            errmsg = errmsg + "Operators in the same row should depend to the same degree on the same output variables";
            logval = false;
            return
        end

        % Check that the domains match
        dom_ij = op_ij.dom;
        if ~isequal(dom.in(idcs2_in,:),dom_ij.in)
            errmsg = errmsg + "Identical input variables should have identical domains.";
            logval = false;
            return
        end
        if ~isequal(dom.out(idcs2_out,:),dom_ij.out)
            errmsg = errmsg + "Identical output variables should have identical domains.";
            logval = false;
            return
        end
    end
end

if logval
    errmsg = "Tensopmat object is properly specified";
end

end
