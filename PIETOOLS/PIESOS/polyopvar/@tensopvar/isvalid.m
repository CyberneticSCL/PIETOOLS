function [logval,errmsg] = isvalid(Cop)
% [LOGVAL,ERRMSG] = ISVALID(COP) checks whether the 'tensopvar' object COP
% is properly specified, i.e. whether the fields are all properly specified
% and of appropriate dimensions.
%
% INPUTS
% - Cop:    'tensopvar' object of which to check whether the different
%           fields are of appropriate type and consistent dimensions;
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
errmsg = "tensopvar object is not properly specified: ";
Cop.chkval = false;

% Check that all fields are of appropriate type
if ~isa(Cop.ops,'cell')
    errmsg = errmsg + "Field 'ops' must be specified as a m x d cell.";
    logval = false;
elseif ~isa(Cop.vars,'struct') || ~isfield(Cop.vars,'in') || ~isfield(Cop.vars,'out')
    errmsg = errmsg + "Field 'vars' must be specified as a struct with fields 'in' and 'out.";
    logval = false;
elseif ~isa(Cop.dom,'struct') || ~isfield(Cop.dom,'in') || ~isfield(Cop.dom,'out')
    errmsg = errmsg + "Field 'dom' must be specified as a struct with fields 'in' and 'out'.";
    logval = false;
elseif ~isa(Cop.deps,'struct') || ~isfield(Cop.deps,'in') || ~isfield(Cop.deps,'out')
    errmsg = errmsg + "Field 'deps' must be specified as a struct with fields 'in' and 'out'.";
    logval = false;
elseif ~isa(Cop.otype,'double') || ~all(size(Cop.otype)==[1,2])
    errmsg = errmsg + "Field 'otype' must be specified as a 1x2 array of type 'double'.";
    logval = false;
end
if ~logval
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
elseif  numel(unique(vars.in))~=N || numel(unique(vars.out))~=M
    errmsg = errmsg + "Input and output variable names must be unique.";
    logval = false;
end
if ~logval
    return
end

% Check that the domains are properly specified
dom = Cop.dom;
if ~isa(dom.in,'double') || ~isa(dom.out,'double')
    errmsg = errmsg + "Input and output domains must be specified as arrays of type 'double'.";
    logval = false;
elseif ~all(size(dom.in)==[N,2]) || ~all(size(dom.out)==[M,2])
    errmsg = errmsg + "Number of input and output domains must match number of input and output variables.";
    logval = false;
elseif any(dom.in(:,1)>=dom.in(:,2)) || any(dom.out(:,1)>=dom.out(:,2))
    errmsg = errmsg + "Lower boundary of domain cannot exceed upper boundary.";
    logval = false;
end
if ~logval
    return
end


% Check that the matrix dimensions match
dims = Cop.dims;
otype = Cop.otype;
if any(isnan(dims))
    errmsg = errmsg + "Matrix dimensions of the different terms must match.";
    logval = false;
elseif otype(1)==1 && max(dims(:,1))~=min(dims(:,1))
    errmsg = errmsg + "Row dimensions of the different factors must match for elementwise product.";
    logval = false;
elseif otype(2)==1 && max(dims(:,2))~=min(dims(:,2))
    errmsg = errmsg + "Column dimensions of the different factors must match for elementwise product.";
    logval = false;
end
if ~logval
    return
end

% Check that the different factors 
for k=1:numel(Cop.ops)
    op_k = Cop.ops{k};
    if ~isa(op_k,'sopvar') && ~isa(op_k,'tensopvar') && ~isa(op_k,'intop')
        errmsg = errmsg + "Factors in the tensor-PI operator must be specified as objects of type 'sopvar'.";
        logval = false;
        return
    end
    % Check that the variables and domains match
    vars_k = op_k.vars;     dom_k = op_k.dom;
    [shared_vars_in,ia_in,ib_in] = intersect(vars_k.in,vars.in,'stable');
    if numel(shared_vars_in)~=numel(vars_k.in)
        errmsg = errmsg + "Input variables of tensor-PI operator must include input variables of all factors.";
        logval = false;
        return
    elseif any(any(dom_k.in(ia_in,:)~=dom.in(ib_in,:)))
        errmsg = errmsg + "Domains of inputs variables of tensor-PI operator must match domains of input variables of all factors.";
        logval = false;
        return
    end
    [shared_vars_out,ia_out,ib_out] = intersect(vars_k.out,vars.out,'stable');
    if numel(shared_vars_out)~=numel(vars_k.out)
        errmsg = errmsg + "Output variables of tensor-PI operator must include output variables of all factors.";
        logval = false;
        return
    elseif any(any(dom_k.out(ia_out,:)~=dom.out(ib_out,:)))
        errmsg = errmsg + "Domains of output variables of tensor-PI operator must match domains of output variables of all factors.";
        logval = false;
        return
    end
end


end