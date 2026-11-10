function [logval,errmsg] = isvalid(Kop)
% [LOGVAL,ERRMSG] = ISVALID(KOP) checks whether the 'intop' object KOP is
% properly specified, i.e. whether the fields are all properly specified
% and of appropriate dimensions.
%
% INPUTS
% - Kop:    'intop' object of which to check whether the different fields
%           are of appropriate type and consistent dimensions;
%
% OUTPUTS
% - logval: boolean value specifying whether (true) or not (false) the
%           object is properly specified;
% - errmsg: 'string' used to declare an error message indicating what is
%           wrong with the object if logval = false;
%
% NOTES
% - The names Kop.Zd of the decision variables are not checked for
%   uniqueness or type.
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
% DJ, 10/22/2026: Check that the intervals of integration are ordered;

% Initialize
logval = true;  % innocent until proven guilty
errmsg = "Intop object is not properly specified: ";
Kop.chkval = false;

% Check that all fields are of appropriate type
if ~isa(Kop.params,'cell')
    errmsg = errmsg + "Field 'params' must be specified as a q1 x ... x qN cell.";
    logval = false;
elseif ~isa(Kop.intmat,'cell')
    errmsg = errmsg + "Field 'intmat' must be specified as a 1 x N cell, for N variables.";
    logval = false;
elseif ~isstruct(Kop.vars) || ~isfield(Kop.vars,'in') || ~isfield(Kop.vars,'out')
    errmsg = errmsg + "Field 'vars' must be specified as a struct with fields 'in' and 'out'.";
    logval = false;
elseif ~isstruct(Kop.dom) || ~isfield(Kop.dom,'in') || ~isfield(Kop.dom,'out')
    errmsg = errmsg + "Field 'dom' must be specified as a struct with fields 'in' and 'out'.";
    logval = false;
elseif ~isstruct(Kop.Z) || ~isfield(Kop.Z,'in') || ~isfield(Kop.Z,'out')
    errmsg = errmsg + "Field 'Z' must be specified as a struct with fields 'in' and 'out'.";
    logval = false;
elseif ~isa(Kop.Zd,'cell') || size(Kop.Zd,2)~=1 && ~isempty(Kop.Zd)
    errmsg = errmsg + "Field 'Zd' must be specified as an nd x 1 cellstr.";
    logval = false;
end
if ~logval
    return
end

% Check that the variables are properly specified
vars = Kop.vars;
if ~iscellstr(vars.in) || ~iscellstr(vars.out)
    errmsg = errmsg + "Input and output variables must be specified as cellstr objects.";
    logval = false;
    return
end
N = numel(vars.in);     M = numel(vars.out);
if M~=0
    errmsg = errmsg + "Functional should not have output variables.";
    logval = false;
elseif ~all(size(vars.in)==[1,N]) || ~all(size(vars.out)==[1,M])
    errmsg = errmsg + "Input and output variables must be specified as row arrays.";
    logval = false;
elseif ~isequal(vars.in,unique(vars.in)) || ~isequal(vars.out,unique(vars.out))
    errmsg = errmsg + "Input and output variable names must be unique and sorted.";
    logval = false;
end
if ~logval
    return
end

% Check that the number of variables is consistent
if ndims(Kop.params)>N && any(size(Kop.params,N+1:ndims(Kop.params))~=1)
    errmsg = errmsg + "Number of dimensions of kernel function cannot exceed number of variables.";
    logval = false;
elseif ~all(size(Kop.intmat)==[1,N])
    errmsg = errmsg + "Field 'intmat' must be specified as a 1 x N cell, for N variables.";
    logval = false;
end
if ~logval
    return
end

% Check the monomials
Z = Kop.Z;
if ~isa(Z.in,'cell') || ~all(size(Z.in)==[1,N]) || ~isa(Z.out,'cell') || ~isempty(Z.out)
    errmsg = errmsg + "Z.in must be specified as a 1 x N cell, and Z.out as an empty cell.";
    logval = false;
    return
end
for i=1:N
    Zi = Z.in{i};
    di = size(Kop.intmat{i},2)/2;
    if round(di)~=di
        errmsg = errmsg + "Each element of 'intmat' should be a qi x 2*di array, where di is the degree of the ith variable.";
        logval = false;
        return
    elseif ~isnumeric(Zi) || ~ismatrix(Zi) || size(Zi,2)~=di || size(Zi,1)<1
        errmsg = errmsg + "Z.in{i} must be a nmons x deps(i) array, with at least one monomial.";
        logval = false;
        return
    elseif any(Zi(:)<0) || any(round(Zi(:))~=Zi(:))
        errmsg = errmsg + "Monomial degrees in Z.in must be nonnegative integers.";
        logval = false;
        return
    elseif size(unique(Zi,'rows'),1)~=size(Zi,1)
        errmsg = errmsg + "Monomials in Z.in{i} must be distinct.";
        logval = false;
        return
    end
end

% Check that the parameters are sparse matrices of the appropriate size
nd = numel(Kop.Zd);
[~,~,nZ] = Zdims(Z.in);
sz1 = [];
for k=1:numel(Kop.params)
    param_k = Kop.params{k};
    if ~isa(param_k,'double') || ~issparse(param_k) || ~ismatrix(param_k)
        errmsg = errmsg + "Parameters must be specified as sparse matrices of type 'double'.";
        logval = false;
        return
    end
    if k==1
        sz1 = size(param_k);
    elseif ~isequal(size(param_k),sz1)
        errmsg = errmsg + "Parameters must all be of the same size.";
        logval = false;
        return
    end
end
if ~isempty(sz1) && (mod(sz1(1),nd+1)~=0 || mod(sz1(2),nZ)~=0)
    errmsg = errmsg + "Parameters must have m*(nd+1) rows and n*nZ columns, for nd decision variables, and nZ monomials.";
    logval = false;
    return
end

% Check that the domains of integration are properly specified
nterms = size(Kop.params,1:N);
for i=1:N
    intmat_i = Kop.intmat{i};
    if all(size(intmat_i)==[1,0])
        continue
    end
    % Check the number of terms
    if size(intmat_i,1)~=nterms(i)
        errmsg = errmsg + "Number of terms specified by 'intmat' should match number of kernels specified by 'params'.";
        logval = false;
        return
    end
    % The degree of the variable is an integer, as checked above
    di = size(intmat_i,2)/2;
    % Check the specified intervals are suitable
    if any(any(~ismember(intmat_i,(0:di+1))))
        errmsg = errmsg + "Elements of 'intmat{i}' must be integers ranging from 0 to di+1, where di is the degree of the ith variable.";
        logval = false;
        return
    end
    % Check that at least one variable is integrated over the full domain
    for k=1:size(intmat_i,1)
        lwr_idcs = (find(intmat_i(k,:)==0)+1)./2;
        upr_idcs = find(intmat_i(k,:)==di+1)./2;
        if any(lwr_idcs~=round(lwr_idcs))
            errmsg = errmsg + "Lower boundary cannot be upper limit of domain of integration.";
            logval = false;
            return
        elseif any(upr_idcs~=round(upr_idcs))
            errmsg = errmsg + "Upper boundary cannot be lower limit of domain of integration.";
            logval = false;
            return
        elseif isempty(intersect(lwr_idcs,upr_idcs)) && ~isempty(intmat_i)
            errmsg = errmsg + "At least one of the variables must be integrated over the full domain.";
            logval = false;
            return
        end
    end
    % Check that the intervals do not impose inconsistent ordering
    if ~is_consistent(intmat_i)
        errmsg = errmsg + "The intervals of integration are incompatible.";
        logval = false;
        return
    end
    % Check that the intervals are ordered, see 'intop.canonicalize_intmat'
    [~,~,is_ord] = intop.canonicalize_intmat(intmat_i,[],true);
    if ~all(is_ord)
        errmsg = errmsg + "Intervals of integration are not ordered. ";
        errmsg = errmsg + "Use intop.canonicalize_intmat to fix.";
        logval = false;
        return
    end
end

% Check that the variable domains are appropriate
dom = Kop.dom;
if ~all(size(dom.in)==[N,2]) || ~all(size(dom.out)==[M,2])
    errmsg = errmsg + "Field 'dom.in' must be specified as a N x 2 array of type 'double', for N variables.";
    logval = false;
elseif any(dom.in(:,1)>=dom.in(:,2))
    errmsg = errmsg + "Lower boundary of domains cannot exceed upper boundary.";
    logval = false;
end


end



%%
function logval = is_consistent(intmat)
% LOGVAL = IS_CONSISTENT(INTMAT_I) checks whether any of the domains of
% integration specified by INTMAT_I is unsuitable, in particular whether
% the domain of integration of some variable involves the variable itself.

logval = true;

% Deal with empty case
if isempty(intmat)
    return
end

% Loop over all specified integrals
[nterms,d] = size(intmat);
d = d/2;
for j=1:nterms
    var_list0 = mat2cell(intmat(j,:),1,2*ones(1,d));
    safe_vars = [0,d+1];
    rem_vars = 1:d;
    var_list = var_list0;
    unsafe_idcs = cellfun(@(a,b) ismember(a,b),num2cell(rem_vars),var_list);
    if any(unsafe_idcs)
        logval = false;
        break
    end
    for k=1:d
        var_list = cellfun(@(a) setdiff(a,safe_vars),var_list,'UniformOutput',false);
        safe_idcs = cellfun(@(a) isempty(a),var_list);
        if all(safe_idcs)
            break
        end
        safe_vars = [safe_vars,find(safe_idcs)];
        rem_vars = rem_vars(~safe_idcs);
        var_list = var_list(~safe_idcs);
        var_list = cellfun(@(a) cell2mat(var_list0(a)),var_list,'UniformOutput',false);
        unsafe_idcs = cellfun(@(a,b) ismember(a,b),num2cell(rem_vars),var_list);
        if any(unsafe_idcs)
            break
        end
    end
    if ~logval || ~all(safe_idcs)
        logval = false;
        break
    end
end

end
