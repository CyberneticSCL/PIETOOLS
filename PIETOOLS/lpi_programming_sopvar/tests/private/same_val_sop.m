function [tf,where] = same_val_sop(a,b,where)
% [TF,WHERE] = SAME_VAL_SOP(A,B,WHERE) bitwise equality of two values,
% recursing through structs and cells; 'polynomial' and 'dpvar' objects are
% compared by their stored arrays (their isequal overloads are elementwise,
% so isequal on two programs would not return one logical). WHERE names
% the first difference found, starting from the given prefix.
%
% Initial coding MMP, 09/29/2026. Shared by the tests of this folder.
if nargin<3,    where = 'value';    end
tf = true;
if ~strcmp(class(a),class(b)),  tf = false;     where = [where ' (class)'];    return,     end
if isstruct(a)
    fa = fieldnames(a);     fb = fieldnames(b);
    if ~isequal(fa,fb) || ~isequal(size(a),size(b))
        tf = false;     where = [where ' (fields)'];   return
    end
    for i = 1:numel(a)
        for k = 1:numel(fa)
            [tf,w] = same_val_sop(a(i).(fa{k}),b(i).(fa{k}),[where '.' fa{k}]);
            if ~tf,     where = w;  return,     end
        end
    end
elseif iscell(a)
    if ~isequal(size(a),size(b)),   tf = false;     where = [where ' (size)'];     return,     end
    for i = 1:numel(a)
        [tf,w] = same_val_sop(a{i},b{i},sprintf('%s{%d}',where,i));
        if ~tf,     where = w;  return,     end
    end
elseif isa(a,'polynomial')
    tf = isequal(a.coefficient,b.coefficient) && isequal(a.degmat,b.degmat) && ...
         isequal(a.varname,b.varname) && isequal(a.matdim,b.matdim);
elseif isa(a,'dpvar')
    tf = isequal(a.C,b.C) && isequal(a.degmat,b.degmat) && isequal(a.varname,b.varname) && ...
         isequal(a.dvarname,b.dvarname) && isequal(a.matdim,b.matdim);
else
    tf = isequaln(a,b) && isequal(issparse(a),issparse(b));
end
if ~tf && ~contains(where,'(')
    where = [where ' (value)'];
end
end
