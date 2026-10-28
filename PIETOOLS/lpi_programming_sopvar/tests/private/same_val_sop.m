function [tf,where] = same_val_sop(a,b,where)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TF,WHERE] = SAME_VAL_SOP(A,B,WHERE) tests whether two values are equal,
% recursing through struct arrays and cells: a struct needs the same field
% names, in the same order, and the same size; a cell needs the same size;
% a 'polynomial' object is compared by its coefficient, degmat, varname
% and matdim, a 'dpvar' object by its C, degmat, varname, dvarname and
% matdim, and any other value by isequaln (NaN equal to NaN) and its
% sparsity. WHERE locates the first difference.
%
% INPUT
% - a, b:   values of any class;
% - where:  (optional) char, the name of A, the prefix of the output WHERE
%           (default 'value').
%
% OUTPUT
% - tf:     logical, true when A and B have the same class and are equal as
%           above;
% - where:  char, the input WHERE when TF is true; otherwise the input
%           WHERE followed by the path to the first difference ('.f' for
%           struct field f, '{i}' for cell entry i) and the tag ' (class)',
%           ' (fields)' (field names or struct size), ' (size)' (cell size)
%           or ' (value)'.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - same_val_sop
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
% Initial coding MMP, 09/29/2026. Shared by the tests of this folder.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body; code unchanged. The deleted header text was
%                the reason the isequal overloads are not used, now a step
%                comment.

% Default prefix; the values are equal until a difference is found
if nargin<3,    where = 'value';    end
tf = true;
% Different classes are a difference
if ~strcmp(class(a),class(b)),  tf = false;     where = [where ' (class)'];    return,     end
if isstruct(a)
    % Struct: the same field names, in order, and the same size; then every
    % field of every element
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
    % Cell: the same size, then every entry
    if ~isequal(size(a),size(b)),   tf = false;     where = [where ' (size)'];     return,     end
    for i = 1:numel(a)
        [tf,w] = same_val_sop(a{i},b{i},sprintf('%s{%d}',where,i));
        if ~tf,     where = w;  return,     end
    end
elseif isa(a,'polynomial')
    % 'polynomial' and 'dpvar': the stored arrays, since their isequal
    % overloads are elementwise (isequal on two programs would not return
    % one logical)
    tf = isequal(a.coefficient,b.coefficient) && isequal(a.degmat,b.degmat) && ...
         isequal(a.varname,b.varname) && isequal(a.matdim,b.matdim);
elseif isa(a,'dpvar')
    tf = isequal(a.C,b.C) && isequal(a.degmat,b.degmat) && isequal(a.varname,b.varname) && ...
         isequal(a.dvarname,b.dvarname) && isequal(a.matdim,b.matdim);
else
    % Any other value: isequaln (NaN equal to NaN) and the same sparsity
    tf = isequaln(a,b) && isequal(issparse(a),issparse(b));
end
% Tag a difference found at this level
if ~tf && ~contains(where,'(')
    where = [where ' (value)'];
end
end
