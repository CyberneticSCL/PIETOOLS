function logval = isequal(PDE1,PDE2,varargin)
% LOGVAL = ISEQUAL(PDE1,PDE2,...) checks whether 'pde_struct' objects are
% equal, comparing all stored properties.
%
% INPUT
% - PDE1, PDE2, ...:    objects to compare, at least one of which is a
%                       'pde_struct' object;
%
% OUTPUT
% - logval:     logical scalar, true if all inputs are equal;
%
% NOTES
% 'polynomial' values are equal if their difference is zero. The built-in
% isequal cannot be used: it crashes MATLAB on objects that store
% 'polynomial' values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 11/24/2026: Initial coding;

logval = is_equal(PDE1,PDE2);
for k=1:numel(varargin)
    if ~logval,     return;     end
    logval = is_equal(PDE1,varargin{k});
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function logval = is_equal(A,B)
% LOGVAL = IS_EQUAL(A,B) checks whether A and B are equal, recursing into
% the stored properties of 'pde_struct' objects, and into structs and
% cells.
%
% INPUT
% - A, B:   values of any type;
%
% OUTPUT
% - logval: logical scalar, true if A and B are equal;

if isa(A,'pde_struct') || isa(B,'pde_struct')
    % Compare all properties that are stored, not the dependent ones, with
    % builtin indexing: x(i) and x.vars mean something else for a PDE
    % variable.
    logval = isa(A,'pde_struct') && isa(B,'pde_struct') && isequal(builtin('size',A),builtin('size',B));
    mc = metaclass(A);      props = mc.PropertyList;
    for ii=1:numel(A)*logval
        for jj=1:numel(props)
            if props(jj).Dependent,     continue;   end
            S = substruct('()',{ii},'.',props(jj).Name);
            logval = is_equal(builtin('subsref',A,S),builtin('subsref',B,S));
            if ~logval,     return;     end
        end
    end
elseif isa(A,'polynomial') || isa(B,'polynomial')
    % The isequal of 'polynomial' objects compares elementwise.
    logval = isequal(size(A),size(B));
    if logval
        D = polynomial(A)-polynomial(B);
        logval = nnz(D.coefficient)==0;
    end
elseif isstruct(A) || isstruct(B)
    % Fields may appear in any order.
    logval = isstruct(A) && isstruct(B) && isequal(size(A),size(B)) ...
                && isempty(setxor(fieldnames(A),fieldnames(B)));
    fnames = fieldnames(A);
    for ii=1:numel(A)*logval
        for jj=1:numel(fnames)
            logval = is_equal(A(ii).(fnames{jj}),B(ii).(fnames{jj}));
            if ~logval,     return;     end
        end
    end
elseif iscell(A) || iscell(B)
    logval = iscell(A) && iscell(B) && isequal(size(A),size(B));
    for ii=1:numel(A)*logval
        logval = is_equal(A{ii},B{ii});
        if ~logval,     return;     end
    end
else
    logval = isequal(A,B);
end

end
