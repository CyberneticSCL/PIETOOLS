function [Pcat] = vertcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = vertcat(varargin) takes n inputs and concatentates them vertically,
% provided they satisfy the following criteria:
% 1) All input must be of type 'sopvar' 
% 2) The input dimensions varargin{j}.dim(:,2) of all opvar objects must
%       match;
% 3) The spatial variables varargin{j}.var1 and varargin{j}.var2, as well 
%       as the domain varargin{j}.I of all opvar objects must match;
% 4) Concatenation should make sense within the context of the opvar class,
%       that is, opvars always mapt to from RxL2, never to L2xR. We cannot
%       concatenate e.g. [A;B] for A:R-->L2 and B:R-->R.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - vertcat
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
% Initial AT 01/21/2026

% Deal with single input case
if nargin==1
    Pcat = varargin{1};
    return
end

% Extract the operators
a = varargin{1};    b = varargin{2};

% Currently support some matrix-opvar concatenation
if ~isa(b,'sopvar') || ~isa(a, 'sopvar')
    error("Currently supported only for ndopvar or sopvar");
    % Only supported if a in fact corresponds to a matrix as well
end

% Check that domain and variables match
if any(any(a.dom_1~=b.dom_1))|| any(any(a.dom_2~=b.dom_2)) || any(any(a.dom_3~=b.dom_3))
    error('Operators being concatenated have different intervals');
end
 

% Check that the output dimensions match
if any(a.dims(2)~=b.dims(2))
    error('Cannot concatenate vertically: Input dimensions of sopvar objects do not match')
end  

% Finally, let's actually concatenate
params = a.params;
for ii=1:numel(a.params)
    params{ii} = [a.params{ii}; b.params{ii}]; % concatenation of quadpoly
end

% Pcat.dim = [a.dim(1), a.dim(2)+b.dim(2)];


% it is assumed that vars_S1, S2, S3 are sorted arrays.
vars_S1 = mergeSortedCellstr(a.vars_S1, b.vars_S1); 
vars_S2 = mergeSortedCellstr(a.vars_S2, b.vars_S2);
vars_S3 = mergeSortedCellstr(a.vars_S3, b.vars_S3);
dom_3 = a.dom_3;
dom_2 = a.dom_2;
dom_1 = a.dom_1;
dims = a.dims;
dims(1) = dims(1) + b.dims(1);

Pcat = sopvar(vars_S3,dom_3,dims,params,vars_S1,dom_1,vars_S2,dom_2);


if nargin>2 
    Pcat = horzcat(Pcat, varargin{3:end});
end


end


% ---------------- local helpers ----------------

function u = mergeSortedCellstr(a, b)
% Sorted union of two sorted-unique cellstr arrays.
i = 1; j = 1;
na = numel(a); nb = numel(b);
u = cell(1, na+nb);
k = 0;

while i <= na && j <= nb
    ai = a{i}; bj = b{j};
    if strcmp(ai, bj)
        k=k+1; u{k} = ai; i=i+1; j=j+1;
    elseif lexLessChar(ai, bj)
        k=k+1; u{k} = ai; i=i+1;
    else
        k=k+1; u{k} = bj; j=j+1;
    end
end
while i <= na, k=k+1; u{k} = a{i}; i=i+1; end
while j <= nb, k=k+1; u{k} = b{j}; j=j+1; end

u = u(1:k);
end

function tf = lexLessChar(s1, s2)
% Lexicographic compare without string allocations.
s1 = char(s1); s2 = char(s2);
L1 = numel(s1); L2 = numel(s2);
L  = min(L1, L2);

d = s1(1:L) - s2(1:L);
idx = find(d ~= 0, 1, 'first');
if isempty(idx)
    tf = (L1 < L2);
else
    tf = (d(idx) < 0);
end
end

 