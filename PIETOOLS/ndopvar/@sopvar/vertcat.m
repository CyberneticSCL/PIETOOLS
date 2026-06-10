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
% Update to new sopvar AT 05/18/26
% Correct order ZL and ZR, split off UnionBasisMonomials, DJ 06/08/2026


% Deal with single input case
if nargin==1
    Pcat = varargin{1};
    return
end

% Extract the operators
a = varargin{1};    b = varargin{2};

% Currently support some matrix-opvar concatenation
if ~isa(b,'sopvar') || ~isa(a, 'sopvar')
    error("Currently supported only for sopvar");
    % Only supported if a in fact corresponds to a matrix as well
end

% Check that domain and variables match
if any(any(a.dom.in~=b.dom.in)) || any(any(a.dom.out~=b.dom.out))
    error('Operators being concatenated have different intervals');
end

% Check that domain and variables match
if any(~strcmp(a.vars.in, b.vars.in))|| any(~strcmp(a.vars.out, b.vars.out))
    error('Operators being concatenated have input/output variables');
end
 

% Check that the output dimensions match
if any(a.dims(2)~=b.dims(2))
    error('Cannot concatenate vertically: Input dimensions of sopvar objects do not match')
end  

ZL1 = a.ZL;
ZL2 = b.ZL;
ZR1 = a.ZR;
ZR2 = b.ZR;

[ZR, C1R, C2R] = UnionBasisMonomials(ZR1, ZR2);
[ZL, C1L, C2L] = UnionBasisMonomials(ZL1, ZL2);


C1L = kron(eye(a.dims(1)), C1L);
C2L = kron(eye(b.dims(1)), C2L);
C1R = kron(eye(a.dims(2)), C1R);
C2R = kron(eye(b.dims(2)), C2R);
% Finally, let's actually concatenate
params = a.params;
for ii=1:numel(a.params)
    params{ii} = [C1L'*a.params{ii}*C1R; C2L'*b.params{ii}*C2R]; % concatenation of quadpoly
end

% Pcat.dim = [a.dim(1), a.dim(2)+b.dim(2)];

dims = a.dims;
dims(1) = dims(1) + b.dims(1);

Pcat = sopvar(params, a.vars, ZR, ZL, a.dom, dims);


if nargin>2 
    Pcat = horzcat(Pcat, varargin{3:end});
end


end