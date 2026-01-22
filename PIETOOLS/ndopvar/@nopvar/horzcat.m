function [Pcat] = horzcat(varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = horzcat(varargin) takes n inputs and concatentates them horizontally,
% provided they satisfy the following criteria:
% 1) All inputs must be of type 'nopvar'
% 2) The output dimensions varargin{j}.dim(:,1) of all opvar objects must
%       match;
% 3) The spatial variables varargin{j}.var1 and varargin{j}.var2, as well 
%       as the domain varargin{j}.I of all opvar objects must match;
% 4) Concatenation should make sense within the context of the opvar class,
%       that is, opvars take inputs from RxL2, never from L2xR. We cannot
%       concatenate e.g. [A,B] for A:L2-->R and B:R-->R.
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - horzcat
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
if ~isa(b,'nopvar') || ~isa(a, 'nopvar')
    error("Currently supported only for nopvar");
    % Only supported if a in fact corresponds to a matrix as well
end


% Check that domain and variables match
if any(any(a.dom~=b.dom))|| any(~strcmp(a.vars(:, 1).varname,b.vars(:, 1).varname)) || any(~strcmp(a.vars(:, 2).varname,b.vars(:, 2).varname))
    error('Operators being concatenated have different intervals or different independent variables');
end

% convert to the same degree 
if any(a.deg~=b.deg)
    max_degree = max(a.deg, b.deg);
    Aop_new = change_degree(a, max_degree);
    Bop_new = change_degree(b, max_degree);
    a = Aop_new;
    b = Bop_new;
    % error('Operators must have the same degree.');
end


% Check that the output dimensions match
a.dim = a.dim;  b.dim = b.dim;  % correction to make components have consistent dimensions 8/27-ss
if any(a.dim(1)~=b.dim(1))
    error('Cannot concatenate horizontally: Output dimensions of opvar objects do not match')
end

% Finally, let's actually concatenate
Pcat = nopvar(); % empty nopvar
Pcat.dom = a.dom;
Pcat.deg = a.deg;
Pcat.vars = a.vars;
N = size(a.dom,1);
Pcat.C = cell([3*ones(1,N),1]);
for ii=1:numel(a.C)

    Pcat.C{ii} = [a.C{ii}, b.C{ii}];
end

% Pcat.dim = [a.dim(1), a.dim(2)+b.dim(2)];




% For concatenation of more than two objects, just repeat
if nargin>2 
    Pcat = horzcat(Pcat, varargin{3:end});
end

end
 