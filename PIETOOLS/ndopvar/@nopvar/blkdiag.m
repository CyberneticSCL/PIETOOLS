function Fdiag = blkdiag(varargin)
% F = blkdiag(varargin) generates a block diagonal nopvar object using
% varargs as the blocks on the diagonal
% 
% INPUTS
% - varargin:   'nopvar' class objects.
% 
% OUTPUTS
% - Fdiag:      'nopvar' object with each of the parameters definined as the
%               block diagonal concatenation of the input operators,  
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  PIETOOLS Team
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
%W
% AT, 02/18/2026: Initial coding;

if nargin==1
    Fdiag=varargin{1};
else % sequentially creating a block diagonal matrix
    A = varargin{1}; B = varargin{2};
    
    if isempty(A)
        Fdiag = B;
    elseif isempty(B)
        Fdiag = A;
    else
    
    % Currently support some matrix-opvar concatenation
    if ~isa(B,'nopvar') || ~isa(A, 'nopvar')
        error("Currently supported only for ndopvar or nopvar");
        % Only supported if a in fact corresponds to a matrix as well
    end
    
    % Check that domain and variables match
    if any(any(A.dom~=B.dom))|| any(~strcmp(A.vars(:, 1).varname,B.vars(:, 1).varname)) || any(~strcmp(A.vars(:, 2).varname,B.vars(:, 2).varname))
        error('Operators have different intervals or different independent variables');
    end
    
    % convert to the same degree
    if any(A.deg~=B.deg)
        max_degree = max(A.deg, B.deg);
        Aop_new = change_degree(A, max_degree);
        Bop_new = change_degree(B, max_degree);
        A = Aop_new;
        B = Bop_new;
    end
    

    % Finally, let's actually concatenate
    Fdiag = nopvar(); % empty nopvar
    Fdiag.dom = A.dom;
    Fdiag.deg = A.deg;
    Fdiag.vars = A.vars;
    N = size(A.dom,1);
    Fdiag.C = cell([3*ones(1,N),1]);
    for ii=1:numel(A.C)
        Fdiag.C{ii} = blkdiag( A.C{ii}, B.C{ii} );
    end
    
    if nargin>2 % repeat when there are more than two block elements
        Fdiag = blkdiag(Fdiag,varargin{3:end});
    end
end
end