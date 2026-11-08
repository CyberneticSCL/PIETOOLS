function Fdiag = blkdiag(varargin)
% F = blkdiag(varargin) generates a block diagonal dopvar object using
% varargs as the blocks on the diagonal
% 
% INPUTS
% - varargin:   'opvar' or 'dopvar' class objects.
% 
% OUTPUTS
% - Fdiag:      'dopvar' object with each of the parameters definined as the
%               block diagonal concatenation of the input operators, so 
%               that e.g. if Cop=blkdiag(Aop,Bop), then
%               Cop.P = [Aop.P, 0; 0, Bop.P], 
%               Cop.R.R0 = [Aop.R.R0, 0; 0, Aop.R.R0]
%               Note that this structure may not necessarily represent the 
%               the block-diagonal concatenation of the input operators;
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2025  PIETOOLS Team
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
% DJ, 09/05/2025: Initial coding;

if nargin==1
    Fdiag=varargin{1};
else % sequentially creating a block diagonal matrix
    A = varargin{1}; B = varargin{2};
    
    if isempty(A)
        Fdiag = B;
    elseif isempty(B)
        Fdiag = A;
    else
        if (~isa(A,'opvar') && ~isa(A,'dopvar')) || (~isa(B,'opvar') && ~isa(B,'dopvar'))
            error("Block-diagonal concatenation of dopvar and non-opvar objects is currently not supported.")
        end
        if any(any(A.I~=B.I))
            error('Spatial domain of dopvar objects must match for block-diagonal concatenation');
        elseif any(~isequal(A.var1,B.var1)) || any(~isequal(A.var2,B.var2))
            error('The spatial variables of dopvar objects must match for block-diagonal concatenation');
        end
        O1 = dopvar();
        O1.dim = [A.dim(:,1),B.dim(:,2)];
        O1.var1 = A.var1;   O1.var2 = A.var2;
        O1.I = A.I;
        O2 = dopvar();
        O2.dim = [B.dim(:,1),A.dim(:,2)];
        O2.var1 = A.var1;   O2.var2 = A.var2;
        O2.I = A.I;
        
        Fdiag = [A, O1; O2, B];
    end
    
    if nargin>2 % repeat when there are more than two block elements
        Fdiag = blkdiag(Fdiag,varargin{3:end});
    end
end
end