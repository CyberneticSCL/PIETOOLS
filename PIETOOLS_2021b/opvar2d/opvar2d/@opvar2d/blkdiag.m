function Fdiag = blkdiag(varargin)
% F = blkdiag(varargin) generates a block diagonal opvar2d object using
% varargs as the diag blocks
% 
% INPUTS:
%   varargs: opvar2d, poly, or double object. At least one object should
%            be opvar2d
% 
% OUTPUTS:
% Fdiag: an opvar2d object
%           - If two objects A,B are opvar2d, blkdiag(A,B) will perform
%             concatenation [A.Rij, 0    ]
%                           [0    , B.Rij] on the different parameters Rij
%             defining A and B
%           - If A is opvar2d but B is not, A must have only 1 nonempty
%             element A.Rij, which will be adjusted to [A.Rij, 0]
%                                                      [0    , B];
%           - If B is opvar2d but A is not, B must have only 1 nonempty
%             element B.Rij, which will be adjusted to [A, 0   ]
%                                                      [0, B.Rij];
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 08/27/2021
% 06/24/2022, DJ - Fix for blkdiag with non-opvar2d object

if nargin==1
    Fdiag=varargin{1};
else % sequentially creating a block diagonal matrix
    A = varargin{1}; B = varargin{2};
    
    if isempty(A)
        Fdiag = B;
    elseif isempty(B)
        Fdiag = A;
    else
        if isa(A,'opvar2d') && isa(B,'opvar2d')
            if any(any(A.I~=B.I))
                error('Spatial domain of opvar2d objects must match for block-diagonal concatenation');
            elseif any(~isequal(A.var1,B.var1)) || any(~isequal(A.var2,B.var2))
                error('The spatial variables of opvar2d objects must match for block-diagonal concatenation');
            end
            O1 = opvar2d([],[A.dim(:,1),B.dim(:,2)],A.I,A.var1,A.var2);
            O2 = opvar2d([],[B.dim(:,1),A.dim(:,2)],A.I,A.var1,A.var2);
            
            Fdiag = [A, O1; O2, B];
            
        elseif isa(A,'opvar2d') && (isa(B,'polynomial') || isa(B,'double'))
            if ~(nnz(A.dim(:,1))==1 && nnz(A.dim(:,2))==1)
                error('Block-diagonal concatenation of opvar2d with polynomial or double is only supported if opvar2d maps only from and to a single space')
            end
            if isa(B,'polynomial') && (any(~ismember(B.varname,A.var1.varname)) || any(~ismember(B.varname,A.var2.varname)))
                error('Variables appearing in the objects to concatenate should match')
            end
            % If B is not opvar2d, take block diagonal of A.Rij with B for the
            % (only) nonempty parameter A.Rij of A
            opvar2d Fdiag;
            Fdiag.I = A.I;
            Fdiag.var1 = A.var1;    Fdiag.var2 = A.var2;
            if A.dim(1,1)~=0
                if A.dim(1,2)~=0
                    Fdiag.R00 = blkdiag(A.R00,B);
                elseif A.dim(2,2)~=0
                    Fdiag.R0x = blkdiag(A.R0x,B);
                elseif A.dim(3,2)~=0
                    Fdiag.R0y = blkdiag(A.R0y,B);
                elseif A.dim(4,2)~=0
                    Fdiag.R02 = blkdiag(A.R02,B);
                end
            elseif A.dim(2,1)~=0
                if A.dim(1,2)~=0
                    Fdiag.Rx0 = blkdiag(A.Rx0,B);
                elseif A.dim(2,2)~=0
                    Rxx = cell(3,1);
                    Rxx{1} = blkdiag(A.Rxx{1},B);
                    Rxx{2} = blkdiag(A.Rxx{2},B);
                    Rxx{3} = blkdiag(A.Rxx{3},B);
                    Fdiag.Rxx = Rxx;
                elseif A.dim(3,2)~=0
                    Fdiag.Rxy = blkdiag(A.Rxy,B);
                elseif A.dim(4,2)~=0
                    Rx2 = cell(3,1);
                    Rx2{1} = blkdiag(A.Rx2{1},B);
                    Rx2{2} = blkdiag(A.Rx2{2},B);
                    Rx2{3} = blkdiag(A.Rx2{3},B);
                    Fdiag.Rx2 = Rx2;
                end
            elseif A.dim(3,1)~=0
                if A.dim(1,2)~=0
                    Fdiag.Ry0 = blkdiag(A.Ry0,B);
                elseif A.dim(2,2)~=0
                    Fdiag.Ryx = blkdiag(A.Ryx,B);
                elseif A.dim(3,2)~=0
                    Ryy = cell(1,3);
                    Ryy{1} = blkdiag(A.Ryy{1},B);
                    Ryy{2} = blkdiag(A.Ryy{2},B);
                    Ryy{3} = blkdiag(A.Ryy{3},B);
                    Fdiag.Ryy = Ryy;
                elseif A.dim(4,2)~=0
                    Ry2 = cell(1,3);
                    Ry2{1} = blkdiag(A.Ry2{1},B);
                    Ry2{2} = blkdiag(A.Ry2{2},B);
                    Ry2{3} = blkdiag(A.Ry2{3},B);
                    Fdiag.Ry2 = Ry2;
                end
            elseif A.dim(4,1)~=0
                if A.dim(1,2)~=0
                    Fdiag.R20 = blkdiag(A.R20,B);
                elseif A.dim(2,2)~=0
                    R2x = cell(3,1);
                    R2x{1} = blkdiag(A.R2x{1},B);
                    R2x{2} = blkdiag(A.R2x{2},B);
                    R2x{3} = blkdiag(A.R2x{3},B);
                    Fdiag.R2x = R2x;
                elseif A.dim(3,2)~=0
                    R2y = cell(1,3);
                    R2y{1} = blkdiag(A.R2y{1},B);
                    R2y{2} = blkdiag(A.R2y{2},B);
                    R2y{3} = blkdiag(A.R2y{3},B);
                    Fdiag.R2y = R2y;
                elseif A.dim(4,2)~=0
                    R22 = cell(3,3);
                    for indx=1:numel(R22)
                        R22{indx} = blkdiag(A.R22{indx},B);
                    end
                    Fdiag.R22 = R22;
                end
            end
            Fdiag.dim = Fdiag.dim;
            
            elseif isa(B,'opvar2d')% && (isa(A,'polynomial') || isa(A,'double'))
            if ~(nnz(B.dim(:,1))==1 && nnz(B.dim(:,2))==1)
                error('Block-diagonal concatenation of opvar2d with polynomial or double is only supported if opvar2d maps only from and to a single space')
            end
            if isa(A,'polynomial') && (any(~ismember(A.varname,B.var1.varname)) || any(~ismember(A.varname,B.var2.varname)))
                error('Variables appearing in the objects to concatenate should match')
            end
            % If A is not opvar2d, take block diagonal of B.Rij with A for the
            % (only) nonempty parameter B.Rij of B
            opvar2d Fdiag;
            Fdiag.I = B.I;
            Fdiag.var1 = B.var1;    Fdiag.var2 = B.var2;
            if B.dim(1,1)~=0
                if B.dim(1,2)~=0
                    Fdiag.R00 = blkdiag(A,B.R00);
                elseif B.dim(2,2)~=0
                    Fdiag.R0x = blkdiag(A,B.R0x);
                elseif B.dim(3,2)~=0
                    Fdiag.R0y = blkdiag(A,B.R0y);
                elseif B.dim(4,2)~=0
                    Fdiag.R02 = blkdiag(A,B.R02);
                end
            elseif B.dim(2,1)~=0
                if B.dim(1,2)~=0
                    Fdiag.Rx0 = blkdiag(A,B.Rx0);
                elseif B.dim(2,2)~=0
                    Rxx = cell(3,1);
                    Rxx{1} = blkdiag(A,B.Rxx{1});
                    Rxx{2} = blkdiag(A,B.Rxx{2});
                    Rxx{3} = blkdiag(A,B.Rxx{3});
                    Fdiag.Rxx = Rxx;
                elseif B.dim(3,2)~=0
                    Fdiag.Rxy = blkdiag(A,B.Rxy);
                elseif B.dim(4,2)~=0
                    Rx2 = cell(3,1);
                    Rx2{1} = blkdiag(A,B.Rx2{1});
                    Rx2{2} = blkdiag(A,B.Rx2{2});
                    Rx2{3} = blkdiag(A,B.Rx2{3});
                    Fdiag.Rx2 = Rx2;
                end
            elseif B.dim(3,1)~=0
                if B.dim(1,2)~=0
                    Fdiag.Ry0 = blkdiag(A,B.Ry0);
                elseif B.dim(2,2)~=0
                    Fdiag.Ryx = blkdiag(A,B.Ryx);
                elseif B.dim(3,2)~=0
                    Ryy = cell(1,3);
                    Ryy{1} = blkdiag(A,B.Ryy{1});
                    Ryy{2} = blkdiag(A,B.Ryy{2});
                    Ryy{3} = blkdiag(A,B.Ryy{3});
                    Fdiag.Ryy = Ryy;
                elseif B.dim(4,2)~=0
                    Ry2 = cell(1,3);
                    Ry2{1} = blkdiag(A,B.Ry2{1});
                    Ry2{2} = blkdiag(A,B.Ry2{2});
                    Ry2{3} = blkdiag(A,B.Ry2{3});
                    Fdiag.Ry2 = Ry2;
                end
            elseif B.dim(4,1)~=0
                if B.dim(1,2)~=0
                    Fdiag.R20 = blkdiag(A,B.R20);
                elseif B.dim(2,2)~=0
                    R2x = cell(3,1);
                    R2x{1} = blkdiag(A,B.R2x{1});
                    R2x{2} = blkdiag(A,B.R2x{2});
                    R2x{3} = blkdiag(A,B.R2x{3});
                    Fdiag.R2x = R2x;
                elseif B.dim(3,2)~=0
                    R2y = cell(1,3);
                    R2y{1} = blkdiag(A,B.R2y{1});
                    R2y{2} = blkdiag(A,B.R2y{2});
                    R2y{3} = blkdiag(A,B.R2y{3});
                    Fdiag.R2y = R2y;
                elseif B.dim(4,2)~=0
                    R22 = cell(3,3);
                    for indx=1:numel(R22)
                        R22{indx} = blkdiag(A,B.R22{indx});
                    end
                    Fdiag.R22 = R22;
                end
            end
            Fdiag.dim = Fdiag.dim;
        end
    end
    
    if nargin>2 % repeat when there are more than two block elements
        Fdiag = blkdiag(Fdiag,varargin{3:end});
    end
end
end