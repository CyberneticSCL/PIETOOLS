function Fdiag = blkdiag(varargin)
% F = blkdiag(varargin) generates a block diagonal quadpoly objects using
% varargs as the blocks on the diagonal
% 
% INPUTS
% - varargin:   'quadpoly' class objects.
% 
% OUTPUTS
% - Fdiag:      'quadpoly' object with each of the parameters definined as the
%               block diagonal concatenation of the input quadpoly
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
%
% AT, 02/18/2026: Initial coding;

if nargin==1
    Fdiag=varargin{1};
else % sequentially creating a block diagonal matrix
    Q = cell(size(varargin));
    for k = 1:nargin
        if isa(varargin{k},'quadPoly')
            Q{k} = varargin{k};
        elseif isnumeric(X)
            A = sparse(varargin{k});
            Q{k} = quadPoly(A, {}, {}, size(A), {}, {});
        else
            error('quadPoly:blkdiag:badType', 'blkdiag supports quadPoly and numeric inputs.');
        end
    end


    % ---- Initialize accumulator with first block A obj----
    A = Q{1};
    nsU = A.ns;  ZsU = A.Zs;
    ntU = A.nt;  ZtU = A.Zt;
    
    Cblkdiag = sparse(A.C);
    mblkdiag = A.dim(1); % nrows
    nblkdiag = A.dim(2); % ncolumns

    for k = 2:numel(Q) % loop as in vertcat/horzcat
        B = Q{k};
        mB = B.dim(1);
        nB = B.dim(2);
        
        % --- merge s-side bases: current union with B ---
        [nsNew, ZsNew, ~, ~, mapS_blkd, mapS_B, sIsblkd, sIsB, dsNew] = unionBasis(nsU, ZsU, B.ns, B.Zs);
        
        % --- merge t-side bases: current union with B ---
        [ntNew, ZtNew, ~, ~, mapT_blkd, mapT_B, tIsblkd, tIsB, dtNew] = unionBasis(ntU, ZtU, B.nt, B.Zt);
        
        % If union changed, lift current concatenated coefficients to new union basis
        if ~(sIsblkd && tIsblkd)
            Cblkdiag = liftIndex(Cblkdiag, mblkdiag, nblkdiag, mapS_blkd, mapT_blkd, dsNew, dtNew);
            nsU = nsNew; ZsU = ZsNew; dsU = dsNew;
            ntU = ntNew; ZtU = ZtNew; dtU = dtNew;
        else
            nsU = nsNew; ZsU = ZsNew; dsU = dsNew;
            ntU = ntNew; ZtU = ZtNew; dtU = dtNew;
        end


        % Lift P into current union basis (if needed)
        if sIsB && tIsB
            CB= sparse(B.C);
        else
            CB = liftIndex(B.C, mB, nB, mapS_B, mapT_B, dsU, dtU);
        end

        % blkdiag of C matrices
        Cblkdiag = blkdiag(Cblkdiag, CB);
        mblkdiag = mblkdiag + mB;
        nblkdiag = nblkdiag + nB;

    end
    Fdiag = quadPoly(Cblkdiag, ZsU, ZtU, [mblkdiag, nblkdiag], nsU, ntU);
end
end