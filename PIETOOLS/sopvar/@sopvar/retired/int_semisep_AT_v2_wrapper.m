% RETIRED 2026-08-31. Repacking wrapper over sopvar/int_semisep.m,
% itself a replacement for the duplicated implementation now in
% int_semisep_AT_v1.m. Superseded when int_semisep gained a 'packed'
% layout option, which builds the packed form directly instead of
% building the cell form and repacking it. @sopvar/mtimes_AT now calls
% int_semisep(...,'packed').

function [CMat,ZL,ZR] = int_semisep_AT_v2_wrapper(G,idxbeta,idxalpha,lims,Csize)
% int_semisep_AT
%
% Evaluates the two-sided semiseparable integral of Sec. 6.1 of the sopvar
% document, exactly as INT_SEMISEP does, and returns the result in the
% packed layout that '@sopvar/mtimes_AT' consumes.
%
% The two routines differ only in how the (beta,alpha) axes are stored:
%
%   int_semisep     C_gam_alp_beta{gamma,beta,alpha}, each block
%                   g1*NL by g2*NR
%
%   int_semisep_AT  CMat{gamma}, one matrix per gamma holding every
%                   (beta,alpha) block, of size g1*NL*nbeta by g2*NR*nalpha,
%                   with block (i,j) at rows (1:g1*NL)+g1*NL*(i-1) and
%                   columns (1:g2*NR)+g2*NR*(j-1)
%
% The integration itself is therefore delegated rather than duplicated. See
% INT_SEMISEP for the inputs, the index convention and the Delta tables.
%
% INPUTS
%   G, idxbeta, idxalpha, lims, Csize : as for INT_SEMISEP
%
% OUTPUTS
%   CMat : 3^ns3a by 1 cell array of packed coefficient matrices
%   ZL   : left monomial basis in s3a
%   ZR   : right monomial basis in s3a_dum
%
% See also INT_SEMISEP.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - int_semisep_AT
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
% AT, 2026: Initial coding
% MMP, 08/30/2026: Reduced to a repacking wrapper around 'int_semisep'. The
%                  two routines held the same integration code, so a fix or
%                  an optimization applied to one silently missed the other;
%                  the lookup table that replaced ismember in the key loop
%                  had gone into 'int_semisep' alone.

if nargin<5
    Csize = [];
end

% % % Do the integration once, in the shared routine.
[C_gam_alp_beta,ZL,ZR] = int_semisep(G,idxbeta,idxalpha,lims,Csize);

% % % Repack {gamma,beta,alpha} into one block matrix per gamma.
nbeta  = size(idxbeta,1);
nalpha = size(idxalpha,1);
ngam   = size(C_gam_alp_beta,1);

NG = prod([cellfun(@numel,G.Z),1]);
g1 = size(G.C,1)/NG;
g2 = size(G.C,2);
NL = prod([cellfun(@numel,ZL),1]);
NR = prod([cellfun(@numel,ZR),1]);
nrow = g1*NL;       ncol = g2*NR;

CMat = cell(ngam,1);
for k = 1:ngam
    M = sparse(nrow*nbeta,ncol*nalpha);
    for i = 1:nbeta
        for j = 1:nalpha
            B = C_gam_alp_beta{k,i,j};
            if isempty(B) || nnz(B)==0
                continue
            end
            M((1:nrow)+nrow*(i-1),(1:ncol)+ncol*(j-1)) = B;
        end
    end
    CMat{k} = M;
end

end
