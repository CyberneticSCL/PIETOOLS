function prog = impose_eq_rows(prog,Eq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = IMPOSE_EQ_ROWS(PROG,EQ) appends to the LPI program PROG, for every
% column l of every block C = EQ.Cs{k}, in order, the equation
%
%       C(1,l) + sum_i C(1+i,l)*d_i = 0,        d_i named EQ.Zd{i},
%
% omitting all-zero columns. Consecutive blocks are joined into one 'eq'
% entry of prog.expr while their nonzeros total at most
% numel(prog.decvartable); a block with more nonzeros is an entry of its
% own.
%
% INPUT
% - prog:   'struct', an LPI program whose prog.decvartable holds every
%           name of EQ.Zd;
% - Eq:     'struct' with fields
%             Cs   cell of sparse double blocks, block k (q+1) x n_k with
%                  q = numel(Zd): row 1 the constant terms, row 1+i the
%                  coefficients of Zd{i}, one column per equation;
%             Zd   q x 1 cellstr, the decision variable names;
%             pos  (optional) double, q x 1, pos(i) the first row of
%                  prog.decvartable named Zd{i}; ignored unless it has q
%                  entries.
%
% OUTPUT
% - prog:   'struct', PROG with one 'eq' entry appended to prog.expr per
%           group of joined blocks; PROG unchanged when EQ.Cs is empty.
%
% NOTES
% EQ.pos is not checked against prog.decvartable: a wrong pos puts the
% coefficients on the wrong decision variables.
%
% See also COLLECT_EQ_ROWS, LPI_EQ_SDOPVAR, LPI_EQ_CDOPVAR, SOSEQ.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - impose_eq_rows
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
% Initial coding MMP, 09/30/2026: the impose mode of 'lpi_eq_sdopvar' as
%                  its own routine: its local function 'impose_rows'
%                  (MMP, 09/26/2026), moved verbatim with its markers and
%                  renamed; unstamped help lines are the rename and the
%                  pointer to the layout note, now in 'collect_eq_rows'.
%                  'lpi_eq_sdopvar' reached it for a struct P, a
%                  mode its public help did not show; 'lpi_eq_cdopvar' now
%                  calls it directly. Same soseq calls, same order, same cap.
% MMP, 10/02/2026: Batches with known table rows (Eq.pos) are written by
%                  'lpi_soseq' instead of soseq: same batches, same entries.
%                  In the 3-D heatNd build the 11 soseq calls cost 6.86 s
%                  in getequation, 5.9 s of it matching the 2.7e6 program
%                  names against each batch (measured). The entry above
%                  ("same soseq calls") now holds for the entries, not the
%                  calls.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), step comments in the body; code
%                  unchanged. The deleted header text was the cost model of
%                  the cap with its measurements, the description of the
%                  soseq and lpi_soseq routes, and a pointer to the layout
%                  note in 'collect_eq_rows'.

if isempty(Eq.Cs),  return,     end                                         % MMP, 09/26/2026
% Rows of Eq.Zd known (collect_eq_rows' Eq.pos): no soseq.                  % MMP, 10/02/2026
direct = isfield(Eq,'pos') && numel(Eq.pos)==numel(Eq.Zd);                  % MMP, 10/02/2026
% Cap on the nonzeros of one entry; the nonzeros of each block
cap = numel(prog.decvartable);                                              % MMP, 09/26/2026
nz = cellfun(@nnz,Eq.Cs);                                                   % MMP, 09/26/2026
% Groups of consecutive blocks: from block k0, add the next block while the
% total stays within the cap
k0 = 1;                                                                     % MMP, 09/26/2026
while k0<=numel(Eq.Cs)                                                      % MMP, 09/26/2026
    k1 = k0;    tot = nz(k0);                                               % MMP, 09/26/2026
    while k1<numel(Eq.Cs) && tot+nz(k1+1)<=cap                              % MMP, 09/26/2026
        k1 = k1+1;  tot = tot+nz(k1);                                       % MMP, 09/26/2026
    end                                                                     % MMP, 09/26/2026
    % The columns of the group, joined in order
    M = [Eq.Cs{k0:k1}];                                                     % MMP, 09/26/2026
    % BEGIN MMP, 10/02/2026: the entry soseq would append, written from the
    % triplets of M with the rows collect_eq_rows found ('lpi_soseq'). The
    % else branch is the former soseq path, kept for a list without pos.
    if direct                                                               % MMP, 10/02/2026
        prog = lpi_soseq(prog,M,Eq.Zd,Eq.pos);                              % MMP, 10/02/2026
        M = [];                                                             % MMP, 10/02/2026
    else                                                                    % MMP, 10/02/2026
    % Keep the constant row and the rows of names in use; compress would    % MMP, 09/26/2026
    % drop the others only after O(q) passes over their names.              % MMP, 09/26/2026
    rows = find(any(M,2));                                                  % MMP, 09/26/2026
    rows = [1; rows(rows>1)];                                               % MMP, 09/26/2026
    % The equations as a 1 x size(M,2) variable-free 'dpvar'
    Dk = dpvar(M(rows,:),zeros(1,0),{},Eq.Zd(rows(2:end)-1),[1,size(M,2)]); % MMP, 09/26/2026
    % Free the joined copy before soseq makes its own.                      % MMP, 09/26/2026
    M = [];                                                                 % MMP, 09/26/2026
    prog = soseq(prog,Dk);                                                  % MMP, 09/26/2026
    end                                                                     % MMP, 10/02/2026
    % END MMP, 10/02/2026
    k0 = k1+1;                                                              % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026
