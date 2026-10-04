function prog = impose_eq_rows(prog,Eq)
% PROG = IMPOSE_EQ_ROWS(PROG,EQ) imposes the blocks EQ.Cs{:}, each (q+1) x n_k
% over the names EQ.Zd with row 1 the constant term, joining consecutive    % MMP, 09/26/2026
% blocks into as few soseq calls as the cap below allows, each on one       % MMP, 09/26/2026
% variable-free row dpvar (see the layout note in 'collect_eq_rows'). Columns
% keep their order, so the At and b produced are those of one soseq per     % MMP, 09/26/2026
% block, concatenated as sossolve concatenates prog.expr.                   % MMP, 09/26/2026
% Cap: at most N = numel(prog.decvartable) nonzeros per soseq; a larger     % MMP, 09/26/2026
% block goes alone. Each soseq pays an O(N) scan of the program's names     % MMP, 09/26/2026
% (getequation), which N nonzeros of work amortize; and soseq holds ~9      % MMP, 09/26/2026
% copies of its coefficients (combine, compress, getequation), so each      % MMP, 09/26/2026
% soseq's transient is O(N) bytes. Eq.Cs itself, one copy of all the        % MMP, 09/26/2026
% collected coefficients (O(nnz)), lives until the last soseq. Measured,    % MMP, 09/26/2026
% 2-D container Hinf, q = 2.3e6: one uncapped soseq peaked at 974 MB,       % MMP, 09/26/2026
% capped 447 MB (per parameter 291 MB, stock lpi_eq_2d 540 MB).             % MMP, 09/26/2026
% Each batch is now written by 'lpi_soseq' when EQ.pos gives the table rows % MMP, 10/02/2026
% of EQ.Zd; soseq runs only without it (a name repeated in EQ.Zd, or a      % MMP, 10/02/2026
% hand-built EQ). lpi_soseq scans no names, so the cap no longer amortizes  % MMP, 10/02/2026
% a scan; it still bounds the per-batch transient (the joined M, 16 B per   % MMP, 10/02/2026
% nonzero, and lpi_soseq's ~80 B), and every entry is the one soseq wrote.  % MMP, 10/02/2026
%
% See also COLLECT_EQ_ROWS, LPI_EQ_SDOPVAR, LPI_EQ_CDOPVAR, SOSEQ.

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

if isempty(Eq.Cs),  return,     end                                         % MMP, 09/26/2026
% Rows of Eq.Zd known (collect_eq_rows' Eq.pos): no soseq.                  % MMP, 10/02/2026
direct = isfield(Eq,'pos') && numel(Eq.pos)==numel(Eq.Zd);                  % MMP, 10/02/2026
cap = numel(prog.decvartable);                                              % MMP, 09/26/2026
nz = cellfun(@nnz,Eq.Cs);                                                   % MMP, 09/26/2026
k0 = 1;                                                                     % MMP, 09/26/2026
while k0<=numel(Eq.Cs)                                                      % MMP, 09/26/2026
    k1 = k0;    tot = nz(k0);                                               % MMP, 09/26/2026
    while k1<numel(Eq.Cs) && tot+nz(k1+1)<=cap                              % MMP, 09/26/2026
        k1 = k1+1;  tot = tot+nz(k1);                                       % MMP, 09/26/2026
    end                                                                     % MMP, 09/26/2026
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
    Dk = dpvar(M(rows,:),zeros(1,0),{},Eq.Zd(rows(2:end)-1),[1,size(M,2)]); % MMP, 09/26/2026
    % Free the joined copy before soseq makes its own.                      % MMP, 09/26/2026
    M = [];                                                                 % MMP, 09/26/2026
    prog = soseq(prog,Dk);                                                  % MMP, 09/26/2026
    end                                                                     % MMP, 10/02/2026
    % END MMP, 10/02/2026
    k0 = k1+1;                                                              % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026
