function [Y,info] = lpigetdual_sop(prog,tag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Y,INFO] = LPIGETDUAL_SOP(PROG,TAG) assembles the dual multipliers of the
% container equality P == 0 whose row record is prog.sopeq{TAG} into a
% fixed container Y with the block grid of P. Each row r of that equality
% sets one coefficient of P to zero: the linear position c of the
% coefficient matrix of parameter cell k of block (i,j). The multiplier y_r
% of that row in prog.solinfo.y is placed at position c of parameter cell
% k of block (i,j) of Y, so that sum_r y_r*coeff_r(P) is the sum, over the
% constrained positions, of the coefficient of P times that of Y.
% Positions without a row are 0; blocks without a record entry are [].
%
% INPUT
% - prog:   'struct', an LPI program, solved, with the fields
%             sopeq    cell, one row record per container equality that
%                      LPI_EQ_CDOPVAR wrote on a program carrying this
%                      field; record TAG a 'struct' with fields rows (1 x 2
%                      double, the first and last row of the equality),
%                      entries (struct array, one entry per block with
%                      rows, in row order, with fields i, j (double, the
%                      block), cells (double row, its parameter cells with
%                      rows), cols (cell, per cell the linear positions c
%                      of its rows, in row order) and blk ('struct' with the
%                      fields dims, NL, NR, psize, vars, ZL, ZR and dom of
%                      the block)), M, N (double, the block grid of P) and
%                      symmetric (logical);
%             solinfo  'struct' with the field y, double, the dual
%                      multipliers of the rows of the program, in the order
%                      the rows were written;
% - tag:    'double', the index of the record in prog.sopeq (the second
%           output of LPI_EQ_SOP).
%
% OUTPUT
% - Y:      'copvar' object with M x N blocks: for each record entry, block
%           (i,j) a 'sopvar' object with the dims, vars, ZL, ZR and dom of
%           block (i,j) of P, whose parameter cells listed in cells hold the
%           multipliers at the positions cols and 0 elsewhere, its other
%           cells []; every other block [];
% - info:   'struct' with fields
%             y          double column, the multipliers of the rows of the
%                        equality;
%             rows       1 x 2 double, its first and last row;
%             nrow       double, its number of rows;
%             blocks     double, one row [i,j] per block of Y with an
%                        entry, in row order;
%             norm1      double, sum(abs(y));
%             symmetric  logical, true for an equality imposed with
%                        'symmetric': rows exist for one block of each
%                        pair (i,j), (j,i) and, in a diagonal block, for
%                        one parameter cell of each adjoint pair, and Y
%                        holds those only.
%
% NOTES
% PROG without the record TAG or without a nonempty prog.solinfo.y, and a
% record whose rows exceed numel(prog.solinfo.y), are errors. So is a
% record with no entry in some row or column of the block grid: the
% 'copvar' constructor cannot size that row or column.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpigetdual_sop
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
% Initial coding MMP, 10/08/2026
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body; code unchanged. The deleted header text was a
%                memo reference, the invariance of the dual under the b
%                normalization, the requirement sos_opts.simplify = 0
%                (sossolve returns y in the rows of the unsimplified
%                program), a pointer to a caller, and the claim that under
%                'symmetric' the lower-triangular blocks are written
%                (LPI_EQ_CDOPVAR writes the blocks (i,j) with i <= j and a
%                block below the diagonal only when its mirror is not an
%                'sdopvar').

% The record TAG and the multipliers must exist
if ~isfield(prog,'sopeq') || tag<1 || tag>numel(prog.sopeq)
    error('lpigetdual_sop:tag','No row record %d in the program; was the equality written by lpi_eq_sop?',tag);
end
if ~isfield(prog,'solinfo') || ~isfield(prog.solinfo,'y') || isempty(prog.solinfo.y)
    error('lpigetdual_sop:y','The program carries no dual multipliers (solve it first; sossolve keeps y).');
end
rec = prog.sopeq{tag};
y = full(prog.solinfo.y(:));
nrow = rec.rows(2)-rec.rows(1)+1;
% The recorded rows must lie within y
if rec.rows(2)>numel(y)
    error('lpigetdual_sop:rows',['The record expects rows up to %d but the dual has %d entries: ' ...
          'the program was simplified or rows were added out of order.'],rec.rows(2),numel(y));
end
% Rows in writing order: per record entry (a block), per parameter cell,
% the multipliers placed at the positions of its rows, 0 elsewhere
C = cell(rec.M,rec.N);
r = rec.rows(1);
blocks = zeros(0,2);
for e = 1:numel(rec.entries)
    en = rec.entries(e);    blk = en.blk;
    m = blk.dims(1);    n = blk.dims(2);    NL = blk.NL;   NR = blk.NR;
    params = cell(blk.psize);
    for t = 1:numel(en.cells)
        k = en.cells(t);    cols = en.cols{t};   nk = numel(cols);
        vals = y(r:r+nk-1);     r = r+nk;
        Pk = sparse(m*NL,n*NR);
        Pk(cols) = vals;
        params{k} = Pk;
    end
    % The block: a fixed 'sopvar' on the spaces and bases of block (i,j)
    vio = struct('out',{reshape(blk.vars.out,1,[])},'in',{reshape(blk.vars.in,1,[])});
    C{en.i,en.j} = sopvar(params,vio,blk.ZL,blk.ZR,blk.dom,blk.dims);
    blocks(end+1,:) = [en.i,en.j]; %#ok<AGROW>
end
% Every recorded row consumed
if r-1~=rec.rows(2)
    error('lpigetdual_sop:count','Row count mismatch: consumed %d, recorded %d.',r-rec.rows(1),nrow);
end
% The container and the multipliers of this equality
Y = copvar(C);
info = struct('y',y(rec.rows(1):rec.rows(2)),'rows',rec.rows,'nrow',nrow,'blocks',blocks, ...
              'norm1',sum(abs(y(rec.rows(1):rec.rows(2)))),'symmetric',rec.symmetric);
end
