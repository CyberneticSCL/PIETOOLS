function [Y,info] = lpigetdual_sop(prog,tag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Y,INFO] = LPIGETDUAL_SOP(PROG,TAG) assembles the solver's dual
% multipliers of one container equality into its signed dual kernel.
%
% The container equality written by LPI_EQ_SOP (TAG its index in
% prog.sopeq) consists of rows "coefficient of monomial (a,b) in cell gamma
% of block (i,j) = 0". The SDP dual carries one multiplier per row,
% prog.solinfo.y. Placing each multiplier at the coefficient position its
% row constrains gives, block by block, a FIXED operator Y of the same class
% as the constrained one: the signed kernel of the proof program's dual
% (minimum_cap_primal_dual.md (4)-(5)), with
%
%     sum_rows y_r coeff_r(P) = <P, Y>     for the constrained operator P
%
% in the pairing that counts each written position once. Under
% 'symmetric' only the lower cells and the lower-triangular blocks are
% written, so Y holds those; its symmetrisation is the caller's (see
% PIE_WITNESS_SOP, which reads the input block's kernel).
%
% INPUT
% - prog: the SOLVED program (lpisolve with sos_opts.simplify = 0, so that
%         the rows keep the writer's order and count);
% - tag:  index into prog.sopeq (second output of LPI_EQ_SOP).
%
% OUTPUT
% - Y:    'copvar' over the constrained container's spaces, the written
%         blocks fixed 'sopvar's, the others [];
% - info: y (the multipliers of this equality), rows, nrow, blocks (the
%         (i,j) written), and norm1 = sum |y|.
%
% NOTES
% Requires prog.solinfo.y, which sossolve stores for SeDuMi and MOSEK;
% errors otherwise. The multipliers are those of the b-normalised program
% sossolve solves; the dual is invariant under that scaling, so no rescaling
% is applied. The positions are those collect_eq_rows recorded (the kept
% canonical positions whose column had a nonzero), in writing order; a
% zero column that lpi_soseq dropped produces no row and is not listed.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(prog,'sopeq') || tag<1 || tag>numel(prog.sopeq)
    error('lpigetdual_sop:tag','No row record %d in the program; was the equality written by lpi_eq_sop?',tag);
end
if ~isfield(prog,'solinfo') || ~isfield(prog.solinfo,'y') || isempty(prog.solinfo.y)
    error('lpigetdual_sop:y','The program carries no dual multipliers (solve it first; sossolve keeps y).');
end
rec = prog.sopeq{tag};
y = full(prog.solinfo.y(:));
nrow = rec.rows(2)-rec.rows(1)+1;
if rec.rows(2)>numel(y)
    error('lpigetdual_sop:rows',['The record expects rows up to %d but the dual has %d entries: ' ...
          'the program was simplified or rows were added out of order.'],rec.rows(2),numel(y));
end
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
    vio = struct('out',{reshape(blk.vars.out,1,[])},'in',{reshape(blk.vars.in,1,[])});
    C{en.i,en.j} = sopvar(params,vio,blk.ZL,blk.ZR,blk.dom,blk.dims);
    blocks(end+1,:) = [en.i,en.j]; %#ok<AGROW>
end
if r-1~=rec.rows(2)
    error('lpigetdual_sop:count','Row count mismatch: consumed %d, recorded %d.',r-rec.rows(1),nrow);
end
Y = copvar(C);
info = struct('y',y(rec.rows(1):rec.rows(2)),'rows',rec.rows,'nrow',nrow,'blocks',blocks, ...
              'norm1',sum(abs(y(rec.rows(1):rec.rows(2)))),'symmetric',rec.symmetric);
end
