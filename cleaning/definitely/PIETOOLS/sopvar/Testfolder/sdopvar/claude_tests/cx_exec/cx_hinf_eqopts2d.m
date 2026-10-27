function eq_opts = cx_hinf_eqopts2d(Km,eq_opts,ztol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EQ_OPTS = CX_HINF_EQOPTS2D(KM,EQ_OPTS,ZTOL) is 'get_eq_opts_2D' read off a
% square 2-D container KM instead of an opvar2d: it switches on the slack
% 'exclude' entries of components that are zero in KM, and the 'sep'
% entries whose integral pairs are equal, with the stock tests verbatim
% (executives/utility_functions/get_eq_opts_2D.m):
%
%   exclude(1)     R00 zero
%   exclude(8:16)  R22{i,j} zero, idcs [8 9 10 11 13 14 12 15 16] over
%                  the 3 x 3 cells in column-major order
%   sep(3)..(6)    all(diff.C <= ztol) for R22{2,1}-R22{3,1}, R22{1,2}-
%                  R22{1,3}, [R22{2,2}-R22{3,2}, R22{2,3}-R22{3,3}] and
%                  [R22{2,2}-R22{2,3}, R22{3,2}-R22{3,3}], under the same
%                  guards. No abs(), as in the stock routine: a difference
%                  whose coefficients are all negative also sets sep.
%
% The stock Qop merges every R part into one R block and every L_2[s1,s2]
% part into one L_2[s1,s2] block; a component is zero iff it is zero in every
% container block of that space pair, and a difference is <= ztol iff it is
% in every one. The executive runs clean_opvar(Qop,1e-12) first, which has
% no container form; coefficients with |c| <= ztol are treated as zero here
% to the same effect. Container cell (g1,g2) of an L_2[s1,s2] block is the
% parameter over sorted S3 = (s1,s2), i.e. R22{g1,g2} (opvar2d2sopvar's map).
% L_2[s1], L_2[s2] spaces are not handled (error), as in cx_hinf_2dpos.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3,    ztol = 1e-12;   end
exc = eq_opts.exclude;      sep = eq_opts.sep;
nvs = sum(Km.space_out,2);
if any(nvs==1) || numel(Km.vars)>2
    error('cx_hinf_eqopts2d:space','Only R^n and L2[s1,s2] spaces are handled.')
end
iR = find(nvs==0);      iX = find(nvs==2);
if ~isempty(iR) && all_zero(Km,iR,1)
    exc(1) = 1;
end
if ~isempty(iX)
    pidx = [8,9,10,11,13,14,12,15,16];
    for ii = 1:9
        if all_zero(Km,iX,ii),  exc(pidx(ii)) = 1;  end
    end
    c = @(g1,g2) sub2ind([3,3],g1,g2);
    if ~exc(9) && ~exc(10) && diff_le(Km,iX,c(2,1),c(3,1),ztol)
        sep(3) = 1;
    end
    if ~exc(11) && ~exc(12) && diff_le(Km,iX,c(1,2),c(1,3),ztol)
        sep(4) = 1;
    end
    if ((~exc(13) && ~exc(14)) || (~exc(15) && ~exc(16))) && ...
            diff_le(Km,iX,c(2,2),c(3,2),ztol) && diff_le(Km,iX,c(2,3),c(3,3),ztol)
        sep(5) = 1;
    end
    if ((~exc(13) && ~exc(15)) || (~exc(14) && ~exc(16))) && ...
            diff_le(Km,iX,c(2,2),c(2,3),ztol) && diff_le(Km,iX,c(3,2),c(3,3),ztol)
        sep(6) = 1;
    end
end
eq_opts.exclude = exc;      eq_opts.sep = sep;
end

function tf = all_zero(Km,idx,q)
% Parameter cell q zero in every block (i,j), i,j in idx.
tf = true;
for i = idx(:)',    for j = idx(:)'
    [A,Bq] = cellAB(Km.C{i,j},q,1e-12);
    if any(A) || nnz(Bq),   tf = false;     return,     end
end,                end
end

function tf = diff_le(Km,idx,q1,q2,ztol)
% All coefficients of cell q1 minus cell q2 <= ztol in every block. The two
% cells of one block share its ZL/ZR, so their coefficients align entrywise.
tf = true;
for i = idx(:)',    for j = idx(:)'
    [A1,B1] = cellAB(Km.C{i,j},q1,ztol);    [A2,B2] = cellAB(Km.C{i,j},q2,ztol);
    dB = B1 - B2;
    if any(A1-A2 > ztol) || any(nonzeros(dB) > ztol)    % implicit zeros pass
        tf = false;     return
    end
end,                end
end

function [A,Bq] = cellAB(B,q,ztol)
% Constant part A (mn x 1) and decision rows Bq (nd x mn) of parameter cell
% q of block B, cleaned at ztol (clean_opvar); empty for an absent block.
A = zeros(0,1);     Bq = sparse(0,0);
if isempty(B),  return,     end
nL = prod([cellfun(@numel,B.ZL),1]);    nR = prod([cellfun(@numel,B.ZR),1]);
mn = B.dims(1)*nL*B.dims(2)*nR;
if isa(B,'sdopvar')
    A = B.params.A{q};      Bq = B.params.B{q};
    if numel(A)~=mn,        A = zeros(mn,1);    end         % 0 / [] shorthand
    if size(Bq,2)~=mn,      Bq = sparse(numel(B.Zd),mn);    end
else                                                        % fixed sopvar block
    A = B.params{q};        Bq = sparse(0,mn);
    if numel(A)~=mn,        A = zeros(mn,1);    end
end
A = full(A(:));     A(abs(A)<=ztol) = 0;
Bq = Bq.*(abs(Bq)>ztol);
end
