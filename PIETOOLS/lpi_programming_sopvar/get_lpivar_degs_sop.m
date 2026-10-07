function Qdeg = get_lpivar_degs_sop(Rop,Top)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QDEG = GET_LPIVAR_DEGS_SOP(ROP[,TOP]) returns the degrees the Q-form
% executives give 'lpivar' for the free operator Q, read off the positive
% operator ROP, for either family:
%
%   ROP opvar/dopvar/opvar2d/dopvar2d   get_lpivar_degs(Rop,Top), unchanged;
%   ROP copvar/cdopvar (1-D)            the same 1-D rule, read off the
%                                       container (TOP not needed).
%
% Stock 1-D rule (executives/utility_functions/get_lpivar_degs.m):
%   deg1 = max monomial degree of Q1, Q2 and R.R0
%   deg2 = max single-variable degree of R.R1 and R.R2
%   Qdeg = [deg1, deg2, deg2-1]
% Container reading, per block and parameter cell, over the monomials whose
% coefficient (A, or any row of B) is nonzero:
%   block on one shared variable, cell 1 (multiplier, R0)   deg1 <- ZL
%   block on one shared variable, cells 2/3 (R1, R2)        deg2 <- ZL, ZR
%   block without a shared variable (Q1, Q2, P)             deg1 <- ZL, ZR
% The multiplier cell reads ZL only: by the canonical multiplier form it
% has no content beyond ZR degree 0.
%
% INPUT
% - Rop:    positive operator (legacy or container);
% - Top:    (legacy only) the PIE's T operator, as get_lpivar_degs takes.
%
% OUTPUT
% - Qdeg:   1 x 3 [deg1, deg2, deg2-1] (legacy: get_lpivar_degs' output).
%
% NOTES
% A stock degmat may list a monomial whose coefficients are all zero; the
% container reading cannot see one, so its Qdeg is entrywise <= the stock
% value for the same operator, and equal whenever every listed monomial has
% a nonzero coefficient (measured on the shipped presets: see the README).
% deg2 = 0 gives Qdeg(3) = -1, as the stock rule does; lpivar_cdopvar
% rejects it. 1-D containers only: a block over two or more variables is an
% error (the stock 2-D rule returns a struct, built from other degrees).
% Cost: one pass over the nonzero coefficients, O(nnz(B)) per cell with a
% transient sparse logical of nnz(B) entries; the stock rule reads degmat
% and is free of the number of decision variables.
%
% See also GET_LPIVAR_DEGS, LPIVAR_CDOPVAR, POSLPIVAR_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_lpivar_degs_sop
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
% Initial coding MMP, 10/06/2026. Replaces the test-folder
%                cx_stability_lpivar_degs and cx_hinf_qdeg (cx_exec, MMP
%                09/25/2026), two readings of the same rule; the body is
%                cx_stability_lpivar_degs' (it also reads fixed sopvar
%                blocks), restricted to 1-D with cx_hinf_qdeg's guard.

if ~(isa(Rop,'copvar') || isa(Rop,'cdopvar'))
    Qdeg = get_lpivar_degs(Rop,Top);
    return
end
deg1 = 0;   deg2 = 0;
for b = 1:numel(Rop.C)
    P = Rop.C{b};
    if isempty(P),  continue,   end
    if numel(P.vars.out)>1 || numel(P.vars.in)>1
        error('get_lpivar_degs_sop:dim',['The 1-D rule reads blocks over at '...
              'most one variable; this block has %d out and %d in.'],...
              numel(P.vars.out),numel(P.vars.in))
    end
    shared = ~isempty(intersect(P.vars.out,P.vars.in));    % L2[s] -> L2[s]
    NL = prod([cellfun(@numel,P.ZL),1]);    NR = prod([cellfun(@numel,P.ZR),1]);
    nrow = P.dims(1)*NL;
    if isa(P,'sdopvar'),    ncell = numel(P.params.A);  else,   ncell = numel(P.params);    end
    for k = 1:ncell                                 % 3 cells if shared, else 1
        if isa(P,'sdopvar')
            A = P.params.A{k};  B = P.params.B{k};
            nz = [];
            if ~isempty(A),     nz = find(A(:)~=0);   end
            if ~isempty(B),     nz = union(nz,find(any(B~=0,1))');  end    % O(nnz(B))
        else
            nz = find(P.params{k}(:)~=0);           % fixed sopvar block
        end
        if isempty(nz),     continue,   end
        r = mod(nz-1,nrow)+1;       c = floor((nz-1)/nrow)+1;
        eL = maxexp(P.ZL,mod(r-1,NL)+1);            % row = (component, monomial inner)
        eR = maxexp(P.ZR,mod(c-1,NR)+1);
        if shared && k==1,      deg1 = max(deg1,eL);            % R0
        elseif shared,          deg2 = max([deg2,eL,eR]);       % R1, R2
        else,                   deg1 = max([deg1,eL,eR]);       % Q1, Q2, P
        end
    end
end
Qdeg = [deg1,deg2,deg2-1];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = maxexp(Z,idx)
% Largest exponent over the monomials IDX of a basis over at most one
% variable (Z = {} or {exponents}).
e = 0;
if ~isempty(Z),     e = max([e; reshape(Z{1}(idx),[],1)]);     end
end
