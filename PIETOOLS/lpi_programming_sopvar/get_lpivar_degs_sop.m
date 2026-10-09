function Qdeg = get_lpivar_degs_sop(Rop,Top)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QDEG = GET_LPIVAR_DEGS_SOP(ROP[,TOP]) returns the degrees the Q-form
% executives give 'lpivar' for the free operator Q, read off the positive
% operator ROP, for either family:
%
%   ROP opvar/dopvar/opvar2d/dopvar2d   get_lpivar_degs(Rop,Top), unchanged;
%   ROP copvar/cdopvar (1-D)            the same 1-D rule, read off the
%                                       container (TOP not needed);
%   ROP copvar/cdopvar (2-D, N-D)       the per-role degree caps of         % MMP, 10/08/2026
%                                       lpivar_cdopvar, read off the        % MMP, 10/08/2026
%                                       container: see NOTES.               % MMP, 10/08/2026
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
% rejects it. A container with a block over two or more variables returns  % MMP, 10/08/2026
% the struct lpivar_cdopvar takes: 'mult' the largest left degree of a      % MMP, 10/08/2026
% nonzero coefficient in a multiplier direction, 'int' the largest [left,   % MMP, 10/08/2026
% right] degrees in an integral direction, 'out' / 'in' the largest degree  % MMP, 10/08/2026
% in an output-only / input-only variable (the cx_hinf_qdeg2d reading of    % MMP, 10/08/2026
% 10/06/2026; the stock 2-D rule builds its struct from other degrees).     % MMP, 10/08/2026
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
% MMP, 10/08/2026: N-D containers no longer error: a block over two or more
%                  variables sends the whole container to the per-role
%                  reading of the test-folder cx_hinf_qdeg2d (MMP,
%                  09/26/2026), copied below as qdeg_nd, so that the 2-D
%                  container executives (executives_sopvar) size their free
%                  operator from one routine. The 1-D rule is unchanged.

if ~(isa(Rop,'copvar') || isa(Rop,'cdopvar'))
    Qdeg = get_lpivar_degs(Rop,Top);
    return
end
for b = 1:numel(Rop.C)                                                      % MMP, 10/08/2026
    P = Rop.C{b};                                                           % MMP, 10/08/2026
    if ~isempty(P) && (numel(P.vars.out)>1 || numel(P.vars.in)>1)           % MMP, 10/08/2026
        Qdeg = qdeg_nd(Rop);    return                                      % MMP, 10/08/2026
    end                                                                     % MMP, 10/08/2026
end                                                                         % MMP, 10/08/2026
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = qdeg_nd(Rm)                                                  % MMP, 10/08/2026
% Per-role caps of lpivar_cdopvar read off an N-D container: for every
% nonzero coefficient of every parameter cell, the exponents of its left and
% right monomials in each shared variable go to 'mult' (gamma_k = 1, left
% only) or 'int' (gamma_k = 2, 3); in an output-only variable to 'out', in
% an input-only one to 'in'. The test-folder cx_hinf_qdeg2d reading.
deg = struct('mult',0,'int',[0 0],'out',0,'in',0);
for i = 1:size(Rm.C,1)
    for j = 1:size(Rm.C,2)
        B = Rm.C{i,j};
        if isempty(B),  continue,   end
        vo = B.vars.out;    vi = B.vars.in;
        S3 = intersect(vo,vi);      n3 = numel(S3);
        EL = kron_exps(B.ZL);       ER = kron_exps(B.ZR);
        nL = size(EL,1);    nR = size(ER,1);
        m = B.dims(1)*nL;   n = B.dims(2)*nR;
        [~,pL] = ismember(S3,vo);   [~,pR] = ismember(S3,vi);
        s2 = ~ismember(vo,S3);      s1 = ~ismember(vi,S3);
        if isa(B,'sdopvar'),    ncell = numel(B.params.A);  else,   ncell = numel(B.params);    end
        for q = 1:ncell
            nz = false(m*n,1);
            if isa(B,'sdopvar')
                A = B.params.A{q};
                if numel(A)==m*n,   nz = nz | full(A(:)~=0);    end
                Bq = B.params.B{q};
                if size(Bq,2)==m*n, nz = nz | full(any(Bq~=0,1)).';    end
            else
                A = B.params{q};
                if numel(A)==m*n,   nz = nz | full(A(:)~=0);    end
            end
            [r,c] = find(reshape(nz,m,n));
            if isempty(r),  continue,   end
            eL = EL(mod(r-1,nL)+1,:);   eR = ER(mod(c-1,nR)+1,:);
            g = cell(1,max(n3,1));      [g{:}] = ind2sub(3*ones(1,max(n3,2)),q);
            for k = 1:n3
                if g{k}==1
                    deg.mult = max(deg.mult,max(eL(:,pL(k))));
                else
                    deg.int = max(deg.int,[max(eL(:,pL(k))), max(eR(:,pR(k)))]);
                end
            end
            if any(s2),     deg.out = max(deg.out,max(max(eL(:,s2))));  end
            if any(s1),     deg.in  = max(deg.in, max(max(eR(:,s1))));  end
        end
    end
end
end


function E = kron_exps(Z)                                                   % MMP, 10/08/2026
% Exponent table of the tensor monomial basis Z{1} x Z{2} x ..., first
% factor outer, one row per basis element.
E = zeros(1,0);
for k = 1:numel(Z)
    z = Z{k}(:);
    E = [kron(E,ones(numel(z),1)), repmat(z,size(E,1),1)];
end
end
