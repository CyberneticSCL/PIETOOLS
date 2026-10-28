function Qdeg = get_lpivar_degs_sop(Rop,Top)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QDEG = GET_LPIVAR_DEGS_SOP(ROP,TOP) returns the degrees of the free
% operator Q of an LPI, for 'lpivar' or LPIVAR_CDOPVAR, read off the
% positive operator ROP.
%
% For a legacy ROP ('opvar', 'dopvar', 'opvar2d', 'dopvar2d') it returns
% get_lpivar_degs(Rop,Top): in 1-D, Qdeg = [deg1, deg2, deg2-1] with deg1
% the largest degree of a monomial of Q1, Q2 and R.R0 and deg2 the largest
% degree in either variable of a monomial of R.R1 and R.R2; in 2-D, the
% largest degrees of the parameters of Rop in each variable.
%
% For a container ROP ('copvar', 'cdopvar'), each parameter cell of a
% block has the coefficient A + B'*d (d the decision variables; A for a
% 'sopvar' block) between the monomial bases ZL of the output variables
% and ZR of the input variables, and only the monomials of a nonzero entry
% of A or a nonzero column of B are read:
%
%   every block over at most one variable on each side: Qdeg = [deg1,
%     deg2, deg2-1], deg1 the largest exponent of ZL in the multiplier cell
%     of a block on one shared variable (R0) and of ZL and ZR in a block
%     without shared variable (P, Q1, Q2), deg2 the largest exponent of ZL
%     and ZR in the integral cells of a block on one shared variable (R1,
%     R2);
%   a block over two or more variables on a side: mult the largest ZL
%     exponent of a shared variable in a multiplier direction, int the
%     largest [ZL, ZR] exponents of a shared variable in an integral
%     direction, out the largest ZL exponent of a variable of the output
%     side only, and in the largest ZR exponent of a variable of the input
%     side only.
%
% INPUT
% - Rop:  'opvar', 'dopvar', 'opvar2d' or 'dopvar2d' object, or 'copvar'
%         or 'cdopvar' object with 'sopvar' and 'sdopvar' blocks;
% - Top:  'opvar', 'dopvar', 'opvar2d' or 'dopvar2d' object, the T operator
%         of the PIE, of the dimension of Rop: required for a legacy Rop,
%         not used for a container.
%
% OUTPUT
% - Qdeg: 1x3 double [deg1, deg2, deg2-1] in 1-D (last entry -1 when
%         deg2 = 0); for a legacy 2-D Rop, a 'struct' with fields dx (3x1
%         cell), dy (1x3 cell) and d2 (3x3 cell) of degree arrays; for a
%         container with a block over two or more variables, a 'struct'
%         with fields mult (double), int (1x2 double), out and in (double).
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
% MMP, 10/16/2026: moved from lpi_programming_sopvar/ to executives_sopvar/
%                  utility_functions/, as the stock get_lpivar_degs is in
%                  executives/utility_functions; the README it cites is
%                  executives_sopvar's. Comment only otherwise.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), step comments in the body and headers
%                  of the local functions; code unchanged. The deleted
%                  header text was the comparison with the stock rule, a
%                  measurement and a README reference, a cost note and
%                  history.

% Legacy operator: the stock rule
if ~(isa(Rop,'copvar') || isa(Rop,'cdopvar'))
    Qdeg = get_lpivar_degs(Rop,Top);
    return
end
% A block over two or more variables on a side: the per-role caps of the
% N-D reading
for b = 1:numel(Rop.C)                                                      % MMP, 10/08/2026
    P = Rop.C{b};                                                           % MMP, 10/08/2026
    if ~isempty(P) && (numel(P.vars.out)>1 || numel(P.vars.in)>1)           % MMP, 10/08/2026
        Qdeg = qdeg_nd(Rop);    return                                      % MMP, 10/08/2026
    end                                                                     % MMP, 10/08/2026
end                                                                         % MMP, 10/08/2026
% 1-D reading: deg1 and deg2 over the nonzero blocks
deg1 = 0;   deg2 = 0;
for b = 1:numel(Rop.C)
    P = Rop.C{b};
    if isempty(P),  continue,   end
    if numel(P.vars.out)>1 || numel(P.vars.in)>1
        error('get_lpivar_degs_sop:dim',['The 1-D rule reads blocks over at '...
              'most one variable; this block has %d out and %d in.'],...
              numel(P.vars.out),numel(P.vars.in))
    end
    % Shared variable, basis sizes, rows of a coefficient matrix, cells
    shared = ~isempty(intersect(P.vars.out,P.vars.in));    % L2[s] -> L2[s]
    NL = prod([cellfun(@numel,P.ZL),1]);    NR = prod([cellfun(@numel,P.ZR),1]);
    nrow = P.dims(1)*NL;
    if isa(P,'sdopvar'),    ncell = numel(P.params.A);  else,   ncell = numel(P.params);    end
    for k = 1:ncell                                 % 3 cells if shared, else 1
        % Linear positions of the nonzero coefficients of cell k
        if isa(P,'sdopvar')
            A = P.params.A{k};  B = P.params.B{k};
            nz = [];
            if ~isempty(A),     nz = find(A(:)~=0);   end
            if ~isempty(B),     nz = union(nz,find(any(B~=0,1))');  end    % O(nnz(B))
        else
            nz = find(P.params{k}(:)~=0);           % fixed sopvar block
        end
        if isempty(nz),     continue,   end
        % Largest left and right exponents over those positions
        r = mod(nz-1,nrow)+1;       c = floor((nz-1)/nrow)+1;
        eL = maxexp(P.ZL,mod(r-1,NL)+1);            % row = (component, monomial inner)
        eR = maxexp(P.ZR,mod(c-1,NR)+1);
        % The multiplier cell reads ZL only: by the canonical multiplier
        % form it has no content beyond ZR degree 0
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
% E = MAXEXP(Z,IDX) returns the largest exponent over the monomials IDX of
% a basis over at most one variable, 0 for a basis over none.
% INPUT  Z:   cell, {} or {exponents} with exponents a double column;
%        idx: double, the indices of the monomials.
% OUTPUT e:   double.
e = 0;
if ~isempty(Z),     e = max([e; reshape(Z{1}(idx),[],1)]);     end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = qdeg_nd(Rm)                                                  % MMP, 10/08/2026
% DEG = QDEG_ND(RM) returns the per-role degree caps of LPIVAR_CDOPVAR read
% off the container RM: for every nonzero coefficient of every parameter
% cell, the exponents of its left and right monomials in each shared
% variable go to 'mult' (gamma_k = 1, left only) or 'int' (gamma_k = 2, 3);
% in an output-only variable to 'out', in an input-only one to 'in'.
% INPUT  Rm:  'copvar' or 'cdopvar' object.
% OUTPUT deg: 'struct' with fields mult (double), int (1x2 double, [left
%             right]), out and in (double).
deg = struct('mult',0,'int',[0 0],'out',0,'in',0);
for i = 1:size(Rm.C,1)
    for j = 1:size(Rm.C,2)
        B = Rm.C{i,j};
        if isempty(B),  continue,   end
        % Shared, output-only and input-only variables; exponent tables of
        % the bases; size of a coefficient matrix
        vo = B.vars.out;    vi = B.vars.in;
        S3 = intersect(vo,vi);      n3 = numel(S3);
        EL = kron_exps(B.ZL);       ER = kron_exps(B.ZR);
        nL = size(EL,1);    nR = size(ER,1);
        m = B.dims(1)*nL;   n = B.dims(2)*nR;
        [~,pL] = ismember(S3,vo);   [~,pR] = ismember(S3,vi);
        s2 = ~ismember(vo,S3);      s1 = ~ismember(vi,S3);
        if isa(B,'sdopvar'),    ncell = numel(B.params.A);  else,   ncell = numel(B.params);    end
        for q = 1:ncell
            % Nonzero entries of the coefficient of cell q: of A, or a
            % column of B with a nonzero entry
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
            % Exponents of their left and right monomials
            [r,c] = find(reshape(nz,m,n));
            if isempty(r),  continue,   end
            eL = EL(mod(r-1,nL)+1,:);   eR = ER(mod(c-1,nR)+1,:);
            % Multi-index of cell q over the shared variables (1
            % multiplier, 2 or 3 integral), and the caps of each
            g = cell(1,max(n3,1));      [g{:}] = ind2sub(3*ones(1,max(n3,2)),q);
            for k = 1:n3
                if g{k}==1
                    deg.mult = max(deg.mult,max(eL(:,pL(k))));
                else
                    deg.int = max(deg.int,[max(eL(:,pL(k))), max(eR(:,pR(k)))]);
                end
            end
            % Caps of the output-only and input-only variables
            if any(s2),     deg.out = max(deg.out,max(max(eL(:,s2))));  end
            if any(s1),     deg.in  = max(deg.in, max(max(eR(:,s1))));  end
        end
    end
end
end


function E = kron_exps(Z)                                                   % MMP, 10/08/2026
% E = KRON_EXPS(Z) returns the exponent table of the tensor monomial basis
% Z{1} x Z{2} x ..., first factor outer, one row per basis element.
% INPUT  Z: cell of double columns, the exponents of each variable.
% OUTPUT E: double, one row per basis element and one column per variable
%           (1 x 0 for an empty Z).
E = zeros(1,0);
for k = 1:numel(Z)
    z = Z{k}(:);
    E = [kron(E,ones(numel(z),1)), repmat(z,size(E,1),1)];
end
end
