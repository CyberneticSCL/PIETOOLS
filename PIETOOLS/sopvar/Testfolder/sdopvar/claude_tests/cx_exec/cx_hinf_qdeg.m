function Qdeg = cx_hinf_qdeg(Rm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QDEG = CX_HINF_QDEG(RM) is 'get_lpivar_degs(Rop,Top)' (1-D branch,
% executives/utility_functions/get_lpivar_degs.m) read off a 1-D container
% RM instead of a dopvar, so the Q-form executives size Q the same way:
%
%   deg1 = max degree of any monomial in Q1, Q2, R.R0      (P is ignored)
%   deg2 = max degree in any single variable of R.R1, R.R2
%   Qdeg = [deg1, deg2, deg2-1]
%
% get_lpivar_degs has no container branch (it errors on anything but
% opvar/dopvar/opvar2d/dopvar2d), hence this. Block roles: an L_2 x L_2
% block's gamma cell 1 is R0 and cells 2,3 are R1, R2; an R x L_2 or
% L_2 x R block is Q1 or Q2. A monomial is present when its coefficient
% A + B'd is structurally nonzero, which is what a dpvar degmat lists.
% Kernel rows are (I_q kron ZL)' C (I_p kron ZR): row index = component
% outer, ZL monomial inner (sdopvar.m header), likewise columns and ZR.
% Cost O(nnz) of each block's coefficients; nothing in q beyond that.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deg1 = 0;   deg2 = 0;
for i = 1:size(Rm.C,1)
    for j = 1:size(Rm.C,2)
        B = Rm.C{i,j};
        if isempty(B),  continue,   end
        if numel(B.vars.out)>1 || numel(B.vars.in)>1
            error('cx_hinf_qdeg:dim','1-D containers only.')
        end
        shared = ~isempty(B.vars.out) && ~isempty(B.vars.in);   % L_2 x L_2
        for q = 1:numel(B.params.A)
            [dl,dr] = cell_degs(B,q);
            if isempty(dl),     continue,   end                 % zero cell
            if shared && q>1,   deg2 = max([deg2,dl,dr]);       % R1, R2
            else,               deg1 = max([deg1,dl,dr]);       % R0, Q1, Q2, P
            end
        end
    end
end
Qdeg = [deg1, deg2, deg2-1];    % -1: T'*Q adds one degree (stock comment)
end

function [dl,dr] = cell_degs(B,q)
% Max ZL and ZR monomial degree over the structurally nonzero coefficients
% of parameter cell q; empty when the cell is zero.
nL = prod([cellfun(@numel,B.ZL),1]);    nR = prod([cellfun(@numel,B.ZR),1]);
m = B.dims(1)*nL;   n = B.dims(2)*nR;
A = B.params.A{q};  nz = false(m*n,1);
if numel(A)==m*n,   nz = nz | full(A(:)~=0);    end
if isfield(B.params,'B') && numel(B.params.B)>=q
    Bq = B.params.B{q};
    if size(Bq,2)==m*n,     nz = nz | full(any(Bq~=0,1)).';     end
end
[r,c] = find(reshape(nz,m,n));
if isempty(r),  dl = [];    dr = [];    return,     end
eL = 0;     if ~isempty(B.ZL),  eL = B.ZL{1}(:);    end
eR = 0;     if ~isempty(B.ZR),  eR = B.ZR{1}(:);    end
dl = max(eL(mod(r-1,nL)+1));    dr = max(eR(mod(c-1,nR)+1));
end
