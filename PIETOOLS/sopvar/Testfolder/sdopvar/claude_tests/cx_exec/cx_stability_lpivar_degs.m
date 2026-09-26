function Qdeg = cx_stability_lpivar_degs(Pm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QDEG = CX_STABILITY_LPIVAR_DEGS(PM) is the container counterpart of the
% 1-D stock  Qdeg = get_lpivar_degs(Pop,Top)  (executives/utility_functions/
% get_lpivar_degs.m:60-81), which has no container overload (it errors at
% line 115 for anything but opvar/dopvar/opvar2d/dopvar2d). Used by the two
% Q-form executives to size Qop = lpivar(prog,Top.dim,Qdeg).
%
% Stock rule, read off the dopvar Pop (Top is not used):
%   deg1 = max degree of Q1, Q2, R.R0  (the multiplier / R^n-coupling terms)
%   deg2 = max single-variable exponent over R.R1, R.R2  (the kernels)
%   Qdeg = [deg1, deg2, deg2-1]
% Container reading, per block of PM and per gamma cell, over the monomials
% with a NONZERO coefficient (A or any row of B):
%   shared variable, multiplier cell (canonical: ZR degree 0)   -> deg1 <- ZL
%   shared variable, integral cell                              -> deg2 <- ZL,ZR
%   no shared variable (R^n <-> L2 coupling, Q1/Q2)             -> deg1 <- ZL,ZR
% A stock degmat may list a monomial whose coefficients are all zero; this
% reader cannot see such a monomial, so equality with the stock Qdeg is a
% measured fact per plant (cx_stability_check), not a guarantee.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deg1 = 0;   deg2 = 0;
for b = 1:numel(Pm.C)
    P = Pm.C{b};
    if isempty(P),  continue,   end
    n3 = numel(intersect(P.vars.out,P.vars.in));    % shared (S3) variables
    NL = prod([cellfun(@numel,P.ZL),1]);    NR = prod([cellfun(@numel,P.ZR),1]);
    nrow = P.dims(1)*NL;
    if isa(P,'sdopvar'),    ncell = numel(P.params.A);  else,   ncell = numel(P.params);    end
    for k = 1:ncell                                 % gamma cells, 3^n3 of them
        if isa(P,'sdopvar')
            A = P.params.A{k};  B = P.params.B{k};
            nz = [];
            if ~isempty(A),     nz = find(A(:)~=0);   end
            if ~isempty(B),     nz = union(nz,find(any(B~=0,1))');  end    % O(nnz(B))
        else
            nz = find(P.params{k}(:)~=0);           % fixed sopvar block: C itself
        end
        if isempty(nz),     continue,   end
        r = mod(nz-1,nrow)+1;       c = floor((nz-1)/nrow)+1;
        eL = maxexp(P.ZL,mod(r-1,NL)+1);            % row = (matrix row, monomial inner)
        eR = maxexp(P.ZR,mod(c-1,NR)+1);
        gam = 1;                                    % multi-index of cell k
        if n3>0,    g = cell(1,n3);     [g{:}] = ind2sub([3*ones(1,n3),1],k);   gam = [g{:}];   end
        if n3>0 && all(gam==1)                      % multiplier: R0
            deg1 = max(deg1,eL);
        elseif n3>0                                 % integral kernel: R1/R2
            deg2 = max([deg2,eL,eR]);
        else                                        % R^n coupling: Q1/Q2
            deg1 = max([deg1,eL,eR]);
        end
    end
end
Qdeg = [deg1,deg2,deg2-1];                          % get_lpivar_degs.m:81
end

function e = maxexp(Z,idx)
% Largest single-variable exponent over the monomials IDX of the Kronecker
% basis Z{1} (x) ... (x) Z{N} (last variable fastest, sdopvar.m header).
if isempty(Z),  e = 0;  return, end
n = cellfun(@numel,Z);  e = 0;
sub = cell(1,numel(n));
[sub{end:-1:1}] = ind2sub([fliplr(n),1],idx(:));  % trailing 1: ind2sub needs 2 dims
for i = 1:numel(n),     e = max([e;Z{i}(sub{i})]);    end
end
