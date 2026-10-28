function [deg,terms,info] = get_lift_degs(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DEG,TERMS,INFO] = GET_LIFT_DEGS(P,OPTS) returns the degrees and the
% Positivstellensatz terms of a positive container operator Pop for the
% constraint P = Pop, Pop >= 0, with P self-adjoint, by the lift/weight
% rule of the separated form (sopvar_lift_notes Sec. 7). In each registry
% direction d, with i and j the exponents of s_d and s'_d in a monomial of
% a kernel of P (the upper cells, adjoints of the lower ones, not read),
%
%   Dmin_d = the largest min(i,j) of the lower cells;
%   Dmax_d = the largest max(i,j) of the lower cells and the largest
%            exponent in d of a block between two L2 spaces of which one
%            lacks d;
%   Mdeg_d = the largest i of the multiplier cells (-1 if none);
%   J_d    = the largest exponent in d of a block between an R^q space and
%            an L2 space (-1 if none).
%
% A positive operator of lift degree D_d and weight degree w_d reaches the
% lower kernels with min(i,j) <= D_d and max(i,j) <= D_d + 2 w_d + 1, the
% multiplier degree 2 w_d and, between R^q and L2, the degree
% D_d + w_d + wR + 1; reaching the support is necessary for P = Pop, not
% sufficient. The rule sets D_d = Dmin_d + OPTS.dD and w_d to the smallest
% weight degree that reaches the support plus OPTS.dw,
%
%   w_d = max(ceil((Dmax_d-D_d-1)/2), ceil(Mdeg_d/2), J_d-D_d-wR-1, 0)
%         + OPTS.dw.
%
% INPUT
% - P:     'cdopvar' or 'copvar' object, square (the same output and input
%          spaces), or an 'sdopvar' or 'sopvar' object from a space to
%          itself (vars.out = vars.in);
% - opts:  (optional) 'struct' with the optional fields
%            dD, dw        double, added to Dmin_d and to the weight rule
%                          (default 1, 1);
%            wR            double, the weight degree of an R^q space
%                          (default 0, the identity basis);
%            psatz         the terms: 'auto' (default; in one variable the
%                          codes [0 1] with offsets [0 1]; in nv ~= 1
%                          variables the codes [0, 3:2nv+2], the plain term
%                          and the 2nv faces, the faces at offset 1 when
%                          every w_d is at least 2, else at offset 0),
%                          'none' (code 0), 'product' (codes [0 1], offsets
%                          [0 1]), 'faces' (codes [0, 3:2nv+2], offsets 0),
%                          or a double row of term codes (0 plain, 1 box
%                          product, 2d+1 / 2d+2 lower / upper face of
%                          direction d, 2nv+2+d box quadratic of direction
%                          d);
%            psatz_offset  double, with a row of codes: the offset of each
%                          term, or one for all (default 0);
%            like          'cdopvar' or 'copvar' object, or a cell of them,
%                          with variables among those of P: their pooled
%                          support is read instead of that of P, and DEG is
%                          laid out over the spaces of P.
%
% OUTPUT
% - deg:   1 x M cell, M the number of spaces of P (1 for a block): for an
%          L2 space a 'struct' with fields int (1 x nv double, w) and mult
%          (1 x nv double, D in the directions of the space, 0 in the
%          others); for an R^q space a 'struct' with field int (1 x nv
%          double, wR);
% - terms: 'struct' with fields codes and offsets, double rows of one
%          length, the term codes and the amounts by which each lowers the
%          weight degree;
% - info:  'struct' with fields Dmin, Dmax, Mdeg, J, D, w (1 x nv double),
%          wR, dD, dw, nv, M (double), vars (1 x nv cellstr, the registry)
%          and isR (1 x M logical, true for an R^q space).
%
% NOTES
% The rule bounds the support of a fixed P. For P linear in a decision
% operator, such as A'*Q + Q'*A with Q free, the support of P is that of
% every value of the decision operator; to size from the positive operator
% of the LPI instead, pass it as OPTS.like.
%
% See also POSCOPVAR_DIRECT, LPI_INEQ_SOP, EQ_OPTS_SOPVAR, GET_LPIVAR_DEGS.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_lift_degs
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
% Initial coding MMP, 10/08/2026. The sizing routine of the separated form,
%                from the measured degree map (memory degree-map-proof-to-
%                code) and the heat benchmark (heatNd_tailor).
% MMP, 10/08/2026: N-D 'auto' terms: the 2N faces at w-1 instead of w when
%                every weight is at least 2. Measured on the 2-D heat
%                benchmark (R (1,1), Q (2,2): 12.336731 at nx 339209
%                against 1.04e6 for the faces at w) and on the 2-D coercive
%                KYP slack (the same gain at 200425 against 557845
%                variables); the Markov-Lukacs tensor set lost both tests
%                (sopvar_lift_notes Sec. 10). A face term at weight 0
%                certified nothing in 1-D DD, hence the floor at w = 1.
% MMP, 10/21/2026: opts.like may be a cell of operators, their supports
%                pooled: LPI_INEQ_SOP passes every positive operator of
%                'poscopvar' its target depends on.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body and headers
%                of the local functions; code unchanged. The deleted
%                header text was measurement results (the heat and KYP
%                benchmarks), the cost note and a second notes reference.
%                The help text of the 10/21/2026 entry is in the INPUT
%                entry of opts.like, without markers.

% Options; the defaults of dD, dw and wR
if nargin<2 || isempty(opts),   opts = struct();    end
dD = 1;     if isfield(opts,'dD') && ~isempty(opts.dD),     dD = opts.dD;   end
dw = 1;     if isfield(opts,'dw') && ~isempty(opts.dw),     dw = opts.dw;   end
wR = 0;     if isfield(opts,'wR') && ~isempty(opts.wR),     wR = opts.wR;   end

% % % Blocks and spaces
if isa(P,'sdopvar') || isa(P,'sopvar')
    % One block: its sorted variables, one space
    vars = union(P.vars.out,P.vars.in);     vars = reshape(sort(vars),1,[]);
    blocks = {P};   bk = 1;  bkp = 1;
    own = {reshape(find(ismember(vars,P.vars.out)),1,[])};
    if ~isequal(sort(P.vars.out),sort(P.vars.in))
        error("A single block must map a space to itself (vars.out = vars.in).")
    end
elseif isa(P,'cdopvar') || isa(P,'copvar')
    if ~isequal(P.space_out,P.space_in) || ~isequal(P.dim_out(:),P.dim_in(:))
        error("P must be square: the same spaces on both sides.")
    end
    % A container: its nonempty blocks with their space indices, and the
    % variables of each space
    vars = reshape(P.vars,1,[]);
    [Mo,Mi] = size(P.C);
    blocks = {};    bk = [];    bkp = [];
    for k = 1:Mo
        for kp = 1:Mi
            if ~isempty(P.C{k,kp}),  blocks{end+1} = P.C{k,kp};  bk(end+1) = k;  bkp(end+1) = kp;  end %#ok<AGROW>
        end
    end
    own = cell(1,Mo);
    for k = 1:Mo,   own{k} = reshape(find(P.space_out(k,:)),1,[]);   end
else
    error("P should be a 'cdopvar', 'copvar', 'sdopvar' or 'sopvar'.")
end
% Registry size, number of spaces, the R^q spaces
nv = numel(vars);   M = numel(own);
isR = cellfun(@isempty,own);
% % % opts.like: the support of another container, laid out over P's spaces
isRs = isR;
if isfield(opts,'like') && ~isempty(opts.like)
% % % BEGIN MMP, 10/21/2026. opts.like may be a cell of operators (the positive
% % % operators LPI_INEQ_SOP finds in its target): their blocks are pooled, so
% % % the support is the union, with each operator's space indices offset
% % % past the previous ones' for the R-space flags isRs. Replaced (initial
% % % coding 10/08/2026): the same lines for one operator, L = opts.like, its
% % % space indices k, kp and isRs = its R-space flags.
    Ls = opts.like;     if ~iscell(Ls),     Ls = {Ls};  end
    blocks = {};    bk = [];    bkp = [];   isRs = false(1,0);
    for j = 1:numel(Ls)
        L = Ls{j};      o = numel(isRs);
        if ~(isa(L,'cdopvar') || isa(L,'copvar'))
            error("opts.like should be a 'cdopvar' or 'copvar', or a cell of them.")
        end
        if ~all(ismember(reshape(L.vars,1,[]),vars))
            error("The variables of opts.like must be among those of P.")
        end
        [Mo,Mi] = size(L.C);
        for k = 1:Mo
            for kp = 1:Mi
                if ~isempty(L.C{k,kp}),  blocks{end+1} = L.C{k,kp};  bk(end+1) = o+k;  bkp(end+1) = o+kp;  end %#ok<AGROW>
            end
        end
        isRs = [isRs, reshape(~any(L.space_out,2),1,[])]; %#ok<AGROW>
    end
% % % END MMP, 10/21/2026
end

% % % The support, per direction
Dmin = zeros(1,nv);     Dmax = zeros(1,nv);     Mdeg = -ones(1,nv);     J = -ones(1,nv);
for b = 1:numel(blocks)
    % Block b: its variables and their registry indices, its monomial
    % grids, dimensions and shared variables S3
    B = blocks{b};
    vout = reshape(B.vars.out,1,[]);    vin = reshape(B.vars.in,1,[]);
    [~,po] = ismember(vout,vars);       [~,pi_] = ismember(vin,vars);
    EL = multiindex_grid(B.ZL,'first_slowest');     NL = size(EL,1);
    ER = multiindex_grid(B.ZR,'first_slowest');     NR = size(ER,1);
    m = B.dims(1);  n = B.dims(2);
    S3 = intersect(vout,vin);   n3 = numel(S3);
    [~,p3o] = ismember(S3,vout);    [~,p3i] = ismember(S3,vin);
    [~,d3] = ismember(S3,vars);
    % True for a block between an R^q space and an L2 space
    otherR = (isRs(bk(b)) && ~isRs(bkp(b))) || (~isRs(bk(b)) && isRs(bkp(b)));
    ncell = numel(cell_params(B,1,'count'));
    for q = 1:ncell
        % The monomial pairs of the nonzeros of cell q, and its multi-index
        [li,ei] = cell_support(B,q,m,n,NL,NR);
        if isempty(li),     continue,   end
        gam = gamma_of_cell(q,n3);
        % shared directions: lower cells give (i,j), multiplier cells i
        for t = 1:n3
            d = d3(t);  i = EL(li,p3o(t));   j = ER(ei,p3i(t));
            if gam(t)==1
                Mdeg(d) = max([Mdeg(d); i(:)]);
            elseif gam(t)==2
                Dmin(d) = max([Dmin(d); min(i(:),j(:))]);
                Dmax(d) = max([Dmax(d); max(i(:),j(:))]);
            end
        end
        % one-sided directions: the exponent of the side that has them
        for t = 1:numel(vout)
            if ismember(vout{t},S3),    continue,   end
            d = find(strcmp(vars,vout{t}));  i = EL(li,t);
            if otherR,  J(d) = max([J(d); i(:)]);  else,   Dmax(d) = max([Dmax(d); i(:)]);  end
        end
        for t = 1:numel(vin)
            if ismember(vin{t},S3),     continue,   end
            d = find(strcmp(vars,vin{t}));   j = ER(ei,t);
            if otherR,  J(d) = max([J(d); j(:)]);  else,   Dmax(d) = max([Dmax(d); j(:)]);  end
        end
    end
end

% % % The rule
D = Dmin + dD;
w = max([ceil((Dmax-D-1)/2); ceil(Mdeg/2); J-D-wR-1; zeros(1,nv)],[],1) + dw;
% Per space: w and the lift D in its own directions, or wR on R^q
deg = cell(1,M);
for k = 1:M
    if isR(k)
        deg{k} = struct('int',wR*ones(1,nv));
    else
        mult = zeros(1,nv);     mult(own{k}) = D(own{k});
        deg{k} = struct('int',w,'mult',mult);
    end
end

% % % The terms
ps = 'auto';
if isfield(opts,'psatz') && ~isempty(opts.psatz),   ps = opts.psatz;    end
% A row of codes with its offsets, or a named set of terms
if isnumeric(ps)
    codes = reshape(ps,1,[]);   offs = zeros(size(codes));
    if isfield(opts,'psatz_offset') && ~isempty(opts.psatz_offset)
        offs = reshape(opts.psatz_offset,1,[]);
        if isscalar(offs),  offs = repmat(offs,size(codes));    end
    end
else
    switch lower(char(ps))
        case 'auto'
            if nv==1,   codes = [0 1];  offs = [0 1];
%           else,       codes = [0, 3:2*nv+2];  offs = zeros(size(codes));  % MMP, 10/08/2026 (was)
            else        % faces one below the plain term when every weight allows it, else at w
                codes = [0, 3:2*nv+2];  offs = [0, double(min(w)>=2)*ones(1,2*nv)]; % MMP, 10/08/2026
            end
        case 'none',    codes = 0;              offs = 0;
        case 'product', codes = [0 1];          offs = [0 1];
        case 'faces',   codes = [0, 3:2*nv+2];  offs = zeros(1,2*nv+1);
        otherwise
            error("'psatz' should be 'auto', 'none', 'product', 'faces' or a row of term codes.")
    end
end
% The outputs
terms = struct('codes',codes,'offsets',offs);
info = struct('Dmin',Dmin,'Dmax',Dmax,'Mdeg',Mdeg,'J',J,'D',D,'w',w,'wR',wR,'dD',dD,'dw',dw,...
              'nv',nv,'M',M,'vars',{vars},'isR',isR);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = cell_params(B,~,~)
% C = CELL_PARAMS(B,~,~) returns the parameter cells of the block B, for
% counting them: params.A of an 'sdopvar', params of a 'sopvar'. The
% second and third inputs are not used.
% INPUT  B: 'sdopvar' or 'sopvar' object.
% OUTPUT c: cell array, 3 x ... x 3 over the shared variables of B.
if isa(B,'sdopvar'),    c = B.params.A;     else,   c = B.params;  end
end


function [li,ei] = cell_support(B,q,m,n,NL,NR)
% [LI,EI] = CELL_SUPPORT(B,Q,M,N,NL,NR) returns the left and right monomial
% indices of every nonzero of parameter cell Q of the block B, over the
% constant part and every decision-variable row.
% INPUT  B:      'sdopvar' or 'sopvar' object;
%        q:      double, the linear index of the parameter cell;
%        m, n:   double, the output and input dimensions of B;
%        NL, NR: double, the numbers of monomials of B.ZL and B.ZR.
% OUTPUT li, ei: double columns, the indices (1..NL and 1..NR) of the
%                left and right monomials of the nonzeros; [] for a zero
%                cell.
if isa(B,'sdopvar')
    % Coefficient positions with a nonzero constant or decision row
    A = B.params.A{q};  Bq = B.params.B{q};
    cols = find(A);     cols = cols(:);
    if ~isempty(Bq),    [~,cB] = find(Bq);   cols = unique([cols; cB(:)]);  end
    if isempty(cols),   li = [];  ei = [];   return,     end
    [ridx,cidx] = ind2sub([m*NL,n*NR],cols);
else
    % The nonzeros of the coefficient matrix
    Cq = B.params{q};
    if isempty(Cq) || nnz(Cq)==0,   li = [];  ei = [];   return,     end
    [ridx,cidx] = find(Cq);
end
% Rows are (matrix row, monomial), columns (matrix column, monomial)
li = mod(ridx-1,NL)+1;  ei = mod(cidx-1,NR)+1;
end
