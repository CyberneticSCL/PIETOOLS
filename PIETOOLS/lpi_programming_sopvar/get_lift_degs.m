function [deg,terms,info] = get_lift_degs(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DEG,TERMS,INFO] = GET_LIFT_DEGS(P,OPTS) sizes a positive container
% operator to cancel the self-adjoint operator P, P = Pop with Pop >= 0, by
% the lift/weight rules of the separated form (sopvar_lift_notes Sec. 7):
% per registry direction d the lower-kernel support {s^i s'^j} of P gives
%
%   Dmin_d = max min(i_d,j_d),   Dmax_d = max max(i_d,j_d),
%   Mdeg_d = the degree of the multiplier cells in d,
%
% and a positive operator of lift degree D_d and weight degree w_d reaches
% the lower kernels {min(i,j) <= D_d, max(i,j) <= D_d + 2 w_d + 1} and the
% multiplier degree 2 w_d (Lemma 7.1 there; necessary, not sufficient). The
% lift is set to Dmin_d + OPTS.dD and the weight to
%
%   w_d = max(ceil((Dmax_d - D_d - 1)/2), ceil(Mdeg_d/2), J_d - D_d - wR - 1, 0) + OPTS.dw,
%
% J_d the one-sided reach needed by a block between an R^q space and an
% L2 space (its kernel in that direction has degree at most D_d + w_d + wR
% + 1). Measured on the heat benchmark (heatNd_tailor, 10/08/2026): the
% certificate needs D = Dmin + 1 and, for the binding operator, w one above
% the support rule at every P degree, which the defaults dD = dw = 1 give;
% in 1-D the pair plain term at w plus product term at w-1 is complete and
% half the cone of two faces at w; in 2-D the face terms at w are needed.
%
% INPUT
% - P:     'cdopvar' or 'copvar', square (output spaces = input spaces), or
%          an 'sdopvar' / 'sopvar' block on one space;
% - opts:  (optional) struct with fields
%   dD, dw        added to Dmin_d and to the weight rule (default 1, 1);
%   wR            weight degree of R^q spaces (default 0, the identity
%                 basis of poslpivar);
%   psatz         'auto' (default: 1-D codes [0 1] with offsets [0 1];
%                 N-D the plain term at w and the 2N faces at w-1 when
%                 every weight is at least 2, at w otherwise), 'none',
%                 'product' ([0 1], [0 1]), 'faces' (the 2N faces at w),
%                 or a row of codes;
%   psatz_offset  with a row of codes, per term (default 0);
%   like          a 'cdopvar' or 'copvar' whose support is read in place of
%                 P's, the specification laid out over P's spaces: for a
%                 target linear in a free operator, pass the fixed positive
%                 operator of the LPI (the KYP slack from R; see NOTES);
%                 its variables must be among P's;
%
% OUTPUT
% - deg:   1 x M cell of structs for 'poscopvar_direct' / 'poscopvar':
%          struct('int',w,'mult',D) with w, D 1 x nv rows for an L2 space
%          (mult 0 in the directions the space lacks), struct('int',wR) for
%          an R^q space;
% - terms: struct with fields codes (row) and offsets (row);
% - info:  struct with the per-direction Dmin, Dmax, Mdeg, J, the chosen D
%          and w, nv, M, and the registry vars.
%
% NOTES
% The support is read from every block: in a shared direction from the
% cells lower in that direction (upper cells are their adjoints) and the
% multiplier cells; in a direction only one side has, from the exponents of
% that side, folded into Dmax (both sides carry theta_d to degree w_d) or
% into J when the other side is an R^q space. Cost: one pass over the
% nonzeros of P, nothing dense in the decision variables.
% The rule is a bound for a FIXED target. A target linear in a free
% decision operator, such as the KYP block A'*Q + Q'*A with Q an lpivar,
% carries the support of the whole Q family under composition (lower
% kernel to degree 9, Dmin = 4, for the 1-D plant io1 of cx_plant), not of
% the optimal Q, and the rule then sizes to the family: (D,w) = (5,3) where
% (2,2) is the measured minimum (hinf_tailor_1d, 10/08/2026). Size such a
% slack from the fixed positive operator of the LPI instead, opts.like: on
% R the rule gives (2,2) at both stock levels.
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

if nargin<2 || isempty(opts),   opts = struct();    end
dD = 1;     if isfield(opts,'dD') && ~isempty(opts.dD),     dD = opts.dD;   end
dw = 1;     if isfield(opts,'dw') && ~isempty(opts.dw),     dw = opts.dw;   end
wR = 0;     if isfield(opts,'wR') && ~isempty(opts.wR),     wR = opts.wR;   end

% % % Blocks and spaces
if isa(P,'sdopvar') || isa(P,'sopvar')
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
nv = numel(vars);   M = numel(own);
isR = cellfun(@isempty,own);
% % % opts.like: the support of another container, laid out over P's spaces
isRs = isR;
if isfield(opts,'like') && ~isempty(opts.like)
    L = opts.like;
    if ~(isa(L,'cdopvar') || isa(L,'copvar'))
        error("opts.like should be a 'cdopvar' or 'copvar'.")
    end
    if ~all(ismember(reshape(L.vars,1,[]),vars))
        error("The variables of opts.like must be among those of P.")
    end
    [Mo,Mi] = size(L.C);
    blocks = {};    bk = [];    bkp = [];
    for k = 1:Mo
        for kp = 1:Mi
            if ~isempty(L.C{k,kp}),  blocks{end+1} = L.C{k,kp};  bk(end+1) = k;  bkp(end+1) = kp;  end %#ok<AGROW>
        end
    end
    isRs = reshape(~any(L.space_out,2),1,[]);
end

% % % The support, per direction
Dmin = zeros(1,nv);     Dmax = zeros(1,nv);     Mdeg = -ones(1,nv);     J = -ones(1,nv);
for b = 1:numel(blocks)
    B = blocks{b};
    vout = reshape(B.vars.out,1,[]);    vin = reshape(B.vars.in,1,[]);
    [~,po] = ismember(vout,vars);       [~,pi_] = ismember(vin,vars);
    EL = multiindex_grid(B.ZL,'first_slowest');     NL = size(EL,1);
    ER = multiindex_grid(B.ZR,'first_slowest');     NR = size(ER,1);
    m = B.dims(1);  n = B.dims(2);
    S3 = intersect(vout,vin);   n3 = numel(S3);
    [~,p3o] = ismember(S3,vout);    [~,p3i] = ismember(S3,vin);
    [~,d3] = ismember(S3,vars);
    otherR = (isRs(bk(b)) && ~isRs(bkp(b))) || (~isRs(bk(b)) && isRs(bkp(b)));
    ncell = numel(cell_params(B,1,'count'));
    for q = 1:ncell
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
terms = struct('codes',codes,'offsets',offs);
info = struct('Dmin',Dmin,'Dmax',Dmax,'Mdeg',Mdeg,'J',J,'D',D,'w',w,'wR',wR,'dD',dD,'dw',dw,...
              'nv',nv,'M',M,'vars',{vars},'isR',isR);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = cell_params(B,~,~)
% The parameter cells of a block, as a cell array, for counting.
if isa(B,'sdopvar'),    c = B.params.A;     else,   c = B.params;  end
end


function [li,ei] = cell_support(B,q,m,n,NL,NR)
% Left and right monomial indices of every nonzero of parameter cell q,
% over the constant part and every decision-variable row.
if isa(B,'sdopvar')
    A = B.params.A{q};  Bq = B.params.B{q};
    cols = find(A);     cols = cols(:);
    if ~isempty(Bq),    [~,cB] = find(Bq);   cols = unique([cols; cB(:)]);  end
    if isempty(cols),   li = [];  ei = [];   return,     end
    [ridx,cidx] = ind2sub([m*NL,n*NR],cols);
else
    Cq = B.params{q};
    if isempty(Cq) || nnz(Cq)==0,   li = [];  ei = [];   return,     end
    [ridx,cidx] = find(Cq);
end
li = mod(ridx-1,NL)+1;  ei = mod(cidx-1,NR)+1;
end
