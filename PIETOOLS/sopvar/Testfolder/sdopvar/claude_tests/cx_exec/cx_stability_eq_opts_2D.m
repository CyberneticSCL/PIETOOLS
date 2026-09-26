function eq_opts = cx_stability_eq_opts_2D(Qm,eq_opts,ztol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EQ_OPTS = CX_STABILITY_EQ_OPTS_2D(QM,EQ_OPTS,ZTOL) is the container
% counterpart of  eq_opts = get_eq_opts_2D(clean_opvar(Qop,ztol),eq_opts,ztol)
% (PIETOOLS_stability_2D.m:174,191; get_eq_opts_2D.m:117-155), which prunes
% the negativity slack from the zero structure of Qop. No container version
% exists; eq_opts_sopvar is deliberately weaker (a nested, larger cone).
%
% Scope: QM a 1x1 'cdopvar' on L2[s1,s2] (a PIE whose state is purely 2-D
% distributed). Its one sdopvar block has a 3x3 gamma cell over (s1,s2) that
% is opvar2d's R22{i,j} cell for cell (opvar2d2sopvar copies the 3-way split
% whole), so the stock rules apply cell for cell:
%   R22{ii} zero (every |coefficient| < ztol, which is clean_opvar then ==0)
%       -> exclude(param_idcs(ii)), param_idcs = [8,9,10,11,13,14,12,15,16];
%   sep(3..6) from the SAME one-sided test the stock uses, all(diff <= ztol)
%       with no abs() (get_eq_opts_2D.m:128,135,143,151), kept so that the
%       two sides prune identically.
% EQ_OPTS stays in poslpivar_2d's vocabulary (exclude 1x16, sep 1x6), for
% settings2possopvar to translate.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3 || isempty(ztol),   ztol = 1e-12;   end
if ~isequal(size(Qm.C),[1 1]) || nnz(Qm.space_out)~=2 || numel(Qm.vars)~=2
    error('cx_stability:scope2D',['Only a state purely in L2[s1,s2] is transcribed; '...
          'the R, L2[x], L2[y] rules of get_eq_opts_2D.m:67-115 are not.'])
end
P = Qm.C{1};
exc = eq_opts.exclude(:)';  sep = eq_opts.sep(:)';
A = P.params.A;     B = P.params.B;
cz = @(k) cellnz(A{k},B{k},ztol);                   % true if cell k is nonzero
pidx = [8,9,10,11,13,14,12,15,16];                  % get_eq_opts_2D.m:119
for ii = 1:9
    if ~cz(ii),     exc(pidx(ii)) = 1;  end
end
le = @(i,j) onesided(A{i},A{j},B{i},B{j},ztol);     % all(R{i}-R{j} <= ztol)
if ~exc(9) && ~exc(10) && le(sub2ind([3 3],2,1),sub2ind([3 3],3,1)),  sep(3) = 1;  end
if ~exc(11) && ~exc(12) && le(sub2ind([3 3],1,2),sub2ind([3 3],1,3)), sep(4) = 1;  end
if ((~exc(13) && ~exc(14)) || (~exc(15) && ~exc(16))) ...
        && le(sub2ind([3 3],2,2),sub2ind([3 3],3,2)) && le(sub2ind([3 3],2,3),sub2ind([3 3],3,3))
    sep(5) = 1;
end
if ((~exc(13) && ~exc(15)) || (~exc(14) && ~exc(16))) ...
        && le(sub2ind([3 3],2,2),sub2ind([3 3],2,3)) && le(sub2ind([3 3],3,2),sub2ind([3 3],3,3))
    sep(6) = 1;
end
eq_opts.exclude = exc;  eq_opts.sep = sep;
end

function tf = cellnz(a,b,ztol)
tf = (~isempty(a) && any(abs(a(:))>=ztol)) || (~isempty(b) && any(abs(nonzeros(b))>=ztol));
end

function tf = onesided(ai,aj,bi,bj,ztol)
% Coefficientwise R{i}-R{j} <= ztol over the constant part and every
% decision variable row; ZL/ZR are shared by all cells of one block.
if isempty(ai), ai = 0; end;    if isempty(aj), aj = 0; end
if isempty(bi), bi = 0; end;    if isempty(bj), bj = 0; end
da = ai - aj;   db = bi - bj;
tf = all(da(:)<=ztol) && all(nonzeros(db)<=ztol);
end
