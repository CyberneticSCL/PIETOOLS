function out = settings2possopvar(settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = SETTINGS2POSSOPVAR(SETTINGS) translates the positive operators that
% a 1-D or 2-D LPI settings struct defines (the Lyapunov operator, the
% operator of the equality constraint, and their psatz terms) into the
% degree and option arguments of 'possopvar' on an L2 space: for each basis
% operator of the settings that is kept, its multi-index (opts.include: 1
% multiplier, 2 lower, 3 upper, 4 full integral, per direction) and its
% monomial degree caps (deg).
%   1-D  The bases Z1, Z2, Z3 of the degrees {d1,d2,d3} have the
%        multi-indices 1, 2, 3; with sep set, Z2 has 4 and Z3 is omitted.
%        Z1 has deg.int = deg.joint = d1, deg.mult = 0; Z2 has deg.int,
%        deg.mult, deg.joint = d2(1), d2(2), d2(3); Z3 likewise from d3.
%   2-D  The basis Z2^{ij} on L2[x,y] has the multi-index [i j]; in a
%        direction that sep separates, 2 becomes 4 and the bases with 3
%        there are omitted. Its array of subset caps d2{i,j} is carried to
%        deg.subset over the variables [theta_1,theta_2,s_1,s_2] (s_k
%        capped at 0 in a multiplier direction), with deg.int, deg.mult
%        the singleton caps and deg.joint the cap of all its variables.
% What the translation cannot carry is listed in OUT.report.lossy.
%
% INPUT
% - settings: 'struct' of LPI settings (settings_PIETOOLS_* or
%             settings_PIETOOLS_*_2D):
%   1-D (no field is2D, or is2D = 0): dd1, options1 (Lyapunov operator);
%             dd12, options12, override1 (its psatz term, translated when
%             override1 is 0 or absent); dd2, options2 (equality
%             operator); dd3, options3, override2 (its psatz term, when
%             override2 is 0 or absent and dd3 exists). Each dd is a 'cell'
%             {d1,d2,d3}, d2 and d3 1x3 or scalar, or a 'double' scalar d
%             ({d,[d d d],[d d d]}); each options a 'struct' with fields
%             exclude (1x4: the R^n, multiplier, lower and upper integral
%             bases), sep and psatz, each optional;
%   2-D (is2D = 1): LF_deg, LF_opts (Lyapunov operator); eq_deg, eq_opts
%             (equality operator); LF_use_psatz, LF_deg_psatz,
%             LF_opts_psatz and eq_use_psatz, eq_deg_psatz, eq_opts_psatz
%             (their psatz terms). Each *_deg is a 'struct' with field d2,
%             the 3x3 'cell' whose entry {i,j} is the array of degree caps
%             of Z2^{ij}, one cap per subset of its variables; each *_opts
%             a 'struct' with fields exclude (1x16), sep (1x6) and psatz,
%             each optional; each *_use_psatz a 'double' array, whose k-th
%             nonzero entry selects entry k of the 'cell' arrays
%             *_deg_psatz and *_opts_psatz;
%   both: (optional) eppos, epneg, and use_sosineq or sosineq_on.
%
% OUTPUT
% - out:  'struct' with fields
%           LF, eq        'struct' with fields deg (1 x N 'cell', one
%                         'struct' per kept basis operator in the order of
%                         opts.include, with fields int, mult and joint,
%                         int and mult 1x2 in 2-D, and in 2-D subset
%                         (1x16)) and opts ('struct' with fields include
%                         (N x 1 in 1-D, N x 2 in 2-D), sep ('logical',
%                         1x2 in 2-D) and psatz (1 when the settings'
%                         psatz is nonzero, else 0)), the arguments of
%                         possopvar(prog,dim,vars,dom,X.deg,X.opts);
%           LF_psatz,     1 x K 'cell' of the same, one per psatz term
%           eq_psatz      translated, each with opts.psatz = 1;
%           eppos, epneg  settings.eppos and settings.epneg unchanged ([]
%                         when absent);
%           use_sosineq   settings.use_sosineq, else settings.sosineq_on,
%                         else [];
%           n3            'double', the number of spatial variables, 1 or 2;
%           report        'struct' with fields lines (what was translated,
%                         and the parts with no L2 -> L2 counterpart that
%                         were left out) and lossy (one entry per part not
%                         translated), each a 'cell' of 'char'.
%
% NOTES
% OUT.report.lossy lists: a missing degree field; in 2-D, a psatz term
% with psatz = 2, which is omitted; sep(2+k) ~= sep(4+k), direction k then
% not separated; a missing d2{i,j}, whose basis operator is omitted.
% OUT.report.lines lists the parts left out: in 1-D the R^n basis that
% exclude(1) = 0 keeps; in 2-D the R^n, L2[x] and L2[y] bases that
% exclude(1:7) keeps, and a nonzero sep(1:2). A d2{i,j} array without one
% cap per subset of the variables of Z2^{ij} is an error.
%
% See also POSSOPVAR, POSLPIVAR, POSLPIVAR_2D.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - settings2possopvar
%
% Copyright (C)2026 PIETOOLS Team
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
% Initial coding MMP, 09/21/2026: the positivity builders take degree
%                  specifications that differ in kind, so a settings file
%                  written for 'poslpivar_2d' cannot be handed to
%                  'possopvar'. Translates the parts that map exactly and
%                  reports the parts that cannot, rather than quietly
%                  building a different cone.
% MMP, 09/21/2026: 'block1d' had deg.int and deg.mult SWAPPED, and carried a
%                  comment asserting that 'possopvar's header had the
%                  correspondence backwards. Measured from 'poslpivar'
%                  itself that the header is right and 'block1d' was wrong:
%                  the first slot of d{2} caps the variable that then
%                  becomes the integration variable. 'caps2spec', the 2D
%                  path in this same file, already had it right, so the two
%                  halves of the file disagreed. Latent until now -- nothing
%                  outside this file calls it -- which is why a green test
%                  suite did not catch it.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), the licence block in the toolbox form
%                  (title line added), step comments in the body; code
%                  unchanged. The deleted header text was the history and
%                  rationale of the degree conversion (EXACTNESS) and the
%                  reasons for the lossy cases; its statement that
%                  exclude(1:7) is reported in report.lossy is replaced by
%                  report.lines, where the code reports it.

% SETTINGS: a 'struct', 2-D when its field is2D is set
if ~isa(settings,'struct')
    error("Settings should be given as a 'struct'.")
end
is2D = isfield(settings,'is2D') && ~isempty(settings.is2D) && settings.is2D;

% The output: number of spatial variables, empty report
out = struct();
out.n3 = 1 + is2D;
out.report = struct('lines',{{}},'lossy',{{}});

% The positive operators, by the 2-D or 1-D translation
if is2D
    out = conv_2d(settings,out);
else
    out = conv_1d(settings,out);
end

% eppos and epneg belong to the executive, not to 'possopvar'. They are
% passed through unchanged; the caller must apply eppos to the component its
% state actually has, which for a purely distributed state is the LAST entry.
out.eppos = getfielddef(settings,'eppos',[]);
out.epneg = getfielddef(settings,'epneg',[]);
% use_sosineq (2-D name) or sosineq_on (1-D name), unchanged
out.use_sosineq = getfielddef(settings,'use_sosineq', ...
                    getfielddef(settings,'sosineq_on',[]));

end


%%
function out = conv_2d(s,out)
% OUT = CONV_2D(S,OUT) adds to OUT the translation of the 2-D settings S:
% the Lyapunov operator (LF_deg, LF_opts), the operator of the equality
% constraint (eq_deg, eq_opts) and their psatz terms.
% INPUT  s:   'struct', 2-D LPI settings;
%        out: 'struct', the output being built (field report).
% OUTPUT out: OUT with the fields LF, eq, LF_psatz, eq_psatz added and the
%             report extended.

% Each call returns an updated 'out', so the field must be attached AFTER
% the struct comes back: '[out.LF,out] = f(...)' would set the field on the
% old copy and then discard it.
[blk,out] = block2d(getfielddef(s,'LF_deg',[]), ...
                    getfielddef(s,'LF_opts',[]),'LF',out);
out.LF = blk;
[blk,out] = block2d(getfielddef(s,'eq_deg',[]), ...
                    getfielddef(s,'eq_opts',[]),'eq',out);
out.eq = blk;

% The psatz terms that are switched on
[lst,out] = psatz2d(s,'LF',out);   out.LF_psatz = lst;
[lst,out] = psatz2d(s,'eq',out);   out.eq_psatz = lst;

end


%%
function [lst,out] = psatz2d(s,tag,out)
% [LST,OUT] = PSATZ2D(S,TAG,OUT) translates the psatz terms of the block TAG
% of the 2-D settings S that are switched on: for the k-th nonzero entry of
% S.<TAG>_use_psatz, entry k of S.<TAG>_deg_psatz and S.<TAG>_opts_psatz,
% omitted and reported lossy when its psatz (that of the options, else the
% use entry) is 2.
% INPUT  s:   'struct', 2-D LPI settings;
%        tag: 'char', 'LF' or 'eq';
%        out: 'struct', the output being built (field report).
% OUTPUT lst: 1 x K 'cell' of 'struct' with fields deg and opts, each with
%             opts.psatz = 1;
%        out: OUT with the report extended.

% The nonzero use entries, one term each
lst = {};
use = getfielddef(s,[tag '_use_psatz'],[]);
use = use(use~=0);
for k = 1:numel(use)
    % Entry k of the degrees and options ([] if absent), and its psatz
    dg = getfielddef(s,[tag '_deg_psatz'],{});
    op = getfielddef(s,[tag '_opts_psatz'],{});
    if numel(dg)>=k, dgk = dg{k}; else, dgk = []; end
    if numel(op)>=k, opk = op{k}; else, opk = []; end
    pv = getfielddef(opk,'psatz',use(k));
    if pv==2
        out.report.lossy{end+1} = sprintf(['%s psatz term %d requests ' ...
            'psatz=2, the ball multiplier. possopvar tests options.psatz ' ...
            'for truthiness only and would return the box multiplier of ' ...
            'psatz=1, a different operator. Term DROPPED.'],tag,k);
        continue
    end
    % The translation, with psatz 1
    [b,out] = block2d(dgk,opk,sprintf('%s_psatz%d',tag,k),out);
    b.opts.psatz = 1;
    lst{end+1} = b;                                                     %#ok<AGROW>
end
end


%%
function [blk,out] = block2d(d,o,tag,out)
% [BLK,OUT] = BLOCK2D(D,O,TAG,OUT) translates one 2-D positive-operator
% specification, degrees D and options O, into 'possopvar' degrees and
% options on an L2[x,y] space.
% INPUT  d:   'struct' with field d2, the 3x3 'cell' of cap arrays ([]:
%             reported lossy, BLK empty);
%        o:   'struct' with fields exclude (1x16), sep (1x6), psatz, each
%             optional;
%        tag: 'char', the name of the specification in the report;
%        out: 'struct', the output being built (field report).
% OUTPUT blk: 'struct' with fields deg (1 x N 'cell') and opts (fields
%             include (N x 2), sep (1x2 'logical'), psatz (0 or 1));
%        out: OUT with the report extended.

% No degrees: an empty translation, reported
blk = struct('deg',{{}},'opts',struct());
if isempty(d)
    out.report.lossy{end+1} = sprintf('%s: no degree field found.',tag);
    return
end

% exclude and sep, padded to 16 and 6 entries
exc = getfielddef(o,'exclude',zeros(1,16));
sep = getfielddef(o,'sep',zeros(1,6));
exc = [reshape(exc,1,[]), zeros(1,max(0,16-numel(exc)))];
sep = [reshape(sep,1,[]), zeros(1,max(0,6-numel(sep)))];

% Components 1..7 are the R^n, L2[x] and L2[y] blocks, which 'possopvar'
% cannot build -- it maps L2 -> L2 only. This is NOT a loss for a purely
% distributed state: 'poslpivar_2d' forces excludeL(1:7)=1 itself whenever
% n0=nx=ny=0, so the settings value is moot there. Reported as a note, and
% only as a real gap if the caller does have those components.
if any(~exc(1:7))
    out.report.lines{end+1} = sprintf(['%s leaves components %s of ' ...
        'exclude on. poslpivar_2d forces them off for a purely ' ...
        'distributed state, so nothing is lost there; they are omitted, ' ...
        'and would be a genuine gap only for a state with R^n, L2[x] or ' ...
        'L2[y] components.'],tag,mat2str(find(~exc(1:7))));
end

% sep(3)/sep(5) act on direction 1 and sep(4)/sep(6) on direction 2. A
% direction is separable only when both of its flags agree; the mixed case
% is the one in which poslpivar_2d itself discards cross terms.
sepv = false(1,2);
for k = 1:2
    lo = sep(2+k);   hi = sep(4+k);
    if lo~=hi
        out.report.lossy{end+1} = sprintf(['%s sets sep(%d)=%d and ' ...
            'sep(%d)=%d for direction %d. possopvar separates per ' ...
            'direction, so that mixed regime cannot be requested; ' ...
            'treated as NOT separable.'],tag,2+k,lo,4+k,hi,k);
    else
        sepv(k) = logical(lo);
    end
end
if any(sep(1:2))
    out.report.lines{end+1} = sprintf(['%s: sep(1:2) act on Rxx/Ryy, ' ...
        'empty for an L2[x,y] state; ignored.'],tag);
end

% Entry 7+r of 'exclude' is the component whose multi-index is map(r,:),
% which is also possopvar's alpha. Note 14 and 15: includeL(14) is Z2ba,
% cell R22{3,2}, and includeL(15) is Z2ab, cell R22{2,3}.
map = [1 1; 2 1; 3 1; 1 2; 1 3; 2 2; 3 2; 2 3; 3 3];
incl = zeros(0,2);   degs = {};
d2 = getfielddef(d,'d2',{});
nkeep = 0;   nmono_a = 0;   nmono_b = 0;
for r = 1:size(map,1)
    if exc(7+r), continue, end
    a = map(r,:);
    % With direction k separable, poslpivar keeps only the lower integral
    % and forces the upper one out, so alpha_k = 2 becomes possopvar's
    % full-domain index 4 and alpha_k = 3 is the discarded duplicate.
    skip = false;   an = a;
    for k = 1:2
        if sepv(k)
            if a(k)==3
                skip = true;
            elseif a(k)==2
                an(k) = 4;
            end
        end
    end
    if skip, continue, end
    if size(d2,1)>=a(1) && size(d2,2)>=a(2)
        A = d2{a(1),a(2)};
    else
        A = [];
    end
    if isempty(A)
        out.report.lossy{end+1} = sprintf('%s: d2{%d,%d} is missing.', ...
            tag,a(1),a(2));
        continue
    end
    [spec,na,nb] = caps2spec(A,a);
    nkeep = nkeep+1;   nmono_a = nmono_a+na;   nmono_b = nmono_b+nb;
    incl(end+1,:) = an;                                                 %#ok<AGROW>
    degs{end+1} = spec;                                                 %#ok<AGROW>
end

% The output (psatz 1 for any nonzero psatz); the monomial count in the
% report
blk.deg  = degs;
blk.opts = struct('include',incl,'sep',sepv, ...
                  'psatz',double(getfielddef(o,'psatz',0)~=0));
if nmono_a > 0
    out.report.lines{end+1} = sprintf(['%s: %d blocks, %d monomials, ' ...
        'EXACT -- the whole subset cap array is carried across.'], ...
        tag,nkeep,nmono_a);
end
end


%%
function [spec,na,nb] = caps2spec(A,alpha)
% [SPEC,NA,NB] = CAPS2SPEC(A,ALPHA) returns the 'possopvar' degree
% specification of the 2-D basis operator of multi-index ALPHA from its
% 'build_monoms' cap array A.
% INPUT  A:      'double' array of 2^nvars caps (layout below);
%        alpha:  1x2 'double', the multi-index (1 multiplier, 2 or 3
%                integral, per direction).
% OUTPUT spec:   'struct' with fields int and mult (1x2, the singleton caps
%                of theta_k and s_k), joint (the cap of all the variables)
%                and subset (1x16, the caps over [theta_1,theta_2,s_1,s_2]);
%        na, nb: 'double', the number of monomials the caps admit (nb = na).
%
% The array holds 2^nvars caps over the variable list
%   [ss1, (tt1 if alpha(1)~=1), ss2, (tt2 if alpha(2)~=1)]
% reshaped from 2x2x...x2, so the cap on the subset with bits B set sits at
% linear index 1+sum(2.^(B-1)); a singleton k is at 1+2^(k-1) and the full
% set at 2^nvars. ss_k carries the INTEGRATION variable, possopvar's
% theta_k, and tt_k the OUTPUT variable, possopvar's s_k -- measured, and
% the opposite of what the names suggest.

% Slots of ss_k (p_int) and tt_k (p_mult, 0 in a multiplier direction) in
% the variable list; one cap per subset of the nv variables
v = A(:);
p_int = [0 0];   p_mult = [0 0];
k = 0;
for dir = 1:2
    k = k+1;   p_int(dir) = k;
    if alpha(dir)~=1
        k = k+1;   p_mult(dir) = k;
    end
end
nv = k;
if numel(v) ~= 2^nv
    error("caps2spec: expected %d caps for alpha %s, got %d.", ...
          2^nv,mat2str(alpha),numel(v));
end

% Singleton caps of theta_k (int) and s_k (mult), the cap of all (joint)
spec = struct();
spec.int   = [v(1+2^(p_int(1)-1)), v(1+2^(p_int(2)-1))];
spec.mult  = [capat(v,p_mult(1)), capat(v,p_mult(2))];
spec.joint = v(2^nv);

% The singleton and total caps alone give only a covering superset, since
% build_monoms also prunes with every intermediate subset. Carry the whole
% cap array across so the conversion is exact. possopvar orders its basis
% variables [theta_1,theta_2,s_1,s_2] and always has all four slots, capping
% s_k at 0 in a multiplier direction, whereas poslpivar omits that slot
% entirely -- so the arrays are permuted, not copied.
slotmap = [p_int(1), p_int(2), p_mult(1), p_mult(2)];
spec.subset = zeros(1,16);
for i = 2:16
    mem = find(bitget(i-1,1:4)>0);
    src = slotmap(mem);
    src = src(src>0);
    if isempty(src)
        % Every member is an s_k that poslpivar does not carry, and that
        % possopvar caps at 0 anyway; bounding their sum by 0 is exact.
        spec.subset(i) = 0;
    else
        spec.subset(i) = v(1 + sum(2.^(src-1)));
    end
end

% Monomial count of the basis
na = count_subset(v,nv);
nb = na;    % exact by construction; asserted by the covering test
end


%%
function c = capat(v,p)
% C = CAPAT(V,P) returns the singleton cap of slot P of the cap vector V,
% or 0 for P = 0.
% INPUT  v: 'double' column of caps;
%        p: 'double', the slot (0: none).
% OUTPUT c: 'double'.
% A multiplier direction has no output-variable slot, so degree 0 in s_k --
% which is also what 'possopvar' enforces internally for such a direction.
if p==0
    c = 0;
else
    c = v(1+2^(p-1));
end
end


%%
function n = count_subset(v,nv)
% N = COUNT_SUBSET(V,NV) returns the number of monomials in NV variables
% whose degree in every subset of the variables is at most the cap of that
% subset in V.
% INPUT  v:  'double' column of 2^nv caps, the subset with bits B at index
%            1+sum(2.^(B-1));
%        nv: 'double', the number of variables.
% OUTPUT n:  'double'.
% Enumeration of the exponent grid 0:max(v) in nv variables, the degree
% axis being small; the exponents within every subset cap are kept
mx = max(v);
g = cell(1,nv);
[g{:}] = ndgrid(0:mx);
E = reshape(cat(nv+1,g{:}),[],nv);
keep = true(size(E,1),1);
for i = 2:2^nv
    b = bitget(i-1,1:nv)>0;
    keep = keep & (sum(E(:,b),2) <= v(i));
end
n = sum(keep);
end


%%
function out = conv_1d(s,out)
% OUT = CONV_1D(S,OUT) adds to OUT the translation of the 1-D settings S:
% dd1, options1 (Lyapunov operator); dd12, options12 (its psatz term, when
% override1 is 0 or absent); dd2, options2 (equality operator); dd3,
% options3 (its psatz term, when override2 is 0 or absent and dd3 exists).
% INPUT  s:   'struct', 1-D LPI settings;
%        out: 'struct', the output being built (field report).
% OUTPUT out: OUT with the fields LF, eq, LF_psatz, eq_psatz added and the
%             report extended.

% See the note in 'conv_2d': attach the field after the struct returns.
[blk,out] = block1d(getfielddef(s,'dd1',[]), ...
                    getfielddef(s,'options1',[]),'LF',out);
out.LF = blk;
[blk,out] = block1d(getfielddef(s,'dd2',[]), ...
                    getfielddef(s,'options2',[]),'eq',out);
out.eq = blk;
% The psatz terms unless overridden, each with psatz 1
out.LF_psatz = {};
out.eq_psatz = {};
if ~getfielddef(s,'override1',0)
    [b,out] = block1d(getfielddef(s,'dd12',[]), ...
                      getfielddef(s,'options12',[]),'LF_psatz1',out);
    b.opts.psatz = 1;
    out.LF_psatz = {b};
end
if ~getfielddef(s,'override2',0) && isfield(s,'dd3')
    [b,out] = block1d(getfielddef(s,'dd3',[]), ...
                      getfielddef(s,'options3',[]),'eq_psatz1',out);
    b.opts.psatz = 1;
    out.eq_psatz = {b};
end
end


%%
function [blk,out] = block1d(d,o,tag,out)
% [BLK,OUT] = BLOCK1D(D,O,TAG,OUT) translates one 1-D positive-operator
% specification, degrees D and options O, into 'possopvar' degrees and
% options on an L2[s] space.
% INPUT  d:   'cell' {d1,d2,d3} (d2, d3 1x3 or scalar) or a 'double'
%             scalar d, read as {d,[d d d],[d d d]} ([]: reported lossy,
%             BLK empty);
%        o:   'struct' with fields exclude (1x4), sep, psatz, each optional;
%        tag: 'char', the name of the specification in the report;
%        out: 'struct', the output being built (field report).
% OUTPUT blk: 'struct' with fields deg (1 x N 'cell') and opts (fields
%             include (N x 1), sep ('logical'), psatz (0 or 1));
%        out: OUT with the report extended.
%
% 'poslpivar' takes d = {d1,d2,d3}: d1 the multiplier degree, and d2/d3 the
% lower/upper integral triples, mapping to (deg.int, deg.mult, deg.joint) in
% that order -- which is what 'sopquadvar's header says.
%
% MEASURED from 'poslpivar' itself, correcting the claim that was here      % MMP, 09/21/2026
% before: poslpivar.m:321 builds Z2sth over [var1;var2] with d{2}(1)        % MMP, 09/21/2026
% capping var1, then lines 343-344 substitute var1 -> sss (the integration  % MMP, 09/21/2026
% variable eta) and var2 -> var1 (the spatial variable s). So the FIRST     % MMP, 09/21/2026
% slot becomes theta and the second becomes s, i.e. d{2}(1) is deg.int.     % MMP, 09/21/2026
% 'caps2spec' in this same file already says so, and                        % MMP, 09/21/2026
% 'test_possopvar_vs_poslpivar' passes 5 of 5 using int = d2(1).            % MMP, 09/21/2026

% No degrees: an empty translation, reported
blk = struct('deg',{{}},'opts',struct());
if isempty(d)
    out.report.lossy{end+1} = sprintf('%s: no degree field found.',tag);
    return
end
% A 'double' d to a cell; exclude padded to 4 entries; sep
if ~iscell(d)
    d = {d,[d d d],[d d d]};
end
exc = getfielddef(o,'exclude',zeros(1,4));
exc = [reshape(exc,1,[]), zeros(1,max(0,4-numel(exc)))];
sepv = logical(getfielddef(o,'sep',0));

% exclude(1) is the real-valued block, which an L2 -> L2 operator has no
% counterpart for; 2..4 are multiplier, lower integral and upper integral.
if ~exc(1)
    out.report.lines{end+1} = sprintf(['%s: exclude(1) is the R^n block, ' ...
        'absent from an L2 -> L2 operator.'],tag);
end
% One basis operator per kept multi-index: 1 multiplier, 2 lower integral
% (4 under sep), 3 upper integral (omitted under sep)
incl = zeros(0,1);   degs = {};
for a = 1:3
    if exc(1+a), continue, end
    if sepv && a==3, continue, end      % the upper integral is the duplicate
    an = a;
    if sepv && a==2, an = 4; end
    if a==1
        t = d{1};
        spec = struct('int',t(1),'mult',0,'joint',t(1));
    else
        t = d{a};
        if isscalar(t), t = [t t t]; end
%       spec = struct('int',t(2),'mult',t(1),'joint',t(3));                  % MMP, 09/21/2026 (was)
        spec = struct('int',t(1),'mult',t(2),'joint',t(3));                  % MMP, 09/21/2026
    end
    incl(end+1,1) = an;                                                 %#ok<AGROW>
    degs{end+1} = spec;                                                 %#ok<AGROW>
end
% The output (psatz 1 for any nonzero psatz); the count in the report
blk.deg  = degs;
blk.opts = struct('include',incl,'sep',sepv, ...
                  'psatz',double(getfielddef(o,'psatz',0)~=0));
out.report.lines{end+1} = sprintf('%s: %d blocks translated.',tag,numel(degs));
end


%%
function v = getfielddef(s,f,dflt)
% V = GETFIELDDEF(S,F,DFLT) returns the field F of S, or DFLT when S is not
% a 'struct' or has no nonempty field F.
% INPUT  s:    any value;
%        f:    'char', the field name;
%        dflt: the default.
% OUTPUT v:    the value of the field, or DFLT.
if isa(s,'struct') && isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
