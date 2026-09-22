function out = settings2possopvar(settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = SETTINGS2POSSOPVAR(SETTINGS) translates a PIETOOLS settings
% structure, as produced by 'settings_PIETOOLS_*' or
% 'settings_PIETOOLS_*_2D', into the degree and option arguments that
% 'possopvar' accepts.
%
% INPUT
% - settings: 'struct' from any of the shipped settings files. The 2D files
%             are recognized by the field 'is2D'.
%
% OUTPUT
% - out: 'struct' with one entry per positive operator the settings define:
%     out.LF, out.eq             struct('deg',...,'opts',...), ready to use
%                                as possopvar(prog,dim,vars,dom,X.deg,X.opts);
%     out.LF_psatz, out.eq_psatz 1 x k cell of the same, holding only the
%                                psatz terms the settings ENABLE through
%                                'LF_use_psatz' / 'eq_use_psatz';
%     out.eppos, out.epneg       passed through UNTRANSLATED, since the
%                                executive applies them and 'possopvar' does
%                                not. For a state with only a distributed
%                                component it is eppos(end) that applies,
%                                not eppos(1);
%     out.n3                     number of spatial variables;
%     out.report                 struct with 'lines' (what was translated)
%                                and 'lossy' (one entry per thing that could
%                                not be).
%
% EXACTNESS
% The degree conversion is EXACT. Both builders bound a component's basis by
% one cap per NONEMPTY SUBSET of its variables -- 'poslpivar_2d' through
% 'build_monoms', and 'possopvar' through 'deg.subset', which was added for
% this purpose. The whole cap array is carried across, permuted from
% poslpivar's per-direction ordering [ss1,(tt1),ss2,(tt2)] into possopvar's
% [theta_1,theta_2,s_1,s_2]; possopvar always carries all four slots and caps
% s_k at 0 in a multiplier direction, where poslpivar omits the slot.
%
% Before 'deg.subset' existed only the singleton and total caps could be
% carried, which gave a covering superset at 1.05x to 1.80x the basis across
% the six shipped 2D files. That is no longer the case, and 'out.report'
% states the monomial count rather than an inflation factor.
%
% Degrees are emitted as a per-block cell in the same order as opts.include,
% which is the order 'possopvar' assigns to its 'alpha_list', so every basis
% operator keeps its own caps instead of sharing one global bound.
%
% NOT REPRESENTABLE. Each of these is reported in out.report.lossy rather
% than silently mapped onto something else:
% - psatz = 2. 'poslpivar_2d' builds the ball multiplier for 2 and the box
%   for 1, while 'possopvar' tests options.psatz for truthiness only and so
%   would return the box. That is a different operator and a different SDP
%   row space, both measured, so such terms are DROPPED.
% - Components 1..7 of a 2D 'exclude', the R^n, L2[x] and L2[y] blocks.
%   'possopvar' maps L2 -> L2 only.
% - A 'sep' pattern that separates one direction's upper integral without its
%   lower one. 'possopvar' separates per direction, and that mixed regime is
%   also where 'poslpivar_2d' itself drops cross terms.
%
% See also POSSOPVAR, POSLPIVAR, POSLPIVAR_2D.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if ~isa(settings,'struct')
    error("Settings should be given as a 'struct'.")
end
is2D = isfield(settings,'is2D') && ~isempty(settings.is2D) && settings.is2D;

out = struct();
out.n3 = 1 + is2D;
out.report = struct('lines',{{}},'lossy',{{}});

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
out.use_sosineq = getfielddef(settings,'use_sosineq', ...
                    getfielddef(settings,'sosineq_on',[]));

end


%%
function out = conv_2d(s,out)
% The 2D settings carry LF_deg/LF_opts for the Lyapunov operator,
% eq_deg/eq_opts for the negativity operator, and one psatz variant per
% entry of LF_use_psatz / eq_use_psatz.

% Each call returns an updated 'out', so the field must be attached AFTER
% the struct comes back: '[out.LF,out] = f(...)' would set the field on the
% old copy and then discard it.
[blk,out] = block2d(getfielddef(s,'LF_deg',[]), ...
                    getfielddef(s,'LF_opts',[]),'LF',out);
out.LF = blk;
[blk,out] = block2d(getfielddef(s,'eq_deg',[]), ...
                    getfielddef(s,'eq_opts',[]),'eq',out);
out.eq = blk;

[lst,out] = psatz2d(s,'LF',out);   out.LF_psatz = lst;
[lst,out] = psatz2d(s,'eq',out);   out.eq_psatz = lst;

end


%%
function [lst,out] = psatz2d(s,tag,out)
% Only the psatz terms the settings switch on are translated: the files
% DEFINE two of them and then enable them separately.

lst = {};
use = getfielddef(s,[tag '_use_psatz'],[]);
use = use(use~=0);
for k = 1:numel(use)
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
    [b,out] = block2d(dgk,opk,sprintf('%s_psatz%d',tag,k),out);
    b.opts.psatz = 1;
    lst{end+1} = b;                                                     %#ok<AGROW>
end
end


%%
function [blk,out] = block2d(d,o,tag,out)
% Translate one 2D positive-operator specification.

blk = struct('deg',{{}},'opts',struct());
if isempty(d)
    out.report.lossy{end+1} = sprintf('%s: no degree field found.',tag);
    return
end

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
% Read the singleton and full-subset caps out of a 'build_monoms' cap array.
%
% The array holds 2^nvars caps over the variable list
%   [ss1, (tt1 if alpha(1)~=1), ss2, (tt2 if alpha(2)~=1)]
% reshaped from 2x2x...x2, so the cap on the subset with bits B set sits at
% linear index 1+sum(2.^(B-1)); a singleton k is at 1+2^(k-1) and the full
% set at 2^nvars. ss_k carries the INTEGRATION variable, possopvar's
% theta_k, and tt_k the OUTPUT variable, possopvar's s_k -- measured, and
% the opposite of what the names suggest.

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

na = count_subset(v,nv);
nb = na;    % exact by construction; asserted by the covering test
end


%%
function c = capat(v,p)
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
% Monomials satisfying every subset cap, by enumeration on the degree axis,
% which is small.
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
% The 1D settings carry dd1/options1 for the Lyapunov operator,
% dd12/options12 for its psatz term, and dd2/options2 and dd3/options3 for
% the negativity side.

% See the note in 'conv_2d': attach the field after the struct returns.
[blk,out] = block1d(getfielddef(s,'dd1',[]), ...
                    getfielddef(s,'options1',[]),'LF',out);
out.LF = blk;
[blk,out] = block1d(getfielddef(s,'dd2',[]), ...
                    getfielddef(s,'options2',[]),'eq',out);
out.eq = blk;
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

blk = struct('deg',{{}},'opts',struct());
if isempty(d)
    out.report.lossy{end+1} = sprintf('%s: no degree field found.',tag);
    return
end
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
blk.deg  = degs;
blk.opts = struct('include',incl,'sep',sepv, ...
                  'psatz',double(getfielddef(o,'psatz',0)~=0));
out.report.lines{end+1} = sprintf('%s: %d blocks translated.',tag,numel(degs));
end


%%
function v = getfielddef(s,f,dflt)
if isa(s,'struct') && isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end
