%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_POSMOPVAR_VS_POSLPIVAR_2D cross-checks 'posmopvar' against
% 'poslpivar_2d' on the FOUR spaces R^n x L_2[x] x L_2[y] x L_2[x,y].
%
% This is the configuration 'mopquadvar' was written for, and a far harder
% test than the two-space 1D case in 'test_posmopvar_vs_poslpivar': there are
% 16 basis operators and 256 block pairs, of which 49 share no spatial
% variable, 126 share one and 81 share two. Every combination of the four
% direction classes -- shared, output-only, input-only, and in neither space
% -- occurs, so all three of the one-sided elimination branches are exercised
% together with the semiseparable integral.
%
% The correspondence of basis operators, read off 'poslpivar_2d's own Z*u
% display:
%
%   u0                                    -> space R^n, the identity
%   Zxo(x)*ux(x)                          -> space L_2[x],   alpha_x = 1
%   int_a^x Zxa(x,tt)*ux(tt) dtt          -> space L_2[x],   alpha_x = 2
%   int_x^b Zxb(x,tt)*ux(tt) dtt          -> space L_2[x],   alpha_x = 3
%   ... the mirror three for L_2[y] ...
%   the nine Z2** over L_2[x,y]           -> space L_2[x,y], (alpha_x,alpha_y)
%
% and of degrees: in each of poslpivar_2d's bases the FIRST slot of a
% direction is the variable that becomes the INTEGRATION variable and the
% second becomes the spatial one, so d(1) -> deg.int and d(2) -> deg.mult.
% That is measured in 'poslpivar.m' at lines 321 and 343-344 for the 1D case
% and at 718-722 for the 2D one, and it is the same mapping
% 'test_posmopvar_vs_poslpivar' uses.
%
% HOW THE COMPARISON IS MADE. Both routines parameterize a LINEAR family of
% operators, so each is characterized by the span of that family, and the two
% families are compared by sampling. No correspondence between decision
% variable names, monomial orderings or basis orderings is assumed anywhere.
%
% The reference side is brought into the container by evaluating
% 'poslpivar_2d's 'dopvar2d' at a value of its decision variables and then
% calling 'opvar2d2mopvar'. That converter is independent of everything under
% test and is separately covered by 'Test_mopvar_converters' (46 checks), and
% using it removes the one real hazard in a hand-written comparison: opvar2d
% gives Rx2 and R2x the identical variable signature (s1, s1_dum, s2) with s2
% the INPUT y in one and the OUTPUT y in the other, so 14 of the 36 leaf
% kernels would misalign if the components were read off by name. The
% converter knows the roles; this test does not have to.
%
% The samples are taken at RANDOM decision values rather than at unit
% vectors. A unit vector leaves whole rows of the container unpopulated, and
% 'mopvar' cannot then derive those rows' space and dimension -- measured, it
% errors. Random values give a dense container, and the span of enough random
% samples is the whole family almost surely. Sampling continues until the
% rank has not increased for a run of consecutive samples, so the family is
% demonstrably spanned rather than assumed to be, and the saturation point is
% reported.
%
% A PRECONDITION IS CHECKED FIRST: the two builders must agree on the total
% Gram size and on the decision variable count. Those depend only on the
% monomial bases, so a mismatch localizes the failure to the degree mapping
% above rather than to the elimination code -- which is the failure this test
% is most likely to suffer.
%
% NOT COVERED. All three cases keep 'poslpivar_2d's d2 at degree zero, so the
% subset-cap mechanism is inert: 'build_exponent_grid' takes one cap per
% variable SUBSET, and reproducing poslpivar_2d's default d2 exactly needs
% that array carried across with its variable order permuted. That mapping is
% the remaining gap, and 'mopquadvar' accepts deg.subset only for a space
% holding every registry variable -- which in this configuration is
% L_2[x,y] alone, the only space that needs it.
%
% MMP, 09/21/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(23);

for nm = {'posmopvar','mopquadvar','poslpivar_2d','opvar2d2mopvar'}
    if numel(which(nm{1},'-all'))~=1
        error('test_posmopvar_vs_poslpivar_2d: %s resolves to %d files.',...
              nm{1},numel(which(nm{1},'-all')));
    end
end

w = warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');

% Component names of an 'opvar2d', in the grid order the converter documents.
NAMES = {'R00','R0x','R0y','R02';
         'Rx0','Rxx','Rxy','Rx2';
         'Ry0','Ryx','Ryy','Ry2';
         'R20','R2x','R2y','R22'};

ngrid = 5;              % Chebyshev points per axis, see 'grid_points'
nsat  = 12;             % consecutive samples with no rank gain = saturated

% {label, dims, dom, dx1, dx2, dy1, dy2}
%   dx1/dy1 the multiplier degree, dx2/dy2 the [int,mult,joint] triple.
%
% Each case earns its place, checked by mutation:
%   * swapping the two middle SPACES in the posmopvar call is caught by every
%     case (22 of 36 block/cell groups differ at the first);
%   * swapping deg.int with deg.mult -- the mapping this whole comparison
%     hinges on -- survives the first three and is caught only by the LAST,
%     because degree 0 and the symmetric triple [1,1,2] both make that swap a
%     no-op. An asymmetric case is therefore not optional here.
cases = {
    {'deg 0, dims 1 1 1 1',   [1;1;1;1], [0,1;0,1], 0, [0,0,0], 0, [0,0,0]}
    {'deg 0, dims 2 1 1 1',   [2;1;1;1], [0,1;0,1], 0, [0,0,0], 0, [0,0,0]}
    {'dx,dy deg 1',           [1;1;1;1], [0,1;0,1], 1, [1,1,2], 1, [1,1,2]}
    {'dx,dy asymmetric',      [1;1;1;1], [0,2;-1,1],1, [1,0,1], 0, [0,0,0]}
    };

npass = 0;
for ic = 1:size(cases,1)
    [lbl,dims,dom,dx1,dx2,dy1,dy2] = deal(cases{ic,1}{:});

    % % % Reference: poslpivar_2d over the same four spaces
    d = struct();
    d.dx = {dx1, dx2(:), dx2(:)};
    d.dy = {dy1, dy2(:).', dy2(:).'};
    d.d2 = cell(3,3);
    for i = 1:3
        for j = 1:3
            d.d2{i,j} = zeros(2+2*(i>1),2+2*(j>1));
        end
    end
    opts = struct('psatz',0,'exclude',zeros(1,16),'sep',zeros(1,6));

    prog = lpiprogram(polynomial({'s1','s2'}),polynomial({'s1_dum','s2_dum'}),dom);
    [~,Pop_pl,Qmat] = poslpivar_2d(prog,dims.',d,opts);
    dv_pl = collect_dvars(Pop_pl,NAMES);

    % % % posmopvar at the matching degrees
    spaces = {{},{'s1'},{'s2'},{'s1','s2'}};
    deg = pl2d_degrees(dx1,dx2,dy1,dy2);
    prog2 = lpiprogram(polynomial({'s1','s2'}),polynomial({'s1_dum','s2_dum'}),dom);
    [~,Pop_pm,Qcell,basis_list] = posmopvar(prog2,dims,spaces,dom,deg);

    % % % Precondition: the bases must have the same total size
    G_pm = 0;
    for c = 1:size(basis_list,1)
        G_pm = G_pm + size(Qcell{c,c},1);
    end
    if size(Qmat,1)~=G_pm
        error(['test_posmopvar_vs_poslpivar_2d: the Gram sides differ for case '...
               '''%s'' -- poslpivar_2d %d, posmopvar %d. The degree mapping is '...
               'wrong; the elimination code is not implicated.'],...
               lbl,size(Qmat,1),G_pm);
    end
    if numel(dv_pl)~=numel(Pop_pm.Zd)
        error(['test_posmopvar_vs_poslpivar_2d: the decision variable counts '...
               'differ for case ''%s'' -- poslpivar_2d %d, posmopvar %d, at '...
               'equal Gram size %d.'],lbl,numel(dv_pl),numel(Pop_pm.Zd),G_pm);
    end

    % % % Sample both families to rank saturation
    pts = grid_points(dom,ngrid);
    [Phi_pl,n_pl,seg] = sample_family(@(dval) opvar2d2mopvar(eval_dopvar2d(Pop_pl,dv_pl,dval,NAMES)), ...
                                      numel(dv_pl),[],pts,nsat);
    [Phi_pm,n_pm] = sample_family(@(dval) Pop_pm,numel(Pop_pm.Zd),'dec',pts,nsat);

    % % % Per block and cell first, so that a genuine difference is localized
    % instead of reported as "the families differ". Each group's rows come
    % from one kernel and are therefore comparably scaled, which makes these
    % comparisons better conditioned than the overall one.
    nbad = 0;
    for s = 1:numel(seg)
        r = seg{s}.rows;
        [b1,b2,bb] = compare_spans(Phi_pl(r,:),Phi_pm(r,:));
        if b1~=b2 || bb~=b1
            nbad = nbad+1;
            fprintf(['    block (%d,%d) cell %d differs: rank %d vs %d, '...
                     'joint %d\n'],seg{s}.i,seg{s}.j,seg{s}.q,b1,b2,bb);
        end
    end
    if nbad>0
        error(['test_posmopvar_vs_poslpivar_2d: %d of %d block/cell groups '...
               'differ for case ''%s''.'],nbad,numel(seg),lbl);
    end

    % % % Compare the spans, reporting the gap at the rank threshold rather
    % than the rank integer alone: a spectrum with no gap there would make the
    % rank an arbitrary choice of tolerance.
    [r1,r2,rb,gap] = compare_spans(Phi_pl,Phi_pm);
    if r1~=r2 || rb~=r1
        error(['test_posmopvar_vs_poslpivar_2d: the two families differ for '...
               'case ''%s'': rank(poslpivar_2d) = %d, rank(posmopvar) = %d, '...
               'rank(joint) = %d.'],lbl,r1,r2,rb);
    end
    if gap < 1e3
        error(['test_posmopvar_vs_poslpivar_2d: the singular values have no '...
               'gap at the rank threshold for case ''%s'' (ratio %.3g), so the '...
               'rank comparison is not meaningful and the verdict must not be '...
               'believed either way.'],lbl,gap);
    end

    fprintf(['  passed: %-22s Gram %3d  ndec %4d  family dim %4d  '...
             'samples %3d/%3d  sv gap %.1e\n'],lbl,G_pm,numel(dv_pl),rb,n_pl,n_pm,gap);
    npass = npass+1;
end

warning(w);
fprintf('test_posmopvar_vs_poslpivar_2d passed (%d of %d cases).\n',npass,size(cases,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = pl2d_degrees(dx1,dx2,dy1,dy2)
% 'poslpivar_2d' degrees to 'posmopvar' degrees, over the sorted registry
% (s1,s2) = (x,y). Each space's caps are given over the whole registry and
% mopquadvar uses the entries for that space's own variables, so a direction
% a space does not have is capped at 0 -- which is what keeps the L_2[x] and
% L_2[y] bases free of the foreign integration variable that poslpivar_2d's
% Zxa and Zya do not carry.
%
% d2 is kept at degree zero by the caller, so the L_2[x,y] specs are all
% zero and the subset-cap mechanism is inert. See the NOT COVERED note.

sR  = struct('int',[0,0]);                       % the identity: no theta
sx  = { struct('int',[dx1,0],   'mult',[0,0],       'joint',dx1)
        struct('int',[dx2(1),0],'mult',[dx2(2),0],  'joint',dx2(3))
        struct('int',[dx2(1),0],'mult',[dx2(2),0],  'joint',dx2(3)) };
sy  = { struct('int',[0,dy1],   'mult',[0,0],       'joint',dy1)
        struct('int',[0,dy2(1)],'mult',[0,dy2(2)],  'joint',dy2(3))
        struct('int',[0,dy2(1)],'mult',[0,dy2(2)],  'joint',dy2(3)) };
s2  = repmat({struct('int',[0,0],'mult',[0,0],'joint',0)},1,9);

deg = { sR, sx(:).', sy(:).', s2 };

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi,nsamp,seg] = sample_family(build,ndec,mode,pts,nsat)
% Sample a linear family of containers until its span stops growing.
%
% 'build' returns a container given a decision value vector; for the
% posmopvar side the container is fixed and the value is substituted when its
% kernels are read, which is what mode 'dec' selects.

Phi = [];   nsamp = 0;  last_rank = 0;  flat = 0;  seg = {};
Cfix = [];
if strcmp(mode,'dec'),  Cfix = build([]);   end
while nsamp < ndec && flat < nsat
    nsamp = nsamp+1;
    dval = randn(ndec,1);
    if strcmp(mode,'dec')
        [col,seg] = sample_container(Cfix,dval,pts);
    else
        [col,seg] = sample_container(build(dval),zeros(0,1),pts);
    end
    if isempty(Phi),    Phi = zeros(numel(col),ndec);   end
    if numel(col)~=size(Phi,1)
        error(['test_posmopvar_vs_poslpivar_2d: sample %d has %d entries '...
               'against %d for the first; the block or cell layout is not '...
               'stable across samples.'],nsamp,numel(col),size(Phi,1));
    end
    Phi(:,nsamp) = col;                                                     %#ok<AGROW>
    r = rank_of(Phi(:,1:nsamp));
    if r > last_rank,   last_rank = r;   flat = 0;   else,   flat = flat+1;   end
end
Phi = Phi(:,1:nsamp);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [col,seg] = sample_container(Mop,dval,pts)
% Every kernel of every block, sampled on the grid and stacked. The block and
% cell layout is determined by the spaces alone -- a block's cell count is
% 3^(shared variables) -- so it is the same on both sides, which is what lets
% the two stacked columns be compared entry by entry.

allvars = {'s1','s2','s1_dum','s2_dum'};
M = size(Mop.C,1);
col = [];   seg = {};
for i = 1:M
    for j = 1:M
        B = Mop.C{i,j};
        if isempty(B)
            error(['test_posmopvar_vs_poslpivar_2d: block (%d,%d) is '...
                   'structurally zero; a dense sample was expected.'],i,j);
        end
        Kc = pi_blk_kernels(B,dval);
        S3 = intersect(B.vars.in,B.vars.out);
        if numel(Kc)~=3^numel(S3)
            error(['test_posmopvar_vs_poslpivar_2d: block (%d,%d) has %d '...
                   'parameter cells for %d shared variables.'],...
                   i,j,numel(Kc),numel(S3));
        end
        for q = 1:numel(Kc)
            Kq = Kc{q};
            % In EACH shared direction where gamma_k = 1 the factor
            % delta(s_k - s_k') identifies the two variables, so the kernel is
            % only determined up to moving degree between s_k and s_k'.
            % Evaluating on the diagonal in that direction removes the
            % ambiguity. 'sdopvar' pins it by keeping all the degree on the
            % left (its canonical multiplier form); an 'opvar2d' from
            % 'poslpivar_2d' has no such convention, so without this the two
            % sample differently while representing the same operator.
            % Measured: doing it only when EVERY direction is a multiplier --
            % which is all the 1D test needs, since there is one direction --
            % left each family with 2 dimensions the other lacked at
            % 'dx,dy deg 1', because block (4,4) has cells that are a
            % multiplier in one direction and an integral in the other.
            if numel(S3)>0
                gam = pi_gamma_index(q,numel(S3));
                for t = 1:numel(S3)
                    if gam(t)==1
                        Kq = subs(Kq,polynomial({[S3{t} '_dum']}),polynomial(S3(t)));
                    end
                end
            end
            V = pi_poly_grid(Kq,allvars,pts);
            n0 = numel(col);
            col = [col;V(:)];                                               %#ok<AGROW>
            seg{end+1} = struct('i',i,'j',j,'q',q, ...                      %#ok<AGROW>
                                'rows',(n0+1:numel(col)).');
        end
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pf = eval_dopvar2d(Pd,dvars,dval,NAMES)
% A 'dopvar2d' evaluated at one value of its decision variables, as a fixed
% 'opvar2d'. Four of the sixteen components are cell-valued -- the 3-way
% alpha split of the PI operator in each direction -- so each entry is
% evaluated in turn.

Pf = opvar2d();
Pf.I = Pd.I;        Pf.var1 = Pd.var1;      Pf.var2 = Pd.var2;
Pf.dim = Pd.dim;
for k = 1:numel(NAMES)
    v = Pd.(NAMES{k});
    if iscell(v)
        c = cell(size(v));
        for q = 1:numel(v)
            c{q} = pi_dpvar_eval(v{q},dvars,dval);
        end
        Pf.(NAMES{k}) = c;
    else
        Pf.(NAMES{k}) = pi_dpvar_eval(v,dvars,dval);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dv = collect_dvars(Pd,NAMES)
% Every decision variable appearing anywhere in a 'dopvar2d'. Collecting
% from a subset of the components would silently drop whole generators: the
% Gram block coupling the R^n basis to itself appears only in R00.

dv = {};
for k = 1:numel(NAMES)
    v = Pd.(NAMES{k});
    if iscell(v)
        for q = 1:numel(v)
            if isa(v{q},'dpvar')
                dv = [dv; cellstr(string(v{q}.dvarname(:)))];               %#ok<AGROW>
            end
        end
    elseif isa(v,'dpvar')
        dv = [dv; cellstr(string(v.dvarname(:)))];                          %#ok<AGROW>
    end
end
dv = unique(dv);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = rank_of(A)

A = normalize_cols(A);
r = rank(A,rank_tol(A));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tol = rank_tol(A)
% A hundred times MATLAB's own default for 'rank'.
%
% The obvious choice, a fixed relative factor such as 1e-8 times the number
% of rows, is wrong here and dangerously so: these matrices have tens of
% thousands of rows, one per kernel sample, and that factor then cuts at
% around 5e-3. Measured at 'dx,dy deg 1': the joint spectrum runs down to
% 2.1e-4 and then falls to 5e-15, a clean gap of 4.1e10 at index 155, but the
% 1e-8 rule truncated at 144 and reported the two families as DIFFERING with
% ranks 149 and 150. The row count is not a measure of numerical rank.

s = max(norm(A),realmin);
tol = 1e2*max(size(A))*eps(s);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r1,r2,rboth,gap] = compare_spans(A,B)
% Ranks of the two column spans and of their union, plus the ratio of the
% last kept singular value to the first dropped one in the joint matrix. A
% rank comparison is only as good as that gap.
%
% The ROWS are equilibrated first, using the norms of the JOINT matrix so the
% same scaling is applied to both sides. Without it the comparison is
% dominated by whichever blocks happen to have the largest kernels: the rows
% come from 36 different kernels whose magnitudes differ by orders of
% magnitude, and a direction that is small in that mixture falls below the
% single rank tolerance inconsistently between the two sides. Measured at
% 'dx,dy deg 1': unequilibrated, the two families reported ranks 149 and 149
% with a joint rank of 151 while EVERY individual block and cell compared
% equal -- an artefact of scaling, not a difference of families.

J = [A,B];
rn = vecnorm(J,2,2);        rn(rn<eps) = 1;
A = A./rn;      B = B./rn;

A = normalize_cols(A);      B = normalize_cols(B);
r1 = rank_of(A);            r2 = rank_of(B);
J = normalize_cols([A,B]);
rboth = rank(J,rank_tol(J));
sv = svd(J);
if rboth>=numel(sv)
    gap = inf;
else
    gap = sv(rboth)/max(sv(rboth+1),eps);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = normalize_cols(A)

nrm = vecnorm(A);
nrm(nrm<eps) = 1;
A = A./nrm;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pts = grid_points(dom,ngrid)
% Tensor grid in (s1,s2,s1_dum,s2_dum), at CHEBYSHEV points rather than
% equispaced ones.
%
% The kernels reach a higher degree than the bases: a Volterra integral
% substitutes a spatial variable into a limit and raises the degree at every
% composition. Equispaced sampling of such a polynomial is badly conditioned,
% and the symptom is a sampled matrix with no clean gap at the rank
% threshold, which makes the rank an artefact of the tolerance. Measured:
% with 4 equispaced points per axis the 'dx,dy deg 1' case had a gap ratio of
% 11; Chebyshev points put it far above that at the same cost.

g1 = cheb(dom(1,:),ngrid);
g2 = cheb(dom(2,:),ngrid);
[A,B,C,D] = ndgrid(g1,g2,g1,g2);
pts = [A(:),B(:),C(:),D(:)];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = cheb(lim,n)

k = (1:n).';
g = (lim(1)+lim(2))/2 + (lim(2)-lim(1))/2*cos((2*k-1)*pi/(2*n));

end
