%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_COPVAR_SILENT_FIXES checks the fixes of 09/29/2026 to four silently
% wrong answers of the containers 'copvar'/'cdopvar' and their blocks
% 'sopvar'/'sdopvar':
%
%   (1) P(...) on a container returned the whole container -> error
%       (@copvar/subsref, @cdopvar/subsref);
%   (2) Q(1,2) = P built an object array, all four classes  -> error
%       (subsasgn of each class);
%   (3) P.' on a container returned P unchanged             -> P'
%       (@copvar/transpose, @cdopvar/transpose);
%   (4) A == B was false for one operator stored in two block partitions
%       -> decided on the common refinement ('eq_copvar'),
%   (5) a comma list of block properties ({P.params{:}}, f(P.ZL{:}),
%       [P.params.A{:}]) kept its first element on a block, and raised
%       MATLAB:TooManyOutputs through a container -> the whole list
%       (@sopvar/@sdopvar subsref and numArgumentsFromSubscript),
%   (6) cat(1|2|3,A,B) built object arrays, all four classes -> vertcat,
%       horzcat, error (cat.m of each class),
%   (7) eq with a zero-count block row or column depended on that row's
%       space, so it was not transitive -> no effect ('eq_copvar'),
%   (8) [A;B;C] of sopvar recursed with horzcat -> vertcat
%       (@sopvar/vertcat),
%   (9) P(I,J).<tail> on a sopvar dropped the tail, and P(I,J).f{:} kept
%       one list element on both block classes -> tail applied, list
%       refused (@sopvar/@sdopvar subsref),
%
% and that nothing else moved: dot and brace reads (chained, comma lists,
% 'end', protected properties, error identifiers), dot assignment, and
% sopvar/sdopvar paren slicing.
%
% Every check is against SEMANTICS or against the builtin, never a round
% trip (CLAUDE.md S4):
% - dot/brace reads are compared with builtin('subsref',...) called from
%   this script, which is the indexing the classes had without an overload;
% - a slice B(I,J) of a block is compared with rows/columns of the block's
%   kernels ('pi_blk_kernels');
% - P.' is checked by <P.'x, y> = <x, P y>, both sides by quadrature;
% - every eq verdict between two partitions is checked against an
%   independent evaluation ('check_eq'): both containers are applied to the
%   same random polynomial test functions, one per GLOBAL component, at the
%   same points, one set per space by NAME, by Gauss quadrature straight
%   from the kernel definition in the 'sopvar.m' header ('heatNd_apply', no
%   class method in the loop), decision blocks at random values drawn once
%   per NAME. The second partition is built by block CONCATENATION
%   ([X,Y], [X;Y], blkdiag, opvar2copvar) while 'eq' cuts blocks by
%   SLICING, and the quadrature decides which answer is right, so a pair of
%   mutually inverse index errors in concatenation and slicing cannot pass;
% - every comma list is compared with the same list formed natively from
%   the underlying cell or struct, read by builtin('subsref',...), and the
%   count numArgumentsFromSubscript gives with the builtin count;
% - cat(1|2,...) is compared with [A;B] / [A,B].
%
% Shapes: 1-D R x L2^2[s1] (the reported blkdiag vs opvar2copvar case),
% 2-D R x L2[s] x L2[s,t] (fixed and decision), and 3-variable L2[s,t,u]
% (with an R row in the eq cases).
%
% Initial coding MMP, 09/29/2026; (5)-(9) added the same day.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(20260929);
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
assert(exist('heatNd_apply','file')==2,'heatNd_apply (PIETOOLS_demos/sopvar_demos) is not on the path');
nchk = 0;
NQ = 6;     % Gauss nodes per direction: exact to degree 11, and the
            % integrands here have degree <= 7 per variable (kernels <= 2,
            % test functions 2, one integration)

% ------------------------------------------------------------ operators
pvar s1 s1_dum
Pop = rand_opvar([1 1;2 2],2,s1,s1_dum,[0 1]);     % R^1 x L2^2[s1] -> same
Pc  = opvar2copvar(Pop);
Rop = mat2opvar(3,[1 1;0 0],[s1,s1_dum],[0 1]);    % R^1 -> R^1
Rc  = opvar2copvar(Rop);
sp2 = struct('out',{{ {}, {'s'}, {'s','t'} }},'in',{{ {}, {'s'}, {'s','t'} }});
dm2 = struct('out',[1;2;1],'in',[1;2;1]);
P2  = rand_copvar(sp2,dm2,[0,1],1,0.7);
D2  = rand_cdopvar(sp2,dm2,[0,1],1,3,0.7);
sp3 = struct('out',{{ {'s','t','u'} }},'in',{{ {'s','t','u'} }});
P3  = rand_copvar(sp3,struct('out',2,'in',2),[0,1],1,0.7);
CONT = {Pc,P2,D2,P3};
NAME = {'1-D copvar','2-D copvar','2-D cdopvar','3-var copvar'};
BLK  = {Pc.C{2,2}, D2.C{3,3}};              % a sopvar and an sdopvar block
dmap = draw_dvals({D2});

% =================================================== (1) paren reference
for k = 1:numel(CONT)
    P = CONT{k};    id = [class(P) ':parenIndex'];
    expect_error(@() P(1),id,'P.C{i,j}');
    expect_error(@() P(1,1),id,'P.C{i,j}');
    expect_error(@() P(:),id,'P.C{i,j}');
    expect_error(@() P(1,:),id,'P.C{i,j}');
    expect_error(@() P(:,:),id,'P.C{i,j}');
    expect_error(@() get_chain(P),id,'P.C{i,j}');         % P(1).C
    expect_error(@() arrayfun(@(x) x,P),id,'');         % arrayfun indexes P(1)
    nchk = nchk+7;
end

% ================================ dot / brace reads unchanged (vs builtin)
for k = 1:numel(CONT)
    P = CONT{k};
    bs = @(varargin) builtin('subsref',P,substruct(varargin{:}));
    assert(isequal(P.C,bs('.','C')),'%s: P.C differs from the builtin',NAME{k});
    [M,N] = size(P.C);
    for i = 1:M
        for j = 1:N
            assert(isequal(P.C{i,j},bs('.','C','{}',{i,j})),'%s: P.C{%d,%d}',NAME{k},i,j);
        end
    end
    assert(isequal(P.C{end,end},bs('.','C','{}',{M,N})),'%s: P.C{end,end}',NAME{k});
    assert(isequal(P.vars,bs('.','vars')) && isequal(P.dom,bs('.','dom')) && ...
           isequal(P.space_out,bs('.','space_out')) && isequal(P.space_in,bs('.','space_in')) && ...
           isequal(P.dim_out(end),P.dim_out(numel(P.dim_out))) && isequal(P.dim_in,bs('.','dim_in')),...
        '%s: metadata reads differ',NAME{k});
    % comma lists: every element, in order, not the first one only
    L = {P.C{:}};                                                           %#ok<CCAT1>
    assert(numel(L)==numel(P.C) && isequal(L,reshape(bs('.','C'),1,[])),...
        '%s: {P.C{:}} is not the full comma list',NAME{k});
    L = {P.C{:,end}};                                                       %#ok<CCAT1>
    assert(numel(L)==M && isequal(L,reshape(P.C(:,N),1,[])),'%s: {P.C{:,end}}',NAME{k});
    if numel(P.C)>=2
        [a,b] = P.C{1:2};
        assert(isequal({a,b},reshape(P.C(1:2),1,[])),'%s: [a,b] = P.C{1:2}',NAME{k});
    end
    % chained into a block
    blk = P.C{1,1};
    assert(isequal(P.C{1,1}.dims,blk.dims) && isequal(P.C{1,1}.vars,blk.vars) && ...
           isequal(P.C{1,1}.ZL,blk.ZL),'%s: chained P.C{1,1}.x',NAME{k});
    if isa(blk,'sdopvar')
        assert(isequal(P.C{1,1}.params.A{1},blk.params.A{1}) && ...
               isequal(P.C{1,1}.params.B{end},blk.params.B{numel(blk.params.B)}),'%s: chained params',NAME{k});
        assert(isequal(P.Zd,bs('.','Zd')) && isequal(P.Zd{end},P.Zd{numel(P.Zd)}),'%s: Zd',NAME{k});
    else
        assert(isequal(P.C{1,1}.params{end},blk.params{numel(blk.params)}),'%s: chained params{end}',NAME{k});
    end
    % method by dot, builtin counts, nesting, and the errors the builtin raises
    assert(isequal(P.metadata,metadata(P)),'%s: P.metadata',NAME{k});
    assert(numel(P)==1 && ~isempty(P) && length(P)==1 && isequal(size(P),size(P.C)) && ...
           isequal(builtin('size',P),[1 1]),'%s: numel/isempty/length/size',NAME{k});
    S = struct('x',P);   c = {P};
    assert(isequal(S.x.C{1,1},blk) && isequal(c{1}.C{1,1},blk),'%s: nested in struct/cell',NAME{k});
    expect_error(@() P{1},'MATLAB:cellRefFromNonCell','');
    expect_error(@() P.nosuchfield,'MATLAB:noSuchMethodOrField','nosuchfield');
    % display, and save/load, which set properties without subsasgn
    txt = evalc('disp(P)');
    assert(contains(txt,'dim_out'),'%s: disp(P) does not list the properties',NAME{k});
    f = [tempname '.mat'];  save(f,'P');  L = load(f);  delete(f);
    assert(isequal(L.P,P),'%s: save/load changed the object',NAME{k});
    nchk = nchk+16;
end
% protected property: still refused, with the builtin's identifier and text
e0 = get_err(@() builtin('subsref',Pc,substruct('.','leftCommonBasis')));
e1 = get_err(@() Pc.leftCommonBasis);
assert(strcmp(e0.identifier,'MATLAB:class:GetProhibited') && strcmp(e1.identifier,e0.identifier) && ...
       strcmp(e1.message,e0.message),'protected read: %s / %s',e1.identifier,e1.message);
e1 = get_err(@() set_prop(Pc,'leftCommonBasis',1));
assert(strcmp(e1.identifier,'MATLAB:class:SetProhibited') && contains(e1.message,'leftCommonBasis'),...
    'protected write on copvar: %s',e1.identifier);
e1 = get_err(@() set_prop(Pc.C{2,2},'vars_S1',{}));
assert(strcmp(e1.identifier,'MATLAB:class:SetProhibited') && contains(e1.message,'vars_S1'),...
    'protected write on sopvar: %s',e1.identifier);
nchk = nchk+3;

% ======================================================== dot assignment
X2 = 2*P2.C{2,2};
Q = P2;     Q.C{2,2} = X2;
assert(isequal(Q.C{2,2},X2) && isequal(Q.C([1 3],:),P2.C([1 3],:)) && isequal(Q.C(2,[1 3]),P2.C(2,[1 3])) ...
       && isequal(metadata(Q),metadata(P2)),'Q.C{2,2} = X changed more than block (2,2)');
Q = P2;     Q.dim_out(1) = 5;
assert(isequal(Q.dim_out,[5;P2.dim_out(2:end)]) && isequal(Q.C,P2.C),'Q.dim_out(1) = 5');
Q = P2;     Q.C{2,2}.params{1} = 3*P2.C{2,2}.params{1};
assert(isequal(Q.C{2,2}.params{1},3*P2.C{2,2}.params{1}) && ...
       isequal(Q.C{2,2}.params(2:end),P2.C{2,2}.params(2:end)),'Q.C{2,2}.params{1} = ...');
Q = P2;     Q.C{end+1,1} = [];
assert(isequal(builtin('size',Q.C),[4 3]) && isempty(Q.C{4,1}),'Q.C{end+1,1} = []');
Q = D2;     Q.C{3,3}.params.A{1}(1) = 7;     Q.Zd{1} = 'renamed';
assert(Q.C{3,3}.params.A{1}(1)==7 && strcmp(Q.Zd{1},'renamed') && isequal(Q.Zd(2:end),D2.Zd(2:end)),...
    'cdopvar dot assignment');
T = D2.C{3,3};  T.params.B{1}(1,1) = 4;  T.Zd{1} = 'q';
assert(T.params.B{1}(1,1)==4 && strcmp(T.Zd{1},'q'),'sdopvar dot assignment');
c = {P2};   c{1}.C{2,2} = X2;
S = struct('x',P2);  S.x.dim_in(2) = 7;
assert(isequal(c{1}.C{2,2},X2) && S.x.dim_in(2)==7,'nested dot assignment');
% comma-list assignment, which needs the builtin count of values
Q = Pc;     [Q.C{1,1:2}] = deal(Pc.C{1,2},Pc.C{1,1});
assert(isequal(Q.C{1,1},Pc.C{1,2}) && isequal(Q.C{1,2},Pc.C{1,1}),'[Q.C{1,1:2}] = deal(...)');
Q = D2;     [Q.Zd{1:2}] = deal('za','zb');
assert(isequal(Q.Zd(1:2),{'za';'zb'}) || isequal(Q.Zd(1:2),{'za','zb'}),'[Qd.Zd{1:2}] = deal(...)');
T = BLK{1}; [T.params{1:2}] = deal(BLK{1}.params{2},BLK{1}.params{1});
assert(isequal(T.params{1},BLK{1}.params{2}) && isequal(T.params{2},BLK{1}.params{1}),'[T.params{1:2}] = deal(...)');
T = BLK{2}; [T.Zd{1:2}] = deal('za','zb');
assert(strcmp(T.Zd{1},'za') && strcmp(T.Zd{2},'zb'),'[Tsd.Zd{1:2}] = deal(...)');
nchk = nchk+11;

% ================================================ (2) paren assignment
OBJ = [CONT, BLK];
for k = 1:numel(OBJ)
    P = OBJ{k};     id = [class(P) ':parenAssign'];
    expect_error(@() asg(P,{1,2},P),id,'');
    expect_error(@() asg(P,{2,1},P),id,'');
    expect_error(@() asg(P,{1},P),id,'');
    expect_error(@() asg(P,{':'},P),id,'');
    expect_error(@() asg(P,{1},[]),id,'');              % deletion Q(1) = []
    expect_error(@() asg_new(P),id,'');                 % y(1) = P, y undefined
    % the failed assignment leaves the variable the scalar object it was
    Q = P;
    try,    Q(1,2) = P;     catch,  end                                     %#ok<NASGU>
    assert(builtin('numel',Q)==1 && isequal(Q,P),'%s: a failed Q(1,2) = P changed Q',class(P));
    nchk = nchk+7;
end
% reached through a container: paren assignment into a block
expect_error(@() asg_nested(Pc,Pc.C{2,2}),'sopvar:parenAssign','');
nchk = nchk+1;

% =================================== sopvar/sdopvar paren slicing unchanged
for kb = 1:numel(BLK)
    B = BLK{kb};
    K = kern(B,dmap);
    m = B.dims(1);  n = B.dims(2);
    for sel = {{1,':'},{':',1},{m:-1:1,':'},{':',[n 1]},{1,1},{':',':'}}
        I = sel{1}{1};  J = sel{1}{2};
        Bs = B(I,J);
        Ks = kern(Bs,dmap);
        if strcmp(I,':'),   I = 1:m;    end
        if strcmp(J,':'),   J = 1:n;    end
        assert(isa(Bs,class(B)) && isequal(Bs.dims,[numel(I),numel(J)]),'%s slice dims',class(B));
        for g = 1:numel(K)
            assert(maxcoef(Ks{g} - K{g}(I,J)) <= 1e-12*max(1,maxcoef(K{g})),...
                '%s(%s,%s): kernels are not the rows/columns of the kernels',class(B),mat2str(I),mat2str(J));
        end
        nchk = nchk+2;
    end
end
% chained, through a container: the same slices as on the block itself
assert(isequal(Pc.C{2,2}(1,:),BLK{1}(1,:)) && isequal(D2.C{3,3}(1,1),BLK{2}(1,1)) && ...
       isequal(D2.C{3,3}(1,1).dims,[1 1]),'slice through a container');
nchk = nchk+1;

% =========================================================== (3) transpose
for k = 1:numel(CONT)
    P = CONT{k};
    Pt = P.';
    assert(strcmp(class(Pt),class(P)) && isequal(size(Pt),fliplr(size(P))),'%s: size of P.''',NAME{k});
    assert(Pt==P' && isequal(Pt,ctranspose(P)),'%s: P.'' is not P''',NAME{k});
    assert(~(Pt==P),'%s: the test operator is self-adjoint, so P.'' = P would pass',NAME{k});
    nchk = nchk+3;
end
% <P.'x, y> = <x, P y> by quadrature, on the domain of each component
for k = 1:3
    P = CONT{k};
    gx = test_funs(P,sum(P.dim_out));    gy = test_funs(P,sum(P.dim_in));
    Gin  = gauss_grids(P,'in',NQ);       Gout = gauss_grids(P,'out',NQ);
    lhs = ip(capply(P.',gx,Gin,dmap,NQ),gy,Gin,comp_list(P,'in'));
    rhs = ip(capply(P,gy,Gout,dmap,NQ),gx,Gout,comp_list(P,'out'));
    assert(abs(lhs-rhs) <= 1e-9*max(1,abs(rhs)),'%s: <P.''x,y> = %.15g but <x,Py> = %.15g',NAME{k},lhs,rhs);
    % calibration: with P in place of P.' (the old answer) the identity fails
    bad = ip(capply(P,gx,Gin,dmap,NQ),gy,Gin,comp_list(P,'in'));
    assert(abs(bad-rhs) > 1e-6*max(1,abs(rhs)),'%s: the adjoint check cannot tell P.'' from P',NAME{k});
    nchk = nchk+2;
end

% ================================================== (4) eq across partitions
% (a) the reported case: R^1 x R^1 x L2^2 (3x3 grid) against R^2 x L2^2 (2x2)
A = blkdiag(Rc,Pc);     B = opvar2copvar(blkdiag(Rop,Pop));
assert(isequal(size(A),[3 3]) && isequal(size(B),[2 2]),'reported case: grids');
check_eq(A,B,true,dmap,NQ,'blkdiag(Rc,Pc) vs opvar2copvar(blkdiag(Rop,Pop))');
% a perturbed component: entry (2,2) of B's R^2 block is A's block (2,2)
Bp = B;     Bp.C{1,1}.params{1}(2,2) = Bp.C{1,1}.params{1}(2,2) + 1e-3;
check_eq(A,Bp,false,dmap,NQ,'reported case, R component perturbed');
Bp = B;     Bp.C{2,2}.params{2}(end,end) = Bp.C{2,2}.params{2}(end,end) + 1e-3;
check_eq(A,Bp,false,dmap,NQ,'reported case, L2 kernel perturbed');
% the two R^1 rows and columns of A exchanged: same sequence of spaces,
% another operator
Ap = copvar(A.C([2 1 3],[2 1 3]));
check_eq(Ap,B,false,dmap,NQ,'R^1 rows exchanged');
% a permuted component ORDER, [L2^2; R^2] against [R^2; L2^2]: another
% ordered product space, so unequal
Bperm = copvar(B.C([2 1],[2 1]));
assert(~(A==Bperm) && ~(Bperm==A),'permuted component order compared equal');
Bperm = copvar(B.C([2 1],:));
assert(~(A==Bperm) && ~(Bperm==A),'permuted output order compared equal');
nchk = nchk+6;

% (b) SPLIT: a container on the fine partition. MERGED: adjacent blocks of
% one space concatenated into one block, by rows and by columns, with the
% block methods horzcat/vertcat. One zero block, in an unmerged position.
SH = { 'copvar', {{}, {'s'}, {'s'}, {'s','t'}},   {1,[2 3],4}
       'cdopvar',{{}, {'s'}, {'s'}, {'s','t'}},   {1,[2 3],4}
       'copvar', {{}, {'s','t','u'}, {'s','t','u'}}, {1,[2 3]} };
for k = 1:size(SH,1)
    [cls,sp,grp] = deal(SH{k,:});
    lbl = sprintf('%s over %d spaces',cls,numel(sp));
    n = numel(sp);  occ = true(n);  occ(1,1) = false;     % R -> R, a singleton group
    spc = struct('out',{sp},'in',{sp});    dmc = struct('out',ones(n,1),'in',ones(n,1));
    if strcmp(cls,'copvar'),    S = rand_copvar(spc,dmc,[0,1],1,0.7,occ);
    else,                       S = rand_cdopvar(spc,dmc,[0,1],1,3,0.7,occ);
    end
    dmap = draw_dvals({S},dmap);
    Mg = merge(S,grp);
    assert(~isequal(size(S),size(Mg)) && isempty(Mg.C{1,1}),'%s: merge did not change the grid',lbl);
    check_eq(S,Mg,true,dmap,NQ,[lbl ': split vs merged']);
    % a coefficient of block (3,2), i.e. inside a merged block and not at
    % its first component, perturbed well above tol
    Mp = merge(perturb(S,1e-3),grp);
    check_eq(S,Mp,false,dmap,NQ,[lbl ': merged block perturbed']);
    % TOLERANCE: the same coefficient perturbed by 1e-10. The refined path
    % (S vs merged) and the same-partition path (S vs split) give the same
    % verdicts, false at the default 1e-14 and true at tol = 1e-9.
    St = perturb(S,1e-10);      Mt = merge(St,grp);
    v = [eq(S,Mt), eq(S,Mt,1e-9), eq(S,St), eq(S,St,1e-9), eq(Mt,St)];
    assert(isequal(v,[false true false true true]),'%s: tolerance verdicts %s, want [0 1 0 1 1]',lbl,mat2str(v));
    % zero operators in the two partitions, and the '== 0' path
    assert((0*S)==(0*Mg) && ~(S==0) && ~(0==Mg) && (Mg-Mg)==0 && (0*Mg)==0,'%s: zero tests',lbl);
    nchk = nchk+4;
    if strcmp(cls,'cdopvar')
        % the same decision operator on a REORDERED list, matched by name
        Mr = reverse_dvars(Mg);
        assert(~isequal(Mr.Zd(:),S.Zd(:)),'reversed list equals the original');
        check_eq(S,Mr,true,dmap,NQ,[lbl ': merged on a reversed decision list']);
        nchk = nchk+2;
    end
end

% (c) different ordered sequences of spaces with the same totals, and
% different totals: false, not an error
A1 = rand_copvar(struct('out',{{ {'s'} }},'in',{{ {'s'} }}),struct('out',2,'in',2),[0,1],1,0.7);
A2 = rand_copvar(struct('out',{{ {'t'} }},'in',{{ {'s'} }}),struct('out',2,'in',2),[0,1],1,0.7);
A3 = rand_copvar(struct('out',{{ {'s'},{'s'} }},'in',{{ {'s'} }}),struct('out',[1;1],'in',2),[0,1],1,0.7);
A4 = rand_copvar(struct('out',{{ {'s'},{'t'} }},'in',{{ {'s'} }}),struct('out',[1;1],'in',2),[0,1],1,0.7);
assert(~(A1==A2) && ~(A2==A1),'L2^2[s] and L2^2[t] outputs compared equal');
assert(~(A3==A2) && ~(A4==A1) && ~(A4==A3),'[s,s], [s,t] and [t,t] output sequences compared equal');
assert(~(A1==P3) && ~(A3==P2),'different totals compared equal');
nchk = nchk+3;

% ======================================= (5) comma lists of block properties
% Each list is compared with the list formed natively from the underlying
% cell or struct, read with builtin('subsref',...): {c{:}}, [c{:}], f(c{:})
% and [a,b] = c{1:2}. Every list here has >= 2 elements unless noted, so a
% first-element truncation fails.
f  = @(varargin) varargin;
BL = {Pc.C{2,2}, P2.C{3,3}, P2.C{2,3}, D2.C{3,3}, D2.C{2,2}};
for kb = 1:numel(BL)
    X = BL{kb};     lbl = sprintf('block %d (%s)',kb,class(X));
    bx = @(varargin) builtin('subsref',X,substruct(varargin{:}));
    ZL = bx('.','ZL');   ZR = bx('.','ZR');   V = bx('.','vars');   prm = bx('.','params');
    same_list({X.ZL{:}},{ZL{:}},[lbl ': {X.ZL{:}}']);                       %#ok<CCAT1>
    same_list({X.ZR{:}},{ZR{:}},[lbl ': {X.ZR{:}}']);                       %#ok<CCAT1>
    same_list(f(X.ZL{:}),f(ZL{:}),[lbl ': f(X.ZL{:})']);
    same_list({X.vars.in{:}},{V.in{:}},[lbl ': {X.vars.in{:}}']);           %#ok<CCAT1>
    same_list({X.vars.out{:}},{V.out{:}},[lbl ': {X.vars.out{:}}']);        %#ok<CCAT1>
    if isa(X,'sdopvar')
        A = prm.A;  Bc = prm.B;  Zd = bx('.','Zd');
        assert(numel(A)>=3 && numel(Zd)>=3,'%s: lists too short to detect truncation',lbl);
        same_list({X.params.A{:}},{A{:}},[lbl ': {X.params.A{:}}']);        %#ok<CCAT1>
        same_list([X.params.A{:}],[A{:}],[lbl ': [X.params.A{:}]']);
        same_list(f(X.params.B{:}),f(Bc{:}),[lbl ': f(X.params.B{:})']);
        same_list({X.params.B{2:end}},{Bc{2:end}},[lbl ': {X.params.B{2:end}}']); %#ok<CCAT1>
        same_list({X.params.A{[]}},{A{[]}},[lbl ': {X.params.A{[]}} (empty)']); %#ok<CCAT1>
        same_list({X.Zd{:}},{Zd{:}},[lbl ': {X.Zd{:}}']);                   %#ok<CCAT1>
        same_list(f(X.Zd{:}),f(Zd{:}),[lbl ': f(X.Zd{:})']);
        [a,b] = X.Zd{1:2};      same_list({a,b},{Zd{1:2}},[lbl ': [a,b] = X.Zd{1:2}']); %#ok<CCAT1>
    else
        assert(numel(prm)>=3,'%s: params too short to detect truncation',lbl);
        same_list({X.params{:}},{prm{:}},[lbl ': {X.params{:}}']);          %#ok<CCAT1>
        same_list([X.params{:}],[prm{:}],[lbl ': [X.params{:}]']);
        same_list(f(X.params{:}),f(prm{:}),[lbl ': f(X.params{:})']);
        same_list({X.params{1:2}},{prm{1:2}},[lbl ': {X.params{1:2}}']);    %#ok<CCAT1>
        same_list({X.params{[true false true]}},{prm{[true false true]}},[lbl ': logical index']); %#ok<CCAT1>
        same_list({X.params{[]}},{prm{[]}},[lbl ': {X.params{[]}} (empty)']); %#ok<CCAT1>
        [a,b] = X.params{1:2};  same_list({a,b},{prm{1:2}},[lbl ': [a,b] = X.params{1:2}']); %#ok<CCAT1>
    end
    % the count MATLAB passes as nargout: equal to the builtin count, for
    % the chains 'numArgumentsFromSubscript' counts itself and for the rest
    ch = { {'.','dims'}, {'.','params'}, {'.','ZL','{}',{1}}, {'.','ZL','{}',{':'}}, ...
           {'.','vars'}, {'.','vars','.','in'}, {'.','vars','.','in','{}',{1}}, ...
           {'.','vars','.','in','{}',{':'}}, {'.','dom','.','in','()',{1,':'}}, {'.','dims','()',{1}} };
    if isa(X,'sdopvar')
        ch = [ch, { {'.','params','.','A','{}',{1}}, {'.','params','.','B','{}',{':'}}, ...
                    {'.','params','.','A','{}',{2:3}}, {'.','Zd','{}',{1}}, {'.','Zd','{}',{':'}}, ...
                    {'.','Zd','()',{1:2}}, {'.','params','.','A','{}',{1},'()',{1}} }];     %#ok<AGROW>
    else
        ch = [ch, { {'.','params','{}',{1}}, {'.','params','{}',{':'}}, {'.','params','{}',{1:2}}, ...
                    {'.','params','{}',{[]}}, {'.','params','()',{1}}, {'.','params','{}',{1},'()',{1,1}} }]; %#ok<AGROW>
    end
    for ctx = [matlab.indexing.IndexingContext.Expression, matlab.indexing.IndexingContext.Statement, ...
               matlab.indexing.IndexingContext.Assignment]
        for kc = 1:numel(ch)
            sc = substruct(ch{kc}{:});
            n1 = numArgumentsFromSubscript(X,sc,ctx);
            n0 = builtin('numArgumentsFromSubscript',X,sc,ctx);
            assert(n1==n0,'%s: count of chain %d is %d, builtin %d (%s)',lbl,kc,n1,n0,char(ctx));
        end
    end
    nchk = nchk + 5 + 8 + 3*numel(ch);
end
% through containers, cells and structs
CT = {Pc,2,2; P2,3,3; P2,2,3; D2,3,3; D2,2,2};
for kc = 1:size(CT,1)
    [P,i,j] = deal(CT{kc,:});   lbl = sprintf('%s.C{%d,%d}',class(P),i,j);
    blk = builtin('subsref',P,substruct('.','C','{}',{i,j}));
    bb  = @(varargin) builtin('subsref',blk,substruct(varargin{:}));
    ZL = bb('.','ZL');  V = bb('.','vars');  prm = bb('.','params');
    same_list(f(P.C{i,j}.ZL{:}),f(ZL{:}),[lbl ': f(P.C{i,j}.ZL{:})']);
    same_list({P.C{i,j}.vars.in{:}},{V.in{:}},[lbl ': {P.C{i,j}.vars.in{:}}']); %#ok<CCAT1>
    if isa(blk,'sdopvar')
        Zd = bb('.','Zd');
        same_list({P.C{i,j}.Zd{:}},{Zd{:}},[lbl ': {P.C{i,j}.Zd{:}}']);     %#ok<CCAT1>
        same_list(f(P.C{i,j}.params.B{:}),f(prm.B{:}),[lbl ': f(P.C{i,j}.params.B{:})']);
        same_list({P.C{i,j}.params.A{:}},{prm.A{:}},[lbl ': {P.C{i,j}.params.A{:}}']); %#ok<CCAT1>
        same_list([P.C{i,j}.params.A{:}],[prm.A{:}],[lbl ': [P.C{i,j}.params.A{:}]']);
    else
        same_list({P.C{i,j}.params{:}},{prm{:}},[lbl ': {P.C{i,j}.params{:}}']); %#ok<CCAT1>
        same_list([P.C{i,j}.params{:}],[prm{:}],[lbl ': [P.C{i,j}.params{:}]']);
        same_list(f(P.C{i,j}.params{:}),f(prm{:}),[lbl ': f(P.C{i,j}.params{:})']);
        same_list({P.C{i,j}.params{1:2}},{prm{1:2}},[lbl ': {P.C{i,j}.params{1:2}}']); %#ok<CCAT1>
    end
    c = {P};    Sx = struct('P',P,'b',blk);     cb = {blk};
    same_list({c{1}.C{i,j}.ZL{:}},{ZL{:}},[lbl ': in a cell']);             %#ok<CCAT1>
    same_list({Sx.P.C{i,j}.ZL{:}},{ZL{:}},[lbl ': in a struct']);           %#ok<CCAT1>
    same_list({Sx.b.ZL{:}},{ZL{:}},[lbl ': block in a struct']);            %#ok<CCAT1>
    same_list({cb{1}.ZL{:}},{ZL{:}},[lbl ': block in a cell']);             %#ok<CCAT1>
    nchk = nchk+10;
end
% the containers' own lists
Zd = builtin('subsref',D2,substruct('.','Zd'));     vv = builtin('subsref',P2,substruct('.','vars'));
same_list({D2.Zd{:}},{Zd{:}},'{D2.Zd{:}}');                                 %#ok<CCAT1>
same_list({P2.vars{:}},{vv{:}},'{P2.vars{:}}');                             %#ok<CCAT1>
nchk = nchk+2;

% ============================================================== (6) cat
% cat(1,...) and cat(2,...) are [A;B] and [A,B], the overloaded operator
% concatenations; any other DIM is an error, never an object array.
CA = {Pc, P2, D2, BLK{1}, BLK{2}, P2.C{2,2}, D2.C{2,2}};
for k = 1:numel(CA)
    X = CA{k};  cls = class(X);
    V1 = cat(1,X,X);    V2 = cat(2,X,X);
    assert(builtin('numel',V1)==1 && builtin('numel',V2)==1,'%s: cat built an object array',cls);
    assert(isequal(V1,[X;X]) && isequal(V2,[X,X]),'%s: cat(1|2,X,X) differs from [X;X] / [X,X]',cls);
    % 3 operands: the same outcome as [X;X;X], value or error (the kernels
    % of [A;B;C] are checked in (8)).
    assert(isequal(outcome(@() cat(1,X,X,X)),outcome(@() [X;X;X])) && ...
           isequal(outcome(@() cat(2,X,X,X)),outcome(@() [X,X,X])) && ...
           isequal(cat(1,X),vertcat(X)),'%s: cat with 3 / 1 operands',cls);
    expect_error(@() cat(3,X,X),[cls ':catDim'],'DIM = 1');
    expect_error(@() cat(0,X,X),[cls ':catDim'],'DIM = 1');
    expect_error(@() cat(1.5,X,X),[cls ':catDim'],'DIM = 1');
    nchk = nchk+6;
end
% mixed operands dispatch as [A;B] / [A,B] do
assert(isequal(cat(1,P2,D2),[P2;D2]) && isequal(cat(2,P2,D2),[P2,D2]),'cat of copvar and cdopvar');
assert(isequal(cat(2,P2.C{2,2},D2.C{2,2}),[P2.C{2,2},D2.C{2,2}]),'cat of sopvar and sdopvar');
e1 = get_err(@() cat(1,Pc,Pc.C{2,2}));   e2 = get_err(@() vertcat(Pc,Pc.C{2,2}));
assert(~isempty(e1.identifier) && strcmp(e1.identifier,e2.identifier),'cat(1,Pc,block): %s vs %s',e1.identifier,e2.identifier);
nchk = nchk+3;

% ====================================== (7) eq with zero-count components
% A block row or column with zero components holds no component, so it
% must not decide A == B: the three partitions of one operator below
% compare equal pairwise (transitivity), checked by the quadrature.
sz = struct('out',{{ {}, {'s1'} }},'in',{{ {}, {'s1'} }});
rng(3);     Z0 = rand_copvar(sz,struct('out',[0;2],'in',[1;2]),[0,1],1,0.8);    % [R^0; L2^2[s1]]
Zr = copvar(Z0.C(2,:));                                                     % the L2 row alone
sz2 = struct('out',{{ {'s1'}, {'s1'} }},'in',{{ {}, {'s1'} }});
rng(3);     Zb = rand_copvar(sz2,struct('out',[0;2],'in',[1;2]),[0,1],1,0.8);
Zb2 = copvar({Zb.C{1,1},Zb.C{1,2}; Z0.C{2,1},Z0.C{2,2}});                   % [L2^0[s1]; L2^2[s1]]
assert(isequal(size(Z0),size(Zb2)) && isequal(Z0.dim_out,Zb2.dim_out) && ...
       ~isequal(Z0.space_out(1,:),Zb2.space_out(1,:)),'zero-row case: same grid and dims, other space');
check_eq(Z0,Zr,true,dmap,NQ,'[R^0;L2^2] vs [L2^2]');
check_eq(Zb2,Zr,true,dmap,NQ,'[L2^0;L2^2] vs [L2^2]');
check_eq(Z0,Zb2,true,dmap,NQ,'[R^0;L2^2] vs [L2^0;L2^2], same grid');
Zp = Zb2;   Zp.C{2,2}.params{2}(1,1) = Zp.C{2,2}.params{2}(1,1) + 1e-3;
check_eq(Z0,Zp,false,dmap,NQ,'zero-row case, L2 kernel perturbed');
% a zero-count COLUMN: [R^0, L2^2] and [L2^0[s1], L2^2] as inputs
sz3 = struct('out',{{ {}, {'s1'} }},'in',{{ {'s1'}, {'s1'} }});
rng(4);     Zc = rand_copvar(sz,struct('out',[1;2],'in',[0;2]),[0,1],1,0.8);
Zcr = copvar(Zc.C(:,2));
rng(5);     Zq = rand_copvar(sz3,struct('out',[1;2],'in',[0;2]),[0,1],1,0.8);
Zc3 = copvar({Zq.C{1,1},Zc.C{1,2}; Zq.C{2,1},Zc.C{2,2}});
check_eq(Zc,Zcr,true,dmap,NQ,'[R^0,L2^2] vs [L2^2] inputs');
check_eq(Zc3,Zcr,true,dmap,NQ,'[L2^0,L2^2] vs [L2^2] inputs');
check_eq(Zc,Zc3,true,dmap,NQ,'[R^0,L2^2] vs [L2^0,L2^2] inputs, same grid');
% no output components at all: L2^2[s1] -> R^0 is the one zero map, so
% equal iff the input sequences are; one zero row or two
Zt  = rand_copvar(struct('out',{{ {} }},'in',{{ {'s1'} }}),struct('out',0,'in',2),[0,1],1,0.8);
Zt2 = rand_copvar(struct('out',{{ {}, {'s1'} }},'in',{{ {'s1'} }}),struct('out',[0;0],'in',2),[0,1],1,0.8);
Zt3 = rand_copvar(struct('out',{{ {} }},'in',{{ {'t'} }}),struct('out',0,'in',2),[0,1],1,0.8);
assert(Zt==Zt2 && Zt2==Zt && ~(Zt==Zt3) && ~(Zt3==Zt2),'operators into R^0');
nchk = nchk+9;

% ==================================== (8) [A;B;C] of sopvar, 3 operands
% @sopvar/vertcat recursed with horzcat, so [A;B;C] was [[A;B],C] whenever
% C had as many rows as [A;B]. Rows 1+2+3 below make that shape legal. The
% kernels of [A;B;C] must be those of A, B, C stacked, cell by cell.
rng(8);
Q3 = rand_copvar(struct('out',{{ {'s'},{'s'},{'s'} }},'in',{{ {'s'} }}),struct('out',[1;2;3],'in',2),[0,1],1,0.7);
[A,B,C] = deal(Q3.C{1,1},Q3.C{2,1},Q3.C{3,1});
V = [A;B;C];
assert(isequal(V.dims,[6 2]),'[A;B;C] has dims %s, expected [6 2]',mat2str(V.dims));
KV = kern(V,dmap);  KA = kern(A,dmap);  KB = kern(B,dmap);  KC = kern(C,dmap);
for ii = 1:numel(KV)
    assert(maxcoef(KV{ii}-[KA{ii};KB{ii};KC{ii}])<1e-12,'[A;B;C]: kernel %d is not [KA;KB;KC]',ii);
end
H = [A.',B.',C.'];     % horzcat of 3, already correct: the same check
KH = kern(H,dmap);
for ii = 1:numel(KH)
    Kt = kern(A.',dmap);  Kt2 = kern(B.',dmap);  Kt3 = kern(C.',dmap);
    assert(maxcoef(KH{ii}-[Kt{ii},Kt2{ii},Kt3{ii}])<1e-12,'[A'',B'',C'']: kernel %d',ii);
end
nchk = nchk+2;

% ===================== (9) indexing after a block slice, P(I,J).<tail>
% The tail after '()' now continues on sopvar (it was dropped: P(1,1).dims
% returned the slice), and a comma-list brace after a slice errors on both
% block classes (a '()' read has one output, so the list kept its first
% element). A brace with scalar subscripts is one value and stays allowed.
for k = 1:2
    X = BLK{k};     cls = class(X);     Y = X(1,1);
    assert(isequal(X(1,1).dims,Y.dims) && isequal(X(1,1).ZL,Y.ZL) && ...
           isequal(X(1,1).ZL{1},Y.ZL{1}),'%s: X(1,1).<tail> is not Y.<tail> of Y = X(1,1)',cls);
    expect_error(@() {X(1,1).ZL{:}},[cls ':parenThenList'],'Slice first');
    expect_error(@() {X(1,1).ZL{[1 1]}},[cls ':parenThenList'],'Slice first');  % 2 values
    nchk = nchk+3;
end
Y = BLK{2}(1,1);
assert(isequal(BLK{2}(1,1).Zd{2},Y.Zd{2}),'sdopvar: X(1,1).Zd{2}');
expect_error(@() {BLK{2}(1,1).Zd{:}},'sdopvar:parenThenList','Slice first');
% through a container, where it raised MATLAB:unassignedOutputs
e = get_err(@() {D2.C{3,3}(1,1).Zd{:}});
assert(strcmp(e.identifier,'sdopvar:parenThenList'),'D2.C{3,3}(1,1).Zd{:}: %s',e.identifier);
nchk = nchk+3;

fprintf('test_copvar_silent_fixes passed (%d checks).\n',nchk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = get_chain(P),      x = P(1).C;             end
function Q = asg(P,idx,val),    Q = P;  Q(idx{:}) = val;    end
function y = asg_new(P),        y(1) = P;               end
function Q = asg_nested(Q,B),   Q.C{2,2}(1,2) = B;      end
function X = set_prop(X,nm,v),  X.(nm) = v;             end

function e = get_err(f)
e = struct('identifier','','message','');
try,        f();
catch ME,   e = ME;
end
end

function o = outcome(h)
% {value} of h(), or {identifier, message} of the error it raises.
try,        o = {h()};
catch ME,   o = {ME.identifier, ME.message};
end
end

function same_list(x,y,lbl)
% x (through the overloads) must be y (native), class, size and value.
assert(strcmp(class(x),class(y)) && isequal(size(x),size(y)) && isequal(x,y),...
    'test_copvar_silent_fixes: %s: %s%s, native %s%s',lbl,class(x),mat2str(size(x)),class(y),mat2str(size(y)));
end

function expect_error(f,id,msg)
% f must raise error 'id' whose message contains 'msg'.
try
    f();
catch ME
    if ~strcmp(ME.identifier,id)
        error('test_copvar_silent_fixes: expected error ''%s'', got ''%s'': %s',id,ME.identifier,ME.message);
    end
    if ~contains(ME.message,msg)
        error('test_copvar_silent_fixes: error ''%s'' does not mention "%s": %s',id,msg,ME.message);
    end
    return
end
error('test_copvar_silent_fixes: expected error ''%s'', but none was raised.',id);
end


% ------------------------------------------------ partitions of one operator
function Mg = merge(S,grp)
% Blocks of S over the groups grp of block rows (and of block columns, the
% same groups) concatenated into one block each: [X,Y] within a row group,
% then vertically. A group of [] blocks stays []; a partly-[] group is not
% generated here.
C = cell(numel(grp));
for I = 1:numel(grp)
    for J = 1:numel(grp)
        blk = S.C(grp{I},grp{J});
        if all(cellfun(@isempty,blk(:))),  continue,   end
        assert(~any(cellfun(@isempty,blk(:))),'merge: partly zero group');
        rows = cell(size(blk,1),1);
        for r = 1:size(blk,1),     rows{r} = horzcat(blk{r,:});    end
        C{I,J} = vertcat(rows{:});
    end
end
Mg = feval(class(S),C);
end

function S = perturb(S,delta)
% Add delta to the first coefficient of the all-integral parameter of block
% (3,2) of S (no canonical-form constraint applies there).
B = S.C{3,2};
n3 = numel(intersect(B.vars.in,B.vars.out));
g = 1 + (3^n3-1)/2;
if isa(B,'sdopvar')
    assert(numel(B.params.A{g})>1,'perturb: zero-shorthand parameter');
    B.params.A{g}(1) = B.params.A{g}(1) + delta;
else
    assert(numel(B.params{g})>1,'perturb: zero-shorthand parameter');
    B.params{g}(1,1) = B.params{g}(1,1) + delta;
end
S.C{3,2} = B;
end


% ------------------------------------------------------ semantic evaluation
function check_eq(A,B,want,dmap,nq,lbl)
% eq(A,B) and eq(B,A) must both be WANT, and so must the verdict of the
% quadrature evaluation of A x and B x.
v = [A==B, B==A];
assert(isequal(v,[want want]),'%s: A==B, B==A = %s, want %d',lbl,mat2str(v),want);
ca = comp_list(A,'out');    cb = comp_list(B,'out');
ia = comp_list(A,'in');     ib = comp_list(B,'in');
assert(isequal({ca.key},{cb.key}) && isequal({ia.key},{ib.key}),...
    '%s: the component sequences differ, so there is no semantic check',lbl);
g = test_funs(A,numel(ia));
pts = rand_grids(A,'out',3);
ya = capply(A,g,pts,dmap,nq);   yb = capply(B,g,pts,dmap,nq);
sc = max(1,max(cellfun(@(y) max(abs(y)),ya)));
err = max(cellfun(@(u,v) max(abs(u-v)),ya,yb))/sc;
if want
    assert(err <= 1e-9,'%s: eq says equal, quadrature differs by %.3g',lbl,err);
else
    assert(err >= 1e-6,'%s: eq says unequal, quadrature differs by only %.3g',lbl,err);
end
end


function cl = comp_list(P,side)
% One entry per global component of one side of P, in order: its variable
% names (sorted), 'nkey' those names joined, 'key' names and domains.
if strcmp(side,'out'),  msk = P.space_out;  d = P.dim_out;
else,                   msk = P.space_in;   d = P.dim_in;
end
cl = struct('names',{},'nkey',{},'key',{});
for i = 1:size(msk,1)
    [nm,o] = sort(P.vars(msk(i,:)));   dm = P.dom(msk(i,:),:);
    e = struct('names',{nm},'nkey',strjoin(nm,','),'key',[strjoin(nm,',') '|' mat2str(dm(o,:))]);
    cl = [cl, repmat(e,1,d(i))];                                            %#ok<AGROW>
end
end


function g = test_funs(P,n)
% n random polynomial test functions, one per global component c:
% g(c,names,X) = a(c) prod_v (1 + b(c,v) x_v + e(c,v) x_v^2), X with one
% column per name; coefficients drawn once per (component, variable NAME).
a = randn(n,1);
cf = containers.Map();
for v = P.vars,     cf(v{1}) = randn(n,2);  end
g = @(c,names,X) tf_eval(a,cf,c,names,X);
end

function y = tf_eval(a,cf,c,names,X)
y = a(c)*ones(size(X,1),1);
for k = 1:numel(names)
    b = cf(names{k});
    y = y .* (1 + b(c,1)*X(:,k) + b(c,2)*X(:,k).^2);
end
end


function G = rand_grids(P,side,ns)
% ns random points per space of one side of P, by name key; an R space has
% the single empty point.
G = containers.Map();
for e = comp_list(P,side)
    if isKey(G,e.nkey),  continue,   end
    [~,iv] = ismember(e.names,P.vars);
    X = zeros(ns,numel(iv));
    for k = 1:numel(iv),    X(:,k) = P.dom(iv(k),1) + diff(P.dom(iv(k),:))*rand(ns,1);    end
    if isempty(iv),     X = zeros(1,0);     end
    G(e.nkey) = struct('names',{e.names},'X',X,'w',ones(size(X,1),1));
end
end


function G = gauss_grids(P,side,nq)
% Tensor Gauss-Legendre grid and weights per space of one side of P.
G = containers.Map();
[x0,w0] = gl01(nq);
for e = comp_list(P,side)
    if isKey(G,e.nkey),  continue,   end
    [~,iv] = ismember(e.names,P.vars);
    X = zeros(1,0);     w = 1;
    for k = 1:numel(iv)
        a = P.dom(iv(k),1);   L = diff(P.dom(iv(k),:));
        n = size(X,1);
        X = [repmat(X,nq,1), repelem(a+L*x0,n,1)];
        w = repmat(w,nq,1).*repelem(L*w0,n,1);
    end
    G(e.nkey) = struct('names',{e.names},'X',X,'w',w);
end
end


function Y = capply(P,g,G,dmap,nq)
% Y{c} = (P x)_c at the points of G for the output component c's space,
% x = g over P's input components; each block applied by 'heatNd_apply'
% from its kernel definition, a decision block at the values in dmap.
ro = [0;cumsum(P.dim_out(:))];     co = [0;cumsum(P.dim_in(:))];
cl = comp_list(P,'out');
Y = cell(1,ro(end));
for i = 1:size(P.C,1)
    if ro(i+1)==ro(i),  continue,   end     % no components in this row
    pt = G(cl(ro(i)+1).nkey);
    for c = ro(i)+1:ro(i+1),    Y{c} = zeros(size(pt.X,1),1);   end
    for j = 1:size(P.C,2)
        B = P.C{i,j};
        if isempty(B) || co(j+1)==co(j),  continue,   end   % [] or no inputs
        B = fixblock(B,dmap);
        vin = B.vars.in;    cj = co(j)+1:co(j+1);
        f = @(X) cell2mat(arrayfun(@(c) g(c,vin,X),cj,'UniformOutput',false));
        [~,loc] = ismember(B.vars.out,pt.names);
        y = heatNd_apply(B,f,pt.X(:,loc),nq);
        for a = 1:ro(i+1)-ro(i),    Y{ro(i)+a} = Y{ro(i)+a} + y(:,a);   end
    end
end
end


function v = ip(Y,g,G,cl)
% sum_c int Y_c g_c over the component's space, by the weights of G.
v = 0;
for c = 1:numel(Y)
    pt = G(cl(c).nkey);
    v = v + sum(pt.w .* Y{c} .* g(c,pt.names,pt.X));
end
end


function B = fixblock(B,dmap)
% An sdopvar block at the decision values of dmap, as an sopvar with
% C_gam = unvec(A_gam + B_gam'*d); a sopvar block is returned as it is.
if ~isa(B,'sdopvar'),   return,     end
d = zeros(numel(B.Zd),1);
for i = 1:numel(B.Zd),  d(i) = dmap(B.Zd{i});   end
nr = B.dims(1)*prod([cellfun(@numel,B.ZL),1]);
nc = B.dims(2)*prod([cellfun(@numel,B.ZR),1]);
prm = cell(size(B.params.A));
for k = 1:numel(prm)
    v = zeros(nr*nc,1);
    Ak = B.params.A{k};     Bk = B.params.B{k};
    if numel(Ak)==nr*nc,    v = v + full(Ak(:));    end
    if size(Bk,2)==nr*nc,   v = v + full(Bk.'*d);   end
    prm{k} = reshape(v,nr,nc);
end
B = sopvar(prm,B.vars,B.ZL,B.ZR,B.dom,B.dims);
end


function [x,w] = gl01(n)
% Gauss-Legendre nodes and weights on [0,1] (Golub-Welsch).
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     w = 2*V(1,ix)'.^2;
x = (x+1)/2;    w = w/2;
end


% ------------------------------------------------------ kernels of a block
function dmap = draw_dvals(Ps,dmap)
% One random value per decision variable NAME, kept across calls.
if nargin<2,    dmap = containers.Map('KeyType','char','ValueType','double');  end
for k = 1:numel(Ps)
    if ~isa(Ps{k},'cdopvar'),   continue,   end
    z = Ps{k}.Zd;
    for i = 1:numel(z)
        if ~isKey(dmap,z{i}),   dmap(z{i}) = randn;     end
    end
end
end

function K = kern(B,dmap)
% Kernels of one block, a decision block at the drawn values.
if isa(B,'sdopvar')
    z = B.Zd;   dval = zeros(numel(z),1);
    for i = 1:numel(z),     dval(i) = dmap(z{i});   end
    K = pi_blk_kernels(B,dval);
else
    K = pi_blk_kernels(B,[]);
end
end

function m = maxcoef(X)
if isnumeric(X)
    m = full(max(abs(X(:))));
else
    X = polynomial(X);  cf = X.coefficient;
    m = full(max(abs(cf(:))));
end
if isempty(m),  m = 0;  end
end

function P = reverse_dvars(P)
% The same decision operator with its decision list reversed: every block's
% B rows permuted with the names (setdvars), then the container rebuilt.
Zr = flipud(P.Zd(:));
C = P.C;
for ii = 1:numel(C)
    if isa(C{ii},'sdopvar'),    C{ii} = setdvars(C{ii},Zr);     end
end
P = cdopvar(C);
end
