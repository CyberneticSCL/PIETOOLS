%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_COPVAR_BLOCKOPS checks the container block operations horzcat,
% vertcat, blkdiag, uminus, minus and eq, on 'copvar' and 'cdopvar'.
%
% Every check is against SEMANTICS, never a round trip (CLAUDE.md S4). The
% oracle is 'pi_blk_kernels', which evaluates each block's kernels from the
% class's own kernel definition; a decision block is evaluated at random
% values of its decision variables, drawn once per NAME so that containers
% on different lists are compared consistently. The identities used go
% through 'mtimes' and 'plus', which are independent of the routines under
% test and are covered by their own suites:
%
%   uminus    kernels(-P)   = -kernels(P)                      blockwise
%   minus     kernels(A-B)  =  kernels(A) - kernels(B)         blockwise
%   horzcat   [A, B]*[X; Y] =  A*X + B*Y
%   vertcat   [A; B]*X      =  [A*X; B*X], read off block rows
%   blkdiag   blkdiag(A,B)*[X; Y] = [A*X; B*Y]
%   eq        against constructed equal / unequal / zero cases
%
% plus metadata checks by variable NAME, registry merging (including an
% R^n-only operand, whose registry is empty), the decision variable list
% merge, class promotion, and the error paths.
%
% Shapes vary the number of spatial variables: the 1-D PIE layout
% R^2 x L2^3[s], a 2-D layout R x L2[s] x L2[s,t], a 3-D layout
% R x L2[s,t,u], and a cross layout L2[s] x L2[t]. Structurally zero blocks
% are included via the 'occ' argument of the generators.
%
% MMP, 09/25/2026: Initial coding
% MMP, 09/25/2026: Review cases. A decision operand whose sdopvar blocks are
%                  on an EMPTY list, in [D,P0], [P0,D] and D+P0 (such blocks
%                  were left off the container list); an operand whose list
%                  already is the union (left in place) and one holding the
%                  same names permuted (moved); 'verify' on every
%                  concatenation; and the 'setdvars' row-map guard.

rng(20260925);
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
nchk = 0;

SH = { struct('name','1D R^2 x L2^3[s]', ...
              'sp',struct('out',{{ {}, {'s'} }},'in',{{ {}, {'s'} }}), ...
              'dm',struct('out',[2;3],'in',[2;3]),'deg',2)
       struct('name','2D R x L2[s] x L2[s,t]', ...
              'sp',struct('out',{{ {}, {'s'}, {'s','t'} }},'in',{{ {}, {'s'}, {'s','t'} }}), ...
              'dm',struct('out',[1;2;1],'in',[1;2;1]),'deg',1)
       struct('name','3D R x L2[s,t,u]', ...
              'sp',struct('out',{{ {}, {'s','t','u'} }},'in',{{ {}, {'s','t','u'} }}), ...
              'dm',struct('out',[1;1],'in',[1;1]),'deg',1)
       struct('name','cross L2[s] x L2[t]', ...
              'sp',struct('out',{{ {'s'}, {'t'} }},'in',{{ {'s'}, {'t'} }}), ...
              'dm',struct('out',[2;2],'in',[2;2]),'deg',2) };

for ish = 1:numel(SH)
    sh = SH{ish};
    M = numel(sh.sp.out);   N = numel(sh.sp.in);
    % A zero block in every container, but no empty row or column, which
    % the generators cannot represent.
    occ = true(M,N);
    if M>=2 && N>=2,    occ(1,2) = false;   end

    P = rand_copvar(sh.sp,sh.dm,[0,1],sh.deg,0.7,occ);
    Q = rand_copvar(sh.sp,sh.dm,[0,1],sh.deg,0.7);
    D = rand_cdopvar(sh.sp,sh.dm,[0,1],sh.deg,3,0.7,occ);
    E = rand_cdopvar(sh.sp,sh.dm,[0,1],sh.deg,5,0.7);    % a different list
    dmap = draw_dvals({D,E});

    % ---------------------------------------------------------- uminus
    for X = {P,D}
        X = X{1};                                                   %#ok<FXSET>
        R = -X;
        assert(strcmp(class(R),class(X)),'uminus changed the class');
        assert(isequal(metadata(R),metadata(X)),'uminus changed the metadata');
        assert(isequal(cellfun(@isempty,R.C),cellfun(@isempty,X.C)),...
            'uminus changed the zero-block pattern');
        assert(sem_equal(R,X,dmap,-1),'uminus: kernels(-P) ~= -kernels(P)');
        nchk = nchk+4;
    end

    % ---------------------------------------------------------- minus
    pairs = {P,Q; D,E; P,D; D,P};
    for k = 1:size(pairs,1)
        A = pairs{k,1};     B = pairs{k,2};
        R = A - B;
        wantcls = 'copvar';
        if isa(A,'cdopvar') || isa(B,'cdopvar'),   wantcls = 'cdopvar';   end
        assert(strcmp(class(R),wantcls),'minus returned class %s',class(R));
        assert(sem_lincomb(R,{A,B},[1,-1],dmap),'minus: kernels(A-B) ~= kernels(A)-kernels(B)');
        nchk = nchk+2;
    end

    % ---------------------------------------------------------- horzcat
    % [A, B]: A = P (rows sp.out, cols sp.in), B maps a different input
    % space (L2[s] or R^1) into the same rows.
    spB = struct('out',{sh.sp.out},'in',{{ sh.sp.in{end} }});
    dmB = struct('out',sh.dm.out,'in',1);
    Bf = rand_copvar(spB,dmB,[0,1],sh.deg,0.7);
    Bd = rand_cdopvar(spB,dmB,[0,1],sh.deg,4,0.7);
    dmap = draw_dvals({D,E,Bd},dmap);
    % A right factor over P's input spaces, for the disjoint-list checks.
    Z0 = on_registry(rand_copvar(struct('out',{sh.sp.in},'in',{{ {} }}), ...
                     struct('out',sh.dm.in,'in',1),[0,1],sh.deg,0.7),P);
    % Right factors: X from P's input spaces, Y from B's, into one space.
    spX = struct('out',{sh.sp.in},'in',{{ {} }});
    spY = struct('out',{{ sh.sp.in{end} }},'in',{{ {} }});
    X = rand_copvar(spX,struct('out',sh.dm.in,'in',1),[0,1],sh.deg,0.7);
    Y = rand_copvar(spY,struct('out',1,'in',1),[0,1],sh.deg,0.7);
    % 'mtimes' and 'plus' require factors on the IDENTICAL registry, and in
    % the cross layout Y's only space is {t} while B's registry is {s,t}.
    % Put the right factors on the full registry; this changes no operator.
    X = on_registry(X,P);       Y = on_registry(Y,P);
    for AB = {{P,Bf,'copvar'},{D,Bf,'cdopvar'},{P,Bd,'cdopvar'},{D,Bd,'cdopvar'}}
        [A,B,wantcls] = deal(AB{1}{:});
        H = [A, B];
        assert(strcmp(class(H),wantcls),'horzcat returned class %s',class(H));
        assert(builtin('numel',H)==1,'horzcat built an object array');
        assert(isequal(size(H),[M, N+1]),'horzcat grid %s',mat2str(size(H)));
        assert(sem_blocks_equal(H,A,1:M,1:N,dmap),'horzcat moved or altered A''s blocks');
        assert(sem_blocks_equal(H,B,1:M,N+1,dmap),'horzcat moved or altered B''s blocks');
        assert(same_spaces_by_name(H,'in',1:N,A,1:N) && same_spaces_by_name(H,'in',N+1,B,1),...
            'horzcat input spaces wrong');
        assert(same_spaces_by_name(H,'out',1:M,A,1:M),'horzcat output spaces wrong');
        lhs = H*[X; Y];
        rhs = A*X + B*Y;
        assert(sem_equal(lhs,rhs,dmap,1),'horzcat: [A,B]*[X;Y] ~= A*X+B*Y');
        v = verify(H);  assert(v.true,'horzcat result fails verify');       % MMP, 09/25/2026
%       nchk = nchk+8;                                                      % MMP, 09/25/2026 (was)
        nchk = nchk+9;                                                      % MMP, 09/25/2026
    end
    % Decision variable list merged when the lists differ, reused when not.
    H = [D, E];
    assert(isequal(sort(H.Zd(:)),sort(unique([D.Zd(:);E.Zd(:)]))),'horzcat Zd is not the union');
    assert(all_blocks_on_list(H),'horzcat left a block off the container list');
    H = [D, D];
    assert(isequal(H.Zd(:),D.Zd(:)),'horzcat of one list did not reuse it');
    % DISJOINT lists, the case two separate declarations produce: the union
    % is the concatenation, and each block is moved by its operand's map.
    Er = rename_dvars(E,'r_');
    dmap = draw_dvals({Er},dmap);
    H = [D, Er];
    assert(isequal(H.Zd(:),[D.Zd(:);Er.Zd(:)]),'disjoint union is not the concatenation');
    assert(all_blocks_on_list(H),'disjoint horzcat left a block off the list');
    % Against the operands' ORIGINAL blocks, evaluated by name. The identity
    % below builds its right side with 'plus', which remaps through the same
    % 'ChangeDecVar' as the concatenation, so on its own it cannot see a
    % remap error: both sides would carry it (CLAUDE.md S4).
    assert(sem_blocks_equal(H,D,1:M,1:N,dmap) && sem_blocks_equal(H,Er,1:M,N+(1:N),dmap),...
        'disjoint horzcat altered a block');
    assert(sem_equal(H*[Z0;Z0],D*Z0 + Er*Z0,dmap,1),'disjoint horzcat: [D,Er]*[Z;Z] ~= D*Z+Er*Z');
    % A list NOT in union order: E with its names reversed. Against D's
    % list the operand's row map is then unsorted, e.g. [4 5 3 2 1].
    Ev = reverse_dvars(E);
    H = [D, Ev];
    assert(all_blocks_on_list(H),'unsorted-map horzcat left a block off the list');
    assert(sem_blocks_equal(H,D,1:M,1:N,dmap) && sem_blocks_equal(H,Ev,1:M,N+(1:N),dmap),...
        'unsorted-map horzcat altered a block');
    assert(sem_equal(H*[Z0;Z0],D*Z0 + Ev*Z0,dmap,1),'unsorted-map horzcat: [D,Ev]*[Z;Z] ~= D*Z+Ev*Z');
    nchk = nchk+3;
    G = blkdiag(D,Er,D);
    assert(all_blocks_on_list(G) && numel(G.Zd)==numel(D.Zd)+numel(Er.Zd),...
        'three-operand blkdiag with a repeated list');
    assert(sem_blocks_equal(G,D,1:M,1:N,dmap) && sem_blocks_equal(G,Er,M+(1:M),N+(1:N),dmap) ...
        && sem_blocks_equal(G,D,2*M+(1:M),2*N+(1:N),dmap),'three-operand blkdiag blocks');
    nchk = nchk+8;
    % Identity-map skip: R2 = D + Er is on [D.Zd; Er.Zd], which IS the union % MMP, 09/25/2026
    % for [D, R2], so only D's blocks move. R3 = Er + D holds the same names % MMP, 09/25/2026
    % permuted: full length, but not the identity map, so it must move.     % MMP, 09/25/2026
    R2 = D + Er;    R3 = Er + D;                                            % MMP, 09/25/2026
    H = [D, R2];                                                            % MMP, 09/25/2026
    assert(isequal(H.Zd(:),R2.Zd(:)) && all_blocks_on_list(H),'identity-map horzcat list'); % MMP, 09/25/2026
    assert(sem_blocks_equal(H,D,1:M,1:N,dmap) && sem_blocks_equal(H,R2,1:M,N+(1:N),dmap),...
        'identity-map horzcat altered a block');                            % MMP, 09/25/2026
    H = [R2, R3];                                                           % MMP, 09/25/2026
    assert(isequal(H.Zd(:),R2.Zd(:)) && all_blocks_on_list(H),'permuted-list horzcat list'); % MMP, 09/25/2026
    assert(sem_blocks_equal(H,R2,1:M,1:N,dmap) && sem_blocks_equal(H,R3,1:M,N+(1:N),dmap),...
        'permuted-list horzcat altered a block');                           % MMP, 09/25/2026
    % A decision operand with NO decision variables - sdopvar blocks on an  % MMP, 09/25/2026
    % empty list - must still be put on the container's list, from either   % MMP, 09/25/2026
    % side, and in a sum.                                                   % MMP, 09/25/2026
    C0 = P.C;                                                               % MMP, 09/25/2026
    for ii = 1:numel(C0)                                                    % MMP, 09/25/2026
        if ~isempty(C0{ii}),    C0{ii} = sopvar2sdopvar(C0{ii});    end     % MMP, 09/25/2026
    end                                                                     % MMP, 09/25/2026
    P0 = cdopvar(C0);                                                       % MMP, 09/25/2026
    assert(isempty(P0.Zd) && isa(P0.C{end,end},'sdopvar'),'setup: P0 should be decision blocks on no list'); % MMP, 09/25/2026
    H1 = [D, P0];   H2 = [P0, D];                                           % MMP, 09/25/2026
    assert(all_blocks_on_list(H1) && all_blocks_on_list(H2) && ...
        isequal(H1.Zd(:),D.Zd(:)) && isequal(H2.Zd(:),D.Zd(:)),'empty-list operand left off the list'); % MMP, 09/25/2026
    assert(sem_blocks_equal(H1,P0,1:M,N+(1:N),dmap) && sem_blocks_equal(H2,P0,1:M,1:N,dmap) ...
        && sem_blocks_equal(H2,D,1:M,N+(1:N),dmap),'empty-list operand altered a block'); % MMP, 09/25/2026
    v1 = verify(H1);    v2 = verify(H2);                                    % MMP, 09/25/2026
    assert(v1.true && v2.true,'empty-list concatenation fails verify');     % MMP, 09/25/2026
    S0 = D + P0;                                                            % MMP, 09/25/2026
    assert(all_blocks_on_list(S0) && sem_lincomb(S0,{D,P0},[1,1],dmap),'D + P0 with an empty-list P0'); % MMP, 09/25/2026
    nchk = nchk+9;                                                          % MMP, 09/25/2026

    % ---------------------------------------------------------- vertcat
    spV = struct('out',{{ sh.sp.out{end} }},'in',{sh.sp.in});
    dmV = struct('out',1,'in',sh.dm.in);
    Vf = rand_copvar(spV,dmV,[0,1],sh.deg,0.7);
    Vd = rand_cdopvar(spV,dmV,[0,1],sh.deg,2,0.7);
    dmap = draw_dvals({Vd},dmap);
    spZ = struct('out',{sh.sp.in},'in',{{ {} }});
    Z = rand_copvar(spZ,struct('out',sh.dm.in,'in',1),[0,1],sh.deg,0.7);
    Z = on_registry(Z,P);
    for AB = {{P,Vf,'copvar'},{D,Vf,'cdopvar'},{P,Vd,'cdopvar'},{D,Vd,'cdopvar'}}
        [A,B,wantcls] = deal(AB{1}{:});
        V = [A; B];
        assert(strcmp(class(V),wantcls),'vertcat returned class %s',class(V));
        assert(builtin('numel',V)==1,'vertcat built an object array');
        assert(isequal(size(V),[M+1, N]),'vertcat grid %s',mat2str(size(V)));
        assert(sem_blocks_equal(V,A,1:M,1:N,dmap),'vertcat moved or altered A''s blocks');
        assert(sem_blocks_equal(V,B,M+1,1:N,dmap),'vertcat moved or altered B''s blocks');
        assert(same_spaces_by_name(V,'out',1:M,A,1:M) && same_spaces_by_name(V,'out',M+1,B,1),...
            'vertcat output spaces wrong');
        VZ = V*Z;   AZ = A*Z;   BZ = B*Z;
        assert(sem_rows_equal(VZ,1:M,AZ,dmap) && sem_rows_equal(VZ,M+1,BZ,dmap),...
            'vertcat: [A;B]*Z rows ~= A*Z, B*Z');
        v = verify(V);  assert(v.true,'vertcat result fails verify');       % MMP, 09/25/2026
%       nchk = nchk+7;                                                      % MMP, 09/25/2026 (was)
        nchk = nchk+8;                                                      % MMP, 09/25/2026
    end

    % ---------------------------------------------------------- blkdiag
    for AB = {{P,Q,'copvar'},{D,Q,'cdopvar'},{P,E,'cdopvar'},{D,E,'cdopvar'}}
        [A,B,wantcls] = deal(AB{1}{:});
        G = blkdiag(A,B);
        assert(strcmp(class(G),wantcls),'blkdiag returned class %s',class(G));
        assert(isequal(size(G),[2*M, 2*N]),'blkdiag grid %s',mat2str(size(G)));
        assert(all(all(cellfun(@isempty,G.C(1:M,N+1:end)))) && ...
               all(all(cellfun(@isempty,G.C(M+1:end,1:N)))),'blkdiag off-diagonal not zero');
        assert(sem_blocks_equal(G,A,1:M,1:N,dmap) && sem_blocks_equal(G,B,M+(1:M),N+(1:N),dmap),...
            'blkdiag moved or altered a diagonal block');
        XY = [Z; Z];
        GXY = G*XY;
        assert(sem_rows_equal(GXY,1:M,A*Z,dmap) && sem_rows_equal(GXY,M+(1:M),B*Z,dmap),...
            'blkdiag: blkdiag(A,B)*[Z;Z] ~= [A*Z; B*Z]');
        v = verify(G);  assert(v.true,'blkdiag result fails verify');       % MMP, 09/25/2026
%       nchk = nchk+5;                                                      % MMP, 09/25/2026 (was)
        nchk = nchk+6;                                                      % MMP, 09/25/2026
    end

    % ---------------------------------------------------------- eq
    assert(P==P && D==D,'eq: an operator is not equal to itself');
    assert(~(P==-P) && ~(D==-D),'eq: P == -P');
    assert((P-P)==0 && (D-D)==0 && 0==(P-P),'eq: P-P is not zero');
    assert(~(P==0) && ~(0==D),'eq: a nonzero operator equals 0');
    assert((2*P)*0.5==P,'eq: (2P)/2 ~= P');
    assert(cdopvar(P)==P && P==cdopvar(P),'eq: promotion changed the operator');
    % A [] block against an explicit zero block of the same spaces.
    [i0,j0] = find(cellfun(@isempty,P.C),1);
    if ~isempty(i0)
        % An explicit zero block with that row's output space and that
        % column's input space.
        Pz = zero_block_like(P,i0,j0,0);
        Cz = P.C;   Cz{i0,j0} = Pz;
        assert(copvar(Cz,metadata(P))==P,'eq: [] ~= explicit zero block');
        % ...and against a NONZERO block there it must be unequal, from
        % either side: a [] block is a zero test, not a wildcard.
        Cn = P.C;   Cn{i0,j0} = zero_block_like(P,i0,j0,1);
        Pn = copvar(Cn,metadata(P));
        assert(~(P==Pn) && ~(Pn==P),'eq: [] block compared equal to a nonzero block');
        nchk = nchk+2;
    end
    % Perturb one coefficient: unequal, but equal at a coarse tolerance.
    Pp = perturb(P,1e-9);
    assert(~(Pp==P) && eq(Pp,P,1e-6),'eq: tolerance handling');
    Dp = perturb(D,1e-9);
    assert(~(Dp==D) && eq(Dp,D,1e-6),'eq: tolerance handling (decision)');
    % A merged registry against one re-derived from the blocks alone: the
    % same operator, so equal, although reached by different routes.
    Hm = [P, Bf];
    assert(Hm==copvar(Hm.C),'eq: merged-registry container ~= rederived one');
    nchk = nchk+10;
    fprintf('  passed: %-26s (%d checks so far)\n',sh.name,nchk);
end

% -------------------------------------------------------------- registries
% [A, B] where B's registry is empty (R^2 <- R^3): the merged registry is
% A's, and the space masks follow the names.
A = rand_copvar(struct('out',{{ {} }},'in',{{ {'s'} }}),struct('out',2,'in',1),[0,1],2,1);
B = rand_copvar(struct('out',{{ {} }},'in',{{ {} }}),struct('out',2,'in',3),[0,1],2,1);
assert(isempty(B.vars),'setup: B should have an empty registry');
H = [A, B];
assert(isequal(H.vars,{'s'}) && isequal(H.space_in,[true;false]) && isequal(H.dim_in,[1;3]),...
    'registry merge with an empty registry');
% Disjoint registries {s} and {t}.
A = rand_copvar(struct('out',{{ {} }},'in',{{ {'s'} }}),struct('out',1,'in',1),[0,1],2,1);
B = rand_copvar(struct('out',{{ {} }},'in',{{ {'t'} }}),struct('out',1,'in',1),[0,1],2,1);
H = [A, B];
assert(isequal(H.vars,{'s','t'}) && isequal(H.space_in,[true false; false true]),...
    'registry merge of disjoint registries');
V = [A'; B'];                                               % rows L2[s], L2[t]
assert(isequal(V.vars,{'s','t'}) && isequal(V.space_out,[true false; false true]),...
    'vertcat registry merge of disjoint registries');
G = blkdiag(A,B);
assert(isequal(G.vars,{'s','t'}) && isequal(G.space_in,[true false; false true]),...
    'blkdiag registry merge');
% Two ZERO operators with the same grid and dimensions but different input
% spaces, R^1 <- L2[s] and R^1 <- L2[t]: equal as numbers, unequal as
% operators, so eq must compare the spaces and not only the blocks.
assert(~((A-A)==(B-B)),'eq: zero operators on different spaces compared equal');
assert(~(A==copvar(A.C,setfield(metadata(A),'dom',[0 2]))),...
    'eq: same blocks, different domain, compared equal');                  %#ok<SFLD>
nchk = nchk+7;

% -------------------------------------------------------------- blocks
% Two bare blocks are the BLOCK classes' business ([S, S] is @sopvar/horzcat,
% one wider sopvar); a container among the operands makes it a container
% concatenation, whichever side the container is on.
S = P.C{end,end};   % a sopvar block
Sd = D.C{end,end};  % an sdopvar block
Pl = copvar({S});
H = [Pl, S];
assert(isa(H,'copvar') && isequal(size(H),[1 2]),'[copvar, sopvar] should be a 1x2 copvar');
H = [S, Pl];
assert(isa(H,'copvar') && isequal(size(H),[1 2]),'[sopvar, copvar] should be a 1x2 copvar');
H = [Pl, Sd];                                               % reaches copvar/horzcat, reroutes
assert(isa(H,'cdopvar') && isequal(size(H),[1 2]),'[copvar, sdopvar] should reroute to cdopvar');
H = [Sd; Pl];
assert(isa(H,'cdopvar') && isequal(size(H),[2 1]),'[sdopvar; copvar] should reroute to cdopvar');
G = blkdiag(Sd,Pl);
assert(isa(G,'cdopvar') && isequal(size(G),[2 2]),'blkdiag(sdopvar, copvar) should be a cdopvar');
H = [Pl, []];
assert(H==Pl,'[P, []] ~= P');
nchk = nchk+6;

% -------------------------------------------------------------- errors
sp1 = struct('out',{{ {}, {'s'} }},'in',{{ {} }});
A = rand_copvar(sp1,struct('out',[2;3],'in',1),[0,1],1,1);
B = rand_copvar(struct('out',{{ {} }},'in',{{ {} }}),struct('out',2,'in',1),[0,1],1,1);
expect_error(@() [A, B],'copvar:horzcatRowMismatch');
B = rand_copvar(struct('out',{{ {}, {'t'} }},'in',{{ {} }}),struct('out',[2;3],'in',1),[0,1],1,1);
expect_error(@() [A, B],'copvar:horzcatSpaceMismatch');
B = rand_copvar(sp1,struct('out',[2;4],'in',1),[0,1],1,1);
expect_error(@() [A, B],'copvar:horzcatDimMismatch');
B = rand_copvar(sp1,struct('out',[2;3],'in',1),[0,2],1,1);
expect_error(@() [A, B],'copvar:domConflict');
expect_error(@() [A', B'],'copvar:domConflict');
Bt = rand_copvar(struct('out',{{ {} }},'in',{{ {'s'} }}),struct('out',2,'in',1),[0,1],1,1);
expect_error(@() [A; Bt],'copvar:vertcatSpaceMismatch');   % R^1 column vs L2[s]
Bt = rand_copvar(struct('out',{{ {} }},'in',{{ {}, {'s'} }}),struct('out',2,'in',[1;1]),[0,1],1,1);
expect_error(@() [A; Bt],'copvar:vertcatColMismatch');     % 1 column vs 2
Bt = rand_copvar(struct('out',{{ {} }},'in',{{ {} }}),struct('out',2,'in',2),[0,1],1,1);
expect_error(@() [A; Bt],'copvar:vertcatDimMismatch');     % R^1 vs R^2
expect_error(@() [A, 3],'copvar:horzcatBadOperand');
expect_error(@() A==3,'eq:badInput');
nchk = nchk+10;
% The row map the container code hands 'setdvars': refused at the wrong     % MMP, 09/25/2026
% length or out of range; a right one moves each row with its name. Bk is   % MMP, 09/25/2026
% a block of the last shape's D.                                            % MMP, 09/25/2026
Bk = D.C{end,end};      n = numel(Bk.Zd);                                   % MMP, 09/25/2026
expect_error(@() setdvars(Bk,Bk.Zd,1:n-1),'ChangeDecVar:badLoc');           % MMP, 09/25/2026
expect_error(@() setdvars(Bk,Bk.Zd,[0,2:n]),'ChangeDecVar:badLoc');         % MMP, 09/25/2026
Zr = [{'x_new'}; flipud(Bk.Zd(:))];     dmap('x_new') = randn;              % MMP, 09/25/2026
Bm = setdvars(Bk,Zr,n+1:-1:2);          % entry k of Bk.Zd is at n+2-k      % MMP, 09/25/2026
Kk = kern(Bk,dmap);                                                         % MMP, 09/25/2026
assert(isequal(Bm.Zd(:),Zr) && kdiff(kern(Bm,dmap),Kk,1) <= 1e-9*kscale(Kk),...
    'setdvars with a map moved a row off its name');                        % MMP, 09/25/2026
nchk = nchk+3;                                                              % MMP, 09/25/2026

fprintf('test_copvar_blockops passed (%d checks).\n',nchk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dmap = draw_dvals(Ps,dmap)
% One random value per decision variable NAME, kept across calls.
if nargin<2,    dmap = containers.Map('KeyType','char','ValueType','double');  end
for k = 1:numel(Ps)
    z = Ps{k}.Zd;
    for i = 1:numel(z)
        if ~isKey(dmap,z{i}),   dmap(z{i}) = randn;     end
    end
end
end


function K = kern(B,dmap)
% Kernels of one block, a decision block at the drawn values.
if isa(B,'sdopvar')
    z = B.Zd;
    dval = zeros(numel(z),1);
    for i = 1:numel(z),     dval(i) = dmap(z{i});   end
    K = pi_blk_kernels(B,dval);
else
    K = pi_blk_kernels(B,[]);
end
end


function d = kdiff(K1,K2,c)
% max |K1 - c*K2| over every gamma cell; K2 may be {} (a zero block).
d = 0;
for g = 1:numel(K1)
    if isempty(K2),     X = K1{g};
    else,               X = K1{g} - c*K2{g};
    end
    d = max(d,maxcoef(X));
end
end


function m = maxcoef(X)
if isnumeric(X)
    m = full(max(abs(X(:))));
else
    X = polynomial(X);
    cf = X.coefficient;
    m = full(max(abs(cf(:))));
end
if isempty(m),  m = 0;  end
end


function s = kscale(K)
s = 0;
for g = 1:numel(K),     s = max(s,maxcoef(K{g}));    end
s = max(s,1);
end


function tf = sem_equal(R,X,dmap,c)
% Every block of R equals c times the same block of X, as kernels.
tf = isequal(size(R),size(X));
if ~tf,     return,     end
for ii = 1:numel(R.C)
    r = R.C{ii};    x = X.C{ii};
    if isempty(r) && isempty(x),    continue,   end
    if isempty(r)
        Kx = kern(x,dmap);      ok = kdiff(Kx,{},1) <= 1e-9*kscale(Kx);
    elseif isempty(x)
        Kr = kern(r,dmap);      ok = kdiff(Kr,{},1) <= 1e-9*kscale(Kr);
    else
        Kr = kern(r,dmap);      Kx = kern(x,dmap);
        ok = numel(Kr)==numel(Kx) && kdiff(Kr,Kx,c) <= 1e-9*max(kscale(Kr),kscale(Kx));
    end
    if ~ok,     tf = false;     return,     end
end
end


function tf = sem_lincomb(R,Xs,cs,dmap)
% Every block of R equals sum_k cs(k)*Xs{k} blockwise, as kernels.
tf = true;
for ii = 1:numel(R.C)
    Kr = {};
    if ~isempty(R.C{ii}),   Kr = kern(R.C{ii},dmap);    end
    Ks = cell(1,numel(Xs));
    for k = 1:numel(Xs)
        if ~isempty(Xs{k}.C{ii}),   Ks{k} = kern(Xs{k}.C{ii},dmap);    end
    end
    ng = max([numel(Kr),cellfun(@numel,Ks)]);
    sc = 1;
    for g = 1:ng
        acc = 0;
        for k = 1:numel(Xs)
            if ~isempty(Ks{k}),     acc = acc + cs(k)*Ks{k}{g};     sc = max(sc,maxcoef(Ks{k}{g}));   end
        end
        if ~isempty(Kr),    acc = Kr{g} - acc;      else,   acc = -acc;     end
        if maxcoef(acc) > 1e-9*sc,  tf = false;     return,     end
    end
end
end


function tf = sem_blocks_equal(H,A,rows,cols,dmap)
% H.C(rows,cols) holds A's blocks, semantically.
tf = true;
for a = 1:numel(rows)
    for b = 1:numel(cols)
        h = H.C{rows(a),cols(b)};   x = A.C{a,b};
        if isempty(h) ~= isempty(x),    tf = false;     return,     end
        if isempty(h),  continue,   end
        Kh = kern(h,dmap);      Kx = kern(x,dmap);
        if numel(Kh)~=numel(Kx) || kdiff(Kh,Kx,1) > 1e-9*kscale(Kx)
            tf = false;     return
        end
    end
end
end


function tf = sem_rows_equal(H,rows,A,dmap)
% Block rows 'rows' of H equal the rows of A, semantically.
tf = size(H.C,2)==size(A.C,2) && numel(rows)==size(A.C,1);
if ~tf,     return,     end
for a = 1:numel(rows)
    for b = 1:size(A.C,2)
        h = H.C{rows(a),b};     x = A.C{a,b};
        if isempty(h) && isempty(x),    continue,   end
        if isempty(h),  Kx = kern(x,dmap);  ok = kdiff(Kx,{},1)<=1e-9*kscale(Kx);
        elseif isempty(x),  Kh = kern(h,dmap);  ok = kdiff(Kh,{},1)<=1e-9*kscale(Kh);
        else
            Kh = kern(h,dmap);  Kx = kern(x,dmap);
            ok = numel(Kh)==numel(Kx) && kdiff(Kh,Kx,1) <= 1e-9*max(kscale(Kh),kscale(Kx));
        end
        if ~ok,     tf = false;     return,     end
    end
end
end


function tf = same_spaces_by_name(H,side,hidx,A,aidx)
% Spaces of H at hidx equal those of A at aidx, compared by variable name.
if strcmp(side,'in')
    mh = H.space_in(hidx,:);    ma = A.space_in(aidx,:);
    dh = H.dim_in(hidx);        da = A.dim_in(aidx);
else
    mh = H.space_out(hidx,:);   ma = A.space_out(aidx,:);
    dh = H.dim_out(hidx);       da = A.dim_out(aidx);
end
tf = isequal(dh(:),da(:));
for k = 1:numel(hidx)
    tf = tf && isequal(sort(H.vars(mh(k,:))),sort(A.vars(ma(k,:))));
end
end


function X = on_registry(X,R)
% The same operator X, with its metadata restated over R's variable
% registry, a superset of X's. Blocks are untouched.
meta = metadata(X);
[~,loc] = ismember(X.vars,R.vars);
so = false(size(X.C,1),numel(R.vars));  so(:,loc) = X.space_out;
si = false(size(X.C,2),numel(R.vars));  si(:,loc) = X.space_in;
meta.vars = R.vars;     meta.dom = R.dom;
meta.space_out = so;    meta.space_in = si;
if isa(X,'cdopvar'),    X = cdopvar(X.C,meta);
else,                   X = copvar(X.C,meta);
end
end


function P = rename_dvars(P,prefix)
% The same operator in renamed decision variables: prefix every name, in the
% container list and in each block, consistently, so the lists stay shared.
Zd = strcat(prefix,P.Zd(:));
for ii = 1:numel(P.C)
    if isa(P.C{ii},'sdopvar'),  P.C{ii}.Zd = Zd;   end
end
P.Zd = Zd;
end


function P = reverse_dvars(P)
% Relabel the decision variables by reversing the NAME list, in the
% container and every block together: row r of each B now belongs to the
% name that was at position end+1-r. A different, consistently labelled
% operator, whose list is out of order against any lexically built union.
Zd = flipud(P.Zd(:));
for ii = 1:numel(P.C)
    if isa(P.C{ii},'sdopvar'),  P.C{ii}.Zd = Zd;   end
end
P.Zd = Zd;
end


function tf = all_blocks_on_list(P)
tf = true;
for ii = 1:numel(P.C)
    if isa(P.C{ii},'sdopvar') && ~isequal(P.C{ii}.Zd(:),P.Zd(:))
        tf = false;     return
    end
end
end


function Z = zero_block_like(P,i0,j0,scale)
% A 'sopvar' block from row i0's output space to column j0's input space,
% scaled by 'scale': 0 gives a zero block, 1 a random nonzero one.
vo = P.vars(P.space_out(i0,:));     vi = P.vars(P.space_in(j0,:));
[~,io] = ismember(vo,P.vars);       [~,ii] = ismember(vi,P.vars);
vars = struct('out',{vo},'in',{vi});
dom  = struct('out',P.dom(io,:),'in',P.dom(ii,:));
degs = struct('out',ones(1,numel(vo)),'in',ones(1,numel(vi)));
Z = scale*rand_sopvar([P.dim_out(i0),P.dim_in(j0)],vars,dom,degs,1);
end


function X = perturb(X,delta)
% Add delta to one stored coefficient of the last populated block.
k = find(~cellfun(@isempty,X.C),1,'last');
B = X.C{k};
if isa(B,'sdopvar')
    ip = find(cellfun(@(a) nnz(a)>0,B.params.A),1);
    if isempty(ip)
        ip = find(cellfun(@(b) nnz(b)>0,B.params.B),1);
        [r,c] = find(B.params.B{ip},1);
        B.params.B{ip}(r,c) = B.params.B{ip}(r,c) + delta;
    else
        l = find(B.params.A{ip},1);
        B.params.A{ip}(l) = B.params.A{ip}(l) + delta;
    end
else
    ip = find(cellfun(@(a) nnz(a)>0,B.params),1);
    l = find(B.params{ip},1);
    B.params{ip}(l) = B.params{ip}(l) + delta;
end
X.C{k} = B;
end


function expect_error(f,id)
try
    f();
catch ME
    if ~strcmp(ME.identifier,id)
        error('test_copvar_blockops: expected error ''%s'', got ''%s'': %s',id,ME.identifier,ME.message);
    end
    return
end
error('test_copvar_blockops: expected error ''%s'', but none was raised.',id);
end
