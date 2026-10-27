function test_container_plus_batch()
% TEST_CONTAINER_PLUS_BATCH checks the container n-ary sum of 10/01/2026,
% '@copvar/plus_batch' (the chain) and '@cdopvar/plus_batch' (one list
% reconciliation), against the chain P1 + P2 + ... + PN, tested against the kernel
% definition elsewhere (test_copvar_blockops, test_copvar_silent_fixes).
% The two must agree BIT FOR BIT - every block's parameters, decision list,
% bases and the container metadata - because the n-ary sum is the chain's
% reconciliation done once, with the blocks added in the same order.
%
% Operands: positive containers from 'poscopvar' (each with its own
% decision variables, as the heatNd face terms), the same operand repeated
% (identical lists), mixed with fixed 'copvar' operands, on a 1-D grid
% R x L2[s] and a 2-D grid L2[s] x L2[s,t], N = 2..5; and errors for a
% grid or space mismatch, as '+' raises them.
%
% Initial coding MMP, 10/01/2026

warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
rng(20261001);
nchk = 0;
cases = { struct('spaces',{{ {}, {'s'} }},'dims',[1;2],'dom',[0 1],'deg',1), ...
          struct('spaces',{{ {'s'}, {'s','t'} }},'dims',[1;1],'dom',[0 1; -1 1],'deg',1) };
for c = 1:numel(cases)
    S = cases{c};
    vars = unique([S.spaces{:}]);
    prog = lpiprogram_sop(vars,S.dom);
    P = cell(1,5);
    for k = 1:5
        [prog,P{k}] = poscopvar(prog,S.dims,S.spaces,S.dom,S.deg);
    end
    F = rand_copvar(struct('out',{S.spaces},'in',{S.spaces},'vars',{vars}), ...
                    struct('out',S.dims,'in',S.dims), ...
                    S.dom,1,0.6);
    sets = { P(1:2), P(1:3), P(1:5), [P(1),P(2),P(1)], [P(3),{F},P(4)], [{F},P(1:2)], ...
             [P(1),P(1)], {F,F} };
    for k = 1:numel(sets)
        ops = sets{k};
        Cb = plus_batch(ops{:});
        Cc = ops{1};
        for t = 2:numel(ops),   Cc = Cc + ops{t};   end
        lbl = sprintf('case %d set %d (%d operands)',c,k,numel(ops));
        assert(strcmp(class(Cb),class(Cc)),'%s: class %s vs %s',lbl,class(Cb),class(Cc));
        same_container(Cb,Cc,lbl);
        nchk = nchk+1;
    end
    % one operand returns itself
    assert(isequal(plus_batch(P{1}),P{1}),'case %d: one operand',c);
    nchk = nchk+1;
end
% errors, as '+' raises them
S = cases{1};
prog = lpiprogram_sop(unique([S.spaces{:}]),S.dom);
[prog,A] = poscopvar(prog,S.dims,S.spaces,S.dom,1);
[~,B] = poscopvar(prog,1,{{'s'}},S.dom,1);              % a 1 x 1 grid
assert(strcmp(errid(@() plus_batch(A,A,B)),errid(@() A+B)),'grid mismatch id');
assert(strcmp(errid(@() plus_batch(A,3)),'plus_batch:badInput'),'numeric summand');
nchk = nchk+2;
fprintf('test_container_plus_batch passed (%d checks).\n',nchk);
end


function same_container(A,B,lbl)
% Bit-for-bit equality of two containers: metadata, grid, every block.
for f = {'vars','dom','space_out','space_in','dim_out','dim_in'}
    assert(isequal(A.(f{1}),B.(f{1})),'%s: %s differs',lbl,f{1});
end
if isa(A,'cdopvar'),    assert(isequal(A.Zd(:),B.Zd(:)),'%s: Zd differs',lbl);   end
assert(isequal(size(A.C),size(B.C)),'%s: grid differs',lbl);
for ii = 1:numel(A.C)
    a = A.C{ii};    b = B.C{ii};
    assert(isequal(class(a),class(b)),'%s: block %d class',lbl,ii);
    if isempty(a),  continue,   end
    assert(isequal(a.vars,b.vars) && isequal(a.dom,b.dom) && isequal(a.dims,b.dims) ...
           && isequal(a.ZL,b.ZL) && isequal(a.ZR,b.ZR),'%s: block %d metadata',lbl,ii);
    if isa(a,'sdopvar')
        assert(isequal(a.Zd(:),b.Zd(:)) && isequal(a.params.A,b.params.A) ...
               && isequal(a.params.B,b.params.B),'%s: block %d parameters',lbl,ii);
    else
        assert(isequal(a.params,b.params),'%s: block %d parameters',lbl,ii);
    end
end
end

function id = errid(f)
id = '';
try,        f();
catch ME,   id = ME.identifier;
end
end
