% Test_tensopmat_is_id_op.m -- tests @tensopmat/is_id_op.m.
%
% Expected values are derived from the kind of each block (identity, zero,
% scaled, ...), never from is_id_op itself. For blocks without integral
% terms, the result is also checked semantically: apply_op(Cop,X)
% must equal X for random X exactly when Cop is flagged identity. Empty
% blocks are zero operators. Swept over block counts (1-4), block
% dimensions, spatial variables per diagonal block (0 or 1, mixed), block
% types ('tensopvar', 'intop') and tolerances.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Requires PIESOS_2D and PIETOOLS on the MATLAB path.

warning('off','sopvar:noncanonicalMultiplier');
nchk = 0;
va = struct('in',{{'s1'}},'out',{{'s1'}});  da = struct('in',[0 1],'out',[0 1]);
Ia = make_id(1,va,da);     I2a = make_id(2,va,da);
Z = @(P) 0*P;
Tv = @(P) tensopvar(P);

%% ---------------------------------------------------------- hand-made cases
% single block
nchk = nchk+1;  chk(is_id_op(tensopmat(Tv(Ia))),true,"1x1 id");
nchk = nchk+1;  chk(is_id_op(tensopmat(Tv(I2a))),true,"1x1 id p=2");
nchk = nchk+1;  chk(is_id_op(tensopmat(Tv(2*Ia))),false,"1x1 2*id");
nchk = nchk+1;  chk(is_id_op(tensopmat(Tv(Z(Ia)))),false,"1x1 zero");
nchk = nchk+1;  chk(is_id_op(tensopmat(intop(eye(2)))),true,"1x1 intop eye(2)");
nchk = nchk+1;  chk(is_id_op(tensopmat(intop(2*eye(2)))),false,"1x1 intop 2*eye(2)");
% block diagonal, empty off-diagonal blocks (zero operators)
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];[],Tv(I2a)})),true,"diag(id,id2), empty off-diag");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];[],Tv(2*I2a)})),false,"diag(id,2*id2)");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];[],Tv(Z(I2a))})),false,"diag(id,0)");
% zero off-diagonal blocks, explicit
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),Tv(Z(Ia));Tv(Z(Ia)),Tv(Ia)})),true,"diag(id,id), explicit zero off-diag");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];Tv(Z(Ia)),Tv(Ia)})),true,"id with one explicit zero below");
% nonzero off-diagonal block
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),Tv(Ia);[],Tv(Ia)})),false,"upper triangular id");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];Tv(Ia),Tv(Ia)})),false,"lower triangular id");
% empty diagonal blocks: permutation, and partially empty diagonal
nchk = nchk+1;  chk(is_id_op(tensopmat({[],Tv(Ia);Tv(Ia),[]})),false,"permutation");
nchk = nchk+1;  chk(is_id_op(tensopmat({[],Tv(Ia);Tv(Ia),Tv(Ia)})),false,"empty (1,1) block");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),Tv(Ia);Tv(Ia),[]})),false,"empty (2,2) block");
% non-square block structure
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),Tv(Ia)})),false,"1x2");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia);Tv(Ia)})),false,"2x1");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),Tv(Ia),[];[],Tv(Ia),Tv(Ia)})),false,"2x3");
% square block structure, rectangular block (dims.in ~= dims.out)
Prect = sopvar({sparse(1,2),sparse(1,2),sparse(1,2)},va,{0},{0},da,[1 2]);
Prect.params{1} = sparse([1 0]);
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];[],Tv(Prect)})),false,"rectangular diagonal block");
% mixed tensopvar / intop diagonal; functional below the diagonal
Kz = intop({polynomial(zeros(2,1))},{[0 2]},{'s1'},[0 1]);
Kn = intop({polynomial(ones(2,1))},{[0 2]},{'s1'},[0 1]);
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];[],intop(eye(2))})),true,"diag(tensopvar id, intop id)");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];[],intop(2*eye(2))})),false,"diag(tensopvar id, intop 2*id)");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];Kz,intop(eye(2))})),true,"zero functional in off-diag");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];Kn,intop(eye(2))})),false,"nonzero functional in off-diag");
% tolerance
tiny = 1e-9;
Ct = tensopmat({Tv(Ia),Tv(tiny*Ia);[],Tv(Ia)});
nchk = nchk+1;  chk(is_id_op(Ct),false,"off-diag 1e-9 default tol");
nchk = nchk+1;  chk(is_id_op(Ct,1e-6),true,"off-diag 1e-9 ztol=1e-6");
Ct = tensopmat({Tv(Ia+tiny*Ia),[];[],Tv(Ia)});
nchk = nchk+1;  chk(is_id_op(Ct),false,"diag id+1e-9 default tol");
nchk = nchk+1;  chk(is_id_op(Ct,1e-6),true,"diag id+1e-9 ztol=1e-6");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia+1e-14*Ia),[];[],Tv(Ia)})),true,"diag id+1e-14 default tol");
% diagonal blocks acting on different variables, and the same variable twice
vb = struct('in',{{'s2'}},'out',{{'s2'}});  Ib = make_id(1,vb,da);
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[];[],Tv(Ib)})),true,"diag(id_s1,id_s2)");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[],[];[],Tv(Ib),[];[],[],intop(eye(3))})),true,"diag(id_s1,id_s2,intop eye 3)");
nchk = nchk+1;  chk(is_id_op(tensopmat({Tv(Ia),[],[];[],Tv(Ib),[];[],[],intop(diag([1 1 0])+0)})),false,"diag(id_s1,id_s2,intop singular)");
% tensor-product block
Tab = sopvar2tensopvar({Ia,Ib},[2,2]);
nchk = nchk+1;  chk(is_id_op(tensopmat({Tab,[];[],Tv(Ia)})),true,"diag(id(x)id, id)");
nchk = nchk+1;  chk(is_id_op(tensopmat({sopvar2tensopvar({Ia,2*Ib},[2,2]),[];[],Tv(Ia)})),false,"diag(id(x)2id, id)");

%% ------------------------------------------------------ randomized sweep
% Every block is p x p on s1 (identity, zero, scaled, tiny, or id+integral),
% or empty. Expected value from the block kinds alone.
rng(1009);
kinds = {'id','zero','scaled','tiny','int','empty'};
nexp_true = 0;  nsem = 0;
for trial=1:200
    nb = randi(4);                  % number of row/column blocks
    p = randi(2);                   % common block dimension
    kd = repmat({'empty'},nb,nb);
    for i=1:nb
        for j=1:nb
            if i==j
                if rand<0.7,  kd{i,j} = 'id';  else,  kd{i,j} = kinds{randi(5)};  end
            else
                if rand<0.5,  kd{i,j} = 'empty';  else,  kd{i,j} = kinds{randi(6)};  end
            end
        end
    end
    % initialize requires a nonempty block in every row and column
    ok = all(any(~strcmp(kd,'empty'),1)) && all(any(~strcmp(kd,'empty'),2));
    if ~ok,  continue,  end
    ops = cell(nb,nb);
    for i=1:nb
        for j=1:nb
            if ~strcmp(kd{i,j},'empty')
                ops{i,j} = Tv(make_block(kd{i,j},p));
            end
        end
    end
    Cop = tensopmat(ops);
    for ztol = [1e-12,1e-6]
        eff = kd;
        eff(strcmp(kd,'tiny')) = {'scaled'};
        if ztol>1e-9,  eff(strcmp(kd,'tiny')) = {'zero'};  end
        isdiag = logical(eye(nb));
        expected = all(strcmp(eff(isdiag),'id')) && all(ismember(eff(~isdiag),{'empty','zero'}));
        nexp_true = nexp_true+expected;
        nchk = nchk+1;
        chk(is_id_op(Cop,ztol),expected,"random trial="+trial+" ztol="+ztol+" kinds="+strjoin(reshape(kd',1,[]),','));
    end
    % semantic oracle: action on random input, for operators without integrals
    if ~any(strcmp(kd(:),'int')) && ~any(strcmp(kd(:),'tiny'))
        X = cell(1,nb);
        for j=1:nb
            X{j} = rand_poly([p,1],polynomial({'s1'}),2,20);
        end
        Y = apply_op(Cop,X);
        err = cleanpoly(Y-vertcat(X{:}),1e-10);
        acts_id = max(abs(full(double(err.C(:)))))<1e-8;
        nsem = nsem+1;
        nchk = nchk+1;
        chk(is_id_op(Cop),acts_id,"semantic trial="+trial+" kinds="+strjoin(reshape(kd',1,[]),','));
    end
end
fprintf('  randomized checks: %d expected identity (of %d tensopmat/ztol pairs), %d semantic comparisons\n',nexp_true,2*nsem,nsem);
assert(nexp_true>10 && nsem>10,"randomized sweep too weak to be meaningful");

fprintf('Test_tensopmat_is_id_op passed (%d checks).\n',nchk);


%% -------------------------------------------------------------- helpers
function chk(val,expected,name)
if ~isequal(logical(val),logical(expected))
    error("Test_tensopmat_is_id_op failed: %s (got %d, expected %d).",name,val,expected);
end
end

function Iop = make_id(p,vars,dom)
% Canonical identity operator I_p on L2^p[dom]: multiplier in every variable.
N = numel(vars.in);
params = [{speye(p)},repmat({sparse(p,p)},1,3^N-1)];
Iop = sopvar(params,vars,repmat({0},1,N),repmat({0},1,N),dom,[p p]);
end

function P = make_block(kind,p)
% p x p operator on s1 in [0,1] of the given kind
vars = struct('in',{{'s1'}},'out',{{'s1'}});    dom = struct('in',[0 1],'out',[0 1]);
Iop = make_id(p,vars,dom);
switch kind
    case 'id',      P = Iop;
    case 'zero',    P = 0*Iop;
    case 'scaled',  P = 2*Iop;
    case 'tiny',    P = 1e-9*Iop;
    case 'int',     P = Iop;    P.params{2} = speye(p);
end
end
