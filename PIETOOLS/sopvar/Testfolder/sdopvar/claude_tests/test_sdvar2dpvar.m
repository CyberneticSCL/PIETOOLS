%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_SDVAR2DPVAR checks 'sdvar2dpvar' (sopvar/Testfolder/converters)
% against the definition in its header,
%
%   D(s,t) = (Im o ZL(s))' unvec(A + B'd) (In o ZR(t)),
%
% ZL(s) = s1.^ZL{1} o ... o sM.^ZL{M} (first variable slowest), ZR(t)
% likewise; unvec to (m*nZL) x (n*nZR), column-major, component outer on
% both axes. Layout taken from the sdopvar class header and confirmed by
% the ANCHOR check, not assumed.
%
% Per case, against the semantics and never a round trip (CLAUDE.md S4):
%   ANCHOR  the numeric reference 'ref_eval' equals the kernel of an
%           'sdopvar' built from the same (A,B), as evaluated by
%           'pi_blk_kernels', on bases of unequal size (a ZL/ZR swap or
%           a transposed vec would show);
%   HAND    the 2x2 case from the defect report, value worked by hand;
%   VALUE   the output dpvar, expanded by the class's own 'dpvar2poly' and
%           evaluated by 'subs' at random (d,s,t), equals 'ref_eval',
%           computed directly from A, B and the degree lists. Neither
%           'sdvar2dpvar' nor 'dpvar2sdvar' enters the reference.
% Sweep: m,n in {1,2,3} x 0/1/2 output x 0/1/2 input variables x 1, 3 or 6
% decision variables (243 cases); degree lists random subsets of
% 0:3 (gaps and lists without a constant included); sparse A and B, every
% decision variable present. Plus cases with no decision variables.
%
% Defects this catches (fixed 09/26/2026): the constant row of matrix row r
% was placed at (r-1)*nZd+1, not (r-1)*(nZd+1)+1, wrong for m >= 2; and a
% single decision variable with two or more nonzeros crashed, since 'find'
% on a row B returns rows.
%
% Initial coding MMP, 09/26/2026

tol = 1e-10;        % relative, Frobenius
nsamp = 3;          % random (d,s,t) per case
nchk = 0;           fails = {};

%% ANCHOR: the reference uses the sdopvar class's own layout
ANCH = {  % m, n, ZL, ZR, nZd
    2, 3, {[0;1;3],2},   {[1;2]},        4
    3, 2, {[0;2]},       {[1;2;3],[0;1]}, 2
    2, 2, {[1;3],[0;2]}, {[0;1;2],3},    5
};
for ia = 1:size(ANCH,1)
    rng(500+ia);
    [m,n,ZL,ZR,nZd] = deal(ANCH{ia,:});
    [P,vars,Zd] = rand_coeffs(m,n,ZL,ZR,nZd);
    dom = struct('in',repmat([0,1],numel(ZR),1),'out',repmat([0,1],numel(ZL),1));
    Pobj = sdopvar(struct('A',{{P.A}},'B',{{P.B}}),vars,Zd,ZL,ZR,dom,[m,n]);
    d = randn(nZd,1);
    [~,loc] = ismember(Pobj.Zd(:),Zd(:));
    Kc = pi_blk_kernels(Pobj,d(loc));
    assert(numel(Kc)==1,'anchor %d: no shared variables, so one kernel',ia);
    s = rand_pt(numel(ZL));     t = rand_pt(numel(ZR));
    K = eval_poly(Kc{1},[vars.out(:);strcat(vars.in(:),'_dum')],[s;t]);
    ref = ref_eval(P,[m,n],ZL,ZR,d,s,t);
    err = norm(K-ref,'fro')/max(1,norm(ref,'fro'));
    assert(err<tol,'anchor %d: reference disagrees with pi_blk_kernels (%.2e)',ia,err);
    nchk = nchk+1;
end
fprintf('ANCHOR: reference matches pi_blk_kernels on %d cases.\n',size(ANCH,1));

%% HAND: 2x2, no monomials, d1 at (1,1), d3 at (2,1), d2 at (1,2)
P = struct('A',sparse([1;0;0;4]),'B',sparse([1 2 3],[1 3 2],1,3,4));
vars = struct('out',{cell(1,0)},'in',{cell(1,0)});
Zd = {'d1';'d2';'d3'};
try
    D = sdvar2dpvar(P,[2,2],vars,cell(1,0),cell(1,0),Zd);
    got = eval_poly(dpvar2poly(D),Zd,[10;20;30]);
    msg = '';
    if ~isequal(got,[11 20; 30 4])
        msg = ['got ',mat2str(got),', expected [11 20;30 4]'];
    end
catch ME
    msg = ['error: ',ME.message];
end
nchk = nchk+1;
if ~isempty(msg),   fails{end+1} = ['HAND: ',msg];  end

%% VALUE: sweep over dimensions, variable counts and degrees
ncase = 0;
for m = 1:3
for n = 1:3
for M = 0:2
for N = 0:2
for nZd = [1,3,6]
    ncase = ncase+1;
    rng(1000+ncase);
    ZL = arrayfun(@(~)rand_degs(),1:M,'un',0);
    ZR = arrayfun(@(~)rand_degs(),1:N,'un',0);
    [P,vars,Zd] = rand_coeffs(m,n,ZL,ZR,nZd);
    msg = check_case(P,[m,n],vars,ZL,ZR,Zd,nsamp,tol);
    nchk = nchk+nsamp;
    if ~isempty(msg)
        fails{end+1} = sprintf('m=%d n=%d M=%d N=%d nZd=%d nnzB=%d: %s',...
                               m,n,M,N,nZd,nnz(P.B),msg);
    end
end
end
end
end
end

% No decision variables: B is 0 x ncoeffs, D is a fixed matrix
NOD = [1 1 0 0; 2 3 1 2; 3 1 2 0; 2 2 2 2];
for k = 1:size(NOD,1)
    rng(2000+k);
    [m,n,M,N] = deal(NOD(k,1),NOD(k,2),NOD(k,3),NOD(k,4));
    ZL = arrayfun(@(~)rand_degs(),1:M,'un',0);
    ZR = arrayfun(@(~)rand_degs(),1:N,'un',0);
    [P,vars,Zd] = rand_coeffs(m,n,ZL,ZR,0);
    msg = check_case(P,[m,n],vars,ZL,ZR,Zd,nsamp,tol);
    nchk = nchk+nsamp;
    if ~isempty(msg)
        fails{end+1} = sprintf('m=%d n=%d M=%d N=%d nZd=0: %s',m,n,M,N,msg);
    end
end

ntot = 1+ncase+size(NOD,1);
for k = 1:numel(fails),     fprintf('FAIL  %s\n',fails{k});     end
assert(isempty(fails),'test_sdvar2dpvar: %d of %d cases FAILED.',numel(fails),ntot);
fprintf('\ntest_sdvar2dpvar passed (%d cases, %d checks).\n',ntot,nchk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = check_case(P,dims,vars,ZL,ZR,Zd,nsamp,tol)
% Empty msg on success. Errors are reported, not thrown, so one crashing
% case does not hide the others.
msg = '';
try
    D = sdvar2dpvar(P,dims,vars,ZL,ZR,Zd);
catch ME
    msg = ['error: ',ME.message];   return
end
if ~isa(D,'dpvar') || ~isequal(size(D),dims)
    msg = 'output is not a dpvar of the requested size';    return
end
Dp = dpvar2poly(D);
names = [Zd(:);vars.out(:);vars.in(:)];
for k = 1:nsamp
    d = randn(numel(Zd),1);
    s = rand_pt(numel(ZL));     t = rand_pt(numel(ZR));
    ref = ref_eval(P,dims,ZL,ZR,d,s,t);
    got = eval_poly(Dp,names,[d;s;t]);
    err = norm(got-ref,'fro')/max(1,norm(ref,'fro'));
    if err>tol
        msg = sprintf('sample %d: relative error %.2e',k,err);  return
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ref = ref_eval(P,dims,ZL,ZR,d,s,t)
% The header definition, numerically: Kronecker monomial vectors at (s,t),
% first variable slowest, and a column-major unvec with the matrix component
% outer on both axes.
zL = 1;
for i = 1:numel(ZL),    zL = kron(zL,s(i).^ZL{i}(:));    end
zR = 1;
for i = 1:numel(ZR),    zR = kron(zR,t(i).^ZR{i}(:));    end
Cm = reshape(full(P.A + P.B.'*d),dims(1)*numel(zL),dims(2)*numel(zR));
ref = kron(eye(dims(1)),zL).'*Cm*kron(eye(dims(2)),zR);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = eval_poly(Pp,names,vals)
% Substitute by name; values as a COLUMN (a row is read as several points).
if isempty(Pp.varname)
    V = double(Pp);     return
end
[tf,loc] = ismember(Pp.varname(:),names(:));
assert(all(tf),'polynomial has a variable with no value');
V = double(subs(Pp,Pp.varname(:),vals(loc)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P,vars,Zd] = rand_coeffs(m,n,ZL,ZR,nZd)
% Sparse random A and B for the given bases. Output variables s1..sM,
% input variables t1..tN (disjoint, as the dpvar varname needs), decision
% variables d1..dnZd, each given at least one nonzero in B.
M = numel(ZL);      N = numel(ZR);
nc = m*n*prod([cellfun(@numel,ZL),1])*prod([cellfun(@numel,ZR),1]);
A = sprandn(nc,1,0.5);
B = sprandn(nZd,nc,0.3) + sparse(1:nZd,randi(nc,1,nZd),randn(1,nZd),nZd,nc);
P = struct('A',A,'B',B);
vars = struct('out',{arrayfun(@(i)sprintf('s%d',i),1:M,'un',0)},...
              'in', {arrayfun(@(i)sprintf('t%d',i),1:N,'un',0)});
Zd = arrayfun(@(i)sprintf('d%d',i),(1:nZd)','un',0);
end

function degs = rand_degs()
% Random nonempty subset of 0:3, as a sorted column.
degs = sort(randperm(4,randi(4))-1).';
end

function x = rand_pt(k)
% k x 1 point with |x_i| in [0.3,0.8] or [1.2,1.7], away from 0 and +-1 so
% distinct degrees give distinct values.
x = (0.3+0.5*rand(k,1)+0.9*(rand(k,1)>0.5)).*sign(randn(k,1));
end
