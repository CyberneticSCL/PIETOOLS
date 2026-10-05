%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_DPVAR2SDVAR checks 'dpvar2sdvar' (sopvar/Testfolder/converters,
% despite its own header's "moved to sopvar/misc" note) against the
% definition in its header,
%
%   D = (Im o ZL{1}(s1) o ... o ZL{M}(sM))^T unvec(P.A + P.B'*d)
%           (In o ZR{1}(t1) o ... o ZR{N}(tN)),
%
% i.e. the inverse direction of sdvar2dpvar (sopvar/Testfolder/converters).
%
% Per case, against the semantics and never a round trip (CLAUDE.md S4):
% the dpvar object D is built directly from the 'dpvar' constructor (never
% from 'sdvar2dpvar'), with a random NATIVE monomial set that need not be a
% full ZL x ZR grid, exercising the basis-building in 'dpvar2sdvar'. The
% ground truth is D's own 'dpvar2poly' expansion, evaluated by 'subs' at
% random (d,s,t); the value under test reconstructs the same point from
% dpvar2sdvar's (P,ZL,ZR) via the header formula above. Neither
% 'dpvar2sdvar' nor 'sdvar2dpvar' enters the ground truth.
%
% Sweep: m,n in {1,2,3} x 0/1/2 output x 0/1/2 input variables x 0, 1 or 3
% decision variables (243 cases); native monomial sets are random (possibly
% sparse/non-grid) subsets of degrees 0:3, every decision variable forced
% to have at least one nonzero. Plus: the nargin==1 branch (all variables
% treated as output); the 'variable in neither vars.in nor vars.out' error
% path; and that P.dvarname is a column even when D.dvarname is a row (see
% dpvar2sdvar.m header, 'sdopvar' expects a column).
%
% Initial coding DJ, 10/05/2026
%

tol = 1e-9;        % relative, Frobenius
nsamp = 4;         % random (d,s,t) per case
nchk = 0;          fails = {};

%% NARGIN1: all variables of D are implicitly output variables
degmat = [1 0; 0 1; 2 1];
C = sparse([1 2 0; 0.5 0 3]);           % rows: [const; d1], m=1, nd=1
D = dpvar(C,degmat,{'x';'y'},{'d1'},[1,1]);
msg = '';
try
    [P,ZL,ZR] = dpvar2sdvar(D);
    if ~isempty(ZR)
        msg = 'ZR should be empty with nargin==1';
    else
        Dp = dpvar2poly(D);
        pt = [1.3,-0.7];    dval = 2.1;
        ref = double(subs(Dp,{'d1';'x';'y'},[dval;pt(1);pt(2)]));
        zL = 1;     for i=1:numel(ZL), zL = kron(zL,pt(i).^ZL{i}(:)); end
        got = full(P.A + P.B'*dval)'*zL;
        if abs(got-ref)>tol
            msg = sprintf('value mismatch: got %.6g, ref %.6g',got,ref);
        end
    end
catch ME
    msg = ['error: ',ME.message];
end
nchk = nchk+1;
if ~isempty(msg),   fails{end+1} = ['NARGIN1: ',msg];  end

%% ERRCHK: a variable in neither vars.in nor vars.out must error
msg = '';
try
    dpvar2sdvar(D,struct('out',{{'x'}},'in',{cell(1,0)}));
    msg = 'expected an error, got none';
catch ME
    if ~contains(ME.message,'neither')
        msg = ['wrong error: ',ME.message];
    end
end
nchk = nchk+1;
if ~isempty(msg),   fails{end+1} = ['ERRCHK: ',msg];  end

%% ORIENT: P.dvarname comes back a column even if D.dvarname is a row         % DJ, 10/05/2026
%  ('dpvar' leaves dvarname's orientation unconstrained; 'sdopvar' and its    % DJ, 10/05/2026
%  machinery -- sync_basis, verify_copvar_meta -- expect a column, see        % DJ, 10/05/2026
%  dpvar2sdvar.m header). Needs >=2 decision variables: a 1x1 cell is both    % DJ, 10/05/2026
%  a row and a column, so it can't tell the orientations apart.               % DJ, 10/05/2026
msg = '';                                                                     % DJ, 10/05/2026
degmat2 = [1 0; 0 1];                                                         % DJ, 10/05/2026
C2 = sparse([1 2; 0.5 0; 0 0.3]);          % rows: [const;d1;d2], m=1, nd=2   % DJ, 10/05/2026
Drow = dpvar(C2,degmat2,{'x';'y'},{'d1','d2'},[1,1]);    % dvarname is a row  % DJ, 10/05/2026
try                                                                           % DJ, 10/05/2026
    Pr = dpvar2sdvar(Drow);                                                  % DJ, 10/05/2026
    if ~iscolumn(Pr.dvarname)                                                % DJ, 10/05/2026
        msg = sprintf('P.dvarname is %dx%d, expected a column',...           % DJ, 10/05/2026
                       size(Pr.dvarname,1),size(Pr.dvarname,2));             % DJ, 10/05/2026
    end                                                                      % DJ, 10/05/2026
catch ME                                                                      % DJ, 10/05/2026
    msg = ['error: ',ME.message];                                            % DJ, 10/05/2026
end                                                                           % DJ, 10/05/2026
nchk = nchk+1;                                                                % DJ, 10/05/2026
if ~isempty(msg),   fails{end+1} = ['ORIENT: ',msg];  end                    % DJ, 10/05/2026

%% VALUE: sweep over dimensions, variable counts and decision variables
ncase = 0;
for m = 1:3
for n = 1:3
for M = 0:2
for N = 0:2
for nd = [0,1,3]
    ncase = ncase+1;
    rng(4000+ncase);
    msg = check_case(m,n,M,N,nd,nsamp,tol);
    nchk = nchk+nsamp;
    if ~isempty(msg)
        fails{end+1} = sprintf('m=%d n=%d M=%d N=%d nd=%d: %s',m,n,M,N,nd,msg);
    end
end
end
end
end
end

ntot = 3+ncase;                                                              % DJ, 10/05/2026 (was 2; added ORIENT)
for k = 1:numel(fails),     fprintf('FAIL  %s\n',fails{k});     end
assert(isempty(fails),'test_dpvar2sdvar: %d of %d cases FAILED.',numel(fails),ntot);
fprintf('\ntest_dpvar2sdvar passed (%d cases, %d checks).\n',ntot,nchk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msg = check_case(m,n,M,N,nd,nsamp,tol)
% Empty msg on success. Errors are reported, not thrown, so one crashing
% case does not hide the others.
msg = '';
varnames_out = arrayfun(@(i)sprintf('s%d',i),1:M,'un',0);
varnames_in  = arrayfun(@(i)sprintf('t%d',i),1:N,'un',0);
varname = [varnames_out(:); varnames_in(:)];
dvarname = arrayfun(@(i)sprintf('d%d',i),(1:nd)','un',0);

% Random native monomial set: NOT necessarily a full ZL x ZR grid, so the
% basis-building (ZL_maps/ZR_maps) in dpvar2sdvar is actually exercised.
nZ = randi([1,5]);
degmat = randi([0,3],nZ,M+N);
degmat = unique(degmat,'rows');
nZ = size(degmat,1);
if nZ==0,   degmat = zeros(0,M+N);     end

C = sprandn((nd+1)*m, nZ*n, 0.6);
for k = 1:nd    % every decision variable gets at least one nonzero
    rr = (randi(m)-1)*(nd+1) + 1 + k;
    cc = randi(max(nZ*n,1));
    C(rr,cc) = randn();
end

try
    D = dpvar(C,degmat,varname,dvarname,[m,n]);
catch ME
    msg = ['dpvar construction error: ',ME.message];    return
end

vars = struct('out',{varnames_out},'in',{varnames_in});
try
    [P,ZL,ZR] = dpvar2sdvar(D,vars);
catch ME
    msg = ['error: ',ME.message];   return
end

Dp = dpvar2poly(D);             % ground truth: D's own semantics
names = [dvarname(:); varnames_out(:); varnames_in(:)];
for k = 1:nsamp
    d = randn(nd,1);
    s = rand_pt(M);     t = rand_pt(N);
    ref = eval_poly(Dp,names,[d;s;t]);
    got = kernel_eval(P,ZL,ZR,[m,n],d,s,t);
    err = norm(got-ref,'fro')/max(1,norm(ref,'fro'));
    if err>tol
        msg = sprintf('sample %d: relative error %.2e',k,err);    return
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = kernel_eval(P,ZL,ZR,dims,d,s,t)
% The dpvar2sdvar header definition, numerically: Kronecker monomial
% vectors at (s,t), first variable slowest, column-major unvec.
zL = 1;     for i=1:numel(ZL),  zL = kron(zL,s(i).^ZL{i}(:));  end
zR = 1;     for i=1:numel(ZR),  zR = kron(zR,t(i).^ZR{i}(:));  end
m = dims(1);    n = dims(2);
Cm = reshape(full(P.A + P.B'*d), m*numel(zL), n*numel(zR));
V = kron(eye(m),zL)'*Cm*kron(eye(n),zR);
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

function x = rand_pt(k)
% k x 1 point with |x_i| in [0.3,0.8] or [1.2,1.7], away from 0 and +-1 so
% distinct degrees give distinct values.
x = (0.3+0.5*rand(k,1)+0.9*(rand(k,1)>0.5)).*sign(randn(k,1));
end
