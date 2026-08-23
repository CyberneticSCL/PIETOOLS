%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_POSSOPVAR_INTEGRATION checks that the operator returned by
% 'possopvar' works with the routines already implemented for the sdopvar
% class, rather than only with the code that builds it.
%
% Every check is made at the level of the operator: the pairing <y,Px> is
% evaluated in closed form from the stored kernels (PI_SDOPVAR_KERNELS and
% PI_KERNEL_PAIRING), so a routine passes only if it transforms the operator
% as intended, independently of how the coefficients happen to be stored.
%
% Routines exercised:
%   plus                (P+Q, and accumulation against a scalar multiple)
%   mtimes  (scalar)    (c*P and P*c)
%   mtimes  (matrix)    (M*P and P*M, checked via <y,(M*P)x> = <M'y,Px>)
%   ctranspose          (P', and P'' = P)
%   transpose           (P.', which forwards to ctranspose)
%   eq                  (P==P, P~=2*P, and comparison against zero)
%   sdopvar2ndopvar     (conversion, and round trip via ndopvar2sdopvar)
%
% The script also reports which routines a caller might reasonably expect
% but which the class does not currently define.
%
% MP, 08/22/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(4);
tol = 1e-8;

cases = {
    {'1D, m=1, deg 1', 1, {'s1'},      [0,1],     1}
    {'1D, m=2, deg 1', 2, {'s1'},      [0,1],     1}
    {'2D, m=1, deg 1', 1, {'s1','s2'}, [0,1;0,1], 1}
    };

nfail = 0;
for ic = 1:numel(cases)
    [lbl,m,vars,dom] = deal(cases{ic}{1:4});
    deg = cases{ic}{5};
    n3 = numel(vars);
    if size(dom,1)==1 && n3~=1
        dom = repmat(dom,n3,1);
    end
    fprintf('\n=== %s ===\n',lbl);

    prog = lpiprogram(polynomial(vars(:)),[],dom);
    [prog,P] = possopvar(prog,m,vars,dom,deg);

    nd = numel(P.Zd);
    dv = randn(nd,1);
    xf = rand_polyvec(vars,m);
    yf = rand_polyvec(vars,m);

    base = pair(P,xf,yf,dv,vars,dom);

    % % % plus
    nfail = nfail + report('plus: <y,(P+P)x> = 2<y,Px>',...
        relerr(pair(P+P,xf,yf,dv,vars,dom),2*base),tol);

    % % % mtimes with a scalar, on either side
    nfail = nfail + report('mtimes: <y,(3*P)x> = 3<y,Px>',...
        relerr(pair(3*P,xf,yf,dv,vars,dom),3*base),tol);
    nfail = nfail + report('mtimes: <y,(P*3)x> = 3<y,Px>',...
        relerr(pair(P*3,xf,yf,dv,vars,dom),3*base),tol);

    % % % plus combined with mtimes
    nfail = nfail + report('plus/mtimes: (P+P) == 2*P',...
        double(~eq(P+P,2*P,tol)),0.5);

    % % % mtimes with a matrix, checked against moving it to the test function
    M = randn(m,m);
    nfail = nfail + report('mtimes: <y,(M*P)x> = <M''y,Px>',...
        relerr(pair(M*P,xf,yf,dv,vars,dom),pair(P,xf,M.'*yf,dv,vars,dom)),tol);
    nfail = nfail + report('mtimes: <y,(P*M)x> = <y,P(Mx)>',...
        relerr(pair(P*M,xf,yf,dv,vars,dom),pair(P,M*xf,yf,dv,vars,dom)),tol);

    % % % ctranspose and transpose
    nfail = nfail + report('ctranspose: <y,P''x> = <x,Py>',...
        relerr(pair(P',xf,yf,dv,vars,dom),pair(P,yf,xf,dv,vars,dom)),tol);
    nfail = nfail + report('ctranspose: <y,P''''x> = <y,Px>',...
        relerr(pair((P')',xf,yf,dv,vars,dom),base),tol);
    nfail = nfail + report('transpose: <y,P.''x> = <y,P''x>',...
        relerr(pair(P.',xf,yf,dv,vars,dom),pair(P',xf,yf,dv,vars,dom)),tol);

    % % % eq
    nfail = nfail + report('eq: P == P',double(~eq(P,P,tol)),0.5);
    nfail = nfail + report('eq: P ~= 2*P',double(eq(P,2*P,tol)),0.5);
    nfail = nfail + report('eq: (P-scaled) 0*P == 0',double(~eq(0*P,0,tol)),0.5);

    % % % Known limitation, recorded here so that a future fix shows up as a
    % % % newly passing check rather than as a silent change.
    if eq(P,P',tol)
        fprintf('  NOTE   eq(P,P'') now succeeds; the multiplier-cell normalization issue appears fixed.\n');
    else
        fprintf('  KNOWN  eq(P,P'') returns false even though <y,P''x> = <x,Py> holds,\n');
        fprintf('         because ctranspose does not renormalize multiplier cells.\n');
    end

    % % % sdopvar2ndopvar, and the round trip back.
    % The round trip is reported but does not count as a failure of
    % 'possopvar': see the note printed at the end of this script.
    try
        Pn = sdopvar2ndopvar(P);
        fprintf('  ok     sdopvar2ndopvar returned a %s\n',class(Pn));
        Pback = ndopvar2sdopvar(Pn);
        errA = 0;   errB = 0;
        [~,loc] = ismember(cellstr(string(Pback.Zd(:))),cellstr(string(P.Zd(:))));
        for k = 1:numel(P.params.A)
            errA = max(errA,max(abs(full(P.params.A{k}-Pback.params.A{k}))));
            B1 = full(P.params.B{k});
            B2 = full(Pback.params.B{k});   B2 = B2(loc,:);
            errB = max(errB,max(abs(B1(:)-B2(:))));
        end
        nfail = nfail + report('round trip preserves the constant term A',errA,tol);
        errop = relerr(pair(Pback,xf,yf,dv_for(Pback,P,dv),vars,dom),base);
        if errop<tol && errB<tol
            fprintf('  ok     round trip preserves <y,Px> and B          (%.1e)\n',errop);
        else
            fprintf('  KNOWN  round trip does NOT preserve B (%.1e) or <y,Px> (%.1e)\n',errB,errop);
        end
    catch ME
        fprintf('  FAIL   sdopvar/ndopvar conversion: %s\n',ME.message);
        nfail = nfail+1;
    end
end

% % % Report routines that a caller might expect but that are not defined.
% Note that a name resolving without error does not mean it does the right
% thing: with no overload, '[P,P]' silently builds a 1x2 object array rather
% than concatenating the operators, and 'size(P)' reports the array size
% rather than the operator dimensions. A non-square operator is used here so
% that the two cannot be confused.
fprintf('\n=== routines not defined for sdopvar ===\n');
prog = lpiprogram(polynomial({'s1'}),[],[0,1]);
[prog,P] = possopvar(prog,2,{'s1'},[0,1],1);
probes = { {'minus (P-Q)',      @() P-P}
           {'uminus (-P)',      @() -P}
           {'size(P)',          @() assert_eq(size(P),P.dims,'size(P)')}
           {'horzcat ([P,P])',  @() assert_eq(numel([P,P]),1,'numel([P,P])')}
           {'vertcat ([P;P])',  @() assert_eq(numel([P;P]),1,'numel([P;P])')}
           {'blkdiag(P,P)',     @() blkdiag(P,P)}
           {'P*P (composition)',@() P*P} };
for k = 1:numel(probes)
    try
        probes{k}{2}();
        fprintf('  available : %s\n',probes{k}{1});
    catch ME
        fprintf('  MISSING   : %-20s (%s)\n',probes{k}{1},strtrim(ME.message));
    end
end

fprintf('\n=== known pre-existing issues, not caused by possopvar ===\n');
fprintf('  1. eq(P,P'') is false for a self-adjoint P with a nonzero multiplier\n');
fprintf('     part, because ctranspose does not restore the canonical form in\n');
fprintf('     which a multiplier cell carries no ZR degree.\n');
fprintf('  2. sdopvar2ndopvar/ndopvar2sdopvar do not round trip the decision\n');
fprintf('     coefficient B whenever there is more than one decision variable.\n');
fprintf('     This reproduces on a hand-built sdopvar in the documented layout\n');
fprintf('     vec(C) = A + B''*d, so it is independent of possopvar. The constant\n');
fprintf('     term A does round trip. @sdopvar/test_script does not catch it\n');
fprintf('     because it only goes ndopvar -> sdopvar -> ndopvar, applying the\n');
fprintf('     same convention in both directions.\n');

fprintf('\n');
if nfail>0
    error('test_possopvar_integration: %d check(s) failed.',nfail);
end
fprintf('possopvar integration test passed.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = pair(P,xf,yf,dv,vars,dom)
% Closed-form evaluation of <y,Px> from the kernels stored in P.

val = pi_kernel_pairing(pi_sdopvar_kernels(P,dv),xf,yf,vars,dom);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dv_new = dv_for(Pnew,Pold,dv)
% Re-express a decision variable assignment in the ordering used by Pnew.

old = cellstr(string(Pold.Zd(:)));
new = cellstr(string(Pnew.Zd(:)));
[tf,loc] = ismember(new,old);
if ~all(tf)
    error('Converted operator introduced unknown decision variables.');
end
dv_new = dv(loc);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = relerr(a,b)

e = abs(a-b)/max(1,abs(b));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nf = report(lbl,err,tol)

if err<tol
    fprintf('  ok     %-46s (%.1e)\n',lbl,err);
    nf = 0;
else
    fprintf('  FAIL   %-46s (%.1e)\n',lbl,err);
    nf = 1;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function assert_eq(got,want,lbl)

if ~isequal(reshape(got,1,[]),reshape(want,1,[]))
    error('%s = %s, expected %s',lbl,mat2str(got),mat2str(want));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = rand_polyvec(vars,m)
% Random m x 1 polynomial vector of degree at most 1 in each variable.

n3 = numel(vars);
if n3==0
    w = polynomial(randn(m,1));
    return
end
grids = cell(1,n3);
for i = 1:n3
    grids{i} = (0:1)';
end
subsc = cell(1,n3);
[subsc{:}] = ndgrid(grids{:});
E = zeros(numel(subsc{1}),n3);
for i = 1:n3
    E(:,i) = subsc{i}(:);
end
T = size(E,1);
w = polynomial(randn(T,m),E,vars(:),[m,1]);

end
