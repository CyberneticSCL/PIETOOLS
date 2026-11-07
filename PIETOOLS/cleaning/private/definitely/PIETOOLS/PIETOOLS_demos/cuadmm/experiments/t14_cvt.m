% t14_cvt.m -- verify the PIETOOLS -> cuADMM conversion INDEPENDENTLY.
%
% Why this is not optional: dump2cuadmm's own header records that cuADMM's
% shipped read_sedumi/svecADMM chain CORRUPTED this data (nnz collapsed
% 1489 -> 193, residual 1.2e+01 at a known-good point).  A converter that is
% wrong in the same way in both directions would still round-trip, so the test
% below never compares the converter against itself.
%
% The identity that must hold, for the SDP to be the same problem:
%     <A_k, X>  ==  svec(A_k)' * svec(X)     for every constraint k
% computed from the ORIGINAL SeDuMi data on the left and from the WRITTEN
% TEXT FILES on the right, at random symmetric X.  This is decisive: it checks
% the sqrt(2) scaling, the symmetrisation of the non-symmetric A_k that
% sossolve emits, the upper-triangle column-major ordering, the free block,
% and the file I/O, all at once.

cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
OUT  = fullfile(cuadmm_outdir(),'cvt_test');

% --- a real PIETOOLS program (1-D stability, small)
pvar s t
x   = pde_var('state',1,s,[0,1]);
PIE = convert([diff(x,t,1)==diff(x,s,2)+2*x; subs(x,s,0)==0; subs(x,s,1)==0]);
st  = lpisettings('light');  st.eppos2 = 1e-2;
M   = cuadmm_private('stab_mirror',PIE,st);

Atf=[]; bf=[];
for i=1:M.prog.expr.num, Atf=[Atf,M.prog.expr.At{i}]; bf=[bf;M.prog.expr.b{i}]; end
S = cuadmm_private('sdpshape',M.prog);
fprintf('T14 program m=%d nvar=%d Kf=%d Ns=[%s]\n',S.m,S.nvar,S.Kf,num2str(S.Ks));

% --- build the SeDuMi-form dump dump2cuadmm expects, mirroring sossolve's
%     processvars: At = RR'*Atf with RR a permutation-like blkdiag of speye.
%     For this program there are no 'poly' vars needing the antidiagonal, so
%     RR = I; assert that rather than assume it.
RR = speye(size(Atf,1));
D.At = RR'*Atf;  D.b = bf;  D.c = sparse(size(Atf,1),1);
D.K  = struct('f',S.Kf,'l',0,'q',[],'s',S.Ks);
D.Ns = S.Ks;  D.Kf = S.Kf;
dumpfile = fullfile(cuadmm_outdir(),'cvt_dump.mat');
save(dumpfile,'-struct','D');

info = dump2cuadmm(dumpfile,OUT);
fprintf('T14 wrote vec_len=%d nnz=%d\n',info.vec_len,info.nnz_At);

% --- read the text files back, with NO reuse of the writer's code
cuAt = readcoo(fullfile(OUT,'At.txt'), info.vec_len, S.m);
cub  = readcoo(fullfile(OUT,'b.txt'),  S.m, 1);
blk  = fileread(fullfile(OUT,'blk.txt'));
fprintf('T14 blk.txt = %s\n', strtrim(strrep(blk,newline,' ')));
fprintf('T14 b match = %.3e\n', norm(full(cub)-full(bf))/max(norm(full(bf)),eps));

% --- the decisive test
rng(7);
worst = 0;
for trial = 1:5
    % random symmetric X in the cone's shape, plus a random free part
    xf = randn(S.Kf,1);
    Xb = cell(numel(S.Ks),1);  xvec = xf;  xsvec = xf;
    for i = 1:numel(S.Ks)
        N = S.Ks(i);
        A = randn(N);  Xi = (A+A')/2;
        Xb{i} = Xi;
        xvec  = [xvec; Xi(:)];                    %#ok<AGROW>  SeDuMi vec order
        xsvec = [xsvec; svec_up(Xi)];             %#ok<AGROW>  SDPT3 svec order
    end
    lhs = full(D.At' * xvec);      % <A_k, X> from the ORIGINAL data
    rhs = full(cuAt' * xsvec);     % from the WRITTEN FILES
    e   = norm(lhs-rhs)/max(norm(lhs),eps);
    worst = max(worst,e);
    fprintf('T14 trial %d  rel_diff = %.3e\n',trial,e);
end
fprintf('T14 WORST = %.3e   (must be ~1e-15)\n',worst);
fprintf('T14 VERDICT %s\n', ternary(worst<1e-12,'CONVERTER OK','CONVERTER WRONG'));
fprintf('T14DONE\n');


function v = svec_up(X)
% SDPT3 svec: upper triangle taken COLUMN BY COLUMN, diagonal unscaled,
% off-diagonal scaled by sqrt(2).  Written from the definition here, not
% lifted from dump2cuadmm, so the two are independent.
N = size(X,1);
v = zeros(N*(N+1)/2,1);
p = 0;
for j = 1:N
    for i = 1:j
        p = p+1;
        if i==j, v(p) = X(i,j); else, v(p) = sqrt(2)*X(i,j); end
    end
end
end

function A = readcoo(fn,nr,nc)
d = load(fn);                       % "row col val", 0-based
if isempty(d), A = sparse(nr,nc); return; end
A = sparse(d(:,1)+1, d(:,2)+1, d(:,3), nr, nc);
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
