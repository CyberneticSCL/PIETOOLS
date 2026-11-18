% INH_FD  Reference thresholds by finite differences (N = 40, 60) and a
% conversion smoke test of every system. SCRATCH, 11/16/2026.
addpath('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces');
S = pf_setup(false);
D = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/inhom';
addpath(D);
z = @(a,b) 0*a;     o = @(a,b) 1+0*a;
Def = { 'A', @(a,b) 1+a.^2, @(a,b) 1+b,    z,          z, z;
        'B', o,             o,             z,          z, @(a,b) -3*a.*b;
        'C', @(a,b) 1+a,    o,             @(a,b) b,   z, z;
        'M2', o,            @(a,b) 1+b.^2, z,          z, z};
T = struct();
for i = 1:size(Def,1)
    for N = [40 60]
        mu = inh_fdmax(Def{i,2:6},N);
        T.(Def{i,1})(N==[40 60]) = -mu;
        fprintf('@@ FD %s N=%d lamstar=%.6f\n',Def{i,1},N,-mu);
    end
end
for N = [40 60 200]
    mu = inh_fdmax(@(s) 1+s,[],@(s) 0*s,[],@(s) 0*s,N);
    fprintf('@@ FD M1 N=%d lamstar=%.6f\n',N,-mu);
    T.M1(find(N==[40 60 200])) = -mu;
end
% Gahlawat (library 7) at lam = 0: x(0)=0, x_s(1)=0 -> ghost node, checked against 4.66
for N = [200 400]
    h = 1/N;    s = (1:N)'*h;   e = ones(N,1);
    D2 = spdiags([e -2*e e],-1:1,N,N);  D2(N,N-1) = 2;  D2 = D2/h^2;
    D1 = spdiags([-e 0*e e],-1:1,N,N);  D1(N,N-1) = 0;  D1 = D1/(2*h);
    a = s.^3-s.^2+2;  b = 3*s.^2-2*s;  c = -0.5*s.^3+1.3*s.^2-1.5*s+0.7;
    L = spdiags(a,0,N,N)*D2 + spdiags(b,0,N,N)*D1 + spdiags(c,0,N,N);
    mu = max(real(eigs(L,4,'largestreal')));
    fprintf('@@ FD G7 N=%d lamstar=%.6f\n',N,-mu);   T.G7(find(N==[200 400])) = -mu;
end
save([D '/results/thresholds.mat'],'T');
for nm = {'A','B','C'}
    t0 = tic;  PIE = inh_sys(nm{1},0);
    fprintf('@@ PIE %s dims T=%s A=%s t=%.1f\n',nm{1},mat2str(size(PIE.T)),mat2str(size(PIE.A)),toc(t0));
end
t0 = tic;  PIE = inh_sys('M',[1 1]);
fprintf('@@ PIE M dims T=%s t=%.1f\n',mat2str(size(PIE.T)),toc(t0));
disp(class(PIE));
disp('## FD DONE');
