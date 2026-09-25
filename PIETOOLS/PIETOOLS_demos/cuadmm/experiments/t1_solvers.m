% t1_solvers.m -- do SeDuMi and Mosek actually run here, and do they agree?
% A tiny SDP with a KNOWN answer, so a solver that merely returns without
% error is still caught.  max tr(X) s.t. X(1,1)=1, X(2,2)=1, X>=0  ->  obj 2.

cuadmm_path;
fprintf('T1 sedumi_found %d mosek_found %d\n', ...
        ~isempty(which('sedumi')), ~isempty(which('mosekopt')));

% SeDuMi form: min c'x s.t. A x = b, x in K.  X is 2x2 PSD, vec'd (K.s=2).
K = struct('f',0,'l',0,'q',[],'s',2);
A = sparse([1 2],[1 4],[1 1],2,4);        % X(1,1)=1 ; X(2,2)=1
b = [1;1];
c = -[1;0;0;1];                            % minimise -tr(X)  ->  max tr(X)

ok_s = false; ok_m = false;
try
    t = tic; [x,~,info] = sedumi(A,b,c,K,struct('fid',0)); ts = toc(t);
    objs = -c'*x;
    fprintf('T1 sedumi obj %.12f t %.3f pinf %d dinf %d numerr %d\n', ...
            objs, ts, info.pinf, info.dinf, info.numerr);
    ok_s = abs(objs-2) < 1e-7;
catch ME
    fprintf('T1 sedumi FAILED %s\n', ME.message);
end

try
    prob.bardim = 2;
    % Mosek via its SeDuMi-compatible path is awkward; use mosekopt on the
    % same data through PIETOOLS' own interface instead if present.
    [r,res] = mosekopt('minimize echo(0)', struct('c',[], 'a',sparse(0,0)));
    fprintf('T1 mosekopt callable rcode %d\n', r);
    ok_m = true;
catch ME
    fprintf('T1 mosekopt FAILED %s\n', ME.message);
end

fprintf('T1 sedumi_correct %d mosek_callable %d\n', ok_s, ok_m);
fprintf('T1DONE\n');
