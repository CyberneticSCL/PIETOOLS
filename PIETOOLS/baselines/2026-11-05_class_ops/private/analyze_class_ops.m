function analyze_class_ops(f,outdir)
% ANALYZE_CLASS_OPS(F,OUTDIR) summarizes the results of BENCH_CLASS_OPS saved
% in the MAT-file F: writes every row to OUTDIR/results.csv and the summary
% tables (markdown) to OUTDIR/tables.md.
% INPUT
% - f:      char, MAT-file with the struct arrays Rc ('const'), Rb ('dvarB'),
%           Rq ('dvarQ'), Rn ('dvarN'), Rs ('dvarS'), Ri ('inv');
% - outdir: char, output folder.
% OUTPUT
% none; the two files. Exponents are least-squares slopes of log t: over
% the (n,d) grid, log t = c + alpha log n + beta log d; over q, the slope
% between the two largest q.

% Initial coding MMP, 11/05/2026: summary of BENCH_CLASS_OPS (README.md).

S = load(f);
T = rows2table(S);
writetable(T,fullfile(outdir,'results.csv'));
fid = fopen(fullfile(outdir,'tables.md'),'w');
c = onCleanup(@() fclose(fid));

% Constant operators and the free-variable dopvar: fits over (n,d)
grid_table(fid,T,'const','opvar','sopvar');
grid_table(fid,T,'dvarB','dopvar','sdopvar');
% q sweeps
q_table(fid,T,'dvarQ','dopvar','sdopvar');
q_table(fid,T,'dvarN','dopvar','sdopvar');
q_table(fid,T,'dvarS','dopvar','sdopvar');
% inverse
Ti = T(strcmp(T.part,'inv'),:);
fprintf(fid,'\n### inv\n\n| n | d | opvar t (s) | sopvar t (s) | ratio | opvar resid | sopvar resid |\n|---|---|---|---|---|---|---|\n');
for k = 1:height(Ti)/2
    o = Ti(2*k-1,:);    s = Ti(2*k,:);
    fprintf(fid,'| %d | %d | %.3g | %.3g | %.2f | %.1e | %.1e |\n',o.n,o.d,o.t,s.t,s.t/o.t,o.agree,s.agree);
end
% storage of the operator itself (the uminus result has the input's size)
fprintf(fid,'\n### storage (bytes of the uminus result = size of the operator)\n\n| part | n | d | q | nnz | legacy | new | new/legacy |\n|---|---|---|---|---|---|---|---|\n');
for p = {'const','dvarB','dvarQ','dvarN','dvarS'}
    Tp = T(strcmp(T.part,p{1}) & strcmp(T.op,'uminus'),:);
    keys = unique(Tp(:,{'n','d','q','nnz'}),'rows','stable');
    for k = 1:height(keys)
        sel = Tp.n==keys.n(k) & Tp.d==keys.d(k) & Tp.q==keys.q(k);
        b = Tp.bytes(sel);
        fprintf(fid,'| %s | %d | %d | %d | %d | %d | %d | %.2f |\n',p{1},keys.n(k),keys.d(k),keys.q(k),keys.nnz(k),b(1),b(2),b(2)/b(1));
    end
end
fprintf('##tables written to %s\n',fullfile(outdir,'tables.md'));
end

function T = rows2table(S)
% All result rows as one table
C = {};
for fn = fieldnames(S)'
    R = S.(fn{1});
    for i = 1:numel(R)
        r = R(i);
        C(end+1,:) = {r.part,r.n,r.d,r.q,r.nnz,r.op,r.cls,r.t,r.nrep,r.bytes,r.err,r.agree}; %#ok<AGROW>
    end
end
T = cell2table(C,'VariableNames',{'part','n','d','q','nnz','op','cls','t','nrep','bytes','err','agree'});
end

function grid_table(fid,T,part,old,new)
% Per operation: time at the largest (n,d), ratio new/old there and its
% median/min/max over the grid, and the exponents alpha (n) and beta (d)
Tp = T(strcmp(T.part,part),:);
if isempty(Tp),     return,     end
ops = unique(Tp.op,'stable');
nmax = max(Tp.n);   dmax = max(Tp.d);
fprintf(fid,'\n### %s: %s vs %s over n x d (largest n = %d, d = %d)\n\n',part,old,new,nmax,dmax);
fprintf(fid,'| operation | %s t (s) | %s t (s) | ratio | ratio median [min, max] | %s alpha, beta | %s alpha, beta |\n|---|---|---|---|---|---|---|\n',old,new,old,new);
for k = 1:numel(ops)
    To = Tp(strcmp(Tp.op,ops{k}) & strcmp(Tp.cls,old),:);
    Tn = Tp(strcmp(Tp.op,ops{k}) & strcmp(Tp.cls,new),:);
    [~,io] = sortrows([To.n To.d]);     [~,in] = sortrows([Tn.n Tn.d]);
    To = To(io,:);  Tn = Tn(in,:);
    rat = Tn.t./To.t;
    big = To.n==nmax & To.d==dmax;
    [ao,bo] = fitnd(To);   [an,bn] = fitnd(Tn);
    fprintf(fid,'| %s | %.3g | %.3g | %.2f | %.2f [%.2f, %.2f] | %.2f, %.2f | %.2f, %.2f |\n', ...
        ops{k},To.t(big),Tn.t(big),rat(big),median(rat),min(rat),max(rat),ao,bo,an,bn);
end
end

function [a,b] = fitnd(Tc)
% least squares log t = c + a log n + b log d
X = [ones(height(Tc),1), log(Tc.n), log(Tc.d)];
w = X\log(Tc.t);    a = w(2);   b = w(3);
end

function q_table(fid,T,part,old,new)
% Per operation: time at the largest q for both classes, the ratio, and the
% slope of log t in log q between the two largest q
Tp = T(strcmp(T.part,part),:);
if isempty(Tp),     return,     end
ops = unique(Tp.op,'stable');
qs = unique(Tp.q);  q1 = qs(end-1);     q2 = qs(end);
nz = unique(Tp.nnz(Tp.q==q2));
fprintf(fid,'\n### %s: %s vs %s over q (n = %d, d = %d; nnz at the largest q = %d)\n\n',part,old,new,Tp.n(1),Tp.d(1),nz(1));
fprintf(fid,'| operation | %s t at q=%d (s) | %s t (s) | ratio | %s q-exponent | %s q-exponent |\n|---|---|---|---|---|---|\n',old,q2,new,old,new);
for k = 1:numel(ops)
    sel = strcmp(Tp.op,ops{k});
    to = Tp.t(sel & strcmp(Tp.cls,old) & Tp.q==q2);    tn = Tp.t(sel & strcmp(Tp.cls,new) & Tp.q==q2);
    to1 = Tp.t(sel & strcmp(Tp.cls,old) & Tp.q==q1);   tn1 = Tp.t(sel & strcmp(Tp.cls,new) & Tp.q==q1);
    fprintf(fid,'| %s | %.3g | %.3g | %.2f | %.2f | %.2f |\n',ops{k},to,tn,tn/to, ...
        log(to/to1)/log(q2/q1),log(tn/tn1)/log(q2/q1));
end
end
