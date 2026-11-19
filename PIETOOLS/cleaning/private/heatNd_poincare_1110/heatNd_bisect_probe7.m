function heatNd_bisect_probe7(P, outfile)
% HEATND_BISECT_PROBE7(P,OUTFILE) rebuilds the family of test_heatNd_poincare section (7)
% (1-D DD, L2, sel full, d = 0, gen 'faces', target 'grad'), as
% build_family of heatNd_poincare.m does, in the export P, with the row
% record of lpi_eq_cdopvar on (prog.sopeq = {}). Prints, per lam, m, the
% equations with no decision variable (zero columns of At) and their b, and
% the record of dropped multiples; saves everything to OUTFILE.
P = strrep(P,'/',filesep);
restoredefaultpath;
old = cd(P);  evalc('pietools_path_update');  cd(tempdir);
warning('off','MATLAB:mpath:privateDirectoriesNotAllowedOnPath');
TP = fullfile(P,'testing_codes','private');
if exist(TP,'dir'),  addpath(genpath(TP));  end
heatNd_path();
assert(startsWith(which('lpi_eq_cdopvar'),P) && startsWith(which('heatNd_poincare'),P));
fprintf('== %s\n', P);
O = heatNd_poincare_ops(heatNd_pie(1,0,{'DD'}));
N = O.N;    vars = O.vars;  dom = O.dom;    I = 1;  mu = sum(O.mu(I));  d = 0;
GG = O.DT{1}'*O.DT{1};      TAs = -1*GG;
Q0 = -1*TAs;    Q1 = -1*(O.T'*O.T);
Iinc = dec2bin(0:2^N-1,N)-'0'+2;
prog = lpiprogram(polynomial(vars(:)),[],dom);
po = struct('include',{{Iinc}});
[prog,G] = poscopvar(prog,1,vars,dom,d,po);
for i = 1:N
    for e = 1:2
        pg = po;    pg.psatz = 2*i+e;
        [prog,Gg] = poscopvar(prog,1,vars,dom,d,pg);
        G = G + Gg;
    end
end
lams = [0.5 0.75 0.9]*mu;
S = struct('lam',{},'At',{},'b',{},'m',{},'rec',{});
for j = 1:3
    p = prog;   p.sopeq = {};
    p = lpi_eq_cdopvar(p,G - (Q0 + lams(j)*Q1),'symmetric');
    D = heatNd_sdp(p);
    S(j) = struct('lam',lams(j),'At',D.At,'b',full(D.b*D.bscl),'m',D.m,'rec',p.sopeq{end});
end
same = isequal(S(1).At,S(2).At) && isequal(S(1).At,S(3).At);
db = (S(2).b-S(1).b)/(lams(2)-lams(1));
berr = norm(S(3).b - (S(1).b + (lams(3)-lams(1))*db))/max(norm(S(3).b),eps);
fprintf('m = %d %d %d, At identical %d, b third-point error %.2e\n', S(1).m, S(2).m, S(3).m, same, berr);
for j = 1:3
    zr = find(full(max(abs(S(j).At),[],1))'<=1e-12);
    E = S(j).rec.entries;
    cells = [];  cols = [];
    for t = 1:numel(E)
        for c = 1:numel(E(t).cells)
            cols = [cols, E(t).cols{c}];  cells = [cells, E(t).cells(c)*ones(1,numel(E(t).cols{c}))]; %#ok<AGROW>
        end
    end
    fprintf('lam = %.6f (%.2f mu): ||b|| %.4e; rows with no variable: %s\n', lams(j), lams(j)/mu, norm(S(j).b), mat2str(zr'));
    for r = zr'
        fprintf('    row %2d = (cell %d, position %d): b = %+.6e\n', r, cells(r), cols(r), S(j).b(r));
    end
    if isfield(E,'dup')
        for t = 1:numel(E)
            dp = E(t).dup;
            for k = 1:numel(dp.col)
                fprintf('    dropped (cell %d, position %d) = %+.6e x written (cell %d, position %d)\n', ...
                        dp.cell(k), dp.col(k), dp.factor(k), dp.rcell(k), dp.rcol(k));
            end
        end
    end
end
save(outfile,'S','lams','mu','same','berr','-v7.3');
cd(old);
end
