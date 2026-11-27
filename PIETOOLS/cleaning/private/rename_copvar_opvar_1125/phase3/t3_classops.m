function t3_classops(root)
% T3_CLASSOPS(ROOT), 11/27/2026: timing spot check of container operations
% on tree ROOT, before (copvar names) or after (opvar names) the rename;
% median of 5 calls after one warm-up call. 1-D: spaces R x L2^2[s]; 2-D:
% R x L2[s1] x L2[s2] x L2[s1,s2]. Operations: positive declaration,
% free declaration, fixed multiplier, composition T'*P*T, sum P + P',
% composition Z*T, and the decision-variable count of each declaration.
root = strrep(root,'/','\');
restoredefaultpath;     addpath(genpath(root));
warning('off','MATLAB:rmpath:DirNotFound');
for f = {'PIETOOLS_demos\cuadmm','PIETOOLS_demos\lowrank_2d_stability','cleaning','orphaned_file'}
    rmpath(genpath(fullfile(root,f{1})));
end
addpath('C:\Program Files\Mosek\11.0\toolbox\r2019b');
if exist('poscopvar','file')      % before the rename
    pos = @poscopvar;   free = @lpivar_cdopvar;     m2o = @mat2copvar_sop;  tag = 'copvar names';
else
    pos = @poslpivar;   free = @lpivar;             m2o = @mat2opvar;       tag = 'opvar names';
end
fprintf('tree %s (%s)\n',root,tag);
cases = {'1-D', {'s'}, [1;2], {{},{'s'}}, [0 1], 2;
         '2-D', {'s1';'s2'}, [1;1;1;1], {{},{'s1'},{'s2'},{'s1','s2'}}, [0 1;0 2], 1};
for c = 1:size(cases,1)
    [lab,v,dm,sp,dom,deg] = cases{c,:};
    prog = lpiprogram(polynomial(v),[],dom);
    [~,P] = pos(prog,dm,sp,dom,deg);    [~,Z] = free(prog,dm,sp,dom,deg);
    T = m2o(eye(sum(dm)),dm,sp,dom);
    ops = {'positive declaration', @() pos(prog,dm,sp,dom,deg);
           'free declaration',     @() free(prog,dm,sp,dom,deg);
           'fixed multiplier',     @() m2o(eye(sum(dm)),dm,sp,dom);
           'T''*P*T',              @() T'*P*T;
           'P + P''',              @() P + P';
           'Z*T',                  @() Z*T};
    for k = 1:size(ops,1)
        f = ops{k,2};   f();    t = zeros(1,5);
        for r = 1:5,    t0 = tic;   f();    t(r) = toc(t0);     end
        fprintf('%s %-22s median %9.2f ms (min %9.2f, max %9.2f)\n',lab,ops{k,1},1e3*median(t),1e3*min(t),1e3*max(t));
    end
    fprintf('%s decision variables: positive %d, free %d\n',lab,numel(P.Zd),numel(Z.Zd));
end
end
