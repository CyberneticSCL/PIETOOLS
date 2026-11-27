function d6_time(root)
% D6_TIME(ROOT) - scratch timing, 11/27/2026. The D6 detection alone for 2
% to 16 spaces, and poslpivar / lpivar end to end in 1-D (2 spaces) and
% 2-D (4 spaces), median of 7 calls after one warm-up call, on tree ROOT.
root = strrep(root,'/','\');
restoredefaultpath;     addpath(genpath(root));
warning('off','MATLAB:rmpath:DirNotFound');
for f = {'PIETOOLS_demos\cuadmm','PIETOOLS_demos\lowrank_2d_stability','cleaning','orphaned_file'}
    rmpath(genpath(fullfile(root,f{1})));
end
addpath('C:\Program Files\Mosek\11.0\toolbox\r2019b');
fprintf('tree %s\n',root);
for nsp = [2 4 8 16]
    spaces = repmat({{'s'}},1,nsp);     dims = ones(nsp,1);     nin = 5; %#ok<NASGU>
    f = @() isnumeric(dims) && ~isempty(dims) && (nin==2 ...
        || (isempty(spaces) && numel(dims)>1) ...
        || (isnumeric(spaces) && ~isempty(spaces)) ...
        || (isstruct(spaces) && ~isfield(spaces,'out')) ...
        || (iscell(spaces) && any(cellfun(@(c) isnumeric(c) && ~isempty(c),spaces(:)))));
    fprintf('detection alone, %2d spaces: %6.2f us\n',nsp,1e6*timeit(f));
end
s = polynomial({'s'});      p1 = lpiprogram(s,[],[0 1]);
p2 = lpiprogram(polynomial({'s1';'s2'}),[],[0 1;0 2]);
sp4 = {{},{'s1'},{'s2'},{'s1','s2'}};
c = {'poslpivar 1-D', @() poslpivar(p1,[1;2],{{},{'s'}},[0 1],2);
     'lpivar 1-D',    @() lpivar(p1,[1;2],{{},{'s'}},[0 1],2);
     'poslpivar 2-D', @() poslpivar(p2,[1;1;1;1],sp4,[0 1;0 2],1);
     'lpivar 2-D',    @() lpivar(p2,[1;1;1;1],sp4,[0 1;0 2],2);
     'mat2opvar 2-D', @() mat2opvar(eye(4),[1;1;1;1],sp4,[0 1;0 2])};
for k = 1:size(c,1)
    f = c{k,2};     f();    t = zeros(1,7);
    for r = 1:7,    t0 = tic;   f();    t(r) = toc(t0);     end
    fprintf('%-14s median %8.2f ms  (min %8.2f, max %8.2f)\n',c{k,1},1e3*median(t),1e3*min(t),1e3*max(t));
end
end
