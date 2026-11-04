function lb_pie_ab(mode)
% LB_PIE_AB(MODE), 10/28/2026: the PIE of every runnable catalog system,
% before and after integrating origin/ndopvar (DJ's pde_struct change sets
% tdiff = NaN for a state without terms).
%
% INPUT
%   mode:   'save'    - build each PIE, write pie_before.mat
%           'compare' - rebuild each PIE and compare it with pie_before.mat,
%                       every stored property bit for bit; then run the two
%                       periodic-converter tests and DEMO09_PBC
%
% OUTPUT (printed)
%   lines 'GIT', 'AB', 'TEST', 'DEMO'; the conversion output goes to the log.

here = fileparts(mfilename('fullpath'));
info = lb_setup();
fprintf('GIT %s\n', info.git);
S = load(fullfile(here,'catalog.mat'));  C = S.C;
C = C([C.run] & ~strcmp({C.exec},'ERROR'));
% One system per library entry: the '/dual' rows share their system
sys = strcat({C.lib},':',{C.id});
[~,iu] = unique(sys,'stable');  C = C(iu);
F = fullfile(here,'pie_before.mat');
if strcmp(mode,'save')
    B = cell(numel(C),2);
    for k = 1:numel(C)
        [B{k,1},B{k,2}] = pie_bytes(C(k));
    end
    keys = strcat({C.lib},':',{C.id});
    save(F,'B','keys');
    fprintf('SAVE %d systems, %d errors\n', numel(C), nnz(~cellfun(@isempty,B(:,2))));
    return
end
L = load(F);
nsame = 0;  ndiff = 0;  nmiss = 0;
for k = 1:numel(C)
    key = [C(k).lib,':',C(k).id];
    j = find(strcmp(L.keys,key),1);
    if isempty(j),  nmiss = nmiss+1;  fprintf('AB %s not in pie_before\n',key);  continue,  end
    [b,e] = pie_bytes(C(k));
    if ~isempty(e) || ~isempty(L.B{j,2})
        same = strcmp(e,L.B{j,2});
        what = sprintf('error before "%s" after "%s"',L.B{j,2},e);
    else
        d = setxor(fieldnames(b),fieldnames(L.B{j,1}));
        f = intersect(fieldnames(b),fieldnames(L.B{j,1}));
        d = [d(:); f(~cellfun(@(p) isequal(b.(p),L.B{j,1}.(p)),f))];
        same = isempty(d);
        what = strjoin(d',' ');
    end
    if same,  nsame = nsame+1;  else,  ndiff = ndiff+1;  fprintf('AB %s DIFF %s\n',key,what);  end
end
fprintf('AB %d same, %d differ, %d missing (of %d systems)\n',nsame,ndiff,nmiss,numel(C));
for tn = {'test_convert_PIETOOLS_PDE_periodic','test_convert_PIETOOLS_PDE_periodic_higher_order'}
    try
        res = feval(tn{1});
        fprintf('TEST %-48s %d of %d pass\n',tn{1},nnz([res.passed]),numel(res));
    catch ME
        fprintf('TEST %-48s %s\n',tn{1},ME.message);
    end
end
t0 = tic;
try,    run_demo(fullfile(info.root,'PIETOOLS_demos','DEMO09_PBC.m'));
        fprintf('DEMO DEMO09_PBC OK (%.0f s)\n',toc(t0));
catch ME,   fprintf('DEMO DEMO09_PBC FAIL %s\n',ME.message);
end
end


function [b,e] = pie_bytes(c)
% Serialized stored properties of the PIE of catalog row C; E the error
% message ('' on success)
b = struct();  e = '';
try
    PIE = lb_build_system(c, [], true);
catch ME
    e = ME.message;  return
end
if isstruct(PIE)
    p = fieldnames(PIE)';
else
    mc = metaclass(PIE);  P = mc.PropertyList';
    p = {P(~[P.Dependent] & ~[P.Constant]).Name};
end
for k = 1:numel(p)
    try,    b.(p{k}) = getByteStreamFromArray(PIE.(p{k}));
    catch,  b.(p{k}) = 'unserializable';
    end
end
end

function run_demo(f)
% Run the script F in this workspace, so that its 'clear' stays here
run(f);
end
