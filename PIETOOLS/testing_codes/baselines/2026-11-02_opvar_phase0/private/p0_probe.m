function p0_probe()
% P0_PROBE: one-off checks before the Phase-0 runs: path resolution in the
% export, the capture shadow (the handle must reach the original), the
% size of the captured operators, the cost of apply_sopvar, the catalog
% builder of the 10/28 harness.
info = p0_path();
fprintf('root %s commit %s MATLAB %s MOSEK %s threads %d\n', info.root, info.commit, info.version, info.mosek, info.threads);
here = fileparts(mfilename('fullpath'));
global P0_CAPTURE
h = @PIESIM_time_integrate;
addpath(fullfile(here,'capture'),'-begin');
w = which('PIESIM_time_integrate','-all');
fprintf('PIESIM_time_integrate resolves first to %s (%d entries)\n', w{1}, numel(w));
fi = functions(h);
fprintf('handle: %s\n', fi.file);
% (1) capture on a 1-D and a 2-D example
for c = [1 10; 2 10]'
    P0_CAPTURE = struct('h',h,'n',0);
    if c(1)==1, [PDE,uinput] = examples_pde_library_PIESIM_1D(c(2));
    else,       [PDE,uinput] = examples_pde_library_PIESIM_2D(c(2));
    end
    opts = struct('N',16,'tf',0.01,'intScheme',1,'Norder',2,'dt',1e-3,'plot','no','ploteig','no','ifexact',false);
    if c(1)==2, opts.N = 8; end
    t0 = tic; [sol,~] = PIESIM(PDE,opts,uinput); t = toc(t0);
    cap = P0_CAPTURE;   cap = rmfield(cap,'h');
    f = [tempname,'.mat'];  save(f,'cap','sol','-v7');  d = dir(f);  delete(f);
    fn = fieldnames(cap.Dop);
    fprintf('%dD ex %d: PIESIM %.2f s, captures %d, Dop fields %d (%s), file %.1f MB\n', c(1), c(2), t, cap.n, numel(fn), strjoin(fn',' '), d.bytes/1e6);
end
% (2) apply_sopvar cost on the operators of a 1-D and a 2-D PIE
for c = [1 39; 2 10]'
    if c(1)==1, [PDE,~] = examples_pde_library_PIESIM_1D(c(2));
    else,       [PDE,~] = examples_pde_library_PIESIM_2D(c(2));
    end
    PIE = convert_PIETOOLS_PDE(PDE);
    t0 = tic; [nb, nops] = deal(0);
    for nm = {'T','Tw','Tu','A','Bw','Bu','Cz','Dzw','Dzu','Cy','Dyw','Dyu'}
        X = PIE.(nm{1});
        if ~isa(X,'copvar') || sum(X.dim_out)==0 || sum(X.dim_in)==0, continue, end
        nops = nops + 1;
        [~,~,Y,nbk] = p0_apply(X, 3, 5, 1);
        nb = nb + nbk;
        assert(all(cellfun(@(y) all(isfinite(y(:))), Y)));
    end
    fprintf('%dD ex %d: apply %d operators, %d blocks, %.2f s\n', c(1), c(2), nops, nb, toc(t0));
end
% (3) the 10/28 catalog builder on the first PDE and the first DDE row
S = load(fullfile(here,'catalog.mat'));  C = S.C;
fprintf('catalog: %d rows, fields %s\n', numel(C), strjoin(fieldnames(C)',' '));
for lib = {'PDE','DDE','NDSDDF'}
    k = find(strcmp({C.lib},lib{1}) & [C.run], 1);
    t0 = tic; [PIE,~,si] = lb_build_system(C(k), [], true);
    fprintf('%s %s: class(T) %s size %s, %.2f s\n', lib{1}, C(k).key, class(PIE.T), si.size, toc(t0));
end
end
