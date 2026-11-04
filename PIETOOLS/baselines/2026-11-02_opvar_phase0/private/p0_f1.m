function p0_f1(pass, cases)
% P0_F1(PASS,CASES) records fixture F1: PIESIM on the examples of the PIESIM
% libraries at fixed options, with the operators PIESIM integrates.
%
% INPUT
% - pass:   double, the pass number; output in ../F1/pass<PASS>/;
% - cases:  (optional) K x 2 double, rows [dimension, example]; default all
%           41 1-D and 32 2-D examples.
%
% OUTPUT (files)
% - <d>D_ex<NN>.mat: label; rec (status, format, seconds, note); PDE and
%   uinput as the library returns them; opts; solution and grid as PIESIM
%   returns them; cap, the inputs of PIESIM_time_integrate (psize, coeff,
%   gridall, Dop: every discretized operator) recorded by the capture shadow;
%   err_exact (1-D: per PDE state, max over the grid of |x - exact| at the
%   final time, as PIESIM_plot_solution forms it; [] otherwise); info (P0_PATH).
% - index.csv: one line per example.
%
% NOTES
% Options: N 16 (1-D) or 8 (2-D), tf 0.1, dt 1e-3, BDF order 2, no plots.
% The capture shadow (capture/PIESIM_time_integrate.m) forwards to the
% original through a handle bound before the shadow enters the path; both
% conditions are asserted.
if nargin<2 || isempty(cases)
    cases = [ones(41,1), (1:41)'; 2*ones(32,1), (1:32)'];
end
info = p0_path();
here = fileparts(mfilename('fullpath'));
out = fullfile(fileparts(here),'F1',sprintf('pass%d',pass));
if ~exist(out,'dir'),   mkdir(out);     end
global P0_CAPTURE
h = @PIESIM_time_integrate;
addpath(fullfile(here,'capture'),'-begin');
w = which('PIESIM_time_integrate','-all');
assert(startsWith(w{1},fullfile(here,'capture')),'p0_f1: capture shadow is not first on the path');
fi = functions(h);
assert(startsWith(fi.file,info.root),'p0_f1: handle not bound to the export');
newf = ~exist(fullfile(out,'index.csv'),'file');
idx = fopen(fullfile(out,'index.csv'),'a');
if newf
    fprintf(idx,'dim,example,format,status,seconds,captures,N,tf,dt,max_abs_final,err_exact_max,note\n');
end
for k = 1:size(cases,1)
    d = cases(k,1);     ex = cases(k,2);
    label = sprintf('%dD_ex%02d', d, ex);
    fprintf('\n@@@ F1 pass %d %s %s\n', pass, label, char(datetime('now')));
    P0_CAPTURE = struct('h',h,'n',0);
    rec = struct('status','ok','format','','seconds',NaN,'note','');
    PDE = [];   uinput = [];    opts = [];  solution = [];  grid = [];  cap = [];   err_exact = [];
    try
        if d==1,    [PDE,uinput] = examples_pde_library_PIESIM_1D(ex);
        else,       [PDE,uinput] = examples_pde_library_PIESIM_2D(ex);
        end
        rec.format = p0_format(PDE);
        opts = struct('N',16,'tf',0.1,'intScheme',1,'Norder',2,'dt',1e-3, ...
                      'plot','no','ploteig','no','ifexact',isfield(uinput,'exact'));
        if d==2,    opts.N = 8;     end
        t0 = tic;
        [solution, grid] = PIESIM(PDE, opts, uinput);
        rec.seconds = toc(t0);
        cap = rmfield(P0_CAPTURE,'h');
        if d==1 && opts.ifexact,    err_exact = p0_err_exact(solution, grid, uinput, cap.psize);    end
    catch ME
        rec.status = 'error';
        rec.note = p0_errstr(ME, info.root);
    end
    close all force
    save(fullfile(out,[label,'.mat']),'label','rec','PDE','uinput','opts','solution','grid', ...
         'cap','err_exact','info','-v7');
    mx = NaN;
    if isstruct(solution) && isfield(solution,'final') && isfield(solution.final,'primary')
        v = cellfun(@(c) max(abs(double(c(:))),[],'omitnan'), solution.final.primary(~cellfun(@isempty,solution.final.primary)));
        if ~isempty(v),     mx = max(v);    end
    end
    ncap = NaN;     if isstruct(cap),   ncap = cap.n;   end
    N = NaN;    tf = NaN;   dt = NaN;
    if isstruct(opts),  N = opts.N;     tf = opts.tf;   dt = opts.dt;   end
    fprintf(idx,'%d,%d,%s,%s,%.3f,%g,%g,%g,%g,%.12g,%.6g,%s\n', d, ex, rec.format, rec.status, ...
            rec.seconds, ncap, N, tf, dt, mx, max([err_exact(:);NaN],[],'omitnan'), strrep(rec.note,',',';'));
end
fclose(idx);
fprintf('\nDONE F1 pass %d\n', pass);
end

function f = p0_format(PDE)
% Input format of a library PDE, as convert_PIETOOLS_PDE dispatches on it.
if isa(PDE,'pde_struct'),                                       f = 'pde_var';
elseif isstruct(PDE) && any(isfield(PDE,{'n0','n1','n2'})),     f = 'batch';
elseif isstruct(PDE) && any(isfield(PDE,{'n','ODE','PDE'})),    f = 'terms_legacy';
elseif isstruct(PDE) && isfield(PDE,'x'),                       f = 'terms_struct';
else,                                                           f = class(PDE);
end
end

function e = p0_err_exact(solution, grid, uinput, psize)
% As PIESIM_plot_solution: exact(n) on the grid at the final time against
% solution.final.primary{2}(:,n), n over the PDE states; NaN where it fails.
syms sx
nx = sum(psize.n);
e = nan(1,nx);
for n = 1:nx
    try
        exs = subs(uinput.exact(n), sx, grid(:)');
        exv = double(subs(exs, solution.tf));
        e(n) = max(abs(solution.final.primary{2}(:,n) - exv(:)));
    catch
    end
end
end
