function p0_f2(pass, which_src)
% P0_F2(PASS,WHICH_SRC) records fixture F2: every library PIE, its operators
% as containers, and their action on random polynomial test functions.
%
% INPUT
% - pass:       double, the pass number; output in ../F2/pass<PASS>/;
% - which_src:  (optional) char, 'all' (default), 'piesim' or 'catalog'.
%
% OUTPUT (files)
% - <label>.mat: label; rec (source, class of PIE.T, number of variables,
%   status, seconds, note); fix, one entry per nonempty operator of the 12
%   (T, Tw, Tu, A, Bw, Bu, Cz, Dzw, Dzu, Cy, Dyw, Dyu) with
%     name, X        the operator, a 'copvar';
%     seed, xin      the seed and the test functions (P0_APPLY, q = 3);
%     pts, Y         points and outputs, Y{i} dim_out(i) x 3 x 5, computed
%                    block by block with apply_sopvar (the sopvar definition);
%     legacy         false; until 11/26/2026, true when the PIE held opvar
%                    operators (legacy formats);
%     Yop, relerr    [] and NaN; until 11/26/2026, for legacy operators,
%                    apply_opvar on the opvar at the same points and
%                    max|Y-Yop|/max(1,max|Y|);
%     dump           []; until 11/26/2026, for legacy operators, the opvar
%                    fields as a plain struct;
%   info (P0_PATH).
% - index.csv: one line per PIE.
%
% NOTES
% Sources: the 41 1-D and 32 2-D PIESIM examples (convert_PIETOOLS_PDE), and
% one row per (library, example) of the 10/28 catalog (lb_build_system,
% default parameters): the PDE, DDE and NDS/DDF libraries.
if nargin<2 || isempty(which_src),  which_src = 'all';  end
info = p0_path();
here = fileparts(mfilename('fullpath'));
out = fullfile(fileparts(here),'F2',sprintf('pass%d',pass));
if ~exist(out,'dir'),   mkdir(out);     end
ops = {'T','Tw','Tu','A','Bw','Bu','Cz','Dzw','Dzu','Cy','Dyw','Dyu'};
src = cell(0,3);
if any(strcmp(which_src,{'all','piesim'}))
    for ex = 1:41,  src(end+1,:) = {sprintf('PIESIM1D_%02d',ex),'PIESIM 1-D',@() pie_piesim(1,ex)};    end %#ok<AGROW>
    for ex = 1:32,  src(end+1,:) = {sprintf('PIESIM2D_%02d',ex),'PIESIM 2-D',@() pie_piesim(2,ex)};    end %#ok<AGROW>
end
if any(strcmp(which_src,{'all','catalog'}))
    S = load(fullfile(here,'catalog.mat'));     C = S.C;
    seen = {};
    for k = 1:numel(C)
        if ~C(k).run || strcmp(C(k).exec,'ERROR'),  continue,   end
        key = [C(k).lib,'_',C(k).id];
        if ismember(key,seen),  continue,   end
        seen{end+1} = key;                                                  %#ok<AGROW>
        m = C(k);
        src(end+1,:) = {['LIB_',regexprep(key,'[^A-Za-z0-9_]','_')],['catalog ',C(k).key],@() lb_build_system(m,[],true)}; %#ok<AGROW>
    end
end
newf = ~exist(fullfile(out,'index.csv'),'file');
idx = fopen(fullfile(out,'index.csv'),'a');
if newf
    fprintf(idx,'label,source,status,class_T,nv,n_ops,n_blocks,seconds,max_abs_Y,max_rel_opvar_vs_sop,note\n');
end
for s = 1:size(src,1)
    label = src{s,1};
    fprintf('\n@@@ F2 pass %d %s %s\n', pass, label, char(datetime('now')));
    rec = struct('source',src{s,2},'status','ok','class_T','','nv',NaN,'seconds',NaN,'note','');
    fix = struct('name',{},'X',{},'seed',{},'xin',{},'pts',{},'Y',{},'nb',{}, ...
                 'legacy',{},'Yop',{},'relerr',{},'dump',{},'err',{});
    t0 = tic;
    try
        PIE = src{s,3}();
        rec.class_T = class(PIE.T);
        for k = 1:numel(ops)
            f = p0_fixop(PIE.(ops{k}), ops{k}, label);
            if ~isempty(f),     fix(end+1) = f;     end                     %#ok<AGROW>
        end
        if ~isempty(fix) && isa(fix(1).X,'copvar'),     rec.nv = numel(fix(1).X.vars);  end
        bad = arrayfun(@(f) ~isempty(f.err), fix);
        if any(bad)
            rec.status = 'op_error';
            rec.note = strjoin(arrayfun(@(f) [f.name,': ',f.err], fix(bad), 'UniformOutput', false), ' ; ');
        end
    catch ME
        rec.status = 'error';
        rec.note = p0_errstr(ME, info.root);
    end
    rec.seconds = toc(t0);
    close all force
    save(fullfile(out,[label,'.mat']),'label','rec','fix','info','-v7');
    my = NaN;   mr = NaN;   nb = 0;
    for f = fix
        if ~isempty(f.Y),   my = max([my, max(cellfun(@(y) max(abs(y(:)),[],'omitnan'), f.Y))],[],'omitnan');  end
        if f.legacy,    mr = max([mr, f.relerr],[],'omitnan');  end
        nb = nb + f.nb;
    end
    fprintf(idx,'%s,%s,%s,%s,%g,%d,%d,%.3f,%.12g,%.3g,%s\n', label, rec.source, rec.status, rec.class_T, ...
            rec.nv, numel(fix), nb, rec.seconds, my, mr, strrep(rec.note,',',';'));
end
fclose(idx);
fprintf('\nDONE F2 pass %d\n', pass);
end

function PIE = pie_piesim(d, ex)
if d==1,    [PDE,~] = examples_pde_library_PIESIM_1D(ex);
else,       [PDE,~] = examples_pde_library_PIESIM_2D(ex);
end
PIE = convert_PIETOOLS_PDE(PDE);
end

function f = p0_fixop(Xin, name, label)
% One operator: [] when absent or without rows or columns.
f = [];
% BEGIN MMP, 11/26/2026: deleted: the branch for opvar and opvar2d inputs,
% converted by op2copvar_sop (retired 11/26; the classes left the path on
% 11/08). Once copvar is renamed opvar, containers would take that branch on
% the rename export only, and F2 would differ for a harness reason. The
% output keeps the field legacy, always false.
legacy = false;
if isa(Xin,'copvar')
    X = Xin;
else
    return                                  % [] or a numeric placeholder
end
% END MMP, 11/26/2026
if sum(X.dim_out)==0 || sum(X.dim_in)==0,   return,     end
seed = p0_seed([label,'|',name]);
f = struct('name',name,'X',X,'seed',seed,'xin',[],'pts',[],'Y',[],'nb',0, ...
           'legacy',legacy,'Yop',[],'relerr',NaN,'dump',[],'err','');
try
    [f.xin, f.pts, f.Y, f.nb] = p0_apply(X, 3, 5, seed);
    % MMP, 11/26/2026: deleted: for an opvar input, its fields as a plain
    % struct (dump) and apply_opvar at the same points (Yop, relerr).
catch ME
    f.err = p0_errstr(ME);
end
end

% MMP, 11/26/2026: deleted: P0_APPLY_OPVAR, apply_opvar on an opvar input at
% the points of the container (its one caller was the branch deleted above).

function s = p0_seed(str)
c = double(str);
s = mod(sum(c.*(1:numel(c))), 2^31-1);
end
