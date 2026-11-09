function p0_f4one(name, pass)
% P0_F4ONE(NAME,PASS) records fixture F4 for one DEMO script: it runs the
% script in a workspace of its own (the script clears it), then saves its
% variables, except the LPI programs, and a summary line.
%
% INPUT
% - name:   char, the DEMO script, e.g. 'DEMO1_stability_and_simulation';
% - pass:   double, the pass number; output in ../F4/pass<PASS>/.
%
% OUTPUT (files)
% - <name>.mat: vars, a struct of the script's variables after it ends
%   (classes kept: double, single, logical, char, string, cell, struct,
%   sym, polynomial, dpvar, copvar, cdopvar, sopvar, pie_struct, pde_struct;
%   names starting with 'prog' and arrays over 50 MB are dropped); rec
%   (status, seconds, note); info (P0_PATH).
% - index.csv: name, status, seconds, and every real scalar variable as
%   name=value.
global P0_DEMO
P0_DEMO = struct('name',name,'pass',pass,'info',p0_path(),'t0',tic,'err',[]);
here = fileparts(mfilename('fullpath'));
P0_DEMO.out = fullfile(fileparts(here),'F4',sprintf('pass%d',pass));
if ~exist(P0_DEMO.out,'dir'),   mkdir(P0_DEMO.out);     end
vars = p0_demo_run(name);
global P0_DEMO                                                              %#ok<REDEFGG>
rec = struct('status','ok','seconds',toc(P0_DEMO.t0),'note','');
if ~isempty(P0_DEMO.err)
    rec.status = 'error';
    rec.note = p0_errstr(P0_DEMO.err, P0_DEMO.info.root);
end
info = P0_DEMO.info;
close all force
save(fullfile(P0_DEMO.out,[name,'.mat']),'vars','rec','info','-v7');
sc = {};
fn = fieldnames(vars);
for k = 1:numel(fn)
    v = vars.(fn{k});
    if (isnumeric(v) || islogical(v)) && isscalar(v) && isreal(v)
        sc{end+1} = sprintf('%s=%.10g', fn{k}, double(v));                  %#ok<AGROW>
    end
end
f = fullfile(P0_DEMO.out,'index.csv');
newf = ~exist(f,'file');
fid = fopen(f,'a');
if newf,    fprintf(fid,'name,status,seconds,scalars,note\n');     end
fprintf(fid,'%s,%s,%.2f,%s,%s\n', name, rec.status, rec.seconds, strjoin(sc,' '), strrep(rec.note,',',';'));
fclose(fid);
end

function vars = p0_demo_run(p0_name__)
% Runs the script here; its 'clear' empties this workspace only.
global P0_DEMO
try
    run(p0_name__);
catch p0_err__
    global P0_DEMO                                                          %#ok<REDEFGG>
    P0_DEMO.err = p0_err__;
end
vars = struct();
p0_w__ = whos;
keep = {'double','single','logical','char','string','cell','struct','sym','polynomial', ...
        'dpvar','copvar','cdopvar','sopvar','sdopvar','pie_struct','pde_struct'};
for p0_k__ = 1:numel(p0_w__)
    p0_n__ = p0_w__(p0_k__).name;
    if startsWith(p0_n__,'p0_') || startsWith(p0_n__,'prog') || strcmp(p0_n__,'P0_DEMO'),   continue,   end
    if ~ismember(p0_w__(p0_k__).class, keep) || p0_w__(p0_k__).bytes > 50e6,  continue,   end
    vars.(p0_n__) = eval(p0_n__);
end
end
