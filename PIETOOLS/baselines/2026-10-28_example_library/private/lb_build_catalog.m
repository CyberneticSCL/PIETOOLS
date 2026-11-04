function C = lb_build_catalog()
% C = LB_BUILD_CATALOG() builds the catalog of the three example
% libraries: one row per (example, executive toggle), written to
% catalog.csv and catalog.mat in this folder. No SDP is solved.
%
% OUTPUT
%   C:  struct array, one row per (example, toggle), fields of lb_meta plus
%       key     - 'LIB:ID' for the primary executive, 'LIB:ID/toggle' for
%                 each further toggle the example sets
%       type    - the toggle ('stability', 'Hinf_gain_dual', ...)
%       exec    - the executive the toggle selects (lb_toggles)
%       dim     - spatial dimension of the PDE (1 for DDE/NDS/DDF)
%       size    - first column of PIE.T.dim ('' for 2-D, not converted)
%       primary - true for the first toggle, which carries the bisection
%
% NOTES
%   - The DDE and NDS/DDF blocks are uncommented into blocks/<id>.m by
%     lb_blocks; each block is then evaluated and converted, so a block
%     extracted wrongly fails here, not in the baseline.
%   - For a bisection row the system is also built at the lower bracket;
%     'setter_effective' is false if the PIE does not change, which
%     detects a parameter override that the example ignores.
%   - 2-D PDEs are initialized but not converted (cost).

info = lb_setup();
M = lb_meta();
libdir = fullfile(info.root,'PIETOOLS_examples');
blockdir = fullfile(info.here,'blocks');
[names, execs] = lb_toggles();
C = struct([]);
for k = 1:numel(M)
    m = M(k);
    base = m;
    base.key = [m.lib,':',m.id];  base.type = '';  base.exec = '';
    base.dim = NaN;  base.size = '';  base.primary = true;
    base.setter_effective = NaN;
    if ~strcmp(m.lib,'PDE')
        lb_blocks(m, libdir, blockdir);
    end
    if ~m.run
        C = lb_cat(C, base);   continue
    end
    fprintf('\n=== catalog: %s\n', base.key);
    try
        [PIE, tg, si] = lb_build_system(m, [], 'auto');
        base.dim = si.dim;  base.size = si.size;
        if m.bisect
            PIE1 = lb_build_system(m, m.hi, 'auto');   % hi: never the default
            % Numeric fingerprints, not isequal on the opvar objects:
            % R2025b crashed (access violation) right after this step.
            f0 = [lb_fp(PIE.A), lb_fp(PIE.T)];   f1 = [lb_fp(PIE1.A), lb_fp(PIE1.T)];
            base.setter_effective = double(numel(f0)~=numel(f1) || ...
                any(abs(f0-f1) > 1e-12*(1+abs(f0))));
        end
    catch ME
        base.exec = 'ERROR';
        base.note = strtrim([base.note,' | catalog build error: ',ME.message]);
        C = lb_cat(C, base);   continue
    end
    on = [];
    for j = 1:numel(names)
        if tg.(names{j})==1, on(end+1) = j; end %#ok<AGROW>
    end
    if isempty(on)
        base.note = strtrim([base.note,' | no executive toggle set']);
        base.run = false;
        C = lb_cat(C, base);   continue
    end
    fam = @(t) erase(t,'_dual');
    for j = 1:numel(on)
        row = base;
        row.type = names{on(j)};  row.exec = execs{on(j)};
        if j>1
            row.key = [base.key,'/',names{on(j)}];
            row.primary = false;   row.bisect = false;
            if contains(row.type,'_dual') && ~isempty(m.docdual)
                row.doc = m.docdual;
            elseif ~strcmp(fam(row.type),fam(names{on(1)}))
                row.doc = '';
            end
            if m.bisect
                row.note = strtrim([row.note,' | single solve at the library default']);
            end
        end
        C = lb_cat(C, row);
    end
end
save(fullfile(info.here,'catalog.mat'),'C');
lb_write_catalog(C, fullfile(info.here,'catalog.csv'));
end

function C = lb_cat(C, row)
if isempty(C), C = row; else, C(end+1) = row; end
end

function lb_write_catalog(C, file)
% One CSV line per catalog row; commas inside fields become ';'.
fid = fopen(file,'w');
fprintf(fid,'library,id,name,type,executive,parameter,documented,bisect,lo,hi,dim,pie_size,run,setter_effective,note\n');
cl = @(s) strrep(strrep(char(string(s)),',',';'),newline,' ');
for k = 1:numel(C)
    c = C(k);
    id = c.key(numel(c.lib)+2:end);
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', c.lib, cl(id), cl(c.name), ...
        cl(c.type), cl(c.exec), cl(c.param), cl(c.doc), lb_yn(c.bisect), num2str(c.lo,6), ...
        num2str(c.hi,6), num2str(c.dim), cl(c.size), lb_yn(c.run), num2str(c.setter_effective), cl(c.note));
end
fclose(fid);
end

function s = lb_yn(b)
if b, s = 'yes'; else, s = 'no'; end
end
