function R = cx_hinf_struct2d(cases,parts,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = CX_HINF_STRUCT2D(CASES,PARTS) compares 2-D executives by STRUCTURE
% ONLY: nothing is solved. CASES is the second output of cx_cases_hinf
% (fields id, exec, setname, solver, plant, gam). PARTS (default
% {'cx','stock'}) selects the container side, the stock side or both, so
% that a container assembly of several minutes can run in its own process
% under the per-run time budget. Per case:
%
%   cx      the container program cx_<exec>(PIE,st,gam), assembly timed;
%   stock0  cx_stock_capture(exec,PIE,st): the executive's objective branch
%           (gain = 0, gam a decision variable), as the shared driver does;
%   stock1  cx_hinf_capture(exec,PIE,st,gam): the executive's own
%           fixed-gain branch, the same LPI the container builds.
%
% Shapes by cx_shape on the UNSOLVED programs (valid here: the 2-D
% executives impose no lpi_ineq, whose slack sossolve would add at solve
% time). 'dup' counts equality rows that repeat an earlier row exactly
% (At column and b), so a row-count difference made only of repeated rows
% is told apart from a real one; repeats do not change the feasible set.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 || isempty(parts),  parts = {'cx','stock'};  end
if nargin<3,    verbose = true;     end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
R = cell(1,numel(cases));
for ic = 1:numel(cases)
    c = cases(ic);
    r = struct('id',c.id,'exec',c.exec,'gam',c.gam,'err','','perr','');
    try
        st = cx_settings(c.setname,c.solver);
        PIE = cx_plant(c.plant{:});
        if ismember('cx',parts)
            t = tic;    pc = feval(['cx_' c.exec],PIE,st,c.gam);   r.t_cx = toc(t);
            r.cx = cx_shape(pc);        r.cx.dup = ndup(pc);    clear pc
        end
        if ismember('stock',parts)
            t = tic;    p0 = cx_stock_capture(c.exec,PIE,st);      r.t_s0 = toc(t);
            r.stock0 = cx_shape(p0);    r.stock0.dup = ndup(p0);    clear p0
            t = tic;    [p1,r.perr] = cx_hinf_capture(c.exec,PIE,st,c.gam);    r.t_s1 = toc(t);
            r.stock1 = cx_shape(p1);    r.stock1.dup = ndup(p1);    clear p1
        end
    catch ME
        r.err = sprintf('[%s] %s (%s:%d)',ME.identifier,ME.message, ...
                        ME.stack(1).name,ME.stack(1).line);
    end
    if verbose,     print_row(r);   end
    R{ic} = r;
end
end

function n = ndup(p)
% Equality rows (columns of [At; b']) repeating an earlier one exactly.
% Two random projections pick candidate pairs in O(nnz); candidates are
% then compared exactly.
M = [];
for i = 1:p.expr.num
    if strcmp(p.expr.type{i},'eq')
        M = [M, [p.expr.At{i}; reshape(full(p.expr.b{i}),1,[])]];          %#ok<AGROW>
    end
end
if isempty(M),  n = 0;  return,     end
rs = RandStream('mt19937ar','Seed',1);
h = [rand(rs,1,size(M,1))*M; rand(rs,1,size(M,1))*M].';
[~,~,g] = unique(round(h*1e12)/1e12,'rows');
n = 0;
for k = find(accumarray(g,1)>1).'
    cols = find(g==k);      keep = cols(1);
    for j = cols(2:end).'
        if any(arrayfun(@(i) isequal(M(:,i),M(:,j)),keep)),     n = n+1;
        else,                                                   keep(end+1) = j;    %#ok<AGROW>
        end
    end
end
end

function print_row(r)
if ~isempty(r.err),     fprintf('  %-32s ERROR %s\n',r.id,r.err);   return,     end
f = @(S) sprintf('ndv %d Kf %d m %d (dup %d) Ks %s',S.ndv,S.Kf,S.m,S.dup,mat2str(S.Ks));
s = sprintf('  %-32s gam %.4g',r.id,r.gam);
if isfield(r,'cx'),     s = [s sprintf(' | cx %s [%.1fs]',f(r.cx),r.t_cx)];     end
if isfield(r,'stock1')
    s = [s sprintf(' | stock(gain) %s [%.1fs] | stock(obj) %s [%.1fs]',f(r.stock1), ...
                   r.t_s1,f(r.stock0),r.t_s0)];
end
if ~isempty(r.perr),    s = [s ' | stock post-capture error: ' r.perr];     end
fprintf('%s\n',s);
end
