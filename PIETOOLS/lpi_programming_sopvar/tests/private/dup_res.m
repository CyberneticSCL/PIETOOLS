function [da,db] = dup_res(prog,ma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DA,DB] = DUP_RES(PROG,MA) returns the residual norms of the equations
% that LPI_EQ_CDOPVAR did not write because each is f times a written one
% (prog.sopeq, entries.dup): sqrt(sum (f*r)^2), r the residual of the
% written row at x = prog.solinfo.RRx; DA over the records that start in
% the first MA rows of prog.expr, DB over the rest.
%
% INPUT
% - prog:   'struct', a solved LPI program with the field sopeq, which
%           LPI_EQ_CDOPVAR fills when it exists (set prog.sopeq = {} before
%           the equalities);
% - ma:     'double', the number of rows of the first equalities; Inf puts
%           every record in DA.
%
% OUTPUT
% - da, db: 'double'.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - dup_res
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 11/02/2026, as a local function of
% e2e_heat3d_stability.
% MMP, 11/07/2026: moved here unchanged, shared by e2e_heat3d_stability,
%                e2e_hinf_gain_1d and e2e_stability_1d, whose checks of the
%                rows against the extracted equalities need the rows not
%                written (since 11/01/2026). The identifier of the assertion
%                was 'heat3d_stability_sop:check', the test's own.

x = prog.solinfo.RRx(:);
% Residual of every row, in the row order of prog.expr (the record's)
res = cell(prog.expr.num,1);
for i = 1:prog.expr.num
    At = prog.expr.At{i};   res{i} = At'*x(1:size(At,1)) - prog.expr.b{i};
end
res = full(vertcat(res{:}));
d2 = [0 0];
for k = 1:numel(prog.sopeq)
    rec = prog.sopeq{k};    r = rec.rows(1);
    for e = 1:numel(rec.entries)
        en = rec.entries(e);
        % First row of each parameter cell of the block, then the next block
        nk = reshape(cellfun(@numel,en.cols),1,[]);     r1 = r + [0, cumsum(nk(1:end-1))];
        r = r + sum(nk);
        % A record written before 11/01/2026 has no field dup: none dropped
        if ~isfield(en,'dup') || isempty(en.dup.rcol),     continue,   end
        % Row of each written equation: its cell's first row plus its place
        dp = en.dup;    row = zeros(size(dp.rcol));
        for t = 1:numel(en.cells)
            s = dp.rcell==en.cells(t);
            [tf,loc] = ismember(dp.rcol(s),en.cols{t});
            row(s) = tf.*(r1(t)+loc-1);
        end
%       assert(all(row>0),'heat3d_stability_sop:check','(2) a dropped equation refers to no written row') % MMP, 11/07/2026 (was)
        assert(all(row>0),'dup_res:row','an equation not written refers to no written row') % MMP, 11/07/2026
        g = 1 + (rec.rows(1)>ma);
        d2(g) = d2(g) + sum((dp.factor.*res(row)).^2);
    end
end
da = sqrt(d2(1));   db = sqrt(d2(2));
end
