function [prog,aux,hist] = lpi_solve_loop_sop(build,st,label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,AUX,HIST] = LPI_SOLVE_LOOP_SOP(BUILD,ST,LABEL) builds and solves an
% LPI, raising the degree of its slacks while the solver status is
% 'unknown' or 'infeasible': for dw = 0, 1, ..., ST.sop.max_raise it builds
% [prog,aux] = BUILD(dw), solves prog by lpisolve with ST.sos_opts and
% classifies the solver exit, and it stops at the first status other than
% 'unknown' and 'infeasible' or after the pass dw = ST.sop.max_raise.
%
% INPUT
% - build:  'function_handle', [prog,aux] = build(dw): for a nonnegative
%           integer dw, the unsolved LPI program ('struct') with the weight
%           degree of every slack raised by dw, and a 'struct' aux whose
%           optional field obj is the objective, a 'dpvar' decision
%           variable ([] for a feasibility LPI);
% - st:     'struct' of LPI settings, completed by EXEC_SOP_SETTINGS; read
%           are sos_opts, sop.max_raise and sop.verbose;
% - label:  (optional) char, the name printed with each pass (default '').
%
% OUTPUT
% - prog:   'struct', the LPI program of the last pass, solved;
% - aux:    'struct', the aux of the last pass with the fields
%             classify  'struct' of the solver exit: status (char,
%                       'optimal', 'inaccurate', 'infeasible', 'unbounded',
%                       'unknown', or 'unsolved' without a solver report),
%                       obj (double, the objective value, NaN without one),
%                       divergent (logical, true for status 'unknown' with
%                       obj NaN or |obj| > divtol), divtol (double, 1e3),
%                       numerr, pinf, dinf, feasratio, cpusec (double, from
%                       the solver report);
%             dw        double, the raise of the last pass;
%             hist      HIST;
%           added;
% - hist:   'struct' array, one entry per pass, with fields dw, status,
%           divergent, obj, numerr, pinf, dinf (as in classify), ndv
%           (double, the number of decision variables,
%           numel(prog.decvartable)), m (double, the number of equality
%           rows of the SDP, the columns of every prog.expr.At{e}), t_build
%           and t_solve (double, seconds).
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_solve_loop_sop
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
% Initial coding MMP, 10/08/2026
% MMP, 10/16/2026: Header only. The earlier text said a raise changes
%                  nothing for a feasibility LPI; the Posopvar session
%                  measured the opposite on the 2-D stability executives
%                  (infeasible at dw 0, certified at dw 1; 2.57 s and
%                  2.48 s against 6.60 s and 5.68 s from the default
%                  dw_P = 2; the 1-D feasibility executives certify at
%                  dw 0 with 2.0-3.7x fewer decision variables at the
%                  same time). Code unchanged.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), the licence block added, step comments
%                  in the body; code unchanged. The deleted header text was
%                  the rationale for raising on each status, a notes
%                  reference and measurement results.

% Default label; complete the container options
if nargin<3,    label = '';     end
st = exec_sop_settings(st);
% Empty history; the first pass at raise 0
hist = struct('dw',{},'status',{},'divergent',{},'obj',{},'numerr',{},'pinf',{},'dinf',{}, ...
              'ndv',{},'m',{},'t_build',{},'t_solve',{});
dw = 0;
for it = 0:st.sop.max_raise
    % Build and solve at raise dw, timing each
    t0 = tic;   [prog,aux] = build(dw);     tb = toc(t0);
    t0 = tic;   prog = lpisolve(prog,st.sos_opts);  ts = toc(t0);
    % Classify the solver exit, reading the objective back when there is one
    obj = [];   if isfield(aux,'obj'),  obj = aux.obj;    end
    cl = lpi_classify_sop(prog,obj);
    % Equality rows of the SDP: the columns of every prog.expr.At{e}
    m = 0;
    if isfield(prog,'expr') && isfield(prog.expr,'At')
        for e = 1:numel(prog.expr.At),  m = m + size(prog.expr.At{e},2);    end
    end
    % Record and print the pass
    hist(end+1) = struct('dw',dw,'status',cl.status,'divergent',cl.divergent,'obj',cl.obj, ...
                         'numerr',cl.numerr,'pinf',cl.pinf,'dinf',cl.dinf, ...
                         'ndv',numel(prog.decvartable),'m',m,'t_build',tb,'t_solve',ts); %#ok<AGROW>
    if st.sop.verbose
        fprintf('  [%s] dw %d: %s%s, objective %s, %d decision variables, %d rows, build %.1f s, solve %.1f s\n', ...
                label,dw,cl.status,tern(cl.divergent,' (divergent)',''),num2str(cl.obj,'%.8g'), ...
                numel(prog.decvartable),m,tb,ts);
    end
    % Stop unless the status is 'unknown' or 'infeasible' and a raise is
    % left; otherwise raise the slack weight by one
%   if ~strcmp(cl.status,'unknown') || it==st.sop.max_raise,   break,  end % MMP, 10/08/2026 (was)
    if ~any(strcmp(cl.status,{'unknown','infeasible'})) || it==st.sop.max_raise,  break,  end % MMP, 10/08/2026
    dw = dw+1;
end
% The last pass, with its classification, its raise and the history
aux.classify = cl;  aux.dw = dw;    aux.hist = hist;
end


function s = tern(c,a,b)
% S = TERN(C,A,B) returns A when C is true and B otherwise.
% INPUT  c: logical; a, b: the two values.
% OUTPUT s: A or B.
if c,   s = a;  else,   s = b;  end
end
