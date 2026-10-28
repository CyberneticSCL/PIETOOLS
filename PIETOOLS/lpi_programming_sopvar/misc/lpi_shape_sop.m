function S = lpi_shape_sop(prog)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = LPI_SHAPE_SOP(PROG) returns the size of the SDP that sossolve passes
% to the solver for the LPI program PROG: the number of decision variables,
% the length of the free block, the sizes of the PSD blocks, the number of
% equality rows and the length of the primal vector.
%
% INPUT
% - prog:   'struct', an LPI program, solved; read are decvartable, var
%           (fields num, type, idx), expr (fields num, At) and, when
%           present, extravar (fields num, idx).
%
% OUTPUT
% - S:      'struct' with fields
%             ndv  double, the number of decision variables,
%                  numel(prog.decvartable);
%             Kf   double, the length of the free block: the
%                  prog.var.idx{1}-1 decision variables declared outside
%                  prog.var and the coefficients of every 'poly' entry of
%                  prog.var;
%             Ks   double row, the side lengths of the PSD blocks, one per
%                  'sos' entry of prog.var and per entry of prog.extravar,
%                  in descending order;
%             m    double, the number of equality rows, the columns of
%                  every prog.expr.At{i};
%             nx   double, the length of the primal vector,
%                  Kf + sum(Ks.^2).
%
% NOTES
% sossolve adds the PSD blocks of the inequality constraints
% (prog.extravar) during the solve, so for an unsolved program Ks omits
% them.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_shape_sop
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
% Initial coding MMP, 10/08/2026 (the count of the test helper cx_shape,
% in the library so that the executives depend on the library only).
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body; code unchanged. No help text deleted; the
%                initial-coding entry moved here from the help text.

% Decision variables; the free block starts with those declared outside
% prog.var
S = struct('ndv',numel(prog.decvartable),'Kf',NaN,'Ks',[],'m',0,'nx',0);
Kf = prog.var.idx{1}-1;
Ks = [];
% Each entry of prog.var: a 'poly' entry joins the free block, an 'sos'
% entry of n^2 coefficients is an n x n PSD block
for i = 1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    switch prog.var.type{i}
        case 'poly',    Kf = Kf + sz;
        case 'sos',     Ks(end+1) = sqrt(sz);                               %#ok<AGROW>
    end
end
% PSD blocks of the inequality constraints, added by sossolve
if isfield(prog,'extravar') && isfield(prog.extravar,'num')
    for i = 1:prog.extravar.num
        Ks(end+1) = sqrt(prog.extravar.idx{i+1}-prog.extravar.idx{i});    %#ok<AGROW>
    end
end
% Equality rows: the columns of every prog.expr.At{i}
m = 0;
for i = 1:prog.expr.num,    m = m + size(prog.expr.At{i},2);    end
% Block sizes in descending order; the primal length
S.Kf = Kf;  S.Ks = sort(round(Ks),'descend');   S.m = m;
S.nx = Kf + sum(S.Ks.^2);
end
