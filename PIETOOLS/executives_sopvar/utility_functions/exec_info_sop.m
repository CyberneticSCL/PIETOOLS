function info = exec_info_sop(prog,aux,st,label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO = EXEC_INFO_SOP(PROG,AUX,ST,LABEL) returns the part of the INFO
% output that every container executive shares: the solver status, the
% degree loop, the slack of the negativity constraint, the size of the SDP
% and, when ST.sop.keepdual is set, the dual kernel of that constraint.
%
% INPUT
% - prog:   'struct', the solved LPI program;
% - aux:    'struct' returned by LPI_SOLVE_LOOP_SOP, with fields classify
%           (the 'struct' of LPI_CLASSIFY_SOP, field status), dw (double,
%           the raise of the last pass), hist (struct array, one entry per
%           pass) and li (the slack info of LPI_INEQ_SOP: fields degrees, a
%           'struct' with D and w, terms and tag);
% - st:     'struct' of LPI settings completed by EXEC_SOP_SETTINGS
%           (fields sop.slack, sop.keepdual, sop.verbose);
% - label:  char, the executive's name for printed messages.
%
% OUTPUT
% - info:   'struct' with fields
%             status      char: 'optimal', 'inaccurate', 'infeasible',
%                         'unbounded' or 'unknown';
%             certified   logical, true when status is 'optimal' or
%                         'inaccurate';
%             dw, hist    the last raise and the loop history;
%             slack       char, 'like' or 'stock';
%             degrees     'struct' with fields D and w, the lift and weight
%                         of the slack ([] for the stock slack);
%             terms       the Positivstellensatz terms of the slack;
%             shape       'struct' of LPI_SHAPE_SOP: ndv (decision
%                         variables), Kf, Ks (PSD block sizes), m (rows),
%                         nx (primal length);
%             dual        'copvar' object, the dual kernel of the
%                         negativity constraint ([] if not assembled);
%             dual_norm1  double, its 1-norm (NaN if not assembled).
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - exec_info_sop
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
% MMP, 10/28/2026: Moved from executives_sopvar/private/ to
%                utility_functions/ with the consolidation of the
%                executives; header to the header protocol, the licence
%                block added, step comments in the body; code unchanged.

% Status of the last pass and the SDP shape
cl = aux.classify;
info = struct('status',cl.status,'certified',any(strcmp(cl.status,{'optimal','inaccurate'})), ...
              'dw',aux.dw,'hist',aux.hist,'slack',st.sop.slack,'degrees',[],'terms',[], ...
              'shape',lpi_shape_sop(prog),'dual',[],'dual_norm1',NaN);
if ~isfield(aux,'li') || ~isstruct(aux.li),    return,     end
% Degrees and terms of the slack
li = aux.li;
if isfield(li,'degrees') && isstruct(li.degrees) && isfield(li.degrees,'D')
    info.degrees = struct('D',li.degrees.D,'w',li.degrees.w);
end
if isfield(li,'terms'),     info.terms = li.terms;     end
% Dual kernel of the tagged constraint, from the solver's multipliers
if st.sop.keepdual && isfield(li,'tag') && ~isempty(li.tag) && li.tag>0
    try
        [Y,di] = lpigetdual_sop(prog,li.tag);
        info.dual = Y;  info.dual_norm1 = di.norm1;
    catch ME
        if st.sop.verbose,  fprintf('  [%s] dual kernel not assembled: %s\n',label,ME.message);  end
    end
end
end
