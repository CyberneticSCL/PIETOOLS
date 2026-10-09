function info = exec_info_sop(prog,aux,st,label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO = EXEC_INFO_SOP(PROG,AUX,ST,LABEL) the part of the INFO output every
% container executive shares, from the solved program and the AUX of
% LPI_SOLVE_LOOP_SOP: status, certified (optimal or inaccurate), dw and the
% loop history, the slack mode, the slack degrees (D, w) and term codes of
% LPI_INEQ_SOP, the SDP shape (LPI_SHAPE_SOP) and, when ST.sop.keepdual is
% set, the dual kernel of the tagged inequality (LPIGETDUAL_SOP) and its
% 1-norm. The builders add their own fields.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cl = aux.classify;
info = struct('status',cl.status,'certified',any(strcmp(cl.status,{'optimal','inaccurate'})), ...
              'dw',aux.dw,'hist',aux.hist,'slack',st.sop.slack,'degrees',[],'terms',[], ...
              'shape',lpi_shape_sop(prog),'dual',[],'dual_norm1',NaN);
if ~isfield(aux,'li') || ~isstruct(aux.li),    return,     end
li = aux.li;
if isfield(li,'degrees') && isstruct(li.degrees) && isfield(li.degrees,'D')
    info.degrees = struct('D',li.degrees.D,'w',li.degrees.w);
end
if isfield(li,'terms'),     info.terms = li.terms;     end
if st.sop.keepdual && isfield(li,'tag') && ~isempty(li.tag) && li.tag>0
    try
        [Y,di] = lpigetdual_sop(prog,li.tag);
        info.dual = Y;  info.dual_norm1 = di.norm1;
    catch ME
        if st.sop.verbose,  fprintf('  [%s] dual kernel not assembled: %s\n',label,ME.message);  end
    end
end
end
