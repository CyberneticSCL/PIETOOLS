function st = exec_sop_settings(st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ST = EXEC_SOP_SETTINGS(ST) completes an LPI settings struct (lpisettings,
% settings_PIETOOLS_*) with the field ST.sop that the container executives
% (executives_sopvar) read. Every existing field of ST.sop is kept; the
% defaults are
%
%   slack      'like' (default) the slack of every negativity constraint is
%              sized by GET_LIFT_DEGS on the storage operator ('like'), with
%              the weight raised by dw_Q in a non-coercive form and by dw_P
%              in a coercive one (measured, lpi_programming_sopvar/README.md:
%              the reader on the storage gives the smallest sufficient slack
%              in 1-D and the clean 2-D solves); 'stock' the stock slack
%              (POSLPIVAR_SETTINGS_SOP 'slack': dd2/options2 + dd3/options3),
%              for comparison;
%   dw_Q       1    weight raise over the reader's rule, non-coercive forms;
%   dw_P       2    the same, coercive forms (T'PA + A'PT);
%   max_raise  2    the degree loop: on MOSEK status UNKNOWN the weight is
%              raised by one and the program rebuilt, at most this often
%              (LPI_SOLVE_LOOP_SOP);
%   margin     0    eps of the e_D margin of LPI_INEQ_SOP (0: none);
%   witness    true evaluate the numerical counterpart of the certificate
%              (PIE_WITNESS_SOP) and report the bracket;
%   N_cheb     24   Chebyshev order of that evaluation;
%   nfreq      240  frequency samples of the gain sweep;
%   keepdual   true record the equality rows (LPI_EQ_SOP) and assemble the
%              dual kernel of the negativity constraint (LPIGETDUAL_SOP).
%              sos_opts.simplify may stay on: sossolve restores the
%              multipliers into the original row coordinates (a removed
%              row gets 0), checked 10/08/2026 (b'y = gamma either way);
%   verbose    true print the stages and the result.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = struct('slack','like','dw_Q',1,'dw_P',2,'max_raise',2,'margin',0,'witness',true, ...
           'N_cheb',24,'nfreq',240,'keepdual',true,'verbose',true);
if ~isstruct(st),   error('exec_sop_settings:st','ST should be a settings struct (lpisettings).');   end
if ~isfield(st,'sop') || isempty(st.sop),   st.sop = struct();     end
fn = fieldnames(d);
for i = 1:numel(fn)
    if ~isfield(st.sop,fn{i}) || isempty(st.sop.(fn{i})),   st.sop.(fn{i}) = d.(fn{i});    end
end
if ~any(strcmp(st.sop.slack,{'like','stock'}))
    error('exec_sop_settings:slack','sop.slack should be ''like'' or ''stock''.');
end
if ~isfield(st,'sos_opts') || isempty(st.sos_opts),    st.sos_opts = struct();    end
if ~isfield(st.sos_opts,'solver') || isempty(st.sos_opts.solver)
    if ~isempty(which('mosekopt')),     st.sos_opts.solver = 'mosek';
    else,                               st.sos_opts.solver = 'sedumi';
    end
end
% if st.sop.keepdual,     st.sos_opts.simplify = 0;   end                   % MMP, 10/08/2026 (was)
end
