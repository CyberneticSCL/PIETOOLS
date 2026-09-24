function st = t1d_set(spec)                                                 % CC, 09/22/2026
% PROVENANCE.  scratchpad/ladder/lset.m verbatim, renamed t1d_set so that the
% banked ladder harness and this suite can be on the MATLAB path at the same
% time without shadowing each other (defect B6: two copies of one name have
% silently tested the wrong tree here before).  Body unmodified.
%                                                                 % CC, 09/22/2026
% lset(spec) -- 1-D LPI settings.  spec is either a stock name
% ('light','heavy','veryheavy',...) or a numeric vector [n1 n2 n3 Dup]
% reproducing the structure of settings_PIETOOLS_heavy with the four degree
% knobs exposed, so m can be pushed smoothly beyond 'veryheavy'.
%   heavy      == [2 1 1 1]
%   veryheavy  == [3 3 3 2]
if ischar(spec) || isstring(spec)
    st = lpisettings(char(spec));
    return
end
n1 = spec(1);  n2 = spec(2);  n3 = spec(3);  Dup = spec(4);
n4 = n2+n3;
st.override1 = 1;
st.options1.sep = 0;        st.options1.exclude = [0 0 0 0];
st.options12.sep = 0;       st.options12.exclude = [0 0 0 0];
st.options12.psatz = 1;
st.dd1  = {n1, [n2 n3 n4], [n2 n3 n4]};
st.dd12 = {n1, [n2 n3 n4], [n2 n3 n4]};
st.ddZ  = [2*n1 2*n2 2*n3];
st.sosineq_on = 0;
st.override2 = 0;
st.options2.psatz = 0;      st.options2.exclude = [0 0 0 0];
st.dd2 = {n1+1, [n2+Dup n3+Dup n4+Dup], [n2+Dup n3+Dup n4+Dup]};
st.options3.psatz = 1;      st.options3.exclude = [0 0 0 0];
st.dd3 = {n1, [n2+Dup-1 n3+Dup-1 n4+Dup-1], [n2+Dup-1 n3+Dup-1 n4+Dup-1]};
% same margins lpisettings installs
st.eppos = 1e-4;   st.eppos2 = 1e-6;   st.epneg = 0;
st.sos_opts.simplify = false;
% Unused in 1-D, kept so the struct is field-for-field what the banked ladder
% produced and a spec can be carried between the two harnesses unchanged.
st.settings_2d = settings_PIETOOLS_heavy_2D();
st.settings_2d.eppos = [1e-4;1e-6;1e-6;1e-6];
st.settings_2d.epneg = 0;
end
