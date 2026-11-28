function st_d2size()
% ST_D2SIZE  Study A2 of PLAN_1118_studies.md: the size (decision
% variables, equality rows m; never solved) of the direct stability form
% 'Q' at the light tier on D2, the 2-D delayed heat equation with its 3-D
% delay state (library case 56) at tau = 0.5 tau*, tau* = 1.034726, for
% one-knob-at-a-time variants around the slack cap rungs: cap_delta(3) 0
% (475,404 variables, m 33,025: certifies only at no delay, D82) and 1
% (2,721,692, m 82,551). Knobs: the triangular cap lpi.ineq.cap_inc, the
% lift raise dD, the degree of Q in three variables deg_Q{3}, the 3-D term
% set. Base as in the Sol jobs: 3-D psatz 'faces1', wcross 0, no raise.
% The size-only LPISOLVE shadow (cleaning/private/settings_1116/d2/shadow)
% prints '@@ LB_SIZE ndv m ineq' and errors before solving; this function
% records '@@ A2' lines. Memory model (D90): MOSEK holds 23-25 bytes per
% m^2, so m <= 45,000 is about 47-51 GB. MMP, 11/18/2026.
% MMP, 11/23/2026: the root R is the checkout that holds this file (was the
%   live repo: a run from an export or a worktree put the live tree on the
%   path; stress test of 11/22/2026).
% R = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';         % MMP, 11/23/2026 (was)
% R: the checkout holding this file, its first ancestor with                % MMP, 11/23/2026
% pietools_path_update.m (the loop ends at the drive root if none);         % MMP, 11/23/2026
% '/' separators as in the literal, for the startsWith check below          % MMP, 11/23/2026
R = fileparts(mfilename('fullpath'));                                       % MMP, 11/23/2026
while ~isfile([R '/pietools_path_update.m']) && ~strcmp(R,fileparts(R))     % MMP, 11/23/2026
    R = fileparts(R);                                                       % MMP, 11/23/2026
end                                                                         % MMP, 11/23/2026
R = strrep(R,'\','/');                                                      % MMP, 11/23/2026
SH = [R '/cleaning/private/settings_1116/d2/shadow'];
t0 = tic;   cd(R);  pietools_path_update;   addpath(SH);
wl = strrep(which('lpisolve'),'\','/');
assert(startsWith(wl,SH),"lpisolve is not the size-only shadow: "+wl);
setappdata(0,'LB_SIZE_ONLY',true);
fprintf('@@ A2 start %s\n',char(datetime('now')));
PDE_t = PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman(0,{'tau=0.5*1.034726;'});
PDE = initialize_PIETOOLS_PDE(PDE_t,true);
PIE = initialize(convert(PDE,{'pie','silent'}));
fprintf('@@ A2 pie %.0f s, dim %d\n',toc(t0),PIE.dim);
% cap_delta(3), cap_inc, dD, deg_Q{3}, psatz{3}
V = {0,   [],    1, 1, 'faces1'      % control: 475,404 / 33,025
     1,   [],    1, 1, 'faces1'      % control: 2,721,692 / 82,551
     1,   'auto',1, 1, 'faces1'
     1,   [],    0, 1, 'faces1'
     1,   [],    1, 0, 'faces1'
     1,   [],    0, 0, 'faces1'
     1,   'auto',0, 1, 'faces1'
     1,   'auto',0, 0, 'faces1'
     2,   [],    0, 0, 'faces1'
     NaN, 'auto',0, 0, 'faces1'
     0,   [],    0, 1, 'faces1'
     1,   [],    1, 1, 'auto'};
for v = 1:size(V,1)
    [cd3,ci,dD,dq,ps] = V{v,:};
    st = lpisettings('light');   st.sos_opts.solver = 'mosek';
    st.lpi.ineq.psatz = {'auto','faces1',ps};   st.lpi.ineq.cap_inc = ci;   st.lpi.ineq.dD = dD;
    st.sop.cap_delta = [NaN NaN cd3];   st.sop.deg_Q{3} = dq;   st.sop.max_raise_Q(3) = 0;
    tag = sprintf('cd=%g capinc=%s dD=%d degQ3=%d psatz3=%s',cd3,mat2str(ci),dD,dq,ps);
    t1 = tic;
    try
        PIETOOLS_stability(PIE,st,'Q');
        fprintf('@@ A2 %s: returned without the size error (%.0f s)\n',tag,toc(t1));
    catch ME
        tok = regexp(ME.message,'LB_SIZE ndv=(\d+) m=(\d+) ineq=(\d+)','tokens','once');
        if isempty(tok)
            fprintf('@@ A2 %s: ERROR %s (%.0f s)\n',tag,strrep(ME.message,newline,' '),toc(t1));
        else
            m = str2double(tok{2});
            fprintf('@@ A2 %s: ndv=%s m=%s ineq=%s, MOSEK at 23-25 B/m^2 %.0f-%.0f GB (%.0f s)\n', ...
                tag,tok{1},tok{2},tok{3},23*m^2/1e9,25*m^2/1e9,toc(t1));
        end
    end
end
fprintf('@@ A2 end %s, %.0f s\n',char(datetime('now')),toc(t0));
end
