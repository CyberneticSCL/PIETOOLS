function probe4d_strip(label)
% PROBE4D_STRIP(LABEL) the 4-variable heat equation of probe4d_pipe with
% the stripped tier: the storage (4 decision variables), then the three
% products of D = T'PA + A'PT one at a time, each timed, so the memory
% watchdog of the runner localizes the wall seen with the light storage.
here = fileparts(mfilename('fullpath'));
diary(fullfile(here,[label '_diary.log']));
root = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
if isempty(which('poscopvar')),  run(fullfile(root,'pietools_path_update.m'));  end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
fprintf('--- %s (%s)\n',label,datestr(now,'yyyy-mm-dd HH:MM:SS'));
n = 4;
nm = arrayfun(@(k) sprintf('s%d',k),1:n,'UniformOutput',false);
s = polynomial(nm(:));
clear stateNameGenerator
x = pde_var(s,repmat([0 1],n,1));
rhs = diff(x,s(1),2);
for k = 2:n,    rhs = rhs + diff(x,s(k),2);    end
rhs = rhs + x;
PDE = diff(x,'t')==rhs;
for k = 1:n,    PDE = [PDE; subs(x,s(k),0)==0; subs(x,s(k),1)==0];    end
PIE = exec_pie_sop(convert_PIETOOLS_PDE(PDE));
Tm = PIE.T;     Am = PIE.A;
t0 = tic;   TT = Tm'*Tm;    fprintf('T''*T (copvar) %.1f s\n',toc(t0));
t0 = tic;   TA = Tm'*Am;    fprintf('T''*A (copvar) %.1f s\n',toc(t0));
st = exec_sop_settings(lpisettings_sop('stripped'));
tl = exec_tools_sop(PIE,st);    st.sop = tl.sop;
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);
t0 = tic;
[prog,Pm] = tl.storage(prog,Tm,'out');
Pm = Pm + tl.eye(1,tl.eppos);
fprintf('stripped storage %.1f s, %d dvars\n',toc(t0),size(prog.decvartable,1));
t0 = tic;   TP = Tm'*Pm;    fprintf('T''*P %.1f s\n',toc(t0));
t0 = tic;   TPA = TP*Am;    fprintf('(T''*P)*A %.1f s\n',toc(t0));
t0 = tic;   Dm = TPA + TPA';    fprintf('D = TPA + TPA'' %.1f s\n',toc(t0));
o = struct('like',Pm,'dw',st.sop.dw_feas,'margin',st.sop.margin,'cap',tl.capfor(Dm));
fprintf('lpi_ineq_sop: dw %d, cap %s\n',o.dw,mat2str(o.cap));
t0 = tic;
try
    [p3,~,li] = lpi_ineq_sop(prog,-Dm,o);
    fprintf('slack built %.1f s, %d dvars, degrees %s\n',toc(t0),size(p3.decvartable,1),mat2str(li.degrees));
catch ME
    fprintf('lpi_ineq_sop ERR %.1f s\n%s\n',toc(t0),getReport(ME,'extended','hyperlinks','off'));
end
fprintf('--- done (%s)\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
diary off
end
