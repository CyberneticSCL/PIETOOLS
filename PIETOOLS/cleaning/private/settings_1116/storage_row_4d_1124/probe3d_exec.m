function probe3d_exec(label)
% PROBE3D_EXEC(LABEL) the 3-variable heat equation of probe4d_pipe:
% (a) the first attempt of PIETOOLS_stability_sop 'P' at light,
% reproduced step by step, with the full report of the error that
% lpi_ineq_sop raised in probe4d_pipe; (b) the executive itself at light.
here = fileparts(mfilename('fullpath'));
diary(fullfile(here,[label '_diary.log']));
root = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
if isempty(which('poscopvar')),  run(fullfile(root,'pietools_path_update.m'));  end
w = strrep(which('poscopvar'),'\','/');
assert(strcmpi(w,[root '/lpi_programming_sopvar/poscopvar.m']),'poscopvar resolves to %s',w);
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
fprintf('--- %s (%s)\n',label,datestr(now,'yyyy-mm-dd HH:MM:SS'));
n = 3;
nm = arrayfun(@(k) sprintf('s%d',k),1:n,'UniformOutput',false);
s = polynomial(nm(:));
clear stateNameGenerator
x = pde_var(s,repmat([0 1],n,1));
rhs = diff(x,s(1),2);
for k = 2:n,    rhs = rhs + diff(x,s(k),2);    end
rhs = rhs + x;
PDE = diff(x,'t')==rhs;
for k = 1:n,    PDE = [PDE; subs(x,s(k),0)==0; subs(x,s(k),1)==0];    end
PIE = convert_PIETOOLS_PDE(PDE);
% (a) the replica, full error report
PIEc = exec_pie_sop(PIE);
st = exec_sop_settings(lpisettings_sop('light'));
tl = exec_tools_sop(PIEc,st);    st.sop = tl.sop;
prog = lpiprogram(PIEc.vars(:,1),PIEc.vars(:,2),PIEc.dom,st);
Tm = PIEc.T;     Am = PIEc.A;
[prog,Pm] = tl.storage(prog,Tm,'out');
Pm = Pm + tl.eye(1,tl.eppos);
t0 = tic;   Dm = Tm'*Pm*Am + Am'*Pm*Tm;    fprintf('(a) D composed %.1f s\n',toc(t0));
o = struct('like',Pm,'dw',st.sop.dw_feas,'margin',st.sop.margin,'cap',tl.capfor(Dm));
t0 = tic;
try
    [p3,~,li] = lpi_ineq_sop(prog,-Dm,o);
    fprintf('(a) slack built %.1f s, %d dvars, degrees %s\n',toc(t0),size(p3.decvartable,1),mat2str(li.degrees));
catch ME
    fprintf('(a) lpi_ineq_sop ERR %.1f s\n%s\n',toc(t0),getReport(ME,'extended','hyperlinks','off'));
end
% (b) the executive
t0 = tic;
try
    [~,~,info] = PIETOOLS_stability_sop(PIE,lpisettings_sop('light'),'P');
    fprintf('(b) executive done %.1f s\n',toc(t0));   disp(info)
catch ME
    fprintf('(b) executive ERR %.1f s\n%s\n',toc(t0),getReport(ME,'extended','hyperlinks','off'));
end
fprintf('--- done (%s)\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
diary off
end
