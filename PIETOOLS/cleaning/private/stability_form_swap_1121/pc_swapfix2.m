function pc_swapfix2(part)
% PC_SWAPFIX2(PART) tests the change of 11/21/2026 to PIETOOLS_Hinf_gain_sop
% and PIETOOLS_H2_norm_sop: every substitution of a form warns and is
% appended to the label; the programs are unchanged.
%
% INPUT
%   part:  'g0' - the rows on the executives before the change;
%          'g1' - the same rows after the change.
% Rows: a substitution case (expect the warning after the change) and a
% case with none (expect no warning), each with gam, the SDP size, the
% status and the last warning. results_SWAP2_<PART>.csv in this folder.
% MMP, 11/21/2026.
here = fileparts(mfilename('fullpath'));
lb_setup();   cd(here);
addpath(fullfile(here,'gate'));
for a = {'LB_ST_HOOK','LB_FORM_MAP','LB_ROW_KEY'}
    if isappdata(0,a{1}),  rmappdata(0,a{1});  end
end
fprintf('## PC_SWAPFIX2 %s %s\n## %s\n## %s\n',part,char(datetime('now')), ...
    which('PIETOOLS_Hinf_gain_sop'),which('PIETOOLS_H2_norm_sop'));
of = fullfile(here,['results_SWAP2_',part,'.csv']);
fid = fopen(of,'w');
fprintf(fid,'system,exec,form,status,exst,gam,ndec,m,sec,warning\n');  fclose(fid);
% The systems: periodic heat (PIE.K), a boundary disturbance (Tw ~= 0),
% 2-D heat (two variables), Dirichlet heat (no substitution)
P = struct();
pvar t s s1 s2
x = pde_var(s,[0,1]);   w = pde_var('in');   z = pde_var('out');
P.periodic = convert(initialize([diff(x,t)==diff(x,s,2)-x+w; z==int(x,s,[0,1]);
    subs(x,s,0)==subs(x,s,1); subs(diff(x,s),s,0)==subs(diff(x,s),s,1)]),'pie');
x = pde_var(s,[0,1]);   w = pde_var('in');   z = pde_var('out');
P.boundary = convert(initialize([diff(x,t)==diff(x,s,2)-x; z==int(x,s,[0,1]);
    subs(x,s,0)==w; subs(x,s,1)==0]),'pie');
x = pde_var(s,[0,1]);   w = pde_var('in');   z = pde_var('out');
P.dirichlet = convert(initialize([diff(x,t)==diff(x,s,2)-x+w; z==int(x,s,[0,1]);
    subs(x,s,0)==0; subs(x,s,1)==0]),'pie');
x = pde_var(1,[s1;s2],[0,1;0,1]);   w = pde_var('in');   z = pde_var('out');
P.heat2d = convert(initialize([diff(x,t)==diff(x,s1,2)+diff(x,s2,2)+w; z==int(x,[s1;s2],[0,1;0,1]);
    subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0]),'pie');
% Rows {system, stock executive, form}
J = {'periodic','PIETOOLS_Hinf_gain','Qd';  'periodic','PIETOOLS_Hinf_gain','Pd';
     'periodic','PIETOOLS_H2_norm_c','c';   'periodic','PIETOOLS_H2_norm_c','cco';
     'boundary','PIETOOLS_Hinf_gain','Q';
     'heat2d','PIETOOLS_Hinf_gain','Qd';
     'dirichlet','PIETOOLS_Hinf_gain','Q';  'dirichlet','PIETOOLS_Hinf_gain','Qd';
     'dirichlet','PIETOOLS_H2_norm_c','c'};
for i = 1:size(J,1)
    [sy,ex,form] = J{i,:};
    fprintf('\n##### SWAP2 %s %s %s  %s\n',sy,ex,form,char(datetime('now')));
    lastwarn('');   t0 = tic;   r = struct('status','error');
    try
        S = lpisettings_sop('light');
        base = 'Q';  tgt = 'PIETOOLS_Hinf_gain_sop';
        if contains(ex,'H2'),  base = 'c';  tgt = 'PIETOOLS_H2_norm_sop';  end
        setappdata(0,'LB_FORM_MAP',struct(base,form,'only',tgt));
        r = lb_exec(['sop:',ex], P.(sy), S);
    catch ME
        r.status = 'error';   r.note = ME.message;
    end
    [msg,id] = lastwarn;
    if isfield(r,'note'),  msg = ['ERROR ',r.note];  id = '';  end
    fid = fopen(of,'a');
    fprintf(fid,'%s,%s,%s,%s,%s,%.8g,%g,%g,%.1f,"[%s] %s"\n',sy,ex,form,f(r,'status'),f(r,'exst'), ...
        g(r,'gam'),g(r,'ndec'),g(r,'m'),toc(t0),id,strrep(strrep(msg,'"',''''),newline,' '));
    fclose(fid);
    fprintf('##### SWAP2 %s %s %s: %s gam %.8g ndec %g m %g; warning [%s] %s\n',sy,ex,form, ...
        f(r,'status'),g(r,'gam'),g(r,'ndec'),g(r,'m'),id,msg);
end
fprintf('\n## PC_SWAPFIX2 %s DONE %s\n',part,char(datetime('now')));
end

function s = f(r,k),  s = '';   if isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end
