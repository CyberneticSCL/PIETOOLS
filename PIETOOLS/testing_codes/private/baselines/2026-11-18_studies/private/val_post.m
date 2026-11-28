function val_post(part)
% VAL_POST(PART) checks post_verify and the opts.post path of ex_search
% before the settings rerun (cleaning/private/settings_rerun_1121/PLAN.md):
%   'a' - library case 36 (PIETOOLS_PDE_Ex_2D_Telegraph_Eq), 'Q' light dD 0
%         wcross 1, refused on 11/21 by the measure gate and passed by pc_post
%         (stage 2 measure 1.121, eta 1.0e-8): post_verify must pass it; the
%         same storage against the PIE shifted by 10 must fail;
%         library case 40 (PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE), certified
%         directly: post_verify must pass it, and fail it shifted by 10;
%         library case 50 (PIETOOLS_PDDE_Ex), dual H2 'c' stripped: refused
%         and failed by pc_post: post_verify must fail it.
%   'b' - ex_search with opts.post on case 40 from its example st, 4 builds.
% results_VALPOST_<PART>.csv; ex_search logs exs_VALPOST*.csv.
% MMP, 11/21/2026.
here = fileparts(mfilename('fullpath'));
lb_setup();   cd(here);
addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
for a = {'LB_ST_HOOK','LB_FORM_MAP','LB_ROW_KEY','LB_LAST'}
    if isappdata(0,a{1}),  rmappdata(0,a{1});  end
end
fprintf('## VAL_POST %s %s\n## %s\n## %s\n',part,char(datetime('now')),which('PIETOOLS_stability'),which('sossolve'));
if strcmp(part,'b')
    sf = fullfile(here,'valpost_starts.csv');
    fid = fopen(sf,'w');  fprintf(fid,'key,kind,form,tier,mods,best\n');
    fprintf(fid,'PDE:40,stab,Q,light,examplest,1\n');  fclose(fid);
    ex_search(sf,'VALPOST',struct('post',true,'maxbuild',4,'stopmin',15,'statedir','_valpost'));
    return
end
of = fullfile(here,['results_VALPOST_',part,'.csv']);
fid = fopen(of,'w');  fprintf(fid,'row,s1_status,s1_fr,s1_exst,s1_eta,s1_ndec,delta,ok,s2_status,s2_fr,s2_exst,s2_eta,s2_ndec,note\n');  fclose(fid);
Cat = load(fullfile(here,'catalog.mat'));  Cat = Cat.C;
J = {'R36','PDE:36','PIETOOLS_PIE2PDEstability','Q','light','st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1',[0 10];
     'R62','PDE:62','PIETOOLS_PIE2PDEstability','P','light','',[0 10];
     'R40','PDE:40','PIETOOLS_PIE2PDEstability','Q','light','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 0',[0 10];
     'H50','PDE:50/H2_norm_dual','PIETOOLS_H2_norm_c','c','stripped','st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = ''auto''',0};
for i = 1:size(J,1)
    [id,key,ex,form,tier,mods,deltas] = J{i,:};
    fprintf('\n##### VAL %s %s %s %s [%s]  %s\n',id,key,form,tier,mods,char(datetime('now')));
    c = Cat(find(strcmp({Cat.key},key),1));
    PIE = lb_build_system(c,[],true);
    S = lpisettings(tier,0,'','mosek');
    eval([strrep(mods,'st.','S.'),';']);
    base = 'Q';  tgt = 'PIETOOLS_stability';
    if strcmp(form,'c'),  base = 'c';  tgt = 'PIETOOLS_H2_norm';  end
    setappdata(0,'LB_FORM_MAP',struct(base,form,'only',tgt));
    setappdata(0,'LB_KEEP_LAST',true);
    r1 = lb_exec(['sop:',ex],PIE,S);
    setappdata(0,'LB_KEEP_LAST',false);
    L = getappdata(0,'LB_LAST');
    for d = deltas
        note = '';  ok = false;  r2 = struct();
        try
            [ok,r2,pre] = post_verify(L.PIE,L.S,form,L.prog,L.info,d);
            note = pre;
        catch ME
            note = ['error: ',ME.message];
        end
        fid = fopen(of,'a');
        fprintf(fid,'%s,%s,%.4g,%s,%.3g,%d,%g,%d,%s,%.4g,%s,%.3g,%g,"%s"\n',id,f(r1,'status'),g(r1,'fr'),f(r1,'exst'), ...
            g(r1,'eta'),g(r1,'ndec'),d,ok,f(r2,'status'),g(r2,'fr'),f(r2,'exst'),g(r2,'eta'),g(r2,'ndec'), ...
            strrep(strrep(note,',',';'),'"',''''));
        fclose(fid);
        fprintf('##### VAL %s delta %g: stage 1 %s fr %.3g exst %s eta %.2g; stage 2 ok %d %s fr %.3g exst %s eta %.2g\n', ...
            id,d,f(r1,'status'),g(r1,'fr'),f(r1,'exst'),g(r1,'eta'),ok,f(r2,'status'),g(r2,'fr'),f(r2,'exst'),g(r2,'eta'));
    end
    rmappdata(0,'LB_LAST');
end
fprintf('\n## VAL_POST %s DONE %s\n',part,char(datetime('now')));
end

function s = f(r,k),  s = '';   if isstruct(r) && isfield(r,k) && ischar(r.(k)),  s = r.(k);  end,  end
function x = g(r,k),  x = NaN;  if isstruct(r) && isfield(r,k) && isnumeric(r.(k)) && isscalar(r.(k)),  x = double(r.(k));  end,  end
