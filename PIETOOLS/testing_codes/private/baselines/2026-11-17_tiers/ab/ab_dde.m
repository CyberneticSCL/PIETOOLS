function ab_dde(tree)
% AB_DDE(TREE) the light tier at three DDE bisection points where the HEAD
% tree certifies and the working tree returns infeasible (library pass
% 11/17/2026): one line per (variant, case) with the harness verdict and
% the program sizes. Run from a harness folder (lb_setup, catalog.mat);
% TREE 'new' adds the knob variants of the new layout, 'head' the old.
here = pwd;
info = lb_setup();  cd(here);  addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',2e6);   setappdata(0,'LB_M_CAP',45000);
S = load(fullfile(here,'catalog.mat'));  C = S.C;
cases = {'DDE:DDE05',1.5707; 'DDE:DDE08',1.3998; 'DDE:DDE03',1.69101; 'DDE:DDE07',1.35191};
V = {};
V{end+1} = struct('name','V0 light',            'f',@(s) s);
V{end+1} = struct('name','V9 no raise',         'f',@(s) setfield(s,'sop',setfield(s.sop,'max_raise_feas',[0 1 0])));
V{end+1} = struct('name','V4 slack stock',      'f',@(s) setfield(s,'sop',setfield(s.sop,'slack','stock')));
V{end+1} = struct('name','V7 dw_feas 1',        'f',@(s) setfield(s,'sop',setfield(s.sop,'dw_feas',[1 1 1])));
V{end+1} = struct('name','V5 dD 2',             'f',@(s) setfield(s,'lpi',setfield(s.lpi,'ineq',setfield(s.lpi.ineq,'dD',2))));
V{end+1} = struct('name','V6 dw 2',             'f',@(s) setfield(s,'lpi',setfield(s.lpi,'ineq',setfield(s.lpi.ineq,'dw',2))));
V{end+1} = struct('name','V8 var.deg heavy',    'f',@(s) setfield(s,'lpi',setfield(s.lpi,'var',struct('deg',{{[4 2 2],1,1}}))));
if strcmp(tree,'new')
    V{end+1} = struct('name','V1 len_sos [3;2;2]','f',@(s) setfield(s,'lpi',setfield(s.lpi,'pos',setfield(s.lpi.pos,'len_sos',{{[3;2;2],[1 1;2 2;0 0],[1 1 1;2 2 2;0 0 0]}}))));
    V{end+1} = struct('name','V1b len_inc/com 3', 'f',@(s) setfield(s,'lpi',setfield(s.lpi,'pos',setfield(setfield(s.lpi.pos,'len_inc',{{3,[0 0],[0 0 0]}}),'len_com',{{3,[2 2],[2 2 2]}}))));
    V{end+1} = struct('name','V2 psatz none',     'f',@(s) setfield(s,'lpi',setfield(s.lpi,'ineq',setfield(s.lpi.ineq,'psatz',{{'none','faces1','auto'}}))));
    V{end+1} = struct('name','V3 psatz faces',    'f',@(s) setfield(s,'lpi',setfield(s.lpi,'ineq',setfield(s.lpi.ineq,'psatz',{{'faces','faces1','auto'}}))));
    V{end+1} = struct('name','V3b psatz product', 'f',@(s) setfield(s,'lpi',setfield(s.lpi,'ineq',setfield(s.lpi.ineq,'psatz',{{'product','faces1','auto'}}))));
    V{end+1} = struct('name','V10 heavy tier',    'f',@(s) lpisettings_sop('heavy',0,'','mosek'));
else
    V{end+1} = struct('name','V1 mult 2',         'f',@(s) setfield(s,'lpi',setfield(s.lpi,'pos',setfield(s.lpi.pos,'mult',[2 0 0]))));
    V{end+1} = struct('name','V10 heavy tier',    'f',@(s) lpisettings_sop('heavy',0,'','mosek'));
end
fprintf('@@AB git %s\n', info.git);
for v = 1:numel(V)
    for k = 1:size(cases,1)
        c = C(strcmp({C.key},cases{k,1}));   c.exec = ['sop:',c.exec];
        tau = cases{k,2};
        try
            S0 = V{v}.f(lpisettings_sop('light',0,'','mosek'));
            PIE = lb_build_system(c, tau, true);
            r = lb_exec(c.exec, PIE, S0);
            fprintf('@@AB %-20s %-10s tau=%-8g %-10s ndv=%-7d m=%-6d nx=%-7d fr=%6.3f pinf=%d dinf=%d numerr=%d %.1fs\n', ...
                V{v}.name, cases{k,1}, tau, r.status, r.ndec, r.m, r.nx, r.fr, r.pinf, r.dinf, r.numerr, r.sec);
        catch e
            fprintf('@@AB %-20s %-10s tau=%-8g ERROR %s: %s\n', V{v}.name, cases{k,1}, tau, e.identifier, strtok(e.message,newline));
        end
    end
end
fprintf('@@AB done %s\n', char(datetime('now')));
end
