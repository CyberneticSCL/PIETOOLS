function ab_2d(tree)
% AB_2D(TREE) PDE 39 (2D_NS_Antonelli) under the light tier, one solve per
% variant (no raise). Library pass 11/17/2026: HEAD certified at 933,696
% variables, the working tree returned unknown (feasratio 7.14) at 648,900.
% 'new': the working tree, varying the 2-D entries of lpi.var.deg (1, the
% light value; 2, the heavy value) and of the term set ('faces1', the
% light value; 'auto'); 'head': the HEAD tier as is. The storage knobs are
% not varied: the executive's 'Q' form has no storage on the working tree.
here = pwd;
info = lb_setup();  cd(here);  addpath(fullfile(here,'gate'));
setappdata(0,'LB_NDV_CAP',3e6);   setappdata(0,'LB_M_CAP',60000);
S = load(fullfile(here,'catalog.mat'));  C = S.C;
c = C(strcmp({C.key},'PDE:39'));   c.exec = ['sop:',c.exec];
noraise = @(s) setfield(s,'sop',setfield(s.sop,'max_raise_feas',[2 0 0]));
V = {};
if strcmp(tree,'new')
    terms  = @(s,t) setfield(s,'lpi',setfield(s.lpi,'ineq',setfield(s.lpi.ineq,'psatz',{'auto',t,'auto'})));
    vardeg = @(s,d) setfield(s,'lpi',setfield(s.lpi,'var',struct('deg',{{[2 2 2],d,1}})));
    V{end+1} = struct('name','W5 vardeg2 faces1',  'f',@(s) noraise(vardeg(s,2)));
    V{end+1} = struct('name','W2 vardeg1 auto',    'f',@(s) noraise(terms(s,'auto')));
    V{end+1} = struct('name','W6 vardeg2 auto',    'f',@(s) noraise(terms(vardeg(s,2),'auto')));
    V{end+1} = struct('name','W0 vardeg1 faces1',  'f',@(s) noraise(s));
else
    V{end+1} = struct('name','H0 head no raise',   'f',@(s) noraise(s));
end
fprintf('@@AB git %s\n', info.git);
PIE = lb_build_system(c, [], true);
for v = 1:numel(V)
    try
        S0 = V{v}.f(lpisettings_sop('light',0,'','mosek'));
        r = lb_exec(c.exec, PIE, S0);
        fprintf('@@AB %-20s %-8s %-10s ndv=%-8d m=%-6d nx=%-8d fr=%7.3f pinf=%d dinf=%d numerr=%d iter=%d %.1fs\n', ...
            V{v}.name, c.key, r.status, r.ndec, r.m, r.nx, r.fr, r.pinf, r.dinf, r.numerr, r.iter, r.sec);
    catch e
        fprintf('@@AB %-20s %-8s ERROR %s: %s\n', V{v}.name, c.key, e.identifier, strtok(e.message,newline));
    end
end
fprintf('@@AB done %s\n', char(datetime('now')));
end
