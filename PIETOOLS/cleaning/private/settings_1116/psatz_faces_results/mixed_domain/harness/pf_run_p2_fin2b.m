% PF_RUN_P2_FIN2B  Phase 2, the rest of PF_RUN_P2_FIN2 after its fin2 run
% under Rheavy (2-D entry all, d1 2, d2 3, d3 [2 2;3 3;3 3]) was killed at
% 18.8 GB and 30 min of MOSEK (07:45; D46): Rheavy is dropped as a tier
% candidate. Light: fin2 and rd2 at 0.9 under Rlight2 (stock, Rdef, Rlight
% recorded in p2fin2.log). Heavy tier: stock, Rlight, Rlight2. Then the
% 3-D heat benchmark at 0.9 kappa* under stock and Rlight (light). Output
% appended to p2fin2.log with its end marker. SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
pvar s1 s2 t
x = pde_var(1,[s1;s2],[0,1;0,1]);  w = pde_var('in',1);  z = pde_var('out',1);
PDE = [diff(x,t)==diff(x,s1,2)+diff(x,s2,2)+0.5*2*pi^2*x+w; z==int(x,[s1;s2],[0,1;0,1]);
       subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
evalc('PIEg = initialize(convert(PDE,{''pie'',''silent''}));');
clear stateNameGenerator
x9 = pde_var(1,[s1;s2],[0,1;0,1]);
PDE9 = [diff(x9,t)==diff(x9,s1,2)+diff(x9,s2,2)+0.9*2*pi^2*x9;
        subs(x9,s1,0)==0; subs(x9,s1,1)==0; subs(x9,s2,0)==0; subs(x9,s2,1)==0];
evalc('PIE9 = initialize(convert(PDE9,{''pie'',''silent''}));');
out = sprintf('%s/results/p2_fin2b.mat',S);   Out = [];
L = {'light',{'Rlight2'}; 'heavy',{'stock','Rlight','Rlight2'}};
for i = 1:size(L,1)
    tier = L{i,1};
    for sto = L{i,2}
        r = pf_exec_row('fin2',PIEg,'hinfP',sto{1},tier);   r.msg = sprintf('gap %+.3e',r.gam/0.068606888-1);
        fprintf('@@ EXGAP fin2 storage=%s tier=%s gap=%+.3e\n',sto{1},tier,r.gam/0.068606888-1);
        if isempty(Out),    Out = r;  else,   Out(end+1) = r;   end         %#ok<SAGROW>
        r = pf_exec_row('rd2_09',PIE9,'stabP',sto{1},tier);
        Out(end+1) = r;                                                     %#ok<SAGROW>
        save(out,'Out');
    end
end
% 3-D: the heat benchmark at 0.9 kappa* (DD x DN x DN), stability 'P', light
PIE3 = [];
try
    pvar s3
    x3 = pde_var(1,[s1;s2;s3],[0,1;0,1;0,1]);
    PDE3 = [diff(x3,t)==diff(x3,s1,2)+diff(x3,s2,2)+diff(x3,s3,2)+0.9*1.5*pi^2*x3;
            subs(x3,s1,0)==0; subs(x3,s1,1)==0; subs(x3,s2,0)==0; subs(diff(x3,s2),s2,1)==0; subs(x3,s3,0)==0; subs(diff(x3,s3),s3,1)==0];
    evalc('PIE3 = initialize(convert(PDE3,{''pie'',''silent''}));');
catch ME
    fprintf('@@ EX heat3d build ERROR %s\n',strrep(ME.message,newline,' '));
end
if ~isempty(PIE3)
    for sto = {'stock','Rlight'}
        r = pf_exec_row('heat3d_09',PIE3,'stabP',sto{1},'light');
        Out(end+1) = r;                                                     %#ok<SAGROW>
        save(out,'Out');
    end
end
disp('## P2FIN2 DONE');
