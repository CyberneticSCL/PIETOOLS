% INH_RUN_2  Item 2: the variable-coefficient S3 cascade at 0.9 of both
% thresholds; stabP light (stock, Rlight0), then the direct PIE2PDE test by
% term set. SCRATCH, 11/16/2026.
addpath('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces');
S = pf_setup(false);
D = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/inhom';
addpath(D);
L = load([D '/results/thresholds.mat']);
r1 = 0.9*L.T.M1(2);   r2 = 0.9*L.T.M2(2);     % N = 60 (the N = 200 eigs run for M1 did not converge)
fprintf('@@ M r1=%.6f r2=%.6f\n',r1,r2);
PIE = inh_sys('M',[r1 r2]);
Out = [];
for sto = {'stock','Rlight0'}
    r = pf_exec_row('M_0.9',PIE,'stabP',sto{1},'light');
    if isempty(Out),  Out = r;  else,  Out(end+1) = r;  end    %#ok<AGROW>
    save([D '/results/item2_exec.mat'],'Out');
end
[o1,n1] = pf_termsets({'faces0','facesS','box'});
[o2,n2] = pf_termsets3({'product2'});
qd1 = {struct('mult',1,'int',[1 1],'out',1,'in',1)};
for j = 1:numel(o1)
    try,  pf_qdirect_pie('M_0.9',PIE,qd1,0,['M_' n1{j} '_dw0'],1e-4,0,o1{j});
    catch ME,  fprintf('@@ QD ERROR M %s: %s\n',n1{j},strrep(ME.message,newline,' '));  end
end
DW1 = {o1{3},'box'; o2{1},'product2'};
for j = 1:2
    try,  pf_qdirect_pie('M_0.9',PIE,qd1,1,['M_' DW1{j,2} '_dw1'],1e-4,0,DW1{j,1});
    catch ME,  fprintf('@@ QD ERROR M %s: %s\n',DW1{j,2},strrep(ME.message,newline,' '));  end
end
disp('## ITEM2 DONE');
