% INH_RUN_4  Item 4: the direct PIE2PDE test on 2-D system A at 0.9 lam*
% by term set (faces0, facesS at dw 0; box, product2 at dw 1).
% SCRATCH, 11/16/2026.
addpath('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces');
S = pf_setup(false);
D = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/inhom';
addpath(D);
L = load([D '/results/thresholds.mat']);
PIE = inh_sys('A',0.9*L.T.A(2));
[o1,n1] = pf_termsets({'faces0','facesS','box'});
[o2,~] = pf_termsets3({'product2'});
qd1 = {struct('mult',1,'int',[1 1],'out',1,'in',1)};
J = {o1{1},'faces0',0; o1{2},'facesS',0; o1{3},'box',1; o2{1},'product2',1};
for j = 1:size(J,1)
    try,  pf_qdirect_pie('A_0.9',PIE,qd1,J{j,3},sprintf('A_%s_dw%d',J{j,2},J{j,3}),1e-4,0,J{j,1});
    catch ME,  fprintf('@@ QD ERROR A %s: %s\n',J{j,2},strrep(ME.message,newline,' '));  end
end
disp('## ITEM4 DONE');
