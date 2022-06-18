

pvar s theta

prog = sosprogram([s,theta]);

n = [5;10]; %size of problem

d = {4,[5,6,7],[5,6,7]};

I = [0,1];

options.psatz = 0;

%disp('Running polynomial poslpivar')
%tic
%pname = 'Profile_poslpivar';
%profile on
%[prog1,P1] = poslpivar(prog,n,I,d,options);
%profile off
%profsave(profile('info'),pname);
%toc

disp('Running DP poslpivar')
tic
%pname = 'Profile_DPposlpivar';
%profile on
[prog2,P2] = DPposlpivar(prog,n,I,d,options);
%profile off
%profsave(profile('info'),pname);
toc

disp('Running DP poslpivar without subsref')
tic
%pname = 'Profile_DPposlpivar_nosubsref';
%profile on
[prog3,P3] = DPposlpivar_nosubsref(prog,n,I,d,options);
%profile off
%profsave(profile('info'),pname);
toc