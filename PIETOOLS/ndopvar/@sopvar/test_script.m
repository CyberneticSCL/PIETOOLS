clc; clear;
nvars_total = 3;
nvars_in = 2;
nvars_out = 2;
vars_total = cellstr("s"+string(1:nvars_total));
vars_in = cellstr("s"+string(sort(randperm(nvars_total, nvars_in))));
vars_out = cellstr("s"+string(sort(randperm(nvars_total,nvars_out))));

vars_dum_total = cellstr(vars_total+"_dum");
vars_dum_in = cellstr(vars_in+"_dum");
vars_dum_out = cellstr(vars_out+"_dum");

pvars_total = polynomial(speye(nvars_total),speye(nvars_total),vars_total,[nvars_total,1]);
pvars_in = polynomial(speye(nvars_in),speye(nvars_in),vars_in,[nvars_in,1]);
pvars_out = polynomial(speye(nvars_out),speye(nvars_out),vars_out,[nvars_out,1]);

pvars_dum_total = polynomial(speye(nvars_total),speye(nvars_total),vars_dum_total,[nvars_total,1]);
pvars_dum_in = polynomial(speye(nvars_in),speye(nvars_in),vars_dum_in,[nvars_in,1]);
pvars_dum_out = polynomial(speye(nvars_out),speye(nvars_out),vars_dum_out,[nvars_out,1]);


s3 = intersect(vars_in,vars_out);
s1 = setdiff(vars_in,s3);
s2 = setdiff(vars_out,s3);
%%
dim = [1,1];
P = sopvar.randsopvar(s1,s2,s3,dim,1,0.05);
Q = sopvar.randsopvar(s1,s2,s3,dim,1,0.05);
%% Addition
R = P+Q;

%% Transpose
Pt = P';

%% mtimes
s4 = {'s4'};
P = sopvar.randsopvar(s2,s1,s3,dim,1,0.05);
Q = sopvar.randsopvar(s4,s2,s3,dim,1,0.05);
S = P*Q;
