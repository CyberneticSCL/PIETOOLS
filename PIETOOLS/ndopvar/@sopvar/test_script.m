clc; clear;
nvars_total = 5;
nvars_in = 3;
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

%%
dim = [1,2];
P = sopvar.randOpvar(vars_in,vars_out,dim,1,0.05);
Q = sopvar.randOpvar(vars_in,vars_out,dim,1,0.05);
%% Addition
R = P+Q;

%% Transpose
Pt = P';
