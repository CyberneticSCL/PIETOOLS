clc; clear;
nvars_total = 5;
nvars_in = 3;
nvars_out = 2;
vars_total = cellstr("s"+string(1:nvars_total));
vars_in = cellstr("s"+string(sort(randperm(nvars_total, nvars_in))));
vars_out = cellstr("s"+string(sort(randperm(nvars_total,nvars_out))));

vars_dum_total = strrep(vars_total,'s','t');
vars_dum_in = strrep(vars_in,'s','t');
vars_dum_out = strrep(vars_out,'s','t');

pvars_total = polynomial(speye(nvars_total),speye(nvars_total),vars_total,[nvars_total,1]);
pvars_in = polynomial(speye(nvars_in),speye(nvars_in),vars_in,[nvars_in,1]);
pvars_out = polynomial(speye(nvars_out),speye(nvars_out),vars_out,[nvars_out,1]);

pvars_dum_total = polynomial(speye(nvars_total),speye(nvars_total),vars_dum_total,[nvars_total,1]);
pvars_dum_in = polynomial(speye(nvars_in),speye(nvars_in),vars_dum_in,[nvars_in,1]);
pvars_dum_out = polynomial(speye(nvars_out),speye(nvars_out),vars_dum_out,[nvars_out,1]);


s3 = intersect(vars_in,vars_out);
s1 = setdiff(vars_total(randperm(nvars_total, nvars_in)),s3);
s2 = setdiff(vars_out,s3);
s4 = setdiff(vars_total(randperm(nvars_total, nvars_in)),s3);
dim = [1,1];
%%
P = sopvar.randsopvar(s1,s2,s3,dim,1,0.05);
Q = sopvar.randsopvar(s1,s2,s3,dim,1,0.05);
%% Addition
R = P+Q;

%% Transpose
Pt = P';

%% mtimes

P = sopvar.randsopvar(s2,s1,s3,fliplr(dim),1,0.05);
Q = sopvar.randsopvar(s1,s2,s3,dim,1,0.05);
S = P*Q;
