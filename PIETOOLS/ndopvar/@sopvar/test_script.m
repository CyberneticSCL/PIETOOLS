clc; clear;
nvars_total = 1;
nvars_in = 1;
nvars_out = 1;
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

%%
pvar s1 s1_dum t1;
P1 = rand_opvar([0,0;1,1],1,pvar('s1'),pvar('s1_dum'),[0,1]);
P2 = rand_opvar([0,0;1,1],1,pvar('s1'),pvar('s1_dum'),[0,1]);

R01quad = quadPoly.polynomial2quadPoly( ...
    subs(P1.R.R0,s1_dum,t1),{'s1'},{'t1'});
R11quad = quadPoly.polynomial2quadPoly( ...
    subs(P1.R.R1,s1_dum,t1),{'s1'},{'t1'});
R21quad = quadPoly.polynomial2quadPoly( ...
    subs(P1.R.R2,s1_dum,t1),{'s1'},{'t1'});
R02quad = quadPoly.polynomial2quadPoly( ...
    subs(P2.R.R0,s1_dum,t1),{'s1'},{'t1'});
R12quad = quadPoly.polynomial2quadPoly( ...
    subs(P2.R.R1,s1_dum,t1),{'s1'},{'t1'});
R22quad = quadPoly.polynomial2quadPoly( ...
    subs(P2.R.R2,s1_dum,t1),{'s1'},{'t1'});

P1sopvar = P; P2sopvar = Q;
P1sopvar.params = {R01quad,R11quad,R12quad};
P2sopvar.params = {R02quad,R12quad,R22quad};

tic;
Q1 = P1*P2;
tmulopvar = toc;
tic;
Q1sopvar = P1sopvar*P2sopvar;
tmulsopvar = toc;

Q1opvar = sopvar2opvar(Q1sopvar);

res=0;
for i=1:3
    sopvarParam = Q1sopvar.params{i};
    polyparam = quadPoly.quadPoly2polynomial(sopvarParam);
    opvarParam = Q1.R.(['R',num2str(i-1)]);
    diff = polyparam-opvarParam;
    res = max(res,max(diff.C(:)));
end