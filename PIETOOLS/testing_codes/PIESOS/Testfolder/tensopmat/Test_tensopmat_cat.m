% Test_tensopmat_cat.m -- tests blkdiag/vertcat/horzcat for 'tensopmat'
% objects via apply_op, using two operators built with entirely
% independent (disjoint) variable pools: concatenation places no
% requirement on shared spaces for blkdiag, and vertcat/horzcat instead
% build B to share A's exact input/output side, respectively.
%
% Oracle: apply_op of the concatenation vs. the stacked/summed apply_op of
% the parts. Sweeps otype(1), nterms and degree per block, p/d/N/M down to
% 1x1, and a sparser 'dnsty' on some reps; 20 fresh random reps.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - blkdiag/vertcat/horzcat only place existing blocks into a bigger cell
%   array and never touch '.otype', so mixing otype(1) here carries none
%   of the apply_op output renaming issue described in rand_tensopmat.m.

tol = 1e-8;
nreps = 20;

for rep_idx=1:nreps

pA = randi([1,3]);  dA = randi([1,3]);
NA = randi([1,3]);  MA = randi([1,3]);
pB = randi([1,3]);  dB = randi([1,3]);
NB = randi([1,3]);  MB = randi([1,3]);
dom = [0,1];
if rand()<0.3
    dnsty = 0.3+0.4*rand();
else
    dnsty = [];
end

[Aop,infoA] = rand_tensopmat(pA,dA,NA,MA,dom,0.7,dnsty,struct(),'a','b','rand',3,3);
[Bop,infoB] = rand_tensopmat(pB,dB,NB,MB,dom,0.7,dnsty,struct(),'c','e','rand',3,3);

XA = cell(1,dA);
for j=1:dA
    XA{j} = rand_poly([infoA.dims_in(j),1],infoA.vars_in_pool(infoA.sel_in{j}),1,5);
end
XB = cell(1,dB);
for j=1:dB
    XB{j} = rand_poly([infoB.dims_in(j),1],infoB.vars_in_pool(infoB.sel_in{j}),1,5);
end
YA = apply_op(Aop,XA);
YB = apply_op(Bop,XB);

% blkdiag: no shared space required; off-diagonal blocks are zero
Cdiag = blkdiag(Aop,Bop);
Xdiag = [XA,XB];
Ydiag = apply_op(Cdiag,Xdiag);
err_diag = cleanpoly(Ydiag-[YA;YB],1e-10);
if max(max(abs(err_diag.C)))>tol
    error("Test_tensopmat_cat failed: blkdiag(A,B) mismatch.")
end

% vertcat: B2 reuses A's exact input side
infoB2_in = struct('vars_in_pool',{infoA.vars_in_pool},'dom_in_pool',infoA.dom_in_pool,...
                    'sel_in',{infoA.sel_in},'dims_in',infoA.dims_in);
[Bop2,infoB2] = rand_tensopmat(pB,[],[],MB,dom,0.7,dnsty,infoB2_in,[],'f','rand',3,3);

Cvert = [Aop;Bop2];
YB2 = apply_op(Bop2,XA);
Yvert = apply_op(Cvert,XA);
err_vert = cleanpoly(Yvert-[YA;YB2],1e-10);
if max(max(abs(err_vert.C)))>tol
    error("Test_tensopmat_cat failed: vertcat(A,B2) mismatch.")
end

% horzcat: B3 reuses A's exact output side
infoB3_out = struct('vars_out_pool',{infoA.vars_out_pool},'dom_out_pool',infoA.dom_out_pool,...
                     'sel_out',{infoA.sel_out},'dims_out',infoA.dims_out);
[Bop3,infoB3] = rand_tensopmat([],dB,NB,[],dom,0.7,dnsty,infoB3_out,'g',[],'rand',3,3);

XB3 = cell(1,dB);
for j=1:dB
    XB3{j} = rand_poly([infoB3.dims_in(j),1],infoB3.vars_in_pool(infoB3.sel_in{j}),1,5);
end
Chorz = [Aop,Bop3];
YB3 = apply_op(Bop3,XB3);
Yhorz = apply_op(Chorz,[XA,XB3]);
err_horz = cleanpoly(Yhorz-(YA+YB3),1e-10);
if max(max(abs(err_horz.C)))>tol
    error("Test_tensopmat_cat failed: horzcat(A,B3) mismatch.")
end

end

disp('Test_tensopmat_cat passed.')
