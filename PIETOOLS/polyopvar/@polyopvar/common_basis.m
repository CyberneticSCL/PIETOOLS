function [A_out,B_out] = common_basis(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_BASIS(A_IN,B_IN) takes two distributed polynomials
% and expresses them in terms of the same basis of distributed monomials
%
% INPUTS
% - A_in, B_in:     m x n 'polyopvar' objects
%
% OUTPUTS
% - A_out, B_out:   m x n 'polyopvar' objects representing the same
%                   distributed polynomials as A_in, B_in, respectively,
%                   but now expressed in terms of the same basis of
%                   monomials, so that
%                   A_out.varname = B_out.varname;
%                   A_out.degmat = B_out.degmat;
%



% Extract the parameters representing each polynomial
var_A = A_in.varname;       var_B = B_in.varname;
pvar_A = A_in.pvarname;     pvar_B = B_in.pvarname;
dom_A = A_in.dom;           dom_B = B_in.dom;
degs_A = A_in.degmat;       degs_B = B_in.degmat;
vmat_A = A_in.varmat;       vmat_B = B_in.varmat;
C_A = A_in.C;               C_B = B_in.C;

% First, combine the list of state variables
p1 = numel(var_A);      p2 = numel(var_B);
var_full = [var_A;var_B];
[var_new,~,new2old_vdcs] = unique(var_full);   % var_new = var_full(old2new_vdcs);
p = numel(var_full);

% Also combine the list of spatial variables
N1 = numel(pvar_A);     N2 = numel(pvar_B);
pvar_full = [pvar_A; pvar_B];
[pvar_new,old2new_pdcs,new2old_pdcs] = unique(pvar_full);   % pvar_full = pvar_new(new2old_pdcs);
dom_full = [dom_A; dom_B];
dom_new = dom_full(old2new_pdcs,:);
N = numel(pvar_full);

% Build a varmat indicating for each of the new state variables on which
% spatial variable it depends.
varmat = logical(zeros(p,N));
varmat(new2old_vdcs(1:p1),new2old_pdcs(1:N1)) = vmat_A;
varmat(new2old_vdcs(p1+(1:p2)),new2old_pdcs(N1+(1:N2))) = vmat_B;

% Augment the monomial basis to incorporate the full list of variables
degs1 = zeros(size(degs_A,1),p);
degs2 = zeros(size(degs_B,1),p);
degs1(:,new2old_vdcs(1:p1)) = degs_A;
degs2(:,new2old_vdcs(p1+(1:p2))) = degs_B;

% Combine the monomial bases into a unique basis
nZ1 = size(degs1,1);
degs_full = [degs1;degs2];
[Pmat,degs_new] = uniquerows_integerTable(degs_full);   % Pmat*Znew = Z_full
nZ = size(degs_new,1);
[r_idcs,new2old_Cdcs] = find(Pmat);
C1 = cell(1,nZ);            C2 = cell(1,nZ);
is_C1 = r_idcs<=nZ1;        is_C2 = r_idcs>nZ1;
C1(new2old_Cdcs(is_C1)) = C_A(r_idcs(is_C1));
C2(new2old_Cdcs(is_C2)) = C_B(r_idcs(is_C2)-nZ1);

% Build polynomials in terms of the shared basis
A_out = polyopvar();            B_out = polyopvar();
A_out.C = C1;                   B_out.C = C2;
A_out.degmat = degs_new;        B_out.degmat = degs_new;
A_out.varname = var_new;        B_out.varname = var_new;
A_out.pvarname = pvar_new;      B_out.pvarname = pvar_new;
A_out.dom = dom_new;            B_out.dom = dom_new;
A_out.varmat = varmat;          B_out.varmat = varmat;

end