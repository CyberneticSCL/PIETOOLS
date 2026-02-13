function [A_out,B_out] = common_vars(A_in,B_in)
% [A_OUT,B_OUT] = COMMON_VARS(A_IN,B_IN) takes two distributed polynomials
% and expresses them in terms of the same state variables
%
% INPUTS
% - A_in, B_in:     m x n 'polyopvar' objects
%
% OUTPUTS
% - A_out, B_out:   m x n 'polyopvar' objects representing the same
%                   distributed polynomials as A_in, B_in, respectively,
%                   but now expressed in terms of the same variables, so
%                   that
%                   A_out.varname = B_out.varname;
%                   A_out.pvarname = B_out.pvarname;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - common_vars
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 01/22/2026: Initial coding


% Extract the parameters representing each polynomial
var_A = A_in.varname;       var_B = B_in.varname;
pvar_A = A_in.pvarname;     pvar_B = B_in.pvarname;
dom_A = A_in.dom;           dom_B = B_in.dom;
degs_A = A_in.degmat;       degs_B = B_in.degmat;
vmat_A = A_in.varmat;       vmat_B = B_in.varmat;
C_A = A_in.C;               C_B = B_in.C;

% First, combine the list of state variables
p1 = numel(var_A);      p2 = numel(var_B);
if isequal(var_A,var_B)
    var_new = var_A;
    new2old_vdcs = [1:p1,1:p1]';
else
    var_full = [var_A;var_B];
    [var_new,~,new2old_vdcs] = unique(var_full);   % var_full = var_new(new2old_vdcs);
end
p = numel(var_new);

% Also combine the list of spatial variables
N1 = numel(pvar_A);     N2 = numel(pvar_B);
pvar_full = [pvar_A; pvar_B];
[pvar_new,old2new_pdcs,new2old_pdcs] = unique(pvar_full);   % pvar_full = pvar_new(new2old_pdcs);
dom_full = [dom_A; dom_B];
dom_new = dom_full(old2new_pdcs,:);
N = numel(pvar_new);

% Build a varmat indicating for each of the new state variables on which
% spatial variable it depends.
varmat = false(p,N);
varmat(new2old_vdcs(1:p1),new2old_pdcs(1:N1)) = vmat_A;
varmat(new2old_vdcs(p1+(1:p2)),new2old_pdcs(N1+(1:N2))) = vmat_B;

% Reorder tensor products of PI operators in C to match new order of vars
if isempty(degs_A)
    d1 = 0;
else
    d1 = max(max(sum(degs_A,2)));
end
if isempty(degs_B)
    d2 = 0;
else
    d2 = max(max(sum(degs_B,2)));
end
dumvar_cell = cell(1,max(d1,d2));
for lidx=1:numel(C_A.ops)
    [~,cidx] = ind2sub(size(C_A.ops),lidx);
    deg_l = degs_A(cidx,:);
    if sum(deg_l)==1 || isempty(C_A.ops{lidx})
        continue
    end
    % Establish which elements of C1.ops{lidx} correspond to each variable
    deg_idcs = mat2cell([[0,cumsum(deg_l(1:end-1))];cumsum(deg_l)],2,ones(1,numel(deg_l)));
    deg_idcs = cellfun(@(a)(a(1)+1:a(2)),deg_idcs,'UniformOutput',false);
    old_deg_idcs = cell(1,p);
    old_deg_idcs(new2old_vdcs(1:p1)) = deg_idcs;
    new_deg_idcs = cell2mat(old_deg_idcs);

    % Reorder to match the new order of the variables
    if isa(C_A.ops{lidx},'intvar')
        Cnew_l = C_A.ops{lidx};
        Cnew_l.pvarname = C_A.ops{lidx}.pvarname(new_deg_idcs);
        % Reorder columns of omat to match new ordering of variables
        for jj=1:numel(new_deg_idcs)
            Cnew_l.omat(C_A.ops{lidx}.omat==new_deg_idcs(jj)) = jj;
        end
        % Use the same dummy variables for monomials of same degree
        if isempty(dumvar_cell{sum(deg_l)})
            dumvar_cell{sum(deg_l)} = Cnew_l.pvarname;
        else
            Cparams = Cnew_l.params;
            varname_new = cell(size(Cparams.varname));
            for kk=1:numel(Cparams.varname)
                old_var = Cparams.varname(kk);
                var_idx = strcmp(Cnew_l.pvarname,old_var);
                varname_new(kk) = dumvar_cell{sum(deg_l)}(var_idx);
            end
            Cparams.varname = varname_new;
            Cnew_l.params = Cparams;
            Cnew_l.pvarname = dumvar_cell{sum(deg_l)};
        end        
    elseif isa(C_A.ops{lidx},'struct') && isfield(C_A.ops{lidx},'params')
        Cnew_l = C_A.ops{lidx};
        Cnew_l.vars = C_A.ops{lidx}.vars(new_deg_idcs);
        % Reorder columns of omat to match new ordering of variables
        for jj=1:numel(new_deg_idcs)
            Cnew_l.omat(C_A.ops{lidx}.omat==new_deg_idcs(jj)) = jj;
        end
        % Use the same dummy variables for monomials of same degree
        if isempty(dumvar_cell{sum(deg_l)})
            dumvar_cell{sum(deg_l)} = cell(size(Cnew_l.vars));
            for kk=1:numel(Cnew_l.vars)
                dumvar_cell{sum(deg_l)}(kk) = Cnew_l.vars(kk).varname;
            end
        else
            Cparams = Cnew_l.params;
            varname_new = cell(size(Cparams.varname));
            for kk=1:numel(Cparams.varname)
                old_var = polynomial(Cparams.varname(kk));
                var_idx = isequal(Cnew_l.vars,old_var);
                varname_new(kk) = dumvar_cell{sum(deg_l)}(var_idx);
            end
            Cparams.varname = varname_new;
            Cnew_l.params = Cparams;
            Cnew_l.vars = polynomial(dumvar_cell{sum(deg_l)});
        end
    else
        Cnew_l = C_A.ops{lidx}(:,new_deg_idcs);
    end
    C_A.ops{lidx} = Cnew_l;
end
for lidx=1:numel(C_B.ops)
    [~,cidx] = ind2sub(size(C_B.ops),lidx);
    deg_l = degs_B(cidx,:);
    if sum(deg_l)==1 || isempty(C_B.ops{lidx})
        continue
    end
    % Establish which elements of C1.ops{lidx} correspond to each variable
    deg_idcs = mat2cell([[0,cumsum(deg_l(1:end-1))];cumsum(deg_l)],2,ones(1,numel(deg_l)));
    deg_idcs = cellfun(@(a)(a(1)+1:a(2)),deg_idcs,'UniformOutput',false);
    old_deg_idcs = cell(1,p);
    old_deg_idcs(new2old_vdcs(p1+(1:p2))) = deg_idcs;
    new_deg_idcs = cell2mat(old_deg_idcs);
    
    % Reorder to match the new order of the variables
    if isa(C_B.ops{lidx},'intvar')
        Cnew_l = C_B.ops{lidx};
        Cnew_l.pvarname = C_B.ops{lidx}.pvarname(new_deg_idcs);
        % Reorder columns of omat to match new ordering of variables
        for jj=1:numel(new_deg_idcs)
            Cnew_l.omat(C_B.ops{lidx}.omat==new_deg_idcs(jj)) = jj;
        end
        % Use the same dummy variables for monomials of same degree
        if isempty(dumvar_cell{sum(deg_l)})
            dumvar_cell{sum(deg_l)} = Cnew_l.pvarname;
        else
            Cparams = Cnew_l.params;
            varname_new = cell(size(Cparams.varname));
            for kk=1:numel(Cparams.varname)
                old_var = Cparams.varname(kk);
                var_idx = strcmp(Cnew_l.pvarname,old_var);
                varname_new(kk) = dumvar_cell{sum(deg_l)}(var_idx);
            end
            Cparams.varname = varname_new;
            Cnew_l.params = Cparams;
            Cnew_l.pvarname = dumvar_cell{sum(deg_l)};
        end 
    elseif isa(C_B.ops{lidx},'struct') && isfield(C_B.ops{lidx},'params')
        Cnew_l = C_B.ops{lidx};
        Cnew_l.vars = C_B.ops{lidx}.vars(new_deg_idcs);
        % Reorder columns of omat to match new ordering of variables
        for jj=1:numel(new_deg_idcs)
            Cnew_l.omat(C_B.ops{lidx}.omat==new_deg_idcs(jj)) = jj;
        end
        % Use the same dummy variables for monomials of same degree
        if isempty(dumvar_cell{sum(deg_l)})
            dumvar_cell{sum(deg_l)} = cell(size(Cnew_l.vars));
            for kk=1:numel(Cnew_l.vars)
                dumvar_cell{sum(deg_l)}(kk) = Cnew_l.vars(kk).varname;
            end
        else
            Cparams = Cnew_l.params;
            varname_new = cell(size(Cparams.varname));
            for kk=1:numel(Cparams.varname)
                old_var = polynomial(Cparams.varname(kk));
                var_idx = isequal(Cnew_l.vars,old_var);
                varname_new(kk) = dumvar_cell{sum(deg_l)}(var_idx);
            end
            Cparams.varname = varname_new;
            Cnew_l.params = Cparams;
            Cnew_l.vars = polynomial(dumvar_cell{sum(deg_l)});
        end
    else
        Cnew_l = C_B.ops{lidx}(:,new_deg_idcs);
    end
    C_B.ops{lidx} = Cnew_l;
end

% Augment the monomial basis to incorporate the full list of variables
degs1 = zeros(size(degs_A,1),p);
degs2 = zeros(size(degs_B,1),p);
degs1(:,new2old_vdcs(1:p1)) = degs_A;
degs2(:,new2old_vdcs(p1+(1:p2))) = degs_B;

% Build polynomials in terms of the shared variable names
A_out = polyopvar();            B_out = polyopvar();
A_out.C = C_A;                  B_out.C = C_B;
A_out.degmat = degs1;           B_out.degmat = degs2;
A_out.varname = var_new;        B_out.varname = var_new;
A_out.pvarname = pvar_new;      B_out.pvarname = pvar_new;
A_out.dom = dom_new;            B_out.dom = dom_new;
A_out.varmat = varmat;          B_out.varmat = varmat;

end