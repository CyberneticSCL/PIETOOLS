clear
pvar s1 s1_dum
var1 = s1;
var2 = s1_dum;
dom = [0,1];
deg = 2;        % maximal monomial degree in independent variables
m = 2;      n = 2;
nvars = 1;
deg1 = randi(2,[1,nvars])-1;
if ~any(deg1)
    deg1(randi(nvars))=1;
end
deg2 = randi(2,[1,nvars])-1;
if ~any(deg2)
    deg2(randi(nvars))=1;
end
d1 = sum(deg1);
d2 = sum(deg2);
m_idx = randi(d1+d2);   % index for which operator may include a multiplier
%m_idx = 2;

% Declare state variable names
xvarname = mat2cell([repmat('x',[nvars,1]),num2str((1:nvars)')],ones(nvars,1));

% Declare d1 random 2-PI operators
Lops_opvar = cell(1,d1);
Lops = tensopvar();
Lops.ops{1} = cell(1,d1);
for ii=1:d1
    Lop_ii = rand_opvar([0,0;m,1],deg,var1,var2,dom);
    if ii~=m_idx
        % Allow only one operator with a nonzero multiplier
        Lop_ii.R.R0 = zeros([m,1]);
    end
    Lops_opvar{ii} = Lop_ii;
    if d1==1
        Lops.ops{1} = dopvar2ndopvar(Lop_ii);
    else
        Lops.ops{1}{ii} = dopvar2ndopvar(Lop_ii);
    end
end
Lmon = polyopvar();
Lmon.varname = xvarname;
Lmon.pvarname = s1.varname;
Lmon.dom = dom;
Lmon.varmat = ones(nvars,1);
Lmon.degmat = deg1;
Lmon.C = Lops;

Rops_opvar = cell(1,d2);
Rops = tensopvar();
Rops.ops{1} = cell(1,d2);
for ii=1:d2
    Rop_ii = rand_opvar([0,0;n,1],deg,var1,var2,dom);
    if ii+d1~=m_idx
        % Allow only one operator with a nonzero multiplier
        Rop_ii.R.R0 = zeros([n,1]);
    end
    Rops_opvar{ii} = Rop_ii;
    if d2==1
        Rops.ops{1} = dopvar2ndopvar(Rop_ii);
    else
        Rops.ops{1}{ii} = dopvar2ndopvar(Rop_ii);
    end
end
Rmon = Lmon;
Rmon.degmat = deg2;
Rmon.C = Rops;

Pmat = rand([m,n]);

% Convert to linear form
[Kop] = quad2lin_term(Pmat,Lmon,Rmon);
Pvec = Kop.params;
idx_mat = Kop.omat;
Kvar2 = Kop.vars;



% Generate random operators to represent RHS of the PI
nfctrs = [1;1];
Fx = cell(1,size(nfctrs,1));
Fops = cell(1,size(nfctrs,1));
m_idx2 = randi(numel(Fx));      % allow only one operator with multiplier terms
for ii=1:numel(Fx)
    Fx{ii} = Rmon;
    Fx{ii}.degmat = nfctrs(ii,:);
    Fx{ii}.C.ops = cell(1,sum(nfctrs(ii,:)));
    Fops{ii} = cell(1,sum(nfctrs(ii)));
    for jj=1:numel(Fops{ii})
        Fops{ii}{jj} = rand_opvar([0,0;1,1],deg,var1,var2,dom);
        if ~isscalar(Fx{ii}) || ii~=m_idx2
            % Do not allow tensor product of multiplier operators
            Fops{ii}{jj}.R.R0 = 0;
        end
        if sum(nfctrs(ii))==1
            Fx{ii}.C.ops{jj} = dopvar2ndopvar(Fops{ii}{jj});
        else
            Fx{ii}.C.ops{1}{jj} = dopvar2ndopvar(Fops{ii}{jj});
        end
    end
end


% Compute the composition of Kop with the factors Fx
Cop = mtimes_functional(Kop,Fx);



% Generate random functions xi(s)
x_tst = polynomial(zeros(1,nvars));
for ii=1:nvars
    x_tst(ii) = rand_poly([1,1],s1,3);
end
ycell = cell(1,numel(Fx));
zval = polynomial(zeros(1,numel(Fx)));
for ii=1:numel(Fx)
    ycell{ii}=cell(1,numel(Fops{ii}));
    zval(ii) = 1;
    strt_idcs = cumsum(nfctrs(ii,:));
    for jj=1:numel(Fops{ii})
        x_idx = find(jj<=strt_idcs,1,'first');
        ycell{ii}{jj} = apply_opvar(Fops{ii}{jj},x_tst(x_idx));
        zval(ii) = zval(ii).*ycell{ii}{jj};
    end
end

% Apply Kop to the product (Rop{1}*x)*(Rop{2}*x)...
fval_true = apply_functional(Kop.params,zval,ones(1,numel(Fx)),Kop.omat,Kop.vars,Kop.dom);

% Apply Cop directly to x*x
fval_alt = apply_functional(Cop.params,x_tst,sum(nfctrs,1),Cop.omat,Cop.vars,dom);



f_err = norm(double(fval_true-fval_alt))