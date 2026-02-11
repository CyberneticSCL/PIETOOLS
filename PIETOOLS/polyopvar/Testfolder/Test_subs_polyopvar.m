clear
pvar s1 s1_dum
var1 = s1;
var2 = s1_dum;
dom = [0,0.5];
pdeg = 2;        % maximal monomial degree in independent variables
m = 3;      n = 3;
nvars = 1;
nvars2 = 2;
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
m_idx = randi(d1);   % index for which operator may include a multiplier
if d1==1 && d2==1
    m_idx2 = d1+randi(d2);   % index for which operator may include a multiplier
else
    m_idx2 = 0;
end
%m_idx = 0;
%m_idx = 2;

% Declare state variable names
xvarname = mat2cell([repmat('x',[nvars,1]),num2str((1:nvars)')],ones(nvars,1));
xvarname2 = mat2cell([repmat('xf',[nvars2,1]),num2str((1:nvars2)')],ones(nvars2,1));
nshare = 1;
share_idx1 = randi(nvars,[1,nshare]);
share_idx1 = unique(share_idx1);
nshare = numel(share_idx1);
share_idx2 = randi(nvars2,[1,nshare]);
share_idx2 = unique(share_idx2);
nshare = numel(share_idx2);
share_idx1 = share_idx1(1:nshare);
xvarname2(share_idx2) = xvarname(share_idx1);


% Declare d1 random 2-PI operators
Lops_opvar = cell(1,d1);
Lops = tensopvar();
Lops.ops{1} = cell(1,d1);
for ii=1:d1
    Lop_ii = rand_opvar([0,0;m,1],pdeg,var1,var2,dom);
    if ii~=m_idx && ii~=m_idx2
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
    Rop_ii = rand_opvar([0,0;n,1],pdeg,var1,var2,dom);
    if ii+d1~=m_idx && ii+d1~=m_idx2
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
Vx = Lmon;
Vx.degmat = Lmon.degmat + Rmon.degmat;
Vx.C.ops = {Kop};



% Generate random operators to represent RHS of the PI
%nfctrs = [3;1];
nfctrs = randi(3,[d1+d2,nvars2])-1;
nfctrs(sum(nfctrs,2)==0,1) = 1;
nterms = randi(2,[d1+d2,1]);
nterms(sum(nfctrs,2)==1) = 1;
Fx = cell(1,size(nfctrs,1));
Fops = cell(1,size(nfctrs,1));
m_idx3 = randi(numel(Fx));      % allow only one operator with multiplier terms
%m_idx2 = 0;
for ii=1:numel(Fx)
    Fx{ii} = Rmon;
    Fx{ii}.varname = xvarname2;
    Fx{ii}.varmat = ones(nvars,1);
    Fx{ii}.degmat = nfctrs(ii,:);
    Fx{ii}.C.ops = {cell(nterms(ii),sum(nfctrs(ii,:)))};
    Fops{ii} = cell(nterms(ii),sum(nfctrs(ii,:)));
    for jj=1:numel(Fops{ii})
        Fops{ii}{jj} = rand_opvar([0,0;1,1],pdeg,var1,var2,dom);
        if ~isscalar(Fops{ii}) || ii~=m_idx3
            % Do not allow tensor product of multiplier operators
            Fops{ii}{jj}.R.R0 = 0;
        end
        if sum(nfctrs(ii,:))==1
            Fx{ii}.C.ops{jj} = dopvar2ndopvar(Fops{ii}{jj});
        else
            Fx{ii}.C.ops{1}{jj} = dopvar2ndopvar(Fops{ii}{jj});
        end
    end
end


% Substitute the polynomial Fx{i} into the ith factor in the monomial
% defining Vx
subs_idx = randi(d1+d2);
Cfun = subs(Vx,subs_idx,Fx{subs_idx});
Cop = Cfun.C.ops{1};
Cdeg = Cfun.degmat;



% Generate random functions xi(s)
x_tst = polynomial(zeros(1,nvars));
for ii=1:nvars
    x_tst(ii) = rand_poly([1,1],s1,3);
end
x_tst2 = polynomial(zeros(1,nvars2));
for ii=1:nvars2
    x_tst2(ii) = rand_poly([1,1],s1,2);
end
x_tst2(share_idx2) = x_tst(share_idx1);
ycell = cell(1,numel(Fx));
zval = polynomial(zeros(1,numel(Fx)));
for ii=1:numel(Fx)
    if ~ismember(ii,subs_idx)
        zval(:,ii) = x_tst(find(ii<=cumsum(deg1+deg2),1,'first'));
        continue
    end
    ycell{ii}=cell(size(Fops{ii}));
    strt_idcs = cumsum(nfctrs(ii,:));
    for ll=1:size(Fops{ii},1)
        zval_ll = 1;
        for jj=1:size(Fops{ii},2)
            x_idx = find(jj<=strt_idcs,1,'first');
            % if ismember(x_idx,share_idx2)
            %     [~,idx] = find(x_idx==share_idx2,1);
            %     x_idx = share_idx1(idx);
            %     ycell{ii}{ll,jj} = apply_opvar(Fops{ii}{ll,jj},x_tst(x_idx));
            % else
                ycell{ii}{ll,jj} = apply_opvar(Fops{ii}{ll,jj},x_tst2(x_idx));
            % end
            zval_ll = zval_ll.*ycell{ii}{ll,jj};
        end
        zval(1,ii) = zval(1,ii) + zval_ll;
    end
end
x_tstf = [x_tst,x_tst2(setdiff(1:nvars2,share_idx2))];

% Apply Kop to the product (Rop{1}*x)*(Rop{2}*x)...
fval_true = apply_functional(Kop,zval,ones(1,numel(Fx)));

% Apply Cop directly to x*x
fval_alt = apply_functional(Cop,x_tstf,Cdeg);



f_err = norm(double(fval_true-fval_alt))


% int_{a}^{b} int_{a}^{t1} int_{a}^{t2} K(t1,t2,t3)*y(t1)*x(t2)*x(t3) dt3 dt2 dt1
% = int_{a}^{b} int_{a}^{s3} int_{a}^{s1} K(s3,s1,s2)*y(s3)*x(s1)*x(s2) ds2 ds1 ds3
% a <= t3 <= t2 <= t1 <= b
% a <= s2 <= s1 <= s3 <= b