function [d] = degbalance(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [d] = degbalance(P,opts) provides degrees to be selected for a DPposlpivar Q such that
% components of P and Q have similar degrees
% 
% INPUT
%   P: PI dopvar2d class object
%   opts: options describing the type of poslpivar (pure, diag)
% 
% OUTPUT 
%   d: vector of degrees that is the input to sos_posopvar
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - degbalance
%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 04_12_2022

% Process the inputs
if nargin==1
    Zexclude = zeros(1,16);
elseif nargin==2
    if isfield(opts,'exclude')
        Zexclude = opts.exclude;
    else
        Zexclude = zeros(1,16);
    end
end

% Extract the dimensions of the operator P
n = P.dim(:,1);
if any(P.dim(:,1)~=P.dim(:,2))
    warning('The input operator is not square; the returned degrees may not be useful for imposing LPI (in)equality constraints')
end
% If the operator does not map to/from certain states, there is no need to
% implement associated monomials
if n(1)==0
    Zexclude(1)=1;
end
if n(2)==0
    Zexclude(2:4)=1;
end
if n(3)==0
    Zexclude(5:7)=1;
end
if n(4)==0
    Zexclude(8:16)=1;
end

% Extract the maximal degrees of the variables appearing in each monomial
% in the operator
[Pdegs,Pmaxdegs,~] = getdeg(P);


% % % We first impose degrees that retrieve all monomials on diagonal of
% % % the operator (Rxx, Ryy, R22)

% % Start with multiplier terms. Degrees of monomials should be half
% % those observed in the operator. Notice that the maximal degrees
% % observed in the operator multiplier terms must be even for the operator
% % parameter to be SOS.
if ~Zexclude(2)
    dx_o = floor(Pmaxdegs.Rxx{1}(2,1)./2);   % flooring the degree should not change the value
else
    dx_o = -Inf;    % Impose an artificial degree to indicate this monomial does not appear
end
if ~Zexclude(5)
    dy_o = floor(Pmaxdegs.Rxx{1}(1,2)./2);
else
    dy_o = -Inf;    % Impose an artificial degree to indicate this monomial does not appear
end
if ~Zexclude(8)
    d2_oo = floor(Pmaxdegs.R22{1,1}(1:2,1:2)./2);
else
    d2_oo = -Inf*ones(2,2);   % Impose an artificial degree to indicate this monomial does not appear
end
% NOTE: The monomials associated to dx_o, dy_o and d2_oo also contribute to
% other (off-diagonal) parameters in the opvar object. However, if the
% operator is SOS, it is impossible for these other parameters to require
% greater degrees than those suggested by the diagonal multipliers. 
% (That would be equivalent to a 0 value on the diagonal of a PSD matrix
%  --> the associated row and column would be zero as well, so the monomial
%  does not contribute)

% % Now for the integrators in Rxx. We assume symmetry:
% % Rxxa(x,tt) = Rxxb^T(tt,x)
if ~Zexclude(3) || ~Zexclude(4)
dx_a = zeros(3,1);                          % Initialize max degrees as zero
Rxxa_maxdeg_tot = Pmaxdegs.Rxx{2}(4,1);     % Extract maximum total degree in x*tt
dx_a_tot = ceil((Rxxa_maxdeg_tot - 1)/2);   % Total degree for monomials of root. Note -1 to account for integration
dx_b_tot = floor((Rxxa_maxdeg_tot - 1)/2);  % Total degree for monomials of root. Note -1 to account for integration
dx_a_tot = max([0,dx_a_tot]);               % Avoid negative degrees
dx_b_tot = max([0,dx_b_tot]);               % Avoid negative degrees
dx_a(3,1) = dx_a_tot;
dx_b(3,1) = dx_b_tot;
dx_a(1:2,1) = [ceil(dx_a_tot/2);floor(dx_a_tot/2)];     % Allow degrees of variables up to total degree dx_a_tot
dx_b(1:2,1) = [ceil(dx_b_tot/2);floor(dx_b_tot/2)];
end
if Zexclude(3)
    dx_a = [-Inf;-Inf;-Inf];
end
if Zexclude(4)
    dx_b = [-Inf;-Inf;-Inf];
end
dx = {dx_o; dx_a; dx_b};


% % Similarly for the integrators in Ryy, we assume symmetry:
% % Ryya(y,nu) = Ryyb^T(nu,y)
if ~Zexclude(6) || ~Zexclude(7)
dy_a = zeros(1,3);                          % Initialize max degrees as zero
Ryya_maxdeg_tot = Pmaxdegs.Ryy{2}(1,4);     % Extract maximum total degree in y*nu
dy_a_tot = ceil((Ryya_maxdeg_tot - 1)/2);   % Total degree for monomials of root. Note -1 to account for integration
dy_b_tot = floor((Ryya_maxdeg_tot - 1)/2);  % Total degree for monomials of root. Note -1 to account for integration
dy_a_tot = max([0,dy_a_tot]);               % Avoid negative degrees
dy_b_tot = max([0,dy_b_tot]);               % Avoid negative degrees
dy_a(1,3) = dy_a_tot;
dy_b(1,3) = dy_b_tot;
dy_a(1,1:2) = [ceil(dy_a_tot/2),floor(dy_a_tot/2)];     % Allow degrees of variables up to total degree dy_a_tot
dy_b(1,1:2) = [ceil(dy_b_tot/2),floor(dy_b_tot/2)];
end
if Zexclude(6)
    dy_a = [-Inf,-Inf,-Inf];
end
if Zexclude(7)
    dy_b = [-Inf,-Inf,-Inf];
end
dy = {dy_o, dy_a, dy_b};

% % Finally for R22, we have double integrators and multiplier-inegrators
% % R22ao(x,y,tt) = R22bo^T(tt,y,x)
if ~all(Zexclude(8:16))
d2_ao = zeros(4,2);                         % Initialize max degrees as zero
R22ao_maxdeg_tot = Pmaxdegs.R22{2,1}(4,2);  % Extract maximum total degree in x*tt*y
R22ao_maxdeg_xtt = Pmaxdegs.R22{2,1}(4,1);  % Extract maximum total degree in x*tt
d2_ao_tot = ceil((R22ao_maxdeg_tot-1)/2);   % Total degree for monomials of root. Note -1 to account for integration
d2_ao_tot = max([0,d2_ao_tot]);             % Avoid negative degrees
d2_ao_xtt = ceil((R22ao_maxdeg_xtt-1)/2);   % Degree for monomials of root in x*tt. Note -1 to account for integration
d2_ao_xtt = max([0,d2_ao_xtt]);             % Avoid negative degrees
d2_ao(4,1:2) = [d2_ao_xtt,d2_ao_tot];
d2_ao(2:3,1) = [ceil(d2_ao_xtt/2);floor(d2_ao_xtt/2)];  % Allow degrees of variables up to total degree d2_ao_xtt
d2_ao(1:3,2) = floor(Pmaxdegs.R22{2,1}(1,2)/2); % Set maximum degree in y

d2_bo_tot = floor((R22ao_maxdeg_tot-1)/2);  % Total degree for monomials of root. Note -1 to account for integration
d2_bo_tot = max([0,d2_bo_tot]);             % Avoid negative degrees
d2_bo_xtt = floor((R22ao_maxdeg_xtt-1)/2);  % Degree for monomials of root in x*tt. Note -1 to account for integration
d2_bo_xtt = max([0,d2_bo_xtt]);             % Avoid negative degrees
d2_bo(4,1:2) = [d2_bo_xtt,d2_bo_tot];
d2_bo(2:3,1) = [ceil(d2_bo_xtt/2);floor(d2_bo_xtt/2)];  % Allow degrees of variables up to total degree d2_bo_xtt
d2_bo(1:3,2) = floor(Pmaxdegs.R22{2,1}(1,2)/2); % Set maximum degree in y

% % R22oa(x,y,nu) = R22ob^T(x,nu,y)
d2_oa = zeros(2,4);                         % Initialize max degrees as zero
R22oa_maxdeg_tot = Pmaxdegs.R22{1,2}(2,4);  % Extract maximum total degree in x*y*nu
R22oa_maxdeg_ynu = Pmaxdegs.R22{1,2}(1,4);  % Extract maximum total degree in y*nu
d2_oa_tot = ceil((R22oa_maxdeg_tot-1)/2);   % Total degree for monomials of root. Note -1 to account for integration
d2_oa_tot = max([0,d2_oa_tot]);             % Avoid negative degrees
d2_oa_ynu = ceil((R22oa_maxdeg_ynu-1)/2);   % Degree for monomials of root in y*nu. Note -1 to account for integration
d2_oa_ynu = max([0,d2_oa_ynu]);             % Avoid negative degrees
d2_oa(1:2,4) = [d2_oa_ynu,d2_oa_tot];
d2_oa(1,2:3) = [ceil(d2_oa_ynu/2);floor(d2_oa_ynu/2)];  % Allow degrees of variables up to total degree d2_oa_ynu
d2_oa(2,1:3) = floor(Pmaxdegs.R22{1,2}(2,1)/2); % Set maximum degree in x

d2_ob_tot = floor((R22oa_maxdeg_tot-1)/2);  % Total degree for monomials of root. Note -1 to account for integration
d2_ob_tot = max([0,d2_ob_tot]);             % Avoid negative degrees
d2_ob_ynu = floor((R22oa_maxdeg_ynu-1)/2);  % Degree for monomials of root in y*nu. Note -1 to account for integration
d2_ob_ynu = max([0,d2_ob_ynu]);             % Avoid negative degrees
d2_ob(1:2,4) = [d2_ob_ynu,d2_ob_tot];
d2_ob(1,2:3) = [ceil(d2_ob_ynu/2);floor(d2_ob_ynu/2)];  % Allow degrees of variables up to total degree d2_ob_ynu
d2_ob(2,1:3) = floor(Pmaxdegs.R22{1,2}(2,1)/2); % Set maximum degree in x

% % R22aa(x,y,tt,nu) = R22bb^T(tt,nu,x,y)
% % R22ba(x,y,tt,nu) = R22ab^T(tt,nu,x,y)
d2_aa = zeros(4,4);                         % Initialize max degrees as zero
R22aa_maxdeg_xtt = Pmaxdegs.R22{2,2}(4,1);  % Extract maximum total degree in x*tt
R22aa_maxdeg_ynu = Pmaxdegs.R22{2,2}(1,4);  % Extract maximum total degree in y*nu
R22aa_maxdeg_tot = Pmaxdegs.R22{2,2}(4,4);  % Extract maximum total degree in x*y*tt*nu

d2_aa_tot = ceil((R22aa_maxdeg_tot-2)/2);  % Degree for monomials of root in x*y*tt*nu. Note -2 to account for double integration
d2_aa_tot = max([0,d2_aa_tot]);             % Avoid negative degrees
d2_aa(4,4) = d2_aa_tot;

d2_ba_tot = floor((R22aa_maxdeg_tot-2)/2);  % Degree for monomials of root in x*y*tt*nu. Note -2 to account for double integration
d2_ba(4,4) = max([0,d2_ba_tot]);
d2_ab(4,4) = d2_ba_tot;
d2_bb(4,4) = max([0,d2_ba_tot - 1]);

d2_aa(4,2:3) = Pmaxdegs.R22{2,2}(4,2:3)-[d2_aa(4,4),d2_ab(4,4)]-2;   % Degree for monomials in x*tt*y and x*tt*nu
d2_aa(2:3,4) = Pmaxdegs.R22{2,2}(2:3,4)-[d2_aa(4,4);d2_ba(4,4)]-2;   % Degree for monomials in x*tt*y and x*tt*nu
d2_aa(2,2) = Pmaxdegs.R22{2,2}(2,2)-d2_aa(4,4)-2;   % Degree for monomials in x*y
d2_aa(3,2) = Pmaxdegs.R22{3,2}(2,2)-d2_aa(4,4)-2;   % Degree for monomials in tt*y
d2_aa(2,3) = Pmaxdegs.R22{2,3}(2,2)-d2_aa(4,4)-2;   % Degree for monomials in x*nu
d2_aa(3,3) = Pmaxdegs.R22{3,3}(2,2)-d2_aa(4,4)-2;   % Degree for monomials in y*nu

d2_ba([2;3],2:4) = d2_aa([3;2],2:4);    d2_ba(4,2:3) = d2_aa(4,2:3);
d2_ab(2:4,[2,3]) = d2_aa(2:4,[3,2]);    d2_ab(2:3,4) = d2_aa(2:3,4);
d2_bb(2:4,[2,3]) = d2_ba(2:4,[3,2]);    d2_bb(2:3,4) = d2_ba(2:3,4);

d2_aa_xtt = ceil((R22aa_maxdeg_xtt-1)/2);   % Degree for monomials of root in x*tt. Note -1 to account for integration
d2_aa_xtt = max([0,d2_aa_xtt]);             % Avoid negative degrees
d2_aa(4,1) = d2_aa_xtt;
d2_aa(2:3,1) = [ceil(d2_aa_xtt/2);floor(d2_aa_xtt/2)]; % Allow degrees of variables x,tt up to total degree d2_aa_xtt
d2_ab(2:4,1) = d2_aa(2:4,1);

d2_ba_xtt = floor((R22aa_maxdeg_xtt-1)/2);  % Degree for monomials of root in x*tt. Note -1 to account for integration
d2_ba_xtt = max([0,d2_ba_xtt]);             % Avoid negative degrees
d2_ba(4,1) = d2_ba_xtt;
d2_ba(2:3,1) = [ceil(d2_ba_xtt/2);floor(d2_ba_xtt/2)]; % Allow degrees of variables x,tt up to total degree d2_ba_xtt
d2_bb(2:4,1) = d2_ba(2:4,1);

d2_aa_ynu = ceil((R22aa_maxdeg_ynu-1)/2);  % Degree for monomials of root in y*nu. Note -1 to account for integration
d2_aa_ynu = max([0,d2_aa_ynu]);             % Avoid negative degrees
d2_aa(1,4) = d2_aa_ynu;
d2_aa(1,2:3) = [ceil(d2_aa_ynu/2),floor(d2_aa_ynu/2)];  % Allow degrees of variables y,nu up to total degree d2_aa_ynu
d2_ba(1,2:4) = d2_aa(1,2:4);

d2_ab_ynu = floor((R22aa_maxdeg_ynu-1)/2);  % Degree for monomials of root in y*nu. Note -1 to account for integration
d2_ab_ynu = max([0,d2_ab_ynu]);             % Avoid negative degrees
d2_ab(1,4) = d2_ab_ynu;
d2_ab(1,2:3) = [ceil(d2_ab_ynu/2),floor(d2_ab_ynu/2)];  % Allow degrees of variables y,nu up to total degree d2_aa_ynu
d2_bb(1,2:4) = d2_ab(1,2:4);
end
if Zexclude(9)
    d2_ao = -Inf*ones(4,2);
end
if Zexclude(10)
    d2_bo = -Inf*ones(4,2);
end
if Zexclude(11)
    d2_oa = -Inf*ones(2,4);
end
if Zexclude(12)
    d2_ob = -Inf*ones(2,4);
end
if Zexclude(13)
    d2_aa = -Inf*ones(4,4);
end
if Zexclude(14)
    d2_ba = -Inf*ones(4,4);
end
if Zexclude(15)
    d2_ab = -Inf*ones(4,4);
end
if Zexclude(16)
    d2_bb = -Inf*ones(4,4);
end
d2 = {d2_oo, d2_oa, d2_ob; d2_ao, d2_aa, d2_ab; d2_bo, d2_ba, d2_bb};
for j=1:numel(d2)
    d2{j} = reduce_joint_degs(d2{j});
end

dZ2 = {[],[],[];[],[],[];[],[],[]};

ij11 = [2,2;3,2;2,3;3,3];   % Loop over parameters R2aa, R2ba, R2ab, R2bb
for k=1:size(ij11,1)
    i1 = ij11(k,1);     j1 = ij11(k,2);
    i2 = 5 - i1;        j2 = 5 - j1;    % 3->2 and 2->3
    
    Pdmat = Pdegs.R22{i1,j1};
    Pdmat11 = Pdmat((Pdmat(:,1)>Pdmat(:,2) & Pdmat(:,3)>Pdmat(:,4)),:); % x>tt, y>nu
    Pdmat21 = Pdmat((Pdmat(:,1)<Pdmat(:,2) & Pdmat(:,3)>Pdmat(:,4)),:);
    Pdmat12 = Pdmat((Pdmat(:,1)>Pdmat(:,2) & Pdmat(:,3)<Pdmat(:,4)),:);
    Pdmat22 = Pdmat((Pdmat(:,1)<Pdmat(:,2) & Pdmat(:,3)<Pdmat(:,4)),:);
    %Pdmat_r = setdiff(Pdmat,[Pdmat11;Pdmat21;Pdmat12;Pdmat22],'rows');

    temp11 = Pdmat11 - sparse([1,0,1,0]);   % Degrees of product before integration
    temp21 = Pdmat21 - sparse([0,1,1,0]);
    temp12 = Pdmat12 - sparse([1,0,0,1]);
    temp22 = Pdmat22 - sparse([0,1,0,1]);
    
    temp11(:,[1,3]) = max(ceil(temp11(:,[1,3]) - temp11(:,[2,4])),0);   % Degrees of "root" monomials
    temp21(:,[2,3]) = max(ceil(temp21(:,[2,3]) - temp21(:,[1,4])),0);
    temp12(:,[1,4]) = max(ceil(temp12(:,[1,4]) - temp12(:,[2,3])),0);
    temp22(:,[2,4]) = max(ceil(temp22(:,[2,4]) - temp22(:,[1,3])),0);

    dZ2{i1,j1} = unique([dZ2{i1,j1}; temp11],'rows');
    dZ2{i2,j1} = unique([dZ2{i2,j2}; temp21],'rows');
    dZ2{i1,j2} = unique([dZ2{i1,j2}; temp12],'rows');
    dZ2{i2,j2} = unique([dZ2{i2,j2}; temp22],'rows');
    
%     %testdeg11 = [[dZ2{i,j}(:,[1,2])+dZ2{i,j}(:,3,4]),dZ2{i,j}(:,3,4)]
%     nZ = size(dZ2{ii,jj},1);
%     dum_degs = repmat(dZ2{ii,jj}(:,[1,2]),[nZ,1]) + kron(dZ2{ii,jj}(:,[1,2]),ones([nZ,1])); % All possible degrees of eta mu;
%     var_degs = repmat([repmat(dZ2{ii,jj}(:,[3,4]),[nZ,1]),kron(dZ2{ii,jj}(:,[3,4]),ones([nZ,1]))],[4,1]);   % Degrees of x,y,tt,nu without integration
%     var_degs(1:nZ^2,1:2) = var_degs(1:nZ^2,1:2) + dum_degs + [1,1]; % Contribution from integration
%     var_degs(nZ^2+1:2*nZ^2,1) = var_degs(nZ^2+1:2*nZ^2,1) + dum_degs(:,1) + 1; % Contribution from integration just in eta
%     var_degs(2*nZ^2+1:3*nZ^2,2) = var_degs(2*nZ^2+1:3*nZ^2,1) + dum_degs(:,2) + 1; % Contribution from integration just in eta
%     var_degs = unique(var_degs,'rows');
%     
%     
%     % Monomials resulting from constant limits of integrals
%     testdeg11 = [repmat(dZ2{i1,j1}(:,[1,2]),[size(dZ2{i1,j1},1),1]),kron(dZ2{i1,j1}(:,[3,4]),ones([size(dZ2{i1,j1},1),1]))];
%     testdeg21 = [repmat(dZ2{i2,j1}(:,[1,2]),[size(dZ2{i2,j1},1),1]),kron(dZ2{i2,j1}(:,[3,4]),ones([size(dZ2{i2,j1},1),1]))];
%     testdeg12 = [repmat(dZ2{i1,j2}(:,[1,2]),[size(dZ2{i1,j2},1),1]),kron(dZ2{i1,j2}(:,[3,4]),ones([size(dZ2{i1,j2},1),1]))];
%     testdeg22 = [repmat(dZ2{i2,j2}(:,[1,2]),[size(dZ2{i2,j2},1),1]),kron(dZ2{i2,j2}(:,[3,4]),ones([size(dZ2{i2,j2},1),1]))];
%     testdeg = [testdeg11;testdeg21;testdeg12;testdeg22];
    
end

ij11 = [2,2;3,2;2,3;3,3];   % Loop over parameters R2aa, R2ba, R2ab, R2bb
for k=1:size(ij11,1)
    i1 = ij11(k,1);     j1 = ij11(k,2);
    
    for l=1:size(ij11,1)    % Loop over contributions from monomials Z2aa, Z2ba, Z2ab, Z2bb
        i2 = ij11(l,1);         j2 = ij11(l,2);
        dum_degs = dZ2{i2,j2}(:,[1,3]) + dZ2{i2,j2}(:,[1,3]);   % All possible degrees of dummy variables eta, mu
        v1 = 1 + abs(i2-i1);   v2 = 3 + abs(j2-j1);             % Indices of variables to which the dummy variables will contribute
        
        var_degs = repmat({[dZ2{i2,j2}(:,[2,4]),dZ2{i2,j2}(:,[2,4])]},[2,2]);   % Degrees of x,y,tt,nu without integration
        var_degs{2,1}(:,v1) = var_degs{1,1}(:,v1) + dum_degs(:,1) + 1;          % Add contribution from integration just in eta
        var_degs{1,2}(:,v2) = var_degs{1,1}(:,v2) + dum_degs(:,2) + 1;          % Add contribution from integration just in mu
        var_degs{2,2}(:,[v1,v2]) = var_degs{1,1}(:,[v1,v2]) + dum_degs + [1,1]; % Contribution from integration in both eta and mu
    
        isneeded_Z = ismember(var_degs{1,1},Pdegs.R22{i1,j1},'rows') | ...
                     ismember(var_degs{2,1},Pdegs.R22{i1,j1},'rows') | ...
                     ismember(var_degs{1,2},Pdegs.R22{i1,j1},'rows') | ...
                     ismember(var_degs{2,2},Pdegs.R22{i1,j1},'rows');       % Check whether any of the four terms the monomial produces is actually present
             
        dZ2{i2,j2} = dZ2{i2,j2}(isneeded_Z,:);  % Remove monomials that do not contribute
    end
end
    

for k = [5,6,8,9]    
    degs_k = dZ2{k};
    for l = 2:2^4
        indcs = cell(1,4);
        [indcs{:}] = ind2sub(2*ones(1,4),l);                % Indcs associated to position l in maxdegs
        log_indcs = cell2mat(indcs)==2;                     % Logical values indicating which vars contribute to the considered joint degree
        maxdegs_k(l) = max(sum(degs_k(:,log_indcs),2));     % Take maximal joint degree in the considered variables
    end
    d22{k} = maxdegs_k;
    
end

% % % Next, test whether this initial guess for the degrees recovers all
% % % the necessary monomials
% Compute maximal degrees for each parameter in the poslpivar
Dmaxdegs = maxdegs_poslpivar(dx,dy,d2);

% % Increase the degrees (if necessary) until all the necessary monomials
% % are included, starting with the parameters on the diagonal of the operator

% Start with Rxx
if n(2)>=1
max_updates = 10;       % Limit the number of times we will try updating the degree
isgood_dx = zeros(3,1); % Binary array indicating if parameter Rxx{i} has all necessary monomials
isgood_dx(1) = 1;       % We know Rxx{1} to have all necessary monomials
isgood_dx = isgood_dx | Zexclude([2;3;4]);
for k=1:max_updates
    % Check parameters Rxxa and Rxxb
    if ~all(isgood_dx(2:3))
    for i=2:3
        check_indcs = {};   % Indices associated to degrees we need to check are sufficient
        update_indcs = {};  % Indices associated to degrees we need to update
        if any(Dmaxdegs.Rxx{i}(4,1)<Pmaxdegs.Rxx{i}(4,1))
            check_indcs = {([4,1]-[0,1])*[1;4]};
            update_indcs = {[3;1;2]};                   % Update degrees pertaining to any combination of x,tt
        end
        if any(Dmaxdegs.Rxx{i}(2:3,1)<Pmaxdegs.Rxx{i}(2:3,1))
            check_indcs = [check_indcs, ([2,1; 3,1]-[0,1])*[1;4]];
            update_indcs = [update_indcs, [1;2]];       % Update degrees pertaining to x, tt separately
        end
        if isempty(update_indcs)
            isgood_dx(i) = 1;
        else
            for m=1:length(check_indcs)
                c_indx = check_indcs{m};
                for u_indx = update_indcs{m}'
                    if all(Dmaxdegs.Rxx{i}(c_indx)>=Pmaxdegs.Rxx{i}(c_indx))
                        break
                    end
                    dx{i}(u_indx) = dx{i}(u_indx)+1;  % Increase degrees suggested by index
                    Dmaxdegs = maxdegs_poslpivar(dx,dy,d2,Dmaxdegs,'Rxx');
                end
            end
            if all(Dmaxdegs.Rxx{i}(2:4,1)>=Pmaxdegs.Rxx{i}(2:4,1))
                isgood_dx(i) = 1;
            end
        end
    end
    end
    if k==max_updates && ~all(isgood_dx)
        warning(['More than the ',num2str(max_updates),' allowed updates to the degrees for P.Rxx are necessary, something might be wrong...'])
    end
end
end

% Then Ryy
if n(3)>=1
max_updates = 10;       % Limit the number of times we will try updating the degree
isgood_dy = zeros(1,3); % Binary array indicating if parameter Ryy{j} has all necessary monomials
isgood_dy(1) = 1;       % We know Ryy{1} to have all necessary monomials
isgood_dy = isgood_dy | Zexclude([5,6,7]);
for k=1:max_updates
    % Check parameters Ryya and Ryyb
    if ~all(isgood_dy(2:3))
    for j=2:3
        check_indcs = {};   % Indices associated to degrees we need to check are sufficient
        update_indcs = {};  % Indices associated to degrees we need to update
        if any(Dmaxdegs.Ryy{j}(1,4)<Pmaxdegs.Ryy{j}(1,4))
            check_indcs = {([1,4]-[0,1])*[1;4]};
            update_indcs = {[3;1;2]};                   % Update degrees pertaining to any combination of y,nu
        end
        if any(Dmaxdegs.Ryy{j}(1,2:3)<Pmaxdegs.Ryy{j}(1,2:3))
            check_indcs = [check_indcs, ([1,2; 1,3]-[0,1])*[1;4]];
            update_indcs = [update_indcs, [1;2]];       % Update degrees pertaining to y, nu separately
        end
        if isempty(update_indcs)
            isgood_dy(j) = 1;
        else
            for m=1:length(check_indcs)
                c_indx = check_indcs{m};
                for u_indx = update_indcs{m}'
                    if all(Dmaxdegs.Ryy{j}(c_indx)>=Pmaxdegs.Ryy{j}(c_indx))
                        break
                    end
                    dy{j}(u_indx) = dy{j}(u_indx)+1;  % Increase degrees suggested by index
                    Dmaxdegs = maxdegs_poslpivar(dx,dy,d2,Dmaxdegs,'Ryy');
                end
            end
            if all(Dmaxdegs.Ryy{j}(1,2:4)>=Pmaxdegs.Ryy{j}(1,2:4))
                isgood_dy(j) = 1;
            end
        end
    end
    end
    if k==max_updates && ~all(isgood_dy)
        warning(['More than the ',num2str(max_updates),' allowed updates to the degrees for P.Ryy are necessary, something might be wrong...'])
    end
end
end

% And finally R22
if n(4)>=1
max_updates = 10;       % Limit the number of times we will try updating the degree
isgood_d2 = zeros(3,3); % Binary array indicating if parameter R22{i,j} has all necessary monomials
isgood_d2(1,1) = 1;     % We know R22{1,1} to have all necessary monomials
isgood_d2 = isgood_d2 | [Zexclude([8,11,12]);Zexclude([9,13,15]);Zexclude([10,14,16])];
for k=1:max_updates
    % Check parameters R22ao and R22bo
    for i=2:3
    if ~isgood_d2(i,1)
        check_indcs = {};   % Indices associated to degrees we need to check are sufficient
        update_indcs = {};  % Indices associated to degrees we need to update
        if any(Dmaxdegs.R22{i,1}(4,2)<Pmaxdegs.R22{i,1}(4,2))
            check_indcs = {([4,2]-[0,1])*[1;4]};
            update_indcs = {([4,2; 2,2; 3,2; 4,1; 2,1; 3,1]-[0,1])*[1;4]};          % Update degrees pertaining to x,y,tt
        end
        if any(Dmaxdegs.R22{i,1}(2:3,2)<Pmaxdegs.R22{i,1}(2:3,2))
            check_indcs = [check_indcs, ([2,2; 3,2]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([2,2; 3,2; 2,1; 3,1]-[0,1])*[1;4]];     % Update degrees pertaining to pairs of x,tt and y
        end
        if any(Dmaxdegs.R22{i,1}(4,1)<Pmaxdegs.R22{i,1}(4,1))
            check_indcs = [check_indcs, ([4,1]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([4,1; 2,1; 3,1]-[0,1])*[1;4]];          % Update degrees pertaining to x,tt
        end
        if any(Dmaxdegs.R22{i,1}(2:3,1)<Pmaxdegs.R22{i,1}(2:3,1))
            check_indcs = [check_indcs, ([2,1; 3,1]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([2,1; 3,1]-[0,1])*[1;4]];    % Update degrees pertaining to x, tt
        end
        if any(Dmaxdegs.R22{i,1}(1,2)<Pmaxdegs.R22{i,1}(1,2))
            check_indcs = [check_indcs, ([1,2]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([1,2]-[0,1])*[1;4]];         % Update degrees pertaining to y
        end
        if isempty(update_indcs)
            isgood_d2(i,1) = 1;
        else
            for m=1:length(check_indcs)
                c_indx = check_indcs{m};
                for u_indx = update_indcs{m}'
                    if all(Dmaxdegs.R22{i,1}(c_indx)>=Pmaxdegs.R22{i,1}(c_indx))
                        break
                    end
                    d2{i,1}(u_indx) = d2{i,1}(u_indx)+1;  % Increase degrees suggested by index
                    Dmaxdegs = maxdegs_poslpivar(dx,dy,d2,Dmaxdegs,'R22');
                end
            end
            if all(all(Dmaxdegs.R22{i,1}(1:4,1:2)>=Pmaxdegs.R22{i,1}(1:4,1:2)))
                isgood_d2(i,1) = 1;
            end
        end
        % Make sure the joint degrees are not smaller than the individual
        % degrees, nor greater than the sum of the individual degrees
        d2{i,1} = reduce_joint_degs(d2{i,1});
    end
    end
    % Check parameters R22oa and R22ob
    for j=2:3
    if ~isgood_d2(1,j)
        check_indcs = {};   % Indices associated to degrees we need to check are sufficient
        update_indcs = {};  % Indices associated to degrees we need to update
        if any(Dmaxdegs.R22{1,j}(2,4)<Pmaxdegs.R22{1,j}(2,4))
            check_indcs = {([2,4]-[0,1])*[1;4]};
            update_indcs = {([2,4; 2,2; 2,3; 1,4; 1,2; 1,3]-[0,1])*[1;4]};          % Update degrees pertaining to x,y,nu
        end
        if any(Dmaxdegs.R22{1,j}(2,2:3)<Pmaxdegs.R22{1,j}(2,2:3))
            check_indcs = [check_indcs, ([2,2; 2,3]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([2,2; 2,3; 1,2; 1,3]-[0,1])*[1;4]];     % Update degrees pertaining to pairs of x and y,nu
        end
        if any(Dmaxdegs.R22{1,j}(1,4)<Pmaxdegs.R22{1,j}(1,4))
            check_indcs = [check_indcs, ([1,4]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([1,4; 1,2; 1,3]-[0,1])*[1;4]];          % Update degrees pertaining to y,nu
        end
        if any(Dmaxdegs.R22{1,j}(1,2:3)<Pmaxdegs.R22{1,j}(1,2:3))
            check_indcs = [check_indcs, ([1,2; 1,3]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([1,2; 1,3]-[0,1])*[1;4]];    % Update degrees pertaining to y,nu
        end
        if any(Dmaxdegs.R22{1,j}(2,1)<Pmaxdegs.R22{1,j}(2,1))
            check_indcs = [check_indcs, ([2,1]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([2,1]-[0,1])*[1;4]];         % Update degrees pertaining to x
        end
        if isempty(update_indcs)
            isgood_d2(1,j) = 1;
        else
            for m=1:length(check_indcs)
                c_indx = check_indcs{m};
                for u_indx = update_indcs{m}'
                    if all(Dmaxdegs.R22{1,j}(c_indx)>=Pmaxdegs.R22{1,j}(c_indx))
                        break
                    end
                    d2{1,j}(u_indx) = d2{1,j}(u_indx)+1;  % Increase degrees suggested by index
                    Dmaxdegs = maxdegs_poslpivar(dx,dy,d2,Dmaxdegs,'R22');
                end
            end
            if all(all(Dmaxdegs.R22{1,j}(1:2,1:4)>=Pmaxdegs.R22{1,j}(1:2,1:4)))
                isgood_d2(1,j) = 1;
            end
        end
        % Make sure the joint degrees are not smaller than the individual
        % degrees, nor greater than the sum of the individual degrees
        d2{1,j} = reduce_joint_degs(d2{1,j});
    end
    end
    % Check parameters R22aa and R22ba
    for ij=[5,6,8,9]
    if ~isgood_d2(ij)
        check_indcs = {};   % Indices associated to degrees we need to check are sufficient
        update_indcs = {};  % Indices associated to degrees we need to update
        if any(Dmaxdegs.R22{ij}(4,4)<Pmaxdegs.R22{ij}(4,4))
            check_indcs = {([4,4]-[0,1])*[1;4]};
            update_indcs = {([4,4; 4,2; 2,4; 4,3; 3,4; 2,2; 3,3; 3,2; 2,3; 4,1; 1,4; 2,1; 1,2; 3,1; 1,3]-[0,1])*[1;4]}; % Update all degrees pertaining to x,y,tt,nu
        end
        if any(Dmaxdegs.R22{ij}(4,2:3)<Pmaxdegs.R22{ij}(4,2:3))
            check_indcs = [check_indcs,([4,2; 4,3]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([4,2; 4,3; 2,2; 3,3; 3,2; 2,3; 4,1; 2,1; 1,2; 3,1; 1,3]-[0,1])*[1;4]]; % Update all degrees pertaining to (x,tt),y,nu
        end
        if any(Dmaxdegs.R22{ij}(2:3,4)<Pmaxdegs.R22{ij}(2:3,4))
            check_indcs = [check_indcs,([2,4; 3,4]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([2,4; 3,4; 2,2; 3,3; 2,3; 3,2; 1,4; 1,2; 2,1; 1,3; 3,1]-[0,1])*[1;4]]; % Update all degrees pertaining to x,tt,(y,nu)
        end
        if any(any(Dmaxdegs.R22{ij}(2:3,2:3)<Pmaxdegs.R22{ij}(2:3,2:3)))
            check_indcs = [check_indcs,([2,2; 3,3; 3,2; 2,3]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([2,2; 3,3; 3,2; 2,3; 2,1; 1,2; 3,1; 1,3]-[0,1])*[1;4]]; % Update all degrees pertaining to pairs of x,tt,y,nu
        end
        if any(Dmaxdegs.R22{ij}(4,1)<Pmaxdegs.R22{ij}(4,1))
            check_indcs = [check_indcs,([4,1]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([4,1; 2,1; 3,1]-[0,1])*[1;4]]; % Update all degrees pertaining to x,tt
        end
        if any(Dmaxdegs.R22{ij}(1,4)<Pmaxdegs.R22{ij}(1,4))
            check_indcs = [check_indcs,([1,4]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([1,4; 1,2; 1,3]-[0,1])*[1;4]]; % Update all degrees pertaining to y,nu
        end
        if any(Dmaxdegs.R22{ij}(2:3,1)<Pmaxdegs.R22{ij}(2:3,1))
            check_indcs = [check_indcs,([2,1; 3,1]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([2,1; 3,1]-[0,1])*[1;4]]; % Update all degrees pertaining to x, tt
        end
        if any(Dmaxdegs.R22{ij}(1,2:3)<Pmaxdegs.R22{ij}(1,2:3))
            check_indcs = [check_indcs,([1,2; 1,3]-[0,1])*[1;4]];
            update_indcs = [update_indcs, ([1,2; 1,3]-[0,1])*[1;4]]; % Update all degrees pertaining to y, nu
        end
        if isempty(update_indcs)
            isgood_d2(ij) = 1;
            continue
        else
            for m=1:length(check_indcs)
                c_indx = check_indcs{m};
                for u_indx = update_indcs{m}'
                    if all(Dmaxdegs.R22{ij}(c_indx)>=Pmaxdegs.R22{ij}(c_indx))
                        break
                    end
                    d2{ij}(u_indx) = d2{ij}(u_indx)+1;  % Increase degrees suggested by index
                    Dmaxdegs = maxdegs_poslpivar(dx,dy,d2,Dmaxdegs,'R22');
                end
            end
            if all(all(Dmaxdegs.R22{ij}(:,:)>=Pmaxdegs.R22{ij}(:,:)))
                isgood_d2(ij) = 1;
                continue
            end
        end
        % Make sure the joint degrees are not smaller than the individual
        % degrees, nor greater than the sum of the individual degrees
        d2{ij} = reduce_joint_degs(d2{ij});
    end
    end
    if k==max_updates && ~all(all(isgood_d2))
        warning(['More than the ',num2str(max_updates),' allowed updates to the degrees for P.R22 are necessary, something might be wrong...'])
    end
end
end

% % % At this point, the degrees dx, dy, d2 should be such that the
% % % parameters on the diagonal of the positive operator will contain all
% % % the necessary monomials. It remains to check the off-diagonal
% % % parameters.

if n(1)>=1  % Start with parameters mapping to \R^n
if n(2)>=1  % Mapping from L_2[x]
    if Dmaxdegs.R0x(2,1)<Pmaxdegs.R0x(2,1)
        [~,param] = max([dx{2}(3),dx{3}(3)]);       % Check which parameter has greatest total degree
        param = param + 1;
        dx{param}(3) = Pmaxdegs.R0x(2,1)-1;         % Increase total degree to necessary value (accounting for integration)
        d_dif = max(dx{param}(3) - sum(dx{param}(1:2)),0);
        dx{param}(1:2) = dx{param}(1:2) + [ceil(d_dif/2);floor(d_dif/2)];   % Increase degrees of individual vars to match desired joint degree
    end
    Dmaxdegs = maxdegs_poslpivar(dx,dy,d2,Dmaxdegs,'Rxx');
end
if n(3)>=1  % Mapping from L_2[y]
    if Dmaxdegs.R0y(1,2)<Pmaxdegs.R0y(1,2)
        [~,param] = max([dy{2}(3),dy{3}(3)]);       % Check which parameter has greatest total degree
        param = param + 1;
        dy{param}(3) = Pmaxdegs.R0y(1,2)-1;         % Increase total degree to necessary value (accounting for integration)
        d_dif = max(dy{param}(3) - sum(dy{param}(1:2)),0);
        dy{param}(1:2) = dy{param}(1:2) + [ceil(d_dif/2),floor(d_dif/2)];   % Increase degrees of individual vars to match desired joint degree
    end
    Dmaxdegs = maxdegs_poslpivar(dx,dy,d2,Dmaxdegs,'Ryy');
end
if n(4)>=1  % Mapping from L_2[x,y]
    if Dmaxdegs.R02(2,2)<Pmaxdegs.R02(2,2)
        [~,param] = max([d2{2,2}(4,4),d2{3,2}(4,4),-Inf,d2{2,3}(4,4),d2{3,3}(4,4)]);       % Check which parameter has greatest total degree
        param = param + 4;
        d2{param}(4,4) = Pmaxdegs.R02(2,2)-2;         % Increase total degree to necessary value (accounting for integration)
        d2{param}(4,1:3) = max(d2{param}(4,1:3),ceil(d2{param}(4,4)/2));
        d2{param}(1:3,4) = max(d2{param}(1:3,4),ceil(d2{param}(4,4)/2));
        d_difx = max(d2{param}(4,1) - sum(d2{param}(2:3,1)),0);
        d_dify = max(d2{param}(1,4) - sum(d2{param}(1,2:3)),0);
        d2{param}(2:3,1) = d2{param}(2:3,1) + [ceil(d_difx/2);floor(d_difx/2)];   % Increase degrees of individual vars to match desired joint degree
        d2{param}(1,2:3) = d2{param}(1,2:3) + [ceil(d_dify/2),floor(d_dify/2)];   % Increase degrees of individual vars to match desired joint degree
        d2{param}(2:3,2:3) = max(d2{param}(2:3,2:3),[ceil(d2{param}(4,4)/2),floor(d2{param}(4,4)/2);floor(d2{param}(4,4)/2),floor(d2{param}(4,4)/2)]);
    end
    Dmaxdegs = maxdegs_poslpivar(dx,dy,d2,Dmaxdegs,'R22');
end
   
end


% Get rid of artificially imposed -Inf degrees
for j=1:numel(dx)
    for k=1:numel(dx{j})
        dx{j}(k) = max(dx{j}(k),0);
    end
end
for j=1:numel(dy)
    for k=1:numel(dy{j})
        dy{j}(k) = max(dy{j}(k),0);
    end
end
for j=1:numel(d2)
    for k=1:numel(d2{j})
        d2{j}(k) = max(d2{j}(k),0);
    end
    d2{j} = reduce_joint_degs(d2{j});
end
d.dx = dx;  d.dy = dy;  d.d2 = d2;

end



function D = maxdegs_poslpivar(dx,dy,d2,D,params)
% This function compute the maximal degrees of each variable in each
% parameter of the positive opvar2d object constructed using "poslpivar2d"
% with monomial degrees dx,dy,d2. Degrees are collected in a struct D, with
% the same fieldnames as a standard opvar2d object:
% D.R00, D.R0x, D.R0y, D.R02;
% D.Rx0, D.Rxx, D.Rxy, D.Rx2;
% D.Ry0, D.Ryx, D.Ryy, D.Ry2;
% D.R20, D.R2x, D.R2y, D.R22;
% Each field contains a 4x4 array providing the maximal degrees of the
% variables ss1, ss2, tt1 and tt2 as they appear in the associated
% parameter of the positive operator. Maximal joint degrees in
% combinations of variables are also provided, ordered as
%
%                 0 |         ss2 |         tt2 |         ss2*tt2
%          ---------|-------------|-------------|-----------------
%               ss1 |     ss1*ss2 |     ss1*tt2 |     ss1*ss2*tt2
%          ---------|-------------|-------------|-----------------
%               tt1 |     tt1*ss2 |     tt1*tt2 |     tt1*ss2*tt2 
%          ---------|-------------|-------------|-----------------
%           ss1*tt1 | ss1*tt1*ss2 | ss1*tt1*tt2 | ss1*tt1*ss2*tt2  
%
% If a variable does not appear in a particular parameter, the associated
% degrees and joint degrees will be set to zero. For example, R0x is a
% function of only ss1, so the field D.R0x will look like
% D.R0x = [0,0,0,0; d,0,0,0; 0,0,0,0; 0,0,0,0] for some value "d", though
% D.R0x = [0,0,0,0; d,0,0,0; 0,0,0,0; d,0,0,0] is also accepted.
%
% See the function poslpivar2d for more information on what the inputs dx,
% dy, d2 should look like, and how an associated positive operator is
% constructed.

if nargin<3
    error('At least 3 inputs are needed')
elseif nargin==3
    D = struct();
    OO = zeros(4,4);
    D.R00 = OO; D.R0x = OO; D.R0y = OO; D.R02 = OO;
    D.Rxx = {OO;OO;OO};     D.Rxy = OO; D.Rx2 = {OO;OO;OO};
    D.Ryy = {OO,OO,OO};     D.Ry2 = {OO,OO,OO};
    D.R22 = {OO,OO,OO;OO,OO,OO;OO,OO,OO};
    
    comp_params = ones(4,4);
elseif nargin==4
    comp_params = ones(4,4);
else
    if ischar(params)
        comp_params = strcmp({'R00','R0x','R0y','R02';
                              'Rx0','Rxx','Rxy','Rx2';
                              'Ry0','Ryx','Ryy','Ry2';
                              'R20','R2x','R2y','R22'},params);
    elseif iscellstr(params)
        comp_params = zeros(4,4);
        for j=numel(params)
            comp_params = strcmp({'R00','R0x','R0y','R02';
                                  'Rx0','Rxx','Rxy','Rx2';
                                  'Ry0','Ryx','Ryy','Ry2';
                                  'R20','R2x','R2y','R22'},params{j}) | comp_params;
        end
    else
        error('Final input argument must be a cellstr specifying for which parameters to update the degrees')
    end
end

% Make sure the joint degrees are not smaller than the individual
% degrees, nor greater than the sum of the individual degrees
for i=2:3
    d2{i,1}(4,2) = min(d2{i,1}(4,2),d2{i,1}(2,1)+d2{i,1}(3,1)+d2{i,1}(1,2));            % x*y*tt
    d2{i,1}(4,1) = min(d2{i,1}(4,1),min(d2{i,1}(2,1)+d2{i,1}(3,1),d2{i,1}(4,2)));       % x*tt
    d2{i,1}(2:3,2) = min(d2{i,1}(2:3,2),min(d2{i,1}(2:3,1)+d2{i,1}(1,2),d2{i,1}(4,1))); % x*y and tt*y
    d2{i,1}(2:3,1) = min(d2{i,1}(2:3,1),min(d2{i,1}(2:3,2),d2{i,1}(4,1)));      % x and tt
    d2{i,1}(1,2) = min(d2{i,1}(1,2),min(d2{i,1}(2:3,2)));                       % y
end
for ij = [5,6,8,9]
    d2{ij}(4,4) = min(d2{ij}(4,4),d2{ij}(2,1)+d2{ij}(3,1)+d2{ij}(1,2)+d2{ij}(1,3));             % x*y*tt*nu
    d2{ij}(4,2:3) = min(d2{ij}(4,2:3),min(d2{ij}(2,1)+d2{ij}(3,1)+d2{ij}(1,2:3),d2{ij}(4,4)));  % x*tt*y and x*tt*nu
    d2{ij}(2:3,4) = min(d2{ij}(2:3,4),min(d2{ij}(1,2)+d2{ij}(1,3)+d2{ij}(2:3,1),d2{ij}(4,4)));  % x*y*nu and tt*y*nu
    d2{ij}(2:3,2:3) = min(d2{ij}(2:3,2:3), min(d2{ij}(2:3,1)+d2{ij}(1,2:3), min([d2{ij}(4,2:3);d2{ij}(4,2:3)],[d2{ij}(2:3,4),d2{ij}(2:3,4)]))); % x*y, tt*y, x*nu and tt*nu
    d2{ij}(4,1) = min(d2{ij}(4,1),min(d2{ij}(2,1)+d2{ij}(3,1),min(d2{ij}(4,2:3))));     % x*tt
    d2{ij}(1,4) = min(d2{ij}(1,4),min(d2{ij}(1,2)+d2{ij}(1,3),min(d2{ij}(2:3,4))));     % y*nu
    d2{ij}(2:3,1) = min(d2{ij}(2:3,1),min(d2{ij}(4,1),min(d2{ij}(2:3,2:3),[],2)));  % x and tt
    d2{ij}(1,2:3) = min(d2{ij}(1,2:3),min(d2{ij}(1,4),min(d2{ij}(2:3,2:3),[],1)));  % y and nu
end


% For example:
% R0x(x) = P_01*Z_1(x) + int_x^b P_02*Z_2(tt,x) dtt
%           + int_a^x P_03*Z_3(tt,x) dtt
% so maximal degree of R0x(x) in x will be max of max degree of Z_1 in x,
% max joint degree of Z_2 in (tt,x), and max joint degree of Z_3 in (tt,x)
if comp_params(1,2) || comp_params(2,1)
    D.R0x(2,1) = max([dx{1}, dx{2}(3)+1, dx{2}(3)+1]);
end
if comp_params(1,3) || comp_params(3,1)
    D.R0y(1,2) = max([dy{1}, dy{2}(3)+1, dy{2}(3)+1]);     % Similarly for R0y
end
if comp_params(1,4) || comp_params(4,1)
D.R02(2,1) = max([d2{1,1}(2,1), d2{1,2}(2,1), d2{1,3}(2,1),...
                      d2{2,1}(4,1)+1, d2{2,2}(4,1)+1, d2{2,3}(4,1)+1,...
                      d2{3,1}(4,1)+1, d2{3,2}(4,1)+1, d2{3,3}(4,1)+1]);     % And also R02
D.R02(1,2) = max([d2{1,1}(1,2), d2{1,2}(1,4)+1, d2{1,3}(1,4)+1,...
                      d2{2,1}(1,2), d2{2,2}(1,4)+1, d2{2,3}(1,4)+1,...
                      d2{3,1}(1,2), d2{3,2}(1,4)+1, d2{3,3}(1,4)+1]);
D.R02(2,2) = max([d2{1,1}(2,2), d2{1,2}(2,4)+1, d2{1,3}(2,4)+1,...
                  d2{2,1}(4,2)+1, d2{2,2}(4,4)+2, d2{2,3}(4,4)+2,...
                  d2{3,1}(4,2)+1, d2{3,2}(4,4)+2, d2{3,3}(4,4)+2]);
end
              
% % % Next, Rxx, Rxy and Rx2
if comp_params(2,2)
D.Rxx{1}(2,1) = 2*dx{1};
D.Rxx{2}(2,1) = max([dx{1}+dx{2}(1), dx{3}(2), dx{2}(3)+dx{2}(1)+1,...
                        dx{2}(1)+dx{3}(3)+1, dx{3}(2)]);
D.Rxx{2}(3,1) = max([dx{2}(2), dx{1}+dx{3}(1), dx{2}(2),...
                        dx{2}(3)+dx{3}(1)+1, dx{3}(3)+dx{3}(1)+1]);
D.Rxx{2}(4,1) = max([dx{1}+dx{2}(3), dx{1}+dx{3}(3), 2*dx{2}(3)+1,...
                        dx{2}(3)+dx{3}(3)+1, 2*dx{3}(3)+1]);
D.Rxx{3} = D.Rxx{2};
D.Rxx{3}([2;3],:) = D.Rxx{3}([3;2],:);  % Rxxb(x,tt) = Rxxa^T(tt,x)
end

% %
if comp_params(2,3) || comp_params(3,2)
D.Rxy(2,1) = max([dx{1}, dx{2}(3)+1, dx{3}(3)+1]);
D.Rxy(1,2) = max([dy{1}, dy{2}(3)+1, dy{3}(3)+1]);
D.Rxy(2,2) = D.Rxy(2,1) + D.Rxy(1,2);   D.Rxy(4,4) = D.Rxy(2,2);
D.Rxy(2,4) = D.Rxy(2,2);                D.Rxy(4,2) = D.Rxy(2,2);
end

% %
if comp_params(2,4) || comp_params(4,2)
D.Rx2{1}(2,1) = max([dx{1}+d2{1,1}(2,1), dx{1}+d2{1,2}(2,1), dx{1}+d2{1,3}(2,1)]);
D.Rx2{1}(1,2) = max([d2{1,1}(1,2), d2{1,2}(1,4)+1, d2{1,3}(1,4)]);
D.Rx2{1}(2,2) = max([dx{1}+d2{1,1}(2,2), dx{1}+d2{1,2}(2,4)+1, dx{1}+d2{1,3}(2,4)]);
D.Rx2{1}(4,2) = D.Rx2{1}(2,2);  % We can omit this step

%
D.Rx2{2}(2,1) = max([dx{1}+d2{2,1}(2,1), dx{3}(2), dx{2}(3)+d2{2,1}(2,1)+1,...
                        dx{3}(3)+d2{2,1}(2,1)+1, dx{3}(2),...
                     dx{1}+d2{2,2}(2,1), dx{3}(2), dx{2}(3)+d2{2,2}(2,1)+1,...
                        dx{3}(3)+d2{2,2}(2,1)+1, dx{3}(2),...
                     dx{1}+d2{2,3}(2,1), dx{3}(2), dx{2}(3)+d2{2,3}(2,1)+1,...
                        dx{3}(3)+d2{2,3}(2,1)+1, dx{3}(2)]);
D.Rx2{2}(3,1) = max([d2{2,1}(3,1), dx{3}(1)+d2{1,1}(2,1), d2{2,1}(3,1),...
                        dx{3}(1)+d2{2,1}(4,1)+1, dx{3}(1)+d2{3,1}(4,1)+1,...
                     d2{2,2}(3,1), dx{3}(1)+d2{1,2}(2,1), d2{2,2}(3,1),...
                        dx{3}(1)+d2{2,2}(4,1)+1, dx{3}(1)+d2{3,2}(4,1)+1,...
                     d2{2,3}(3,1), dx{3}(1)+d2{1,3}(2,1), d2{2,3}(3,1),...
                        dx{3}(1)+d2{2,3}(4,1)+1, dx{3}(1)+d2{3,3}(4,1)+1]);
D.Rx2{2}(4,1) = max([dx{1}+d2{2,1}(4,1), dx{3}(3)+d2{1,1}(2,1), dx{2}(3)+d2{2,1}(4,1)+1,...
                        dx{3}(3)+d2{2,1}(4,1)+1, dx{3}(3)+d2{3,1}(4,1)+1,...
                     dx{1}+d2{2,2}(4,1), dx{3}(3)+d2{1,2}(2,1), dx{2}(3)+d2{2,2}(4,1)+1,...
                        dx{3}(3)+d2{2,2}(4,1)+1, dx{3}(3)+d2{3,2}(4,1)+1,...
                     dx{1}+d2{2,3}(4,1), dx{3}(3)+d2{1,3}(2,1), dx{2}(3)+d2{2,3}(4,1)+1,...
                        dx{3}(3)+d2{2,3}(4,1)+1, dx{3}(3)+d2{3,3}(4,1)+1]);
D.Rx2{2}(1,2) = max([d2{1,1}(1,2), d2{1,2}(1,4)+1, d2{1,3}(1,4)+1,...
                     d2{2,1}(1,2), d2{2,2}(1,4)+1, d2{2,3}(1,4)+1,...
                     d2{3,1}(1,2), d2{3,2}(1,4)+1, d2{3,3}(1,4)+1]);
D.Rx2{2}(2,2) = max([dx{1}+d2{2,1}(2,2), dx{3}(2)+d2{1,1}(1,2), dx{2}(3)+d2{2,1}(2,2)+1,...
                        dx{3}(3)+d2{2,1}(2,2)+1, dx{3}(2)+d2{3,1}(1,2),...
                     dx{1}+d2{2,2}(2,4)+1, dx{3}(2)+d2{1,2}(1,4)+1, dx{2}(3)+d2{2,2}(2,4)+2,...
                        dx{3}(3)+d2{2,2}(2,4)+2, dx{3}(2)+d2{3,2}(1,4)+1,...
                     dx{1}+d2{2,3}(2,4)+1, dx{3}(2)+d2{1,3}(1,4)+1, dx{2}(3)+d2{2,3}(2,4)+2,...
                        dx{3}(3)+d2{2,3}(2,4)+2, dx{3}(2)+d2{3,3}(1,4)+1]);
D.Rx2{2}(3,2) = max([d2{2,1}(3,2), dx{3}(1)+d2{1,1}(2,2), d2{2,1}(3,2),...
                        dx{3}(1)+d2{2,1}(4,2)+1, dx{3}(1)+d2{3,1}(4,2)+1,...
                     d2{2,2}(3,4)+1, dx{3}(1)+d2{1,2}(2,4)+1, d2{2,2}(3,4)+1,...
                        dx{3}(1)+d2{2,2}(4,4)+2, dx{3}(1)+d2{3,2}(4,4)+2,...
                     d2{2,3}(3,4)+1, dx{3}(1)+d2{1,3}(2,4)+1, d2{2,3}(3,4)+1,...
                        dx{3}(1)+d2{2,3}(4,4)+2, dx{3}(1)+d2{3,3}(4,4)+2]);
D.Rx2{2}(4,2) = max([dx{1}+d2{2,1}(4,2), dx{3}(3)+d2{1,1}(2,2), dx{2}(3)+d2{2,1}(4,2)+1,...
                        dx{3}(3)+d2{2,1}(4,2)+1, dx{3}(3)+d2{3,1}(4,2)+1,...
                     dx{1}+d2{2,2}(4,4)+1, dx{3}(3)+d2{1,2}(2,4)+1, dx{2}(3)+d2{2,2}(4,4)+2,...
                        dx{3}(3)+d2{2,2}(4,4)+2, dx{3}(3)+d2{3,2}(4,4)+2,...
                     dx{1}+d2{2,3}(4,4)+1, dx{3}(3)+d2{1,3}(2,4)+1, dx{2}(3)+d2{2,3}(4,4)+2,...
                        dx{3}(3)+d2{2,3}(4,4)+2, dx{3}(3)+d2{3,3}(4,4)+2]);
D.Rx2{2}(:,4) = D.Rx2{2}(:,2);  % We can omit this step

%
D.Rx2{3}(2,1) = max([dx{1}+d2{3,1}(2,1), dx{2}(2), dx{2}(2),...
                        dx{2}(3)+d2{3,1}(2,1)+1, dx{3}(3)+d2{3,1}(2,1)+1,...
                     dx{1}+d2{3,2}(2,1), dx{2}(2), dx{2}(2),...
                        dx{2}(3)+d2{3,2}(2,1)+1, dx{3}(3)+d2{3,2}(2,1)+1,...   
                     dx{1}+d2{3,3}(2,1), dx{2}(2), dx{2}(2),...
                        dx{2}(3)+d2{3,3}(2,1)+1, dx{3}(3)+d2{3,3}(2,1)+1]);
D.Rx2{3}(3,1) = max([d2{3,1}(3,1), dx{2}(1)+d2{1,1}(2,1), dx{2}(1)+d2{2,1}(4,1)+1,...
                        dx{2}(1)+d2{3,1}(4,1)+1, d2{3,1}(3,1),...
                     d2{3,2}(3,1), dx{2}(1)+d2{1,2}(2,1), dx{2}(1)+d2{2,2}(4,1)+1,...
                        dx{2}(1)+d2{3,2}(4,1)+1, d2{3,2}(3,1),...
                     d2{3,3}(3,1), dx{2}(1)+d2{1,3}(2,1), dx{2}(1)+d2{2,3}(4,1)+1,...
                        dx{2}(1)+d2{3,3}(4,1)+1, d2{3,3}(3,1)]);
D.Rx2{3}(4,1) = max([dx{1}+d2{3,1}(4,1), dx{2}(3)+d2{1,1}(2,1), dx{2}(3)+d2{2,1}(4,1)+1,...
                        dx{2}(3)+d2{3,1}(4,1)+1, dx{3}(3)+d2{3,1}(4,1)+1,...
                     dx{1}+d2{3,2}(4,1), dx{2}(3)+d2{1,2}(2,1), dx{2}(3)+d2{2,2}(4,1)+1,...
                        dx{2}(3)+d2{3,2}(4,1)+1, dx{3}(3)+d2{3,2}(4,1)+1,...   
                     dx{1}+d2{3,3}(4,1), dx{2}(3)+d2{1,3}(2,1), dx{2}(3)+d2{2,3}(4,1)+1,...
                        dx{2}(3)+d2{3,3}(4,1)+1, dx{3}(3)+d2{3,3}(4,1)+1]);
D.Rx2{2}(1,2) = max([d2{1,1}(1,2), d2{1,2}(1,4)+1, d2{1,3}(1,4)+1,...
                     d2{2,1}(1,2), d2{2,2}(1,4)+1, d2{2,3}(1,4)+1,...
                     d2{3,1}(1,2), d2{3,2}(1,4)+1, d2{3,3}(1,4)+1]);                    
D.Rx2{3}(2,2) = max([dx{1}+d2{3,1}(2,2), dx{2}(2)+d2{1,1}(1,2), dx{2}(2)+d2{2,1}(1,2),...
                        dx{2}(3)+d2{3,1}(2,2)+1, dx{3}(3)+d2{3,1}(2,2)+1,...
                     dx{1}+d2{3,2}(2,4)+1, dx{2}(2)+d2{1,2}(1,4)+1, dx{2}(2)+d2{2,2}(1,4)+1,...
                        dx{2}(3)+d2{3,2}(2,4)+2, dx{3}(3)+d2{3,2}(2,4)+2,... 
                     dx{1}+d2{3,3}(2,4)+1, dx{2}(2)+d2{1,3}(1,4)+1, dx{2}(2)+d2{2,3}(1,4)+1,...
                        dx{2}(3)+d2{3,3}(2,4)+2, dx{3}(3)+d2{3,3}(2,4)+2]);
D.Rx2{3}(3,2) = max([d2{3,1}(3,2), dx{2}(1)+d2{1,1}(2,2), dx{2}(1)+d2{2,1}(4,2)+1,...
                        dx{2}(1)+d2{3,1}(4,2)+1, d2{3,2}(3,1),...
                     d2{3,2}(3,4)+1, dx{2}(1)+d2{1,2}(2,4)+1, dx{2}(1)+d2{2,2}(4,4)+2,...
                        dx{2}(1)+d2{3,2}(4,4)+2, d2{3,2}(3,4)+1,...
                     d2{3,3}(3,4)+1, dx{2}(1)+d2{1,3}(2,4)+1, dx{2}(1)+d2{2,3}(4,4)+2,...
                        dx{2}(1)+d2{3,3}(4,4)+2, d2{3,3}(3,4)+1]);
D.Rx2{3}(4,2) = max([dx{1}+d2{3,1}(4,2), dx{2}(3)+d2{1,1}(2,2), dx{2}(3)+d2{2,1}(4,2)+1,...
                        dx{2}(3)+d2{3,1}(4,2)+1, dx{3}(3)+d2{3,1}(4,2)+1,...
                     dx{1}+d2{3,2}(4,4)+1, dx{2}(3)+d2{1,2}(2,4)+1, dx{2}(3)+d2{2,2}(4,4)+2,...
                        dx{2}(3)+d2{3,2}(4,4)+2, dx{3}(3)+d2{3,2}(4,4)+2,...   
                     dx{1}+d2{3,3}(4,4)+1, dx{2}(3)+d2{1,3}(2,4)+1, dx{2}(3)+d2{2,3}(4,4)+2,...
                        dx{2}(3)+d2{3,3}(4,4)+2, dx{3}(3)+d2{3,3}(4,4)+2]);
D.Rx2{3}(:,4) = D.Rx2{3}(:,2);  % We can omit this step
end

% % % Ryy and Ry2
if comp_params(3,3)
D.Ryy{1}(1,2) = 2*dy{1};
D.Ryy{2}(1,2) = max([dy{1}+dy{2}(1), dy{3}(2), dy{2}(3)+dy{2}(1)+1,...
                        dy{2}(1)+dy{3}(3)+1, dy{3}(2)]);
D.Ryy{2}(1,3) = max([dy{2}(2), dy{1}+dy{3}(1), dy{2}(2),...
                        dy{2}(3)+dy{3}(1)+1, dy{3}(3)+dy{3}(1)+1]);
D.Ryy{2}(1,4) = max([dy{1}+dy{2}(3), dy{1}+dy{3}(3), 2*dy{2}(3)+1,...
                        dy{2}(3)+dy{3}(3)+1, 2*dy{3}(3)+1]);
D.Ryy{3} = D.Ryy{2};
D.Ryy{3}(:,[2,3]) = D.Ryy{3}(:,[3,2]);  % Ryyb(y,nu) = Ryya^T(nu,y)   
end

% %
if comp_params(3,4) || comp_params(4,3)
D.Ry2{1}(1,2) = max([dy{1}+d2{1,1}(1,2), dy{1}+d2{2,1}(1,2), dy{1}+d2{3,1}(1,2)]);
D.Ry2{1}(2,1) = max([d2{1,1}(2,1), d2{2,1}(4,1)+1, d2{3,1}(4,1)]);
D.Ry2{1}(2,2) = max([dy{1}+d2{1,1}(2,2), dy{1}+d2{2,1}(4,2)+1, dy{1}+d2{3,1}(4,2)]);
D.Ry2{1}(2,4) = D.Ry2{1}(2,2);  % We can omit this step

%
D.Ry2{2}(1,2) = max([dy{1}+d2{1,2}(1,2), dy{3}(2), dy{2}(3)+d2{1,2}(1,2)+1,...
                        dy{3}(3)+d2{1,2}(1,2)+1, dy{3}(2),...
                     dy{1}+d2{2,2}(1,2), dy{3}(2), dy{2}(3)+d2{2,2}(1,2)+1,...
                        dy{3}(3)+d2{2,2}(1,2)+1, dy{3}(2),...
                     dy{1}+d2{3,2}(1,2), dy{3}(2), dy{2}(3)+d2{3,2}(1,2)+1,...
                        dy{3}(3)+d2{3,2}(1,2)+1, dy{3}(2)]);
D.Ry2{2}(1,3) = max([d2{1,2}(1,3), dy{3}(1)+d2{1,1}(1,2), d2{1,2}(1,3),...
                        dy{3}(1)+d2{1,2}(1,4)+1, dy{3}(1)+d2{1,3}(1,4)+1,...
                     d2{2,2}(1,3), dy{3}(1)+d2{2,1}(1,2), d2{2,2}(1,3),...
                        dy{3}(1)+d2{2,2}(1,4)+1, dy{3}(1)+d2{2,3}(1,4)+1,...
                     d2{3,2}(1,3), dy{3}(1)+d2{3,1}(1,2), d2{3,2}(1,3),...
                        dy{3}(1)+d2{3,2}(1,4)+1, dy{3}(1)+d2{3,3}(1,4)+1]);
D.Ry2{2}(1,4) = max([dy{1}+d2{1,2}(1,4), dy{3}(3)+d2{1,1}(1,2), dy{2}(3)+d2{1,2}(1,4)+1,...
                        dy{3}(3)+d2{1,2}(1,4)+1, dy{3}(3)+d2{1,3}(1,4)+1,...
                     dy{1}+d2{2,2}(1,4), dy{3}(3)+d2{2,1}(1,2), dy{2}(3)+d2{2,2}(1,4)+1,...
                        dy{3}(3)+d2{2,2}(1,4)+1, dy{3}(3)+d2{2,3}(1,4)+1,...
                     dy{1}+d2{3,2}(1,4), dy{3}(3)+d2{3,1}(2,1), dy{2}(3)+d2{3,2}(1,4)+1,...
                        dy{3}(3)+d2{3,2}(1,4)+1, dy{3}(3)+d2{3,3}(1,4)+1]);
D.Ry2{2}(2,1) = max([d2{1,1}(2,1), d2{1,2}(2,1), d2{1,3}(2,1),...
                     d2{2,1}(4,1)+1, d2{2,2}(4,1)+1, d2{2,3}(4,1)+1,...
                     d2{3,1}(4,1)+1, d2{3,2}(4,1)+1, d2{3,3}(4,1)+1]);
D.Ry2{2}(2,2) = max([dy{1}+d2{1,2}(2,2), dy{3}(2)+d2{1,1}(2,1), dy{2}(3)+d2{1,2}(2,2)+1,...
                        dy{3}(3)+d2{1,2}(2,2)+1, dy{3}(2)+d2{1,3}(2,1),...
                     dy{1}+d2{2,2}(4,2)+1, dy{3}(2)+d2{2,1}(4,1)+1, dy{2}(3)+d2{2,2}(4,2)+2,...
                        dy{3}(3)+d2{2,2}(4,2)+2, dy{3}(2)+d2{2,3}(4,1)+1,...
                     dy{1}+d2{3,2}(4,2)+1, dy{3}(2)+d2{3,1}(4,1)+1, dy{2}(3)+d2{3,2}(4,2)+2,...
                        dy{3}(3)+d2{3,2}(4,2)+2, dy{3}(2)+d2{3,3}(4,1)+1]);
D.Ry2{2}(2,3) = max([d2{1,2}(2,3), dy{3}(1)+d2{1,1}(2,2), d2{1,2}(2,3),...
                        dy{3}(1)+d2{1,2}(2,4)+1, dy{3}(1)+d2{1,3}(2,4)+1,...
                     d2{2,2}(4,3)+1, dy{3}(1)+d2{2,1}(4,2)+1, d2{2,2}(4,3)+1,...
                        dy{3}(1)+d2{2,2}(4,4)+2, dy{3}(1)+d2{2,3}(4,4)+2,...
                     d2{3,2}(4,3)+1, dy{3}(1)+d2{3,1}(4,2)+1, d2{3,2}(4,3)+1,...
                        dy{3}(1)+d2{3,2}(4,4)+2, dy{3}(1)+d2{3,3}(4,4)+2]);
D.Ry2{2}(2,4) = max([dy{1}+d2{1,2}(2,4), dy{3}(3)+d2{1,1}(2,2), dy{2}(3)+d2{1,2}(2,4)+1,...
                        dy{3}(3)+d2{1,2}(2,4)+1, dy{3}(3)+d2{1,3}(2,4)+1,...
                     dy{1}+d2{2,2}(4,4)+1, dy{3}(3)+d2{2,1}(4,2)+1, dy{2}(3)+d2{2,2}(4,4)+2,...
                        dy{3}(3)+d2{2,2}(4,4)+2, dy{3}(3)+d2{2,3}(4,4)+2,...
                     dy{1}+d2{3,2}(4,4)+1, dy{3}(3)+d2{3,1}(4,2)+1, dy{2}(3)+d2{3,2}(4,4)+2,...
                        dy{3}(3)+d2{3,2}(4,4)+2, dy{3}(3)+d2{3,3}(4,4)+2]);
D.Ry2{2}(4,:) = D.Ry2{2}(2,:);  % We can omit this step

%
D.Ry2{3}(1,2) = max([dy{1}+d2{1,3}(1,2), dy{2}(2), dy{2}(2),...
                        dy{2}(3)+d2{1,3}(1,2)+1, dy{3}(3)+d2{1,3}(1,2)+1,...
                     dy{1}+d2{2,3}(1,2), dy{2}(2), dy{2}(2),...
                        dy{2}(3)+d2{2,3}(1,2)+1, dy{3}(3)+d2{2,3}(1,2)+1,...   
                     dy{1}+d2{3,3}(1,2), dy{2}(2), dy{2}(2),...
                        dy{2}(3)+d2{3,3}(1,2)+1, dy{3}(3)+d2{3,3}(1,2)+1]);
D.Ry2{3}(1,3) = max([d2{1,3}(1,3), dy{2}(1)+d2{1,1}(1,2), dy{2}(1)+d2{1,2}(1,4)+1,...
                        dy{2}(1)+d2{1,3}(1,4)+1, d2{1,3}(1,3),...
                     d2{2,3}(1,3), dy{2}(1)+d2{2,1}(1,2), dy{2}(1)+d2{2,2}(1,4)+1,...
                        dy{2}(1)+d2{2,3}(1,4)+1, d2{2,3}(1,3),...
                     d2{3,3}(1,3), dy{2}(1)+d2{3,1}(1,2), dy{2}(1)+d2{3,2}(1,4)+1,...
                        dy{2}(1)+d2{3,3}(1,4)+1, d2{3,3}(1,3)]);
D.Ry2{3}(1,4) = max([dy{1}+d2{1,3}(1,4), dy{2}(3)+d2{1,1}(1,2), dy{2}(3)+d2{1,2}(1,4)+1,...
                        dy{2}(3)+d2{1,3}(1,4)+1, dy{3}(3)+d2{1,3}(1,4)+1,...
                     dy{1}+d2{2,3}(1,4), dy{2}(3)+d2{2,1}(1,2), dy{2}(3)+d2{2,2}(1,4)+1,...
                        dy{2}(3)+d2{2,3}(1,4)+1, dy{3}(3)+d2{2,3}(1,4)+1,...   
                     dy{1}+d2{3,3}(1,4), dy{2}(3)+d2{3,1}(1,2), dy{2}(3)+d2{3,2}(1,4)+1,...
                        dy{2}(3)+d2{3,3}(1,4)+1, dy{3}(3)+d2{3,3}(1,4)+1]);
D.Ry2{2}(2,1) = max([d2{1,1}(2,1), d2{1,2}(2,1), d2{1,3}(2,1),...
                     d2{2,1}(4,1)+1, d2{2,2}(4,1)+1, d2{2,3}(4,1)+1,...
                     d2{3,1}(4,1)+1, d2{3,2}(4,1)+1, d2{3,3}(4,1)+1]);                    
D.Ry2{3}(2,2) = max([dy{1}+d2{1,3}(2,2), dy{2}(2)+d2{1,1}(2,1), dy{2}(2)+d2{1,2}(2,1),...
                        dy{2}(3)+d2{1,3}(2,2)+1, dy{3}(3)+d2{1,3}(2,2)+1,...
                     dy{1}+d2{2,3}(4,2)+1, dy{2}(2)+d2{2,1}(4,1)+1, dy{2}(2)+d2{2,2}(4,1)+1,...
                        dy{2}(3)+d2{2,3}(4,2)+2, dy{3}(3)+d2{2,3}(4,2)+2,... 
                     dy{1}+d2{3,3}(4,2)+1, dy{2}(2)+d2{3,1}(4,1)+1, dy{2}(2)+d2{3,2}(4,1)+1,...
                        dy{2}(3)+d2{3,3}(4,2)+2, dy{3}(3)+d2{3,3}(4,2)+2]);
D.Ry2{3}(2,3) = max([d2{1,3}(2,3), dy{2}(1)+d2{1,1}(2,2), dy{2}(1)+d2{1,2}(2,4)+1,...
                        dy{2}(1)+d2{1,3}(2,4)+1, d2{2,3}(1,3),...
                     d2{2,3}(4,3)+1, dy{2}(1)+d2{2,1}(4,2)+1, dy{2}(1)+d2{2,2}(4,4)+2,...
                        dy{2}(1)+d2{2,3}(4,4)+2, d2{2,3}(4,3)+1,...
                     d2{3,3}(4,3)+1, dy{2}(1)+d2{3,1}(4,2)+1, dy{2}(1)+d2{3,2}(4,4)+2,...
                        dy{2}(1)+d2{3,3}(4,4)+2, d2{3,3}(4,3)+1]);
D.Ry2{3}(2,4) = max([dy{1}+d2{1,3}(2,4), dy{2}(3)+d2{1,1}(2,2), dy{2}(3)+d2{1,2}(2,4)+1,...
                        dy{2}(3)+d2{1,3}(2,4)+1, dy{3}(3)+d2{1,3}(2,4)+1,...
                     dy{1}+d2{2,3}(4,4)+1, dy{2}(3)+d2{2,1}(4,2)+1, dy{2}(3)+d2{2,2}(4,4)+2,...
                        dy{2}(3)+d2{2,3}(4,4)+2, dy{3}(3)+d2{2,3}(4,4)+2,...   
                     dy{1}+d2{3,3}(4,4)+1, dy{2}(3)+d2{3,1}(4,2)+1, dy{2}(3)+d2{3,2}(4,4)+2,...
                        dy{2}(3)+d2{3,3}(4,4)+2, dy{3}(3)+d2{3,3}(4,4)+2]);
D.Ry2{3}(4,:) = D.Ry2{3}(2,:);  % We can omit this step
end

% % % Finally, R22
if comp_params(4,4)
D.R22{1,1}(1:2,1:2) = 2*d2{1,1}(1:2,1:2);

%
D.R22{2,1}(2,1) = max([d2{1,1}(2,1)+d2{2,1}(2,1), d2{3,1}(3,1),...
                        d2{2,1}(4,1)+d2{2,1}(2,1)+1, d2{3,1}(4,1)+d2{2,1}(2,1)+1, d2{3,1}(3,1)]);
D.R22{2,1}(3,1) = max([d2{2,1}(3,1), d2{3,1}(2,1)+d2{1,1}(2,1),...
                        d2{2,1}(3,1), d2{3,1}(2,1)+d2{2,1}(4,1)+1, d2{3,1}(2,1)+d2{3,1}(4,1)+1]);
D.R22{2,1}(4,1) = max([d2{1,1}(2,1)+d2{2,1}(4,1), d2{3,1}(4,1)+d2{1,1}(2,1),...
                        2*d2{2,1}(4,1)+1, d2{3,1}(4,1)+d2{2,1}(4,1)+1, 2*d2{3,1}(4,1)+1]);
D.R22{2,1}(1,2) = max([d2{1,1}(1,2)+d2{2,1}(1,2), d2{3,1}(1,2)+d2{1,1}(1,2),...
                        2*d2{2,1}(1,2), d2{3,1}(1,2)+d2{2,1}(1,2), 2*d2{3,1}(1,2)]);
D.R22{2,1}(2,2) = max([d2{1,1}(2,2)+d2{2,1}(2,2), d2{3,1}(3,2)+d2{1,1}(1,2),...
                        d2{2,1}(4,2)+d2{2,1}(2,2)+1, d2{3,1}(4,2)+d2{2,2}(2,2)+1, d2{3,1}(3,2)+d2{3,1}(1,2)]);
D.R22{2,1}(3,2) = max([d2{1,1}(1,2)+d2{2,1}(3,2), d2{3,1}(2,2)+d2{1,1}(2,2),...
                        d2{2,1}(1,2)+d2{2,1}(3,2), d2{3,1}(2,2)+d2{2,1}(4,2)+1, d2{3,1}(2,2)+d2{3,1}(4,2)+1]);
D.R22{2,1}(4,2) = max([d2{1,1}(2,2)+d2{2,1}(4,2), d2{3,1}(4,2)+d2{1,1}(2,2),...
                        2*d2{2,1}(4,2)+1, d2{3,1}(4,2)+d2{2,2}(4,2)+1, 2*d2{3,1}(4,2)]);                  
D.R22{3,1} = D.R22{2,1};
D.R22{3,1}([2;3],:) = D.R22{3,1}([3;2],:);  % R22bo(x,y,tt) = R22ao^T(tt,y,x)

%
D.R22{1,2}(1,2) = max([d2{1,1}(1,2)+d2{1,2}(1,2), d2{1,3}(1,3),...
                        d2{1,2}(1,4)+d2{1,2}(1,2)+1, d2{1,3}(1,4)+d2{1,2}(1,2)+1, d2{1,3}(1,3)]);
D.R22{1,2}(1,3) = max([d2{1,2}(1,3), d2{1,3}(1,2)+d2{1,1}(1,2),...
                        d2{1,2}(1,3), d2{1,3}(1,2)+d2{1,2}(1,4)+1, d2{1,3}(1,2)+d2{1,3}(1,4)+1]);
D.R22{1,2}(1,4) = max([d2{1,1}(1,2)+d2{1,2}(1,4), d2{1,3}(1,4)+d2{1,1}(1,2),...
                        2*d2{1,2}(1,4)+1, d2{1,3}(1,4)+d2{1,2}(1,4)+1, 2*d2{1,3}(1,4)+1]);
D.R22{1,2}(2,1) = max([d2{1,1}(2,1)+d2{1,2}(2,1), d2{1,3}(2,1)+d2{1,1}(2,1),...
                        2*d2{1,2}(2,1), d2{1,3}(2,1)+d2{1,2}(2,1), 2*d2{1,3}(2,1)]);
D.R22{1,2}(2,2) = max([d2{1,1}(2,2)+d2{1,2}(2,2), d2{1,3}(2,3)+d2{1,1}(2,1),...
                        d2{1,2}(2,4)+d2{1,2}(2,2)+1, d2{1,3}(2,4)+d2{2,2}(2,2)+1, d2{1,3}(2,3)+d2{1,3}(2,1)]);
D.R22{1,2}(2,3) = max([d2{1,1}(2,1)+d2{1,2}(2,3), d2{1,3}(2,2)+d2{1,1}(2,2),...
                        d2{1,2}(2,1)+d2{1,2}(2,3), d2{1,3}(2,2)+d2{1,2}(2,4)+1, d2{1,3}(2,2)+d2{1,3}(2,4)+1]);
D.R22{1,2}(2,4) = max([d2{1,1}(2,2)+d2{1,2}(2,4), d2{1,3}(2,4)+d2{1,1}(2,2),...
                        2*d2{1,2}(2,4)+1, d2{1,3}(2,4)+d2{2,2}(2,4)+1, 2*d2{1,3}(2,4)]);                  
D.R22{1,3} = D.R22{1,2};
D.R22{1,3}(:,[2,3]) = D.R22{1,3}(:,[3,2]);  % R22ob(x,y,nu) = R22ob^T(x,nu,y)

%
D.R22{2,2}(2,1) = max([d2{1,1}(2,1)+d2{2,2}(2,1), d2{3,1}(3,1),...
                        d2{2,1}(4,1)+d2{2,2}(2,1)+1, d2{3,1}(4,1)+d2{2,2}(2,1)+1, d2{3,1}(3,1),...
                       d2{1,3}(2,1)+d2{2,1}(2,1), d2{3,3}(3,1),...
                        d2{2,3}(4,1)+d2{2,1}(2,1)+1, d2{3,3}(4,1)+d2{2,1}(2,1)+1, d2{3,3}(3,1),...
                       d2{1,2}(2,1)+d2{2,2}(2,1), d2{3,2}(3,1),...
                        d2{2,2}(4,1)+d2{2,2}(2,1)+1, d2{3,2}(4,1)+d2{2,2}(2,1)+1, d2{3,2}(3,1),...
                       d2{1,3}(2,1)+d2{2,2}(2,1), d2{3,3}(3,1),...
                        d2{2,3}(4,1)+d2{2,2}(2,1)+1, d2{3,3}(4,1)+d2{2,2}(2,1)+1, d2{3,3}(3,1),...
                       d2{1,3}(2,1)+d2{2,3}(2,1), d2{3,3}(3,1),...
                        d2{2,3}(4,1)+d2{2,3}(2,1)+1, d2{3,3}(4,1)+d2{2,3}(2,1)+1, d2{3,3}(3,1)]);
D.R22{2,2}(3,1) = max([d2{2,2}(3,1), d2{3,1}(2,1)+d2{1,2}(2,1),...
                        d2{2,2}(3,1), d2{3,1}(2,1)+d2{2,2}(4,1)+1, d2{3,1}(2,1)+d2{3,2}(4,1)+1,...
                       d2{2,1}(3,1), d2{3,3}(2,1)+d2{1,1}(2,1),...
                        d2{2,1}(3,1), d2{3,3}(2,1)+d2{2,1}(4,1)+1, d2{3,3}(2,1)+d2{3,1}(4,1)+1,...
                       d2{2,2}(3,1), d2{3,2}(2,1)+d2{1,2}(2,1),...
                        d2{2,2}(3,1), d2{3,2}(2,1)+d2{2,2}(4,1)+1, d2{3,2}(2,1)+d2{3,2}(4,1)+1,...
                       d2{2,2}(3,1), d2{3,3}(2,1)+d2{1,2}(2,1),...
                        d2{2,2}(3,1), d2{3,3}(2,1)+d2{2,2}(4,1)+1, d2{3,3}(2,1)+d2{3,2}(4,1)+1,...
                       d2{2,3}(3,1), d2{3,3}(2,1)+d2{1,3}(2,1),...
                        d2{2,3}(3,1), d2{3,3}(2,1)+d2{2,3}(4,1)+1, d2{3,3}(2,1)+d2{3,3}(4,1)+1]);
D.R22{2,2}(4,1) = max([d2{1,1}(2,1)+d2{2,2}(4,1), d2{3,1}(4,1)+d2{1,2}(2,1),...
                        d2{2,1}(4,1)+d2{2,2}(4,1)+1, d2{3,1}(4,1)+d2{2,2}(4,1)+1, d2{3,1}(4,1)+d2{3,2}(4,1)+1,...
                       d2{1,3}(2,1)+d2{2,1}(4,1), d2{3,3}(4,1)+d2{1,1}(2,1),...
                        d2{2,3}(4,1)+d2{2,1}(4,1)+1, d2{3,3}(4,1)+d2{2,1}(4,1)+1, d2{3,3}(4,1)+d2{3,1}(4,1)+1,...
                       d2{1,2}(2,1)+d2{2,2}(4,1), d2{3,2}(4,1)+d2{1,2}(2,1),...
                        d2{2,2}(4,1)+d2{2,2}(4,1)+1, d2{3,2}(4,1)+d2{2,2}(4,1)+1, d2{3,2}(4,1)+d2{3,2}(4,1)+1,...
                       d2{1,3}(2,1)+d2{2,2}(4,1), d2{3,3}(4,1)+d2{1,2}(2,1),...
                        d2{2,3}(4,1)+d2{2,2}(4,1)+1, d2{3,3}(4,1)+d2{2,2}(4,1)+1, d2{3,3}(4,1)+d2{3,2}(4,1)+1,...
                       d2{1,3}(2,1)+d2{2,3}(4,1), d2{3,3}(4,1)+d2{1,3}(2,1),...
                        d2{2,3}(4,1)+d2{2,3}(4,1)+1, d2{3,3}(4,1)+d2{2,3}(4,1)+1, d2{3,3}(4,1)+d2{3,3}(4,1)+1]);
D.R22{2,2}(1,2) = max([d2{1,1}(1,2)+d2{2,2}(1,2), d2{1,3}(1,3),...
                        d2{1,2}(1,4)+d2{2,2}(1,2)+1, d2{1,3}(1,4)+d2{2,2}(1,2)+1, d2{1,3}(1,3),...
                       d2{3,1}(1,2)+d2{1,2}(1,2), d2{3,3}(1,3),...
                        d2{3,2}(1,4)+d2{1,2}(1,2)+1, d2{3,3}(1,4)+d2{1,2}(1,2)+1, d2{3,3}(1,3),...
                       d2{2,1}(1,2)+d2{2,2}(1,2), d2{2,3}(1,3),...
                        d2{2,2}(1,4)+d2{2,2}(1,2)+1, d2{2,3}(1,4)+d2{2,2}(1,2)+1, d2{2,3}(1,3),...
                       d2{3,1}(1,2)+d2{2,2}(1,2), d2{3,3}(1,3),...
                        d2{3,2}(1,4)+d2{2,2}(1,2)+1, d2{3,3}(1,4)+d2{2,2}(1,2)+1, d2{3,3}(1,3),...
                       d2{3,1}(1,2)+d2{3,2}(1,2), d2{3,3}(1,3),...
                        d2{3,2}(1,4)+d2{3,2}(1,2)+1, d2{3,3}(1,4)+d2{3,2}(1,2)+1, d2{3,3}(1,3)]);
D.R22{2,2}(1,3) = max([d2{2,2}(1,3), d2{1,3}(1,2)+d2{2,1}(1,2),...
                        d2{2,2}(1,3), d2{1,3}(1,2)+d2{2,2}(1,4)+1, d2{1,3}(1,2)+d2{2,3}(1,4)+1,...
                       d2{1,2}(1,3), d2{3,3}(1,2)+d2{1,1}(1,2),...
                        d2{1,2}(1,3), d2{3,3}(1,2)+d2{1,2}(1,4)+1, d2{3,3}(1,2)+d2{1,3}(1,4)+1,...
                       d2{2,2}(1,3), d2{2,3}(1,2)+d2{2,1}(1,2),...
                        d2{2,2}(1,3), d2{2,3}(1,2)+d2{2,2}(1,4)+1, d2{2,3}(1,2)+d2{2,3}(1,4)+1,...
                       d2{2,2}(1,3), d2{3,3}(1,2)+d2{2,1}(1,2),...
                        d2{2,2}(1,3), d2{3,3}(1,2)+d2{2,2}(1,4)+1, d2{3,3}(1,2)+d2{2,3}(1,4)+1,...
                       d2{3,2}(1,3), d2{3,3}(1,2)+d2{3,1}(1,2),...
                        d2{3,2}(1,3), d2{3,3}(1,2)+d2{3,2}(1,4)+1, d2{3,3}(1,2)+d2{3,3}(1,4)+1]);
D.R22{2,2}(1,4) = max([d2{1,1}(1,2)+d2{2,2}(1,4), d2{1,3}(1,4)+d2{2,1}(1,2),...
                        d2{1,2}(1,4)+d2{2,2}(1,4)+1, d2{1,3}(1,4)+d2{2,2}(1,4)+1, d2{1,3}(1,4)+d2{2,3}(1,4)+1,...
                       d2{3,1}(1,2)+d2{1,2}(1,4), d2{3,3}(1,4)+d2{1,1}(1,2),...
                        d2{3,2}(1,4)+d2{1,2}(1,4)+1, d2{3,3}(1,4)+d2{1,2}(1,4)+1, d2{3,3}(1,4)+d2{1,3}(1,4)+1,...
                       d2{2,1}(1,2)+d2{2,2}(1,4), d2{2,3}(1,4)+d2{2,1}(1,2),...
                        d2{2,2}(1,4)+d2{2,2}(1,4)+1, d2{2,3}(1,4)+d2{2,2}(1,4)+1, d2{2,3}(1,4)+d2{2,3}(1,4)+1,...
                       d2{3,1}(1,2)+d2{2,2}(1,4), d2{3,3}(1,4)+d2{2,1}(1,2),...
                        d2{3,2}(1,4)+d2{2,2}(1,4)+1, d2{3,3}(1,4)+d2{2,2}(1,4)+1, d2{3,3}(1,4)+d2{2,3}(1,4)+1,...
                       d2{3,1}(1,2)+d2{3,2}(1,4), d2{3,3}(1,4)+d2{3,1}(1,2),...
                        d2{3,2}(1,4)+d2{3,2}(1,4)+1, d2{3,3}(1,4)+d2{3,2}(1,4)+1, d2{3,3}(1,4)+d2{3,3}(1,4)+1]);
D.R22{2,2}(2,2) = max([d2{1,1}(2,2)+d2{2,2}(2,2), d2{3,1}(3,2)+d2{1,2}(1,2),...
                        d2{2,1}(4,2)+d2{2,2}(2,2)+1, d2{3,1}(4,2)+d2{2,2}(2,2)+1, d2{3,1}(3,2)+d2{3,2}(1,2),...
                       d2{1,3}(2,3)+d2{2,1}(2,1), d2{3,3}(3,1),...
                        d2{2,3}(4,3)+d2{2,1}(2,1)+1, d2{3,3}(4,3)+d2{2,1}(2,1)+1, d2{3,3}(3,3),...
                       d2{1,2}(2,4)+d2{2,2}(2,2)+1, d2{3,2}(3,4)+d2{1,2}(1,2)+1,...
                        d2{2,2}(4,4)+d2{2,2}(2,2)+2, d2{3,2}(4,4)+d2{2,2}(2,2)+2, d2{3,2}(3,4)+d2{3,2}(1,2)+1,...
                       d2{1,3}(2,4)+d2{2,2}(2,2)+1, d2{3,3}(3,4)+d2{1,2}(1,2)+1,...
                        d2{2,3}(4,4)+d2{2,2}(2,2)+2, d2{3,3}(4,4)+d2{2,2}(2,2)+2, d2{3,3}(3,4)+d2{3,2}(1,2)+1,...
                       d2{1,3}(2,3)+d2{2,3}(2,1), d2{3,3}(3,3),...
                        d2{2,3}(4,3)+d2{2,3}(2,1)+1, d2{3,3}(4,3)+d2{2,3}(2,1)+1, d2{3,3}(3,3)]);
D.R22{2,2}(3,2) = max([d2{1,1}(1,2)+d2{2,2}(3,2), d2{3,1}(2,2)+d2{1,2}(2,2),...
                        d2{2,1}(1,2)+d2{2,2}(3,2), d2{3,1}(2,2)+d2{2,2}(4,2)+1, d2{3,1}(2,2)+d2{3,2}(4,2)+1,...
                       d2{1,3}(1,3)+d2{2,1}(3,1), d2{3,3}(2,3)+d2{1,1}(2,1),...
                        d2{2,3}(1,3)+d2{2,1}(3,1), d2{3,3}(2,3)+d2{2,1}(4,1)+1, d2{3,3}(2,3)+d2{3,1}(4,1)+1,...
                       d2{1,2}(1,4)+d2{2,2}(3,2)+1, d2{3,2}(2,4)+d2{1,2}(2,2)+1,...
                        d2{2,2}(1,4)+d2{2,2}(3,2)+1, d2{3,2}(2,4)+d2{2,2}(4,2)+2, d2{3,2}(2,4)+d2{3,2}(4,2)+2,...
                       d2{1,3}(1,4)+d2{2,2}(3,2)+1, d2{3,3}(2,4)+d2{1,2}(2,2)+1,...
                        d2{2,3}(1,4)+d2{2,2}(3,2)+1, d2{3,3}(2,4)+d2{2,2}(4,2)+2, d2{3,3}(2,4)+d2{3,2}(4,2)+2,...
                       d2{1,3}(1,3)+d2{2,3}(3,1), d2{3,3}(2,3)+d2{1,3}(2,1),...
                        d2{2,3}(1,3)+d2{2,3}(3,1), d2{3,3}(2,3)+d2{2,3}(4,1)+1, d2{3,3}(2,3)+d2{3,3}(4,1)+1]);
D.R22{2,2}(4,2) = max([d2{1,1}(2,2)+d2{2,2}(3,2), d2{3,1}(4,2)+d2{1,2}(1,2),...
                        d2{2,1}(4,2)+d2{2,2}(4,2)+1, d2{3,1}(4,2)+d2{2,2}(4,2)+1, d2{3,1}(4,2)+d2{3,2}(4,2)+1,...
                       d2{1,3}(2,3)+d2{2,1}(4,1), d2{3,3}(4,1)+d2{1,1}(2,1),...
                        d2{2,3}(4,3)+d2{2,1}(4,1)+1, d2{3,3}(4,3)+d2{2,1}(4,1)+1, d2{3,3}(4,3)+d2{3,1}(4,1),...
                       d2{1,2}(2,4)+d2{2,2}(4,2)+1, d2{3,2}(4,4)+d2{1,2}(2,2)+1,...
                        d2{2,2}(4,4)+d2{2,2}(4,2)+2, d2{3,2}(4,4)+d2{2,2}(4,2)+2, d2{3,2}(4,4)+d2{3,2}(4,2)+2,...
                       d2{1,3}(2,4)+d2{2,2}(4,2)+1, d2{3,3}(4,4)+d2{1,2}(2,2)+1,...
                        d2{2,3}(4,4)+d2{2,2}(4,2)+2, d2{3,3}(4,4)+d2{2,2}(4,2)+2, d2{3,3}(4,4)+d2{3,2}(4,2)+2,...
                       d2{1,3}(2,3)+d2{2,3}(4,1), d2{3,3}(4,3)+d2{1,3}(2,1),...
                        d2{2,3}(4,3)+d2{2,3}(4,1)+1, d2{3,3}(4,3)+d2{2,3}(4,1)+1, d2{3,3}(4,3)+d2{3,3}(4,1)]);
D.R22{2,2}(2,3) = max([d2{1,1}(2,1)+d2{2,2}(2,3), d2{1,3}(2,2)+d2{2,1}(2,2),...
                        d2{1,2}(2,1)+d2{2,2}(2,3), d2{1,3}(2,2)+d2{2,2}(2,4)+1, d2{1,3}(2,2)+d2{2,3}(2,4)+1,...
                       d2{3,1}(3,1)+d2{1,2}(1,3), d2{3,3}(3,2)+d2{1,1}(1,2),...
                        d2{3,2}(3,1)+d2{1,2}(1,3), d2{3,3}(3,2)+d2{1,2}(1,4)+1, d2{3,3}(3,2)+d2{1,3}(1,4)+1,...
                       d2{2,1}(4,1)+d2{2,2}(2,3)+1, d2{2,3}(4,2)+d2{2,1}(2,2)+1,...
                        d2{2,2}(4,1)+d2{2,2}(2,3)+1, d2{2,3}(4,2)+d2{2,2}(2,4)+2, d2{2,3}(4,2)+d2{2,3}(2,4)+2,...
                       d2{3,1}(4,1)+d2{2,2}(2,3)+1, d2{3,3}(4,2)+d2{2,1}(2,2)+1,...
                        d2{3,2}(4,1)+d2{2,2}(2,3)+1, d2{3,3}(4,2)+d2{2,2}(2,4)+2, d2{3,3}(4,2)+d2{2,3}(2,4)+2,...
                       d2{3,1}(3,1)+d2{3,2}(1,3), d2{3,3}(3,2)+d2{3,1}(1,2),...
                        d2{3,2}(3,1)+d2{3,2}(1,3), d2{3,3}(3,2)+d2{3,2}(1,4)+1, d2{3,3}(3,2)+d2{3,3}(1,4)+1]);
D.R22{2,2}(2,4) = max([d2{1,1}(2,2)+d2{2,2}(2,3), d2{1,3}(2,4)+d2{2,1}(2,1),...
                        d2{1,2}(2,4)+d2{2,2}(2,4)+1, d2{1,3}(2,4)+d2{2,2}(2,4)+1, d2{1,3}(2,4)+d2{2,3}(2,4)+1,...
                       d2{3,1}(3,2)+d2{1,2}(1,4), d2{3,3}(1,4)+d2{1,1}(1,2),...
                        d2{3,2}(3,4)+d2{1,2}(1,4)+1, d2{3,3}(3,4)+d2{1,2}(1,4)+1, d2{3,3}(3,4)+d2{1,3}(1,4),...
                       d2{2,1}(4,2)+d2{2,2}(2,4)+1, d2{2,3}(4,4)+d2{2,1}(2,2)+1,...
                        d2{2,2}(4,4)+d2{2,2}(2,4)+2, d2{2,3}(4,4)+d2{2,2}(2,4)+2, d2{2,3}(4,4)+d2{2,3}(2,4)+2,...
                       d2{3,1}(4,2)+d2{2,2}(2,4)+1, d2{3,3}(4,4)+d2{2,1}(2,2)+1,...
                        d2{3,2}(4,4)+d2{2,2}(2,4)+2, d2{3,3}(4,4)+d2{2,2}(2,4)+2, d2{3,3}(4,4)+d2{2,3}(2,4)+2,...
                       d2{3,1}(3,2)+d2{3,2}(1,4), d2{3,3}(3,4)+d2{3,1}(1,2),...
                        d2{3,2}(3,4)+d2{3,2}(1,4)+1, d2{3,3}(3,4)+d2{3,2}(1,4)+1, d2{3,3}(3,4)+d2{3,3}(1,4)]);
D.R22{2,2}(3,3) = max([d2{2,2}(3,3), d2{3,1}(2,1)+d2{1,2}(2,3),...
                        d2{2,2}(3,3), d2{3,1}(2,1)+d2{2,2}(4,3)+1, d2{3,1}(2,1)+d2{3,2}(4,3)+1,...
                       d2{1,3}(1,2)+d2{2,1}(3,2), d2{3,3}(2,2)+d2{1,1}(2,2),...
                        d2{2,3}(1,2)+d2{2,1}(3,2), d2{3,3}(2,2)+d2{2,1}(4,2)+1, d2{3,3}(2,2)+d2{3,1}(4,2)+1,...
                       d2{2,2}(3,3), d2{3,2}(2,1)+d2{1,2}(2,3),...
                        d2{2,2}(3,3), d2{3,2}(2,1)+d2{2,2}(4,3)+1, d2{3,2}(2,1)+d2{3,2}(4,3)+1,...
                       d2{1,3}(1,2)+d2{2,2}(3,4)+1, d2{3,3}(2,2)+d2{1,2}(2,4)+1,...
                        d2{2,3}(1,2)+d2{2,2}(3,4)+1, d2{3,3}(2,2)+d2{2,2}(4,4)+2, d2{3,3}(2,2)+d2{3,2}(4,4)+2,...
                       d2{1,3}(1,2)+d2{2,3}(3,4)+1, d2{3,3}(2,2)+d2{1,3}(2,4)+1,...
                        d2{2,3}(1,2)+d2{2,3}(3,4)+1, d2{3,3}(2,2)+d2{2,3}(4,4)+2, d2{3,3}(2,2)+d2{3,3}(4,4)+2]);
D.R22{2,2}(4,3) = max([d2{1,1}(2,1)+d2{2,2}(4,3), d2{3,1}(4,1)+d2{1,2}(2,3),...
                        d2{2,1}(4,1)+d2{2,2}(4,3)+1, d2{3,1}(4,1)+d2{2,2}(4,3)+1, d2{3,1}(4,1)+d2{3,2}(4,3)+1,...
                       d2{1,3}(2,2)+d2{2,1}(4,2), d2{3,3}(4,2)+d2{1,1}(2,2),...
                        d2{2,3}(4,2)+d2{2,1}(4,2)+1, d2{3,3}(4,2)+d2{2,1}(4,2)+1, d2{3,3}(4,2)+d2{3,1}(4,2)+1,...
                       d2{1,2}(2,1)+d2{2,2}(4,3), d2{3,2}(4,1)+d2{1,2}(2,3),...
                        d2{2,2}(4,1)+d2{2,2}(4,3)+1, d2{3,2}(4,1)+d2{2,2}(4,3)+1, d2{3,2}(4,1)+d2{3,2}(4,3)+1,...
                       d2{1,3}(2,2)+d2{2,2}(4,4)+1, d2{3,3}(4,2)+d2{1,2}(2,4)+1,...
                        d2{2,3}(4,2)+d2{2,2}(4,4)+2, d2{3,3}(4,2)+d2{2,2}(4,4)+2, d2{3,3}(4,2)+d2{3,2}(4,4)+2,...
                       d2{1,3}(2,2)+d2{2,3}(4,4)+1, d2{3,3}(4,2)+d2{1,3}(2,4)+1,...
                        d2{2,3}(4,2)+d2{2,3}(4,4)+2, d2{3,3}(4,2)+d2{2,3}(4,4)+2, d2{3,3}(4,2)+d2{3,3}(4,4)+2]);
D.R22{2,2}(3,4) = max([d2{1,1}(1,2)+d2{2,2}(3,4), d2{1,3}(1,4)+d2{2,1}(3,2),...
                        d2{1,2}(1,4)+d2{2,2}(3,4)+1, d2{1,3}(1,4)+d2{2,2}(3,4)+1, d2{1,3}(1,4)+d2{2,3}(3,4)+1,...
                       d2{3,1}(2,2)+d2{1,2}(2,4), d2{3,3}(2,4)+d2{1,1}(2,2),...
                        d2{3,2}(2,4)+d2{1,2}(2,4)+1, d2{3,3}(2,4)+d2{1,2}(2,4)+1, d2{3,3}(2,4)+d2{1,3}(2,4)+1,...
                       d2{2,1}(1,2)+d2{2,2}(3,4), d2{2,3}(1,4)+d2{2,1}(3,2),...
                        d2{2,2}(1,4)+d2{2,2}(3,4)+1, d2{2,3}(1,4)+d2{2,2}(3,4)+1, d2{2,3}(1,4)+d2{2,3}(3,4)+1,...
                       d2{3,1}(2,2)+d2{2,2}(4,4)+1, d2{3,3}(2,4)+d2{2,1}(4,2)+1,...
                        d2{3,2}(2,4)+d2{2,2}(4,4)+2, d2{3,3}(2,4)+d2{2,2}(4,4)+2, d2{3,3}(2,4)+d2{2,3}(4,4)+2,...
                       d2{3,1}(2,2)+d2{3,2}(4,4)+1, d2{3,3}(2,4)+d2{3,1}(4,2)+1,...
                        d2{3,2}(2,4)+d2{3,2}(4,4)+2, d2{3,3}(2,4)+d2{3,2}(4,4)+2, d2{3,3}(2,4)+d2{3,3}(4,4)+2]);
D.R22{2,2}(4,4) = max([d2{1,1}(2,2)+d2{2,2}(4,4), d2{3,1}(4,2)+d2{1,2}(2,4),...
                        d2{2,1}(4,2)+d2{2,2}(4,4)+1, d2{3,1}(4,2)+d2{2,2}(4,4)+1, d2{3,1}(4,2)+d2{3,2}(4,4)+1,...
                       d2{1,3}(2,4)+d2{2,1}(4,2), d2{3,3}(4,4)+d2{1,1}(2,2),...
                        d2{2,3}(4,4)+d2{2,1}(4,2)+1, d2{3,3}(4,4)+d2{2,1}(4,2)+1, d2{3,3}(4,4)+d2{3,1}(4,2)+1,...
                       d2{1,2}(2,4)+d2{2,2}(4,4)+1, d2{3,2}(4,4)+d2{1,2}(2,4)+1,...
                        d2{2,2}(4,4)+d2{2,2}(4,4)+2, d2{3,2}(4,4)+d2{2,2}(4,4)+2, d2{3,2}(4,4)+d2{3,2}(4,4)+2,...
                       d2{1,3}(2,4)+d2{2,2}(4,4)+1, d2{3,3}(4,4)+d2{1,2}(2,4)+1,...
                        d2{2,3}(4,4)+d2{2,2}(4,4)+2, d2{3,3}(4,4)+d2{2,2}(4,4)+2, d2{3,3}(4,4)+d2{3,2}(4,4)+2,...
                       d2{1,3}(2,4)+d2{2,3}(4,4)+1, d2{3,3}(4,4)+d2{1,3}(2,4)+1,...
                        d2{2,3}(4,4)+d2{2,3}(4,4)+2, d2{3,3}(4,4)+d2{2,3}(4,4)+2, d2{3,3}(4,4)+d2{3,3}(4,4)+2]);
D.R22{3,2} = D.R22{2,2};
D.R22{3,2}([2;3],:) = D.R22{3,2}([3;2],:);
D.R22{2,3} = D.R22{2,2};
D.R22{2,3}(:,[2,3]) = D.R22{2,3}(:,[3,2]);
D.R22{3,3} = D.R22{3,2};
D.R22{3,3}(:,[2,3]) = D.R22{3,3}(:,[3,2]);
end
end


function [maxdegs,n_updates] = reduce_joint_degs(maxdegs)
% This function reduces the joint degrees presented in maxdegs to sensible
% values, making sure that the individual degrees do not exceed the joint
% degrees, nor that the sum of the individual degrees exceed the joint
% degrees.
%
% INPUT
%   maxdegs: An array of 2^nvars elements, for nvars variables. Required to
%            be a (reshaped version of a) 2x2x...x2 array, with each
%            dimension corresponding to a single variable, so that e.g.
%            element (2,1) corresponds to the degree of ss1, element (1,2)
%            corresponds to the degree of ss2, and element (2,2)
%            corresponds to the joint degree in ss1*ss2. Note that the
%            standard degree object used in e.g. poslpivar, taking the form
%
%                 0 |         ss2 |         tt2 |         ss2*tt2
%          ---------|-------------|-------------|-----------------
%               ss1 |     ss1*ss2 |     ss1*tt2 |     ss1*ss2*tt2
%          ---------|-------------|-------------|-----------------
%               tt1 |     tt1*ss2 |     tt1*tt2 |     tt1*ss2*tt2 
%          ---------|-------------|-------------|-----------------
%           ss1*tt1 | ss1*tt1*ss2 | ss1*tt1*tt2 | ss1*tt1*ss2*tt2  
%
%            can be easily reshaped to the 2x2x2x2 form. Also note that the
%            first element must always be zero (anything else is ignored).
%
% OUTPUT
%   maxdegs: An array of the same dimensions as the input, where now the
%            degrees are reduced to such an extent that:
%           - The joint degree in any set of variables (x_1,...,x_n) does
%             not exceed the joint degree in the set
%             (x_1,...,x_{i-1},x_{i+1},...,x_n) for any i\in\{1,...,n\};
%           - The joint degree in any set of variables (x_1,...,x_n) does
%             not exceed the sum of the degrees in x_i and in the set
%             (x_1,...,x_{i-1},x_{i+1},...,x_n) for any i\in\{1,...,n\};
%            Note that degrees are only reduced, never increased.
%   n_updates: Number of cycles of updates performed to obtain the final
%              maxdegs output. If the joint degrees in the input maxdegs
%              are already sensible, n_updates will be zero. Otherwise, it
%              should be just 1, but there may be particular (edge) cases
%              for which more than 1 cycle is necessary. Please share such
%              cases with us. The number of cycles is limited to 5.

% % Check how many variables we're considering
if all(size(maxdegs)==2)
    nvars = ndims(maxdegs);
else
    ndegs = numel(maxdegs);
    nvars = round(log(ndegs)/log(2));
    if ndegs==2^nvars-1
        maxdegs = [0;maxdegs(:)];
    elseif ndegs~=2^nvars
        error('Number of degrees in the input should be 2^nvars for nvars variables')
    end
end
dsize_prod = 2.^(0:nvars-1);    % Vector translating step 1 increase in dimension d to associated increase in linear index

if maxdegs(1)>0 && ~isnan(maxdegs(1))
    warning('The very first element of the maximal degree object does not correspond to any variable; the input value will be ignored')
    maxdegs(1) = 0;
end

% % To be safe, we perform several cycles of updates to the degrees, though
% % 1 should be sufficient (in general)
maxdegs_old = maxdegs;
max_updates = 5;
n_updates = max_updates;
for m=1:max_updates
% % Iterate over each of the (joint) degrees, making sure the joint degree
% % does not decrease upon including an additional variable, nor increase
% % with more than the maximal degree of this variable
for k=2^nvars-1:-1:2
    b_indx = str2num(dec2bin(k-1,nvars)')';  
    v_indx = b_indx(end:-1:1)>0;     % Logical array indicating which variables contribute to this (joint) degree
    for v=find(~v_indx)     % For each variable that does not contribute...
        i_indx = dsize_prod(v) + 1; % Establish the linear index associated to this individual variable      
        j_indx = dsize_prod(v) + k; % Establish the linear index associated to the joint degree with this variable
        
        % Make sure the joint degree without the variable does not exceed that with the variable
        maxdegs(k) = min(maxdegs(k),maxdegs(j_indx));
        
        % Make sure the joint degree with the variable does not exceed that 
        % of the sum of the degree without the variable and that of the 
        % individual variable
        maxdegs(j_indx) = min(maxdegs(j_indx),maxdegs(k)+maxdegs(i_indx)); 
    end
end
if ~any(maxdegs_old(:)-maxdegs(:))
    n_updates = m-1;
    break
elseif m==max_updates
    warning(['More than the allowed ',num2str(max_updates),' updates to the degree structure are necessary in reducing the joint degrees; something might be wrong...'])
else
    maxdegs_old = maxdegs;
end
end

end