% Test_quad_expand.m -- test quad_expand(params,varnames), which expresses
% a cell of parameters in a common set of monomials Z in the variables
% varnames.out and varnames.in, so that
%   params{k} = kron(eye(m),ZL)'*Ccell{k}*kron(eye(n),ZR).
% The parameters are random, of type 'poly', 'quadPoly' or 'dpvar' as set
% by clss ('dpvar' by default). The parameters are reconstructed from
% Ccell and Z, and the summed error is displayed as err_tot; no pass/fail
% check is made.
%
% DJ (CC), 10/15/2026: Initial coding;

%
nparams = 5;
m_arr = randi(3,[1,nparams]);
n_arr = randi(3,[1,nparams]);

% Set parameters
ndvars = 5;   % number of decision variables
M = randi([0,3]);       Mtot = M+randi([0,3]);              
N = randi([0,3]);       Ntot = N+randi([0,3]);
deg = 3;
nZ = 25;

% Set the variables
svars = cell(1,Mtot);
for i=1:Mtot
    svars{i} = ['s',num2str(i)];
end
tvars = cell(1,Ntot);
for i=1:Ntot
    tvars{i} = ['t',num2str(i)];
end
svar_cell = cell(1,nparams);
tvar_cell = cell(1,nparams);
for i=1:nparams
    s_idcs = randperm(Mtot,M);       
    t_idcs = randperm(Ntot,N);
    svar_cell{i} = svars(s_idcs);
    tvar_cell{i} = tvars(t_idcs);
end

pdegs = randi([0,deg],[N+M,1]);

% Generate random arrays of desired dimensions
clss = 'dpvar';
params = cell(1,nparams);
params_poly = cell(1,nparams);
for i=1:nparams
    if strcmp(clss,'double')
        params{i} = rand([m1,n1]);
    elseif strcmp(clss,'poly') || strcmp(clss,'quadPoly')
        A_poly = rand_poly([m_arr(i),n_arr(i)],[svar_cell{i},tvar_cell{i}]',pdegs,nZ);
        if strcmp(clss,'quadPoly')
            params{i} = quadPoly.polynomial2quadPoly(A_poly,svar_cell{i},tvar_cell{i});
            params_poly{i} = A_poly;
        else
            params{i} = A_poly;
        end
    elseif strcmp(clss,'dpvar')
        params{i} = rand_dpvar([m_arr(i),n_arr(i)],ndvars,[svar_cell{i},tvar_cell{i}],pdegs,nZ);
    else
        error("Type not supported.")
    end
end

% Take the product
varnames = struct();
varnames.in = tvars;
varnames.out = svars;
[Ccell,Z] = quad_expand(params,varnames);

% Compute parameters from coefficients
params_new = cell(1,nparams);
ZLfull = 1;
for i=1:Mtot
    ZL_i = polynomial(svars(i)).^Z.out{i};
    ZLfull = kron(ZLfull,ZL_i);
end
ZRfull = 1;
for i=1:Ntot
    ZR_i = polynomial(tvars(i)).^Z.in{i};
    ZRfull = kron(ZRfull,ZR_i);
end
err_tot = 0;
err_cell = cell(1,nparams);
for i=1:nparams
    params_new{i} = kron(eye(m_arr(i)),ZLfull)'*Ccell{i}*kron(eye(n_arr(i)),ZRfull);
    if isa(params{i},'quadPoly')
        err_cell{i} = params_poly{i} - params_new{i};
    else
        err_cell{i} = params{i} - params_new{i};
    end
    err_tot = err_tot + max(max(err_cell{i}.C));
end

disp(err_tot);

