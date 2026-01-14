function P_out = ndopvar2dopvar(P_in)
% P_OUT = NDOPVAR2DOPVAR(P_IN) takes an 'ndopvar' object, P_IN and builds 
% an 'opvar' or 'dopvar' object representing the same PI operator. 
% In particular, the parameters defining this operator will be given by
%   R_{i}(s,t) = (Im o Zd(s))^T (Ik o [1;dvars])^T C{i} (In o Zd(t))
% where the coefficient C{i} are stored in P_IN.C, the degree d in 
% P_IN.deg, the decision variables dvars in P_IN.dvarname,
% and [m,n] = P_IN.dim;

% Check the inputs
if ~isa(P_in,'ndopvar')
    error("Input must be of type 'ndopvar'.")
end

% Establish the number of variables
N = size(P_in.dom,1);
if N>2
    error("Only 1D or 2D PI operators are supported.")
end
m = P_in.dim(1);    n = P_in.dim(2);

% Get the variables
d = P_in.deg;
var1 = P_in.vars(:,1);
var2 = P_in.vars(:,2);

% Declare the monomial bases in each variable
Z1t = 1;
Z2_cell = cell(N,1);
for ii=1:N
    Z1t = kron(Z1t,polynomial(eye(d(ii)+1),(0:d(ii))',var1(ii).varname,[1,d(ii)+1]));
    Z2_cell{ii} = polynomial(eye(d(ii)+1),(0:d(ii))',var2(ii).varname,[d(ii)+1,1]);
end
Z1t = kron(eye(m),Z1t);

% Include the decision variables
dvarname = P_in.dvarname;       q = numel(dvarname);
if q>0
    k = size(Z1t,2);    
    Zdec = dpvar(speye(k*(q+1)),zeros(1,0),{},dvarname,[k,k*(q+1)]);
    Z1t = Z1t*Zdec;
end

% Compute the parameters associated with each coefficient matrix
C = P_in.C;
sz_C = size(C);
sz_C = sz_C(1:N);
for ii=1:numel(C)
    % Determine the index of element ii along each dimension of the cell C
    idcs = cell(1,N);
    [idcs{:}] = ind2sub(sz_C,ii);
    idcs = cell2mat(idcs);
    % If element ii corresponds to an integral, we need to include the
    % monomial basis in the associated dummy variable
    is_int = logical(idcs-1);
    Z2_tmp = Z2_cell;
    Z2_tmp(~is_int) = {1};
    Z2 = 1;
    for jj=1:N
        Z2 = kron(Z2,Z2_tmp{jj});
    end
    % Compute the parameter associated with coefficient matrix C{ii}
    C{ii} = Z1t*C{ii}*kron(eye(n),Z2);
end

% Declare the operator in opvar format
if N==1
    if q==0
        P_out = opvar();
    else
        P_out = dopvar();
    end
    P_out.R.R0 = C{1};
    P_out.R.R1 = C{2};
    P_out.R.R2 = C{3};
else
    if q==0
        P_out = opvar2d();
    else
        P_out = dopvar2d();
    end
    for ii=1:numel(P_out.R22)
        P_out.R22{ii} = C{ii};
    end
end
P_out.var1 = var1;      P_out.var2 = var2;
P_out.I = P_in.dom;

end