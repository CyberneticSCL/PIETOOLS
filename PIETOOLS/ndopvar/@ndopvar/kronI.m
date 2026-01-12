function P_out = kronI(P_in,k)
% P_OUT = KRONI(P_IN,K) returns the 'ndopvar' object representing the
% Kronecker product of the operator defined by 'ndopvar' object P_IN and
% an identity matrix of dimension K x K
%
% INPUTS
% - P_in:   'ndopvar' object representing a PI operator Pop;
% - k:      scalar integer specifying the dimension of the identity
%           matrix;
%
% OUTPUTS
% - P_out:  'ndopvar' object representing the operator (Pop o I_{k});

% if ~isempty(P_in.dvarname)
%     error("Decision variable operators are currently not supported.")
% end


% Construct the permutation matrices (I_{m(d+1)} o P)^T and 
% (I_{n(d+1)} o P) for P defined such that
%   Z_{d}(s) o I_{k} = P*(I_{k} o Z_{d}(s))
deg = P_in.deg;
N = numel(deg);
d = prod(deg+1);
m = P_in.dim(1);   n = P_in.dim(2);

%r_idcs = 1:k*d;
%c_idcs = (0:k-1)'*d + (1:d);
%P = sparse(r_idcs,c_idcs,1,q*d,q*d);
%Pt = sparse(c_idcs,r_idcs,1,k*d,k*d);
%IP_m = spIkron(m,Pt);
%IP_n = spIkron(n,P);

% Declare the (k,d) commutation matrix, so that
%   Z_{d}(s)' o I_{k} = (I_{k} o Z_{d}(s)')P_{k,d}'
Pt = commat(k,d,'transpose');

% In case of decision variables, we get an extra factor
if isempty(P_in.dvarname)
    IP_q = 1;
    IP_m = spIkron(m,Pt);
else
    q = numel(P_in.dvarname);
    Pt2 = commat(k,q+1,'transpose');
    IP_q = spIkron(d,Pt2);
    IP_m1 = kron(Pt,speye(q+1));
    IP_m = spIkron(m,IP_m1*IP_q);
    % [(I_{m} o Zd(s))' (I_{md} o [1;dvars])' o I_{k}]
    %   = (I_{mk} o Zd(s))' (I_{mdk} o [1;dvars])' IP_m
end

% Construct the coefficients representing C_str o I_{q}
C_cell = P_in.C;
D_cell = cell(size(C_cell));
sz_C = size(C_cell);
for ii=1:numel(D_cell)
    % Determine the index of element ii along each dimension of the cell C
    idcs = cell(1,N);
    [idcs{:}] = ind2sub(sz_C,ii);
    idcs = cell2mat(idcs);
    % If element ii corresponds to an integral, we need to account for the
    % monomial basis in the associated dummy variable
    is_int = logical(idcs-1);
    d_tmp = prod(deg(is_int)+1);
    r_idcs = 1:k*d_tmp;
    c_idcs = (0:k-1)'*d_tmp + (1:d_tmp);
    P = sparse(r_idcs,c_idcs,1,k*d_tmp,k*d_tmp);
    IP_n = spIkron(n,P);

    D_cell{ii} = IP_m*kron(C_cell{ii},speye(k))*IP_n;
end
P_out = P_in;
%P_out.dim = k*P_in.dim;
P_out.C = D_cell;

end

