%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_terms.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE_out=convert_PIETOOLS_PDE_terms_internal(PDE)
% convert_PIETOOLS_PDE_terms is the new version of the PDE converter
% file which uses the term-based data structure and performs the following tasks.
% The ODE-PDE is parametrized as
%  [x_t] = [A  Bxw Bxu Bxr] [x]
%  [z]     [Cz Dzw Dzu Dzr] [w]
%  [y]     [Cy Dyw Dyu Dyr] [u]
%  [v]     [Cv Dvw Dvu 0  ] [r]
%
%  [r]   = [    \int_a^b [Crp]           ] [X_D] + [0   Drb] [v  ]
%  [X_t]   [Ap0+\int_a^s Ap1+\int_s^b AP2] [   ]   [Bpv Bpb] [X_b]
%
%  B X_b = \int_a^b Bp X_D + Bv v
% where X_D is all derivatives ordered as X_D = [X_0; X_1...;X_N; dX_1...;dX_N;....;d^(N)X_N]
% and X_b is all boundary values ordered as X_b =
% [X_1(a)...;X_N(a);dX_2(a)...;dX_N(a)...;X_1(b)...;X_N(b);dX_2(b)...;d^(N-1)X_N(b)].


% A Partial Integral Equation is defined by 12 PI operators as
%
% BT1op \dot{w}(t)+BT2op \dot{u}(t)+Top \dot{x}(t)= Aop x(t) + B2op u(t) + B1op w(t)
%                                             z(t)= C1op x(t) + D12op u(t) + D11op w(t)
%                                             y(t)= C2op x(t) + D22op u(t) + D21op w(t)
%
% This script takes a user-defined PDE system in the format outlined in the
% header of solver_PIETOOLS_PDE and converts it to a PIE by defining the 11
% PI operators {BT1op,BT2op,Top,Aop,B1op,B2op,C1op,D11op,D12op,C2op,D21op,D22op}
%
%
% ~ NOTE: The workspace is cleaned at the end of the file, but feel free to
% exclude variables from this cleaning process if they are needed for
% anything.
%
% Perhaps we should create a PIE structure similar to the PDE structure,
% collecting all the necessary/useful variables. I don't know what other
% codes would then also have to be adjusted though.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following script performs error checking operations and sets
% undefined operators to empty opbjects.
PDE=initialize_PIETOOLS_PDE_terms_legacy(PDE);


% Converts ODE-PDE to PIE and defines PI operators
fprintf('\n --- Converting ODE-PDE to PIE --- \n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define auxiliary variables to convert the ODE-PDE to a PIE
X = PDE.dom;

polyvars = PDE.vars;
PIE.vars = PDE.vars;
s1 = polyvars(1);
s1_dum = polyvars(2);
n = PDE.n;
n_pde = n.n_pde;
N = length(n_pde)-1;

if N<1
    PIE_out = 'No distributed state. Cannot convert to PIE';
    return
end

a = X(1);      b = X(2);

np_all_der = sum((1:N+1));  % used only in indexing
np_all_der_full = sum((1:N+1).*n_pde); % size of X_D

no= n.nx; np = sum(n_pde);
nz = n.nz; ny = n.ny;
nw = n.nw; nu = n.nu;
nv = n.nv; nr = n.nr;
nBVs = 2*sum((0:N).*PDE.n.n_pde); % total boundary values, size of X_b

% now gather all parameters from terms format
Ap0=[]; Ap1 = []; Ap2 = []; Bpb = []; % pde dynamics parameters
Crp = []; Drb =[]; % pde to ode signal parameters
B = []; Bp = []; % BC parameters
for j=1:N+1
    TAp0 = []; TAp1 = []; TAp2 = []; % < - temp variables
    for i=1:np_all_der
        tmp = PDE.PDE.A{i+(j-1)*np_all_der};
        tmp2 = PDE.PDE.A{i+(j-1)*np_all_der+(N+1)*np_all_der};
        tmp3 = PDE.PDE.A{i+(j-1)*np_all_der+2*(N+1)*np_all_der};
        TAp0 = [TAp0 tmp.coeff];
        TAp1 = [TAp1 tmp2.coeff];
        TAp2 = [TAp2 tmp3.coeff];
    end
    Ap0 = [Ap0;TAp0];
    Ap1 = [Ap1;TAp1];
    Ap2 = [Ap2;TAp2];
end
for j=1:N+1
    TBpb = [];
    for i=1:2*(np_all_der-N-1)
        tmp = PDE.PDE.Bpb{i+(j-1)*2*(np_all_der-N-1)};
        TBpb = [TBpb tmp.coeff];
    end
    Bpb = [Bpb;TBpb];
end
for i=1:np_all_der
    tmp4 = PDE.PDE.Crp{i}; tmp5 = PDE.BC.Ebp{i};
    Crp = [Crp tmp4.coeff];
    Bp = [Bp tmp5.coeff];
end

for i=1:2*(np_all_der-N-1)
    tmp = PDE.BC.Ebb{i}; tmp2 = PDE.PDE.Drb{i};
    Drb = [Drb tmp2.coeff];
    B = [B tmp.coeff];
end

Bv = PDE.BC.Ebv;


% Next we construct T, U1, U2, Q, T1, Q1 (from the paper)
bigI = eye(np_all_der_full);

% Split bigI into U1 and U2 such that x_D = U2x_bc(s)+ U1 x_f
idx2 = [];
ci=0;
for i=0:N
    idx2 = [idx2, ci+1:ci+n_pde(i+1)];
    ci = ci+ sum(n_pde(i+1:end));
end
U1 = bigI(:,idx2);
U2 = bigI(:,setdiff(1:end,idx2));

T = []; Q = [];
for i =1:N
    Ti=[]; Qi=[]; % <- temp variables
    for k=1:N
        if i>k
            tau{i,k} = zeros(sum(n_pde(i+1:end)),sum(n_pde(k+1:end)));
        else
            tau{i,k} = ((s1-a)^(k-i)/factorial(k-i))...
                *[zeros(sum(n_pde(i+1:k)),sum(n_pde(k+1:end)));eye(sum(n_pde(k+1:end)))];
            Qi = blkdiag(Qi,((s1-s1_dum)^(k-i)/factorial(k-i))*eye(n_pde(k+1)));
        end
        Ti = [Ti,tau{i,k}];
    end
    if i==1
        T1 = Ti;
        Q1 = [zeros(size(Qi,1),n_pde(1)), Qi];
    end
    Q = [Q;[zeros(sum(n_pde(i+1:end)),sum(n_pde(1:i))) Qi]];
    T = [T; Ti];
end

G0 = [eye(n_pde(1)) zeros(n_pde(1),np-n_pde(1)); 
     zeros(np-n_pde(1),n_pde(1)) zeros(np-n_pde(1))];


pvar eta;
% now find E_T (B_T in paper) and Q_T (B_T*B_Q in paper)
B_T = B*[subs(T,s1,a);subs(T,s1,b)]-int(Bp*U2*T,s1,a,b);

try B_T = double(B_T);
catch
    error('Boundary conditions seem to be function of space');
end

% Test for well-posedness of Boundary conditions
if abs(det(B_T))>eps
    BTinv = inv(B_T);
    singularET = 0;
else % do SVD
    %     warning('Defined boundary conditions are rank deficient or have prohibited boundary conditions. Your system is likely ill-posed.');
    %     [U,S,V] = svd(ET);
    singularET = 1;
end
QT = B*[zeros(size(Q));subs(subs(Q,s1,b),s1_dum,s1)] - Bp*U1 - int(subs(subs(Bp*U2*Q,s1,eta),s1_dum,s1),eta,s1,b);
B_Q = BTinv*QT;

if ~singularET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % next we initialize intermdiate opvars that will be used in the script
    % X = Thatop Xhat + Tvop v
    opvar Thatop Tvop;
    Thatop.dim = [0 0; np np]; Thatop.var1 = s1; Thatop.var2 = s1_dum; Thatop.I = X;
    Tvop.dim = [0 nv; np 0]; Tvop.var1 = s1; Tvop.var2 = s1_dum; Tvop.I = X;
    
    Thatop.R.R0 = G0; 
    Thatop.R.R2 = [zeros(n_pde(1),np);T1*subs(B_Q,s1,s1_dum)]; 
    Thatop.R.R1 = Thatop.R.R2+[zeros(n_pde(1),np); Q1];
    
    Tvop.Q2 = [zeros(n_pde(1),np); T1*BTinv*Bv];

    % construct Aop_ODE for ODE subsystem, [x_t; z; y; v] = Aop_ODE*[x; w; u; r]
    opvar Aop_ODE;

    Aop_ODE.dim = [no+nz+ny+nv no+nw+nu+nr; 0 0]; 
    Aop_ODE.var1 = s1; Aop_ODE.var2 = s1_dum; Aop_ODE.I = X;
    Aop_ODE.P = [PDE.ODE.A PDE.ODE.Bxw PDE.ODE.Bxu PDE.ODE.Bxr;
                 PDE.ODE.Cz PDE.ODE.Dzw PDE.ODE.Dzu PDE.ODE.Dzr;
                 PDE.ODE.Cy PDE.ODE.Dyw PDE.ODE.Dyu PDE.ODE.Dyr;
                 PDE.ODE.Cv PDE.ODE.Dvw PDE.ODE.Dvu zeros(nv,nr)];

    % construct Vop_out for ODE output to PDE, v = Vop_out*[x; w; u]
    Vop_out = Aop_ODE(no+nw+nu+1:no+nw+nu+nv, 1:no+nw+nu);

    % now we construct Aop_PDE for PDE subsystem, 
    % [r; X_t] = Aop_PDE*[X_all_der] + Bop_PDE*[v; X_b]
    opvar Aop_PDE Bop_PDE;

    Aop_PDE.dim = [nr 0; np np_all_der_full]; 
    Aop_PDE.var1 = s1; Aop_PDE.var2 = s1_dum; Aop_PDE.I = X;
    Bop_PDE.dim = [nr nv+nBVs; np 0]; 
    Bop_PDE.var1 = s1; Bop_PDE.var2 = s1_dum; Bop_PDE.I = X;

    Aop_PDE.R.R0 = Ap0; Aop_PDE.R.R1 = Ap1; Aop_PDE.R.R2 = Ap2;
    Aop_PDE.Q1 = Crp;

    Bop_PDE.P = [zeros(nr,nv) Drb]; Bop_PDE.Q2 = [PDE.PDE.Bpv Bpb];

    %%% We now construct Top_vbd: [v; X_b; X_D] = Top_vbd[v; X_f] 
    opvar Top_vbd;
    Top_vbd.dim = [nv+nBVs nv;np_all_der_full np]; 
    Top_vbd.var1 = s1; Top_vbd.var2 = s1_dum; Top_vbd.I = X;
    Top_vbd.P =  [eye(nv); [subs(T,s1,a);subs(T,s1,b)]*BTinv*Bv];
    Top_vbd.Q1 = [zeros(nv,np); B_Q; subs(T,s1,b)*B_Q+subs(subs(Q,s1,b),s1_dum,s1)];



    %%% We now construct the PbigP, [r; X_t] = PbigP*[x;w;u;x_f]
    % = Aop_PDE*[X_D] + Bop_PDE*[v; X_b] = [Bop_PDE Aop_PDE]*[v; X_b; X_D]
    % = [Bop_PDE Aop_PDE]*Top_vbd*[v; X_f] = [Bop_PDE Aop_PDE]*Top_vbd*[Vop_out 0; 0 I]*[x; w; u; X_f]
    opvar ZZ1 ZZ2 Inp;
    ZZ1.dim = [nv 0; 0 np]; ZZ2.dim = [0 no+nw+nu; np 0]; 
    Inp.dim = [0 0; np np]; Inp.R.R0 = eye(np);
    PbigP = [Bop_PDE Aop_PDE]*Top_vbd*[Vop_out ZZ1; ZZ2 Inp];
    
    %%% We now construct the PbigO, [x_t;z;y;v] = PbigO*[x;w;u;X_f]
    % = Aop_ODE*[x; w; u; r] 
    % = Aop_ODE*([I 0 0;0 I 0; 0 I 0; 0 0 0]*[x; w; u]
    %      + [0; 0; 0; I]*r)
    % = Aop_ODE([I 0 0 0;0 I 0 0; 0 I 0 0; 0 0 0 0]*[x; w; u; X_f]
    %      +[0; 0; 0; I]*TbigP(1:nr,:)*[x; w; u; X_f])
    
    opvar Ixwuf; Ixwuf.dim = [no+nw+nu+nr no+nw+nu; 0 np];
    Ixwuf.P = [eye(no+nw+nu); zeros(nr,no+nw+nu)];

    PbigO = Aop_ODE*(Ixwuf+...
        [zeros(no,nr); zeros(nw,nr); zeros(nu,nr); eye(nr)]*PbigP(1:nr,1:no+nw+nu+np));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first construct all opvars needed for PIE
    opvar Top Twop Tuop; %LHS opvars
    Top.dim = [no no; np np]; Top.var1 = s1; Top.var2 = s1_dum; Top.I = X;
    Twop.dim = [no nw; np 0]; Twop.var1 = s1; Twop.var2 = s1_dum; Twop.I = X;
    Tuop.dim = [no nu; np 0]; Tuop.var1 = s1; Tuop.var2 = s1_dum; Tuop.I = X;

    opvar Aop B1op B2op C1op C2op D11op D12op D21op D22op; % RHS opvars

    Aop.dim = [no no; np np]; Aop.var1 = s1; Aop.var2 = s1_dum; Aop.I = X;
    B1op.dim = [no nw; np 0]; B1op.var1 = s1; B1op.var2 = s1_dum; B1op.I = X;
    B2op.dim = [no nu; np 0]; B2op.var1 = s1; B2op.var2 = s1_dum; B2op.I = X;
    C1op.dim = [nz no; 0 np]; C1op.var1 = s1; C1op.var2 = s1_dum; C1op.I = X;
    C2op.dim = [ny no; 0 np]; C2op.var1 = s1; C2op.var2 = s1_dum; C2op.I = X;
    D11op.dim = [nz nw; 0 0]; D11op.var1 = s1; D11op.var2 = s1_dum; D11op.I = X;
    D12op.dim = [nz nu; 0 0]; D12op.var1 = s1; D12op.var2 = s1_dum; D12op.I = X;
    D21op.dim = [ny nw; 0 0]; D21op.var1 = s1; D21op.var2 = s1_dum; D21op.I = X;
    D22op.dim = [ny nu; 0 0]; D22op.var1 = s1; D22op.var2 = s1_dum; D22op.I = X;


    %%% Now we have to partition PbigO and PbigP to get the desired pieces
    D11op=PbigO(no+1:no+nz,no+1:no+nw);
    D21op=PbigO((no+nz+1):(no+nz+ny),no+1:no+nw);
    D12op=PbigO(no+1:no+nz,(no+nw+1):(no+nw+nu));
    D22op=PbigO((no+nz+1):(no+nz+ny),(no+nw+1):(no+nw+nu));
    C1op=[PbigO(no+1:no+nz,(1):(no)),PbigO(no+1:no+nz,(no+nw+nu+1):(no+nw+nu+np))];
    C2op=[PbigO((no+nz+1):(no+nz+ny),(1):(no)),PbigO((no+nz+1):(no+nz+ny),(no+nw+nu+1):(no+nw+nu+np))];
    B1op=[PbigO((1):(no),no+1:no+nw);PbigP(nr+1:nr+np,no+1:no+nw)];
    B2op=[PbigO((1):(no),no+nw+1:no+nw+nu);PbigP(nr+1:nr+np,no+nw+1:no+nw+nu)];
    Aop=[PbigO((1):(no),(1):(no))   PbigO((1):(no),(no+nw+nu+1):(no+nw+nu+np));
         PbigP(nr+1:nr+np,1:no)     PbigP(nr+1:nr+np,(no+nw+nu+1):(no+nw+nu+np))];

    % LHS PI operators
    tmp = Tvop*Vop_out;
    Twop.Q2 = tmp(1:np,(no+1):(no+nw)).Q2; Twop.R = tmp(1:np,(no+1):(no+nw)).R;
    Tuop.Q2 = tmp(1:np,(no+nw+1):(no+nw+nu)).Q2; Tuop.R = tmp(1:np,(no+nw+1):(no+nw+nu)).R;
    Top = Thatop;
    Top.P = eye(no); Top.Q2 = tmp(1:np,1:no).Q2;
else
    % do nothing
    % fill the singular case formulae
    PIE_out = 'Singular B_T: A PIE representation of this PDE may not exist';
    return
end


%nx1 = no; nx2 = np;

% DJ 09/29/2021: fill in empty parameters with appropriate zeros
opnames = {'Top','Twop','Tuop','Aop','B1op','B2op','C1op','D11op','D12op','C2op','D21op','D22op'};
for j=1:length(opnames)
    op = opnames{j};
    eval([op,'.dim=',op,'.dim;'])
end

PIE.T = Top; PIE.Tw = Twop; PIE.Tu = Tuop;
PIE.A = Aop; PIE.B1 = B1op; PIE.B2 = B2op;
PIE.C1 = C1op; PIE.D11 = D11op; PIE.D12 = D12op;
PIE.C2 = C2op; PIE.D21 = D21op; PIE.D22 = D22op;
PIE.dom = PDE.dom;

PIE_out=PIE;

% %remove temporary opvars
% clear TbigO TbigP Bop_PDE Bop_PDE_sliced Vop_out Vop_out_extended
% clear Aop_ODE Aop_ODE_sliced Aop_PDE
% clear TBpb Tvop TDop TBVop Tu Tw Ti Thatop T
% clear Ap0 Ap1 Ap2 Bpb Crp Drb Ebb Ebp Ebv
% clear K L1 Q Qi QT ET ET ETinv G0 TAp0 TAp1 TAp2 U1 U2 X singularET
% clear n no np nr nu nv nw nx1 nx2 ny nz n_pde nBVs np_all_der np_all_der_full N
% clear a b i j k idx2 ci bigI
% clear tmp tmp2 tmp3 tmp4 tmp5
% clear polyvars eta  tau
end
