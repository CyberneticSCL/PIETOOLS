function PDE = getPDEstructure_old(PDE_GUI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getPDEstructure.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is exclusively used by the PDE_GUI to convert GUI input PDE
% object to standard terms format of the PDE object.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following inputs must be defined passed:
%
% PDE_GUI - PDE object from the GUI
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications


% added modification to allow repeated terms
PDE.dim = 1;
PDE.dom = [0 1];


PDE.n = PDE_GUI.n;
PDE.dom = [0 1]; % GUI input domain is always 0 to 1
PDE.vars = PDE_GUI.vars;
s = PDE.vars(1); theta = PDE.vars(2);


nr = length(PDE_GUI.ODE.Bxb)+length(PDE_GUI.ODE.Bxi)+length(PDE_GUI.ODE.Dzb)...
    +length(PDE_GUI.ODE.Dzi)+length(PDE_GUI.ODE.Dyb)+length(PDE_GUI.ODE.Dyi);
PDE.n.nr = nr;
PDE.n.nv = PDE.n.nx+PDE.n.nw+PDE.n.nu;


varidx = PDE_GUI.varidx;
Xvaridx = PDE_GUI.Xvaridx;

% storing dimensions for convenience
PDE = initialize_PIETOOLS_PDE_terms_legacy(PDE);

no = PDE.n.nx; nw = PDE.n.nw; nu = PDE.n.nu;
n_pde = PDE.n.n_pde; nz = PDE.n.nz; ny = PDE.n.ny;
N = length(n_pde)-1; %max derivative
nBC = sum(n_pde.*(0:N)); nBV = 2*nBC;
np_all_derivatives = sum(1:N+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PDE.ODE.Cv = [eye(no,no);zeros(nw,no);zeros(nu,no)]; % GUI input v is always [x;w;u]
PDE.ODE.Dvw = [zeros(no,nw);eye(nw);zeros(nu,nw)];
PDE.ODE.Dvu = [zeros(no,nu);zeros(nw,nu);eye(nu)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get ODE parameters out
for i=1:length(PDE_GUI.ODE.A)
    tmp = PDE_GUI.ODE.A{i};
    j = varidx(tmp.Lstate); 
    k = varidx(tmp.Rstate); 
    PDE.ODE.A(j,k) = PDE.ODE.A(j,k)+tmp.coeff;
end
for i=1:length(PDE_GUI.ODE.Bxw)
    tmp = PDE_GUI.ODE.Bxw{i};
    j = varidx(tmp.Lstate); 
    k = varidx(tmp.Rstate);
    PDE.ODE.Bxw(j,k) = PDE.ODE.Bxw(j,k)+tmp.coeff;
end
for i=1:length(PDE_GUI.ODE.Bxu)
    tmp = PDE_GUI.ODE.Bxu{i};
    j = varidx(tmp.Lstate); 
    k = varidx(tmp.Rstate);
    PDE.ODE.Bxu(j,k) = PDE.ODE.Bxu(j,k)+tmp.coeff;
end
for m=1:length(PDE_GUI.ODE.Bxb)
        tmp = PDE_GUI.ODE.Bxb{m};
        k = tmp.D; l = tmp.Rstate(1:end-3);
        delta = str2double(tmp.Rstate(end-1));
        locidx = Xvaridx(l);
        PDE.ODE.Bxr(varidx(tmp.Lstate),m) = 1; 
        % add to r(t)
        locB = delta*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + locidx(1);
        cf = PDE.PDE.Drb{locB}.coeff; cf(m,locidx(2)) = cf(m,locidx(2))+double(tmp.coeff);
        PDE.PDE.Drb{locB} = struct('D',k,'Rstate',locidx(1),'delta',delta,'coeff',cf);
end
M = length(PDE_GUI.ODE.Bxb);
for m=1:length(PDE_GUI.ODE.Bxi)
        tmp = PDE_GUI.ODE.Bxi{m};
        k = tmp.D; l = tmp.Rstate;
        locidx = Xvaridx(l);
        PDE.ODE.Bxr(varidx(tmp.Lstate),M+m)=1;
        % add to r(t)
        loc = sum(N:-1:N-k+1) + locidx(1) + 1; 
        cf = PDE.PDE.Crp{loc}.coeff; cf(M+m,locidx(2)) = cf(M+m,locidx(2))+tmp.coeff;
        PDE.PDE.Crp{loc} = struct('D',k,'Rstate',locidx(1),'coeff',cf);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get reg out parameters
for i=1:length(PDE_GUI.ODE.Cz)
    tmp = PDE_GUI.ODE.Cz{i};
    j = varidx(tmp.Lstate); 
    k = varidx(tmp.Rstate); 
    PDE.ODE.Cz(j,k) = PDE.ODE.Cz(j,k)+tmp.coeff;
end
for i=1:length(PDE_GUI.ODE.Dzw)
    tmp = PDE_GUI.ODE.Dzw{i};
    j = varidx(tmp.Lstate); 
    k = varidx(tmp.Rstate);
    PDE.ODE.Dzw(j,k) = PDE.ODE.Dzw(j,k)+tmp.coeff;
end
for i=1:length(PDE_GUI.ODE.Dzu)
    tmp = PDE_GUI.ODE.Dzu{i};
    j = varidx(tmp.Lstate); 
    k = varidx(tmp.Rstate);
    PDE.ODE.Dzu(j,k) = PDE.ODE.Dzu(j,k)+tmp.coeff;
end
for m=1:length(PDE_GUI.ODE.Dzb)
        tmp = PDE_GUI.ODE.Dzb{m};
        j = str2double(tmp.Rstate(end-1)); k = tmp.D; l = tmp.Rstate(1:end-3);
        locidx = Xvaridx(l);
        loc = M+length(PDE_GUI.ODE.Bxi)+m;
        PDE.ODE.Dzr(varidx(tmp.Lstate),loc) = 1; 
        % add coeff to r(t)
        locB = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + locidx(1);
        cf = PDE.PDE.Drb{locB}.coeff; cf(loc,locidx(2)) = cf(loc,locidx(2))+double(tmp.coeff);
        PDE.PDE.Drb{locB} = struct('D',k,'Rstate',locidx(1),'delta',j,'coeff',cf);
end
M = length(PDE_GUI.ODE.Bxb)+length(PDE_GUI.ODE.Bxi)+length(PDE_GUI.ODE.Dzb);
for m=1:length(PDE_GUI.ODE.Dzi)
        tmp = PDE_GUI.ODE.Dzi{m};
        k = tmp.D; l = tmp.Rstate;
        locidx = Xvaridx(l);
        PDE.ODE.Dzr(varidx(tmp.Lstate),M+m)=1;
        % add to r(t)
        loc = sum(N:-1:N-k+1) + locidx(1) + 1; 
        cf = PDE.PDE.Crp{loc}.coeff; cf(M+m,locidx(2)) = cf(M+m,locidx(2))+tmp.coeff;
        PDE.PDE.Crp{loc} = struct('D',k,'Rstate',locidx(1),'coeff',cf);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get obs out parameters
for i=1:length(PDE_GUI.ODE.Cy)
    tmp = PDE_GUI.ODE.Cy{i};
    j = varidx(tmp.Lstate); 
    k = varidx(tmp.Rstate); 
    PDE.ODE.Cy(j,k) = PDE.ODE.Cy(j,k)+tmp.coeff;
end
for i=1:length(PDE_GUI.ODE.Dyw)
    tmp = PDE_GUI.ODE.Dyw{i};
    j = varidx(tmp.Lstate); 
    k = varidx(tmp.Rstate);
    PDE.ODE.Dyw(j,k) = PDE.ODE.Dyw(j,k)+tmp.coeff;
end
for i=1:length(PDE_GUI.ODE.Dyu)
    tmp = PDE_GUI.ODE.Dyu{i};
    j = varidx(tmp.Lstate); 
    k = varidx(tmp.Rstate);
    PDE.ODE.Dyu(j,k) = PDE.ODE.Dyu(j,k)+tmp.coeff;
end
for m=1:length(PDE_GUI.ODE.Dyb)
        tmp = PDE_GUI.ODE.Dyb{m};
        j = str2double(tmp.Rstate(end-1)); k = tmp.D; l = tmp.Rstate(1:end-3);
        locidx = Xvaridx(l);
        loc = M+length(PDE_GUI.ODE.Dzi)+m; 
        PDE.ODE.Dyr(varidx(tmp.Lstate),loc) = 1; 
        % add coeff to r(t)
        locB = j*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + locidx(1);
        cf = PDE.PDE.Drb{locB}.coeff; cf(loc,locidx(2)) = cf(loc,locidx(2))+double(tmp.coeff);
        PDE.PDE.Drb{locB} = struct('D',k,'Rstate',locidx(1),'delta',j,'coeff',cf);
end
M = M+length(PDE_GUI.ODE.Dzi)+length(PDE_GUI.ODE.Dyb);
for m=1:length(PDE_GUI.ODE.Dyi)
        tmp = PDE_GUI.ODE.Dyi{m};
        k = tmp.D; l = tmp.Rstate;
        locidx = Xvaridx(l);
        PDE.ODE.Dyr(varidx(tmp.Lstate),M+m)=1;
        % add to r(t)
        loc = sum(N:-1:N-k+1) + locidx(1) + 1; 
        cf = PDE.PDE.Crp{loc}.coeff; cf(M+m,locidx(2)) = cf(M+m,locidx(2))+tmp.coeff;
        PDE.PDE.Crp{loc} = struct('D',k,'Rstate',locidx(1),'coeff',cf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get BC parameters
for m=1:length(PDE_GUI.BC.Ebb)
        tmp = PDE_GUI.BC.Ebb{m};
        delta = str2double(tmp.Rstate(end-1)); k = tmp.D; l = tmp.Rstate(1:end-3);
        locidx = Xvaridx(l);
        
        
        loc = delta*(np_all_derivatives-N-1) + sum(N-1:-1:N-k) + locidx(1); 
        PDE.BC.Ebb{loc}.coeff(tmp.Lstate,locidx(2)) = PDE.BC.Ebb{loc}.coeff(tmp.Lstate,locidx(2))+tmp.coeff;
end
for m=1:length(PDE_GUI.BC.Ebp)
        tmp = PDE_GUI.BC.Ebp{m};
        k = tmp.D; l = tmp.Rstate;
        locidx = Xvaridx(l);
        
        loc = sum(N:-1:N-k+1) + locidx(1) + 1;
        PDE.BC.Ebp{loc}.coeff(tmp.Lstate,locidx(2)) = PDE.BC.Ebp{loc}.coeff(tmp.Lstate,locidx(2))+polynomial(tmp.coeff);
end
Ex = PDE.BC.Ebv(:,1:PDE.n.nx); Ew = PDE.BC.Ebv(:,PDE.n.nx+1:PDE.n.nx+PDE.n.nw); 
Eu = PDE.BC.Ebv(:,PDE.n.nx+PDE.n.nw+1:end); 
% zeros(nBC,no);Ew = zeros(nBC,nw);Eu = zeros(nBC,nu);
for m=1:length(PDE_GUI.BC.Ex)
        tmp = PDE_GUI.BC.Ex{m};
        l = tmp.Rstate;
        locidx = varidx(l);
        Ex(tmp.Lstate,locidx) = Ex(tmp.Lstate,locidx)+tmp.coeff;
end
for m=1:length(PDE_GUI.BC.Ew)
        tmp = PDE_GUI.BC.Ew{m};
        l = tmp.Rstate;
        locidx = varidx(l);
        Ew(tmp.Lstate,locidx) = Ew(tmp.Lstate,locidx)+tmp.coeff;
end
for m=1:length(PDE_GUI.BC.Eu)
        tmp = PDE_GUI.BC.Eu{m};
        l = tmp.Rstate;
        locidx = varidx(l);
        Eu(tmp.Lstate,locidx) = Eu(tmp.Lstate,locidx)+tmp.coeff;
end
PDE.BC.Ebv = [Ex Ew Eu];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get PDE parameters
for m=1:length(PDE_GUI.PDE.App)
        tmp = PDE_GUI.PDE.App{m};
        j = tmp.Lstate; k = tmp.Rstate(1:end-3); l = tmp.D; 
        locidxj = Xvaridx(j); locidxk = Xvaridx(k);
        
        loc = (locidxj(1))*np_all_derivatives + sum(N:-1:N-l+1) + locidxk(1) + 1; 
        PDE.PDE.A{loc}.coeff(locidxj(2),locidxk(2)) = PDE.PDE.A{loc}.coeff(locidxj(2),locidxk(2))+tmp.coeff;
end
for m=1:length(PDE_GUI.PDE.Bpi)
        tmp = PDE_GUI.PDE.Bpi{m};
        j = tmp.Lstate; k = tmp.Rstate; l=tmp.D; i = tmp.I;
        if isequal(i(2),s)
            i=1;
        else
            i=2;
        end
        locidxj = Xvaridx(j); locidxk = Xvaridx(k);
        
        loc = (i)*(N+1)*np_all_derivatives +(locidxj(1))*np_all_derivatives + sum(N:-1:N-l+1) + locidxk(1) + 1; 
        PDE.PDE.A{loc}.coeff(locidxj(2),locidxk(2)) = PDE.PDE.A{loc}.coeff(locidxj(2),locidxk(2))+polynomial(tmp.coeff);
end
for m=1:length(PDE_GUI.PDE.Bpb)
        tmp = PDE_GUI.PDE.Bpb{m};
        j = tmp.Lstate; k = tmp.Rstate; l = tmp.D; delta = k(end-1);
        locidxj = Xvaridx(j); locidxk = Xvaridx(k(1:end-3));
        loc=(locidxj(1))*(2*np_all_derivatives-2*(N+1)) + delta*(np_all_derivatives-N-1)+(l*N-(l)*(l+1)/2) + locidxk(1); 
        PDE.PDE.Bpb{loc}.coeff(locidxj(2),locidxk(2)) = PDE.PDE.Bpb{loc}.coeff(locidxj(2),locidxk(2))+polynomial(tmp.coeff);
end
Bpx = PDE.PDE.Bpv(:,1:PDE.n.nx); Bpw = PDE.PDE.Bpv(:,PDE.n.nx+1:PDE.n.nx+PDE.n.nw); 
Bpu = PDE.PDE.Bpv(:,PDE.n.nx+PDE.n.nw+1:end); 
% polynomial(zeros(sum(n_pde),no));Bpw = polynomial(zeros(sum(n_pde),nw));
% Bpu = polynomial(zeros(sum(n_pde),nu));
for m=1:length(PDE_GUI.PDE.Bpx)
        tmp = PDE_GUI.PDE.Bpx{m};
        j = Xvaridx(tmp.Lstate); k = tmp.Rstate; 
        loc = getLocA(n_pde,j(1),j(2),0,0);
        Bpx(loc,varidx(k)) = Bpx(loc,varidx(k))+polynomial(tmp.coeff);
end
for m=1:length(PDE_GUI.PDE.Bpu)
        tmp = PDE_GUI.PDE.Bpu{m};
        j = Xvaridx(tmp.Lstate); k = tmp.Rstate; 
        loc = getLocA(n_pde,j(1),j(2),0,0);
        Bpu(loc,varidx(k)) = Bpu(loc,varidx(k))+polynomial(tmp.coeff);
end
for m=1:length(PDE_GUI.PDE.Bpw)
        tmp = PDE_GUI.PDE.Bpw{m};
        j = Xvaridx(tmp.Lstate); k = tmp.Rstate; 
        loc = getLocA(n_pde,j(1),j(2),0,0);
        Bpw(loc,varidx(k)) = Bpw(loc,varidx(k))+polynomial(tmp.coeff);
end
PDE.PDE.Bpv = [Bpx Bpw Bpu];

for m=1:length(PDE_GUI.PDE.Bpb)
        tmp = PDE_GUI.PDE.Bpb{m};
        j = tmp.Lstate; k = tmp.D; l = tmp.Rstate(1:end-3); delta = tmp.Rstate(end-1);
        locidx = Xvaridx(l); locidxj = Xvaridx(j);
        loc = locidxj(1)*(2*np_all_derivatives-2*(N+1)) + delta*(np_all_derivatives-N-1)+(k*N-(k)*(k+1)/2) + locidx(1);
        PDE.PDE.Bpb{loc}.coeff(locidxj(2),locidx(2)) = PDE.PDE.Bpb{loc}.coeff(locidxj(2),locidx(2))+polynomial(tmp.coeff);
end

end
function loc = getLocA(n_pde, whichXi, whereInXi, D, I)
loc = I*sum(n_pde.*(1:length(n_pde)));
for j=0:D-1
    loc = loc+sum(n_pde(j+1:end));
end
loc = loc+sum(n_pde(D+1:whichXi))+whereInXi;
end