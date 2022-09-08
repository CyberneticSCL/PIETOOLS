function PIE = convert_PIETOOLS_PDE_terms(PDE)
% For now convert the new terms format to old terms format

% extract standard information from the new terms format (domain,
% dimension, vars, etc.)
dom = PDE.dom;

odeID = PDE.x_tab(PDE.x_tab(:,3)==0,1);
pdetab = PDE.x_tab(PDE.x_tab(:,3)~=0,:);
[~,sortidx] = sort(pdetab(:,3),1);
pdetab = pdetab(sortidx,:); % PDE table sorted based on order of differentiability
pdeID = pdetab(sortidx,1);
[~,reorder] = ismember([odeID;pdeID],PDE.x_tab(:,1));
x = PDE.x(reorder); z = PDE.z; y = PDE.y; BC = PDE.BC;
xtab = PDE.x_tab(reorder,:);
odestate = xtab(:,3)==0;
odeID = xtab(odestate);

difforder = pdetab(:,4);
pdesize = pdetab(:,2);
for i=0:max(difforder)
    n.n_pde(i+1) = sum(pdesize(difforder==i));
end
N = length(n.n_pde)-1;
% find total PDE states
np = sum(n.n_pde(1:N+1));
np_all_derivatives = sum(1:N+1);
% find total possible Boundary values
nBVs=2*sum((0:N).*n.n_pde); nBC = nBVs/2;


n.nx = sum(PDE.x_tab(:,2).*odestate);
n.nz = sum(PDE.z_tab(:,2));
n.ny = sum(PDE.y_tab(:,2));
n.nw = sum(PDE.w_tab(:,2));
n.nu = sum(PDE.u_tab(:,2));




% find partitioning
wpartition = [0 cumsum(PDE.w_tab(:,2))]+1;
upartition = [0 cumsum(PDE.u_tab(:,2))]+1;
zpartition = [0 cumsum(PDE.z_tab(:,2))]+1;
ypartition = [0 cumsum(PDE.y_tab(:,2))]+1;
xpartition = [0 cumsum(PDE.x_tab(odestate,2))]+1;
Xpartition = [0 cumsum(pdetab(:,2))]+1;
rpartition = 1;
vpartition = 1;

% calculate nr and nv the interconnection signals
nr = 0; nv = 0;

for i=1:length(z) % look for nr terms in z
    for j=1:length(z{i}.term)
        nrlen = size(z{i}.term{j}.C,1);
        if isfield(z{i}.term{j},'x')&& ismember(z{i}.term{j}.x,pdeID)
            nr=nr+nrlen;
            rpartition(end+1) = rpartition(end)+nrlen;
        end
    end
end
for i=1:length(y) % look for nr terms in y
    for j=1:length(y{i}.term)
        nrlen = size(y{i}.term{j}.C,1);
        if isfield(y{i}.term{j},'x')&& ismember(y{i}.term{j}.x,pdeID)
            nr=nr+nrlen;
            rpartition(end+1) = rpartition(end)+nrlen;
        end
    end
end
for i=1:length(x)
    if odestate(i) % ode dynamics; look for nr terms
        for j=1:length(x{i}.term)
            nrlen = size(x{i}.term{j}.C,1);
            if isfield(x{i}.term{j},'x')&& ismember(x{i}.term{j}.x,pdeID)
                nr=nr+nrlen;
                rpartition(end+1) = rpartition(end)+nrlen;
            end
        end
    else % pde dynamics; look for nv terms
        for j=1:length(x{i}.term)
            nvlen = size(x{i}.term{j}.C,1);
            if isfield(x{i}.term{j},'x')&& ismember(x{i}.term{j}.x,odeID)
                nv=nv+nvlen;
                vpartition(end+1) = vpartition(end)+nvlen;
            elseif isfield(x{i}.term{j},'w')
                nv=nv+nvlen;
                vpartition(end+1) = vpartition(end)+nvlen;
            elseif isfield(x{i}.term{j},'u')
                nv=nv+nvlen;
                vpartition(end+1) = vpartition(end)+nvlen;
            end
        end
    end
end
for i=1:length(BC) % look for nv terms in BC
    for j=1:length(BC{i}.term)
        nvlen = size(BC{i}.term{j}.C,1);
        if isfield(BC{i}.term{j},'x')&& ismember(BC{i}.term{j}.x,odeID)
            nv=nv+nvlen;
            vpartition(end+1) = vpartition(end)+nvlen;
        elseif isfield(BC{i}.term{j},'w')
            nv=nv+nvlen;
            vpartition(end+1) = vpartition(end)+nvlen;
        elseif isfield(BC{i}.term{j},'u')
            nv=nv+nvlen;
            vpartition(end+1) = vpartition(end)+nvlen;
        end
    end
end

n.nr = nr; n.nv = nv;

% set all the parameters to an appropriate structure
PDE_out.n = n;
PDE_out.dom = dom;
PDE_out = initialize_PIETOOLS_PDE_terms_legacy(PDE_out);

% next place the parameters in the correct location
% collect parameters for z
k= 1; % set counter k for rows in r signal
for i=1:length(z)
    rows = zpartition(i):zpartition(i+1)-1;
    for j=1:length(z{i}.term)
        tmpterm = z{i}.term;
        if isfield(tmpterm,'x') % either ode or pde coeff; place in Cz or Dzr
            if ismember(tmpterm.x,odeID) % ode state; place in Cz
                cols = find(tmpterm.x==PDE.x_tab(:,1));
                cols = xpartition(cols):xpartition(cols+1)-1;
                PDE_out.ODE.Cz(rows,cols) = tmpterm.C;
            else % pde state; place in Dzr and PDE.Crp/Drb
                rrows = rpartition(k):rpartition(k+1)-1;
                PDE_out.ODE.Dzr(rows,cols) = tmpterm.C;
                if isfield(tmpterm,'I') %integral
                    l = pdetab((tmpterm.x==pdeID),3); % max derivative
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    rcols = pdetab([pdetab(:,3)==l],:);
                    idLoc = find(tmpterm.x==rcols(:,1));
                    cSum = [0 cumsum(rcols(:,2))]+1;
                    rcols = cSum(idLoc):cSum(idLoc)-1;
                    der = tmpterm.D;
                    loc = sum(N:-1:N-der+1) + l + 1; % find this
                    PDE_out.PDE.Crp{loc}.coeff(rrows,rcols) = tmpterm.C;
                else % boundary value
                    del = poly2double(tmpterm.loc(2)); l = pdetab((tmpterm.x==pdeID),3); % max derivative
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    der = tmpterm.D;
                    loc = del*(np_all_derivatives-N-1) + sum(N-1:-1:N-der) + l;
                    PDE_out.PDE.Drb{loc}.coeff(rrows,rcols) = tmpterm.C;
                end
                k = k+1;
            end
        elseif isfield(tmpterm,'w') % dist term; place in Dzw
            cols = find(tmpterm.w==PDE.w_tab(:,1));
            cols = wpartition(cols):wpartition(cols+1)-1;
            PDE_out.ODE.Dzw(rows,cols) = tmpterm.C;
        else % control term; place in Dzu
            cols = find(tmpterm.u==PDE.u_tab(:,1));
            cols = upartition(cols):upartition(cols+1)-1;
            PDE_out.ODE.Dzu(rows,cols) = tmpterm.C;
        end
    end
end
% collect parameters for y % k is a continuining index that tracks
% placement in r
for i=1:length(y)
    rows = ypartition(i):ypartition(i+1)-1;
    for j=1:length(y{i}.term)
        tmpterm = y{i}.term;
        if isfield(tmpterm,'x') % either ode or pde coeff; place in Cz or Dzr
            if ismember(tmpterm.x,odeID) % ode state; place in Cz
                cols = find(tmpterm.x==PDE.x_tab(:,1));
                cols = xpartition(cols):xpartition(cols+1)-1;
                PDE_out.ODE.Cy(rows,cols) = tmpterm.C;
            else % pde state; place in Dzr
                rrows = rpartition(k):rpartition(k+1)-1;
                PDE_out.ODE.Dyr(rows,cols) = tmpterm.C;
                if isfield(tmpterm,'I') %integral
                    l = pdetab((tmpterm.x==pdeID),3); % max derivative
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    rcols = pdetab([pdetab(:,3)==l],:);
                    idLoc = find(tmpterm.x==rcols(:,1));
                    cSum = [0 cumsum(rcols(:,2))]+1;
                    rcols = cSum(idLoc):cSum(idLoc)-1;
                    der = tmpterm.D;
                    loc = sum(N:-1:N-der+1) + l + 1; % find this
                    PDE_out.PDE.Crp{loc}.coeff(rrows,rcols) = tmpterm.C;
                else % boundary value
                    del = poly2double(tmpterm.loc(2)); l = pdetab((tmpterm.x==pdeID),3); % max derivative
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    der = tmpterm.D;
                    loc = del*(np_all_derivatives-N-1) + sum(N-1:-1:N-der) + l;
                    PDE_out.PDE.Drb{loc}.coeff(rrows,rcols) = tmpterm.C;
                end
                k = k+1;
            end
        elseif isfield(tmpterm,'w') % dist term; place in Dzw
            cols = find(tmpterm.w==PDE.w_tab(:,1));
            cols = wpartition(cols):wpartition(cols+1)-1;
            PDE_out.ODE.Dyw(rows,cols) = tmpterm.C;
        else % control term; place in Dzu
            cols = find(tmpterm.u==PDE.u_tab(:,1));
            cols = upartition(cols):upartition(cols+1)-1;
            PDE_out.ODE.Dyu(rows,cols) = tmpterm.C;
        end
    end
end

% find ode parameters
for i=1:length(x)
    if odestate(i) % only go through ODE dynamics for now
        rows = xpartition(i):xpartition(i+1)-1;
        for j=1:length(x{i}.term)
            tmpterm = x{i}.term;
            if isfield(tmpterm,'x') % either ode or pde coeff; place in Cz or Dzr
                if ismember(tmpterm.x,odeID) % ode state; place in Cz
                    cols = find(tmpterm.x==PDE.x_tab(:,1));
                    cols = xpartition(cols):xpartition(cols+1)-1;
                    PDE_out.ODE.A(rows,cols) = tmpterm.C;
                else % pde state; place in Dzr and PDE.Crp/Drb
                    rrows = rpartition(k):rpartition(k+1)-1;
                    PDE_out.ODE.Bxr(rows,cols) = tmpterm.C;
                    if isfield(tmpterm,'I') %integral
                        l = pdetab((tmpterm.x==pdeID),3); % max derivative
                        if ~isfield(tmpterm,'D')
                            tmpterm.D = 0;
                        end
                        rcols = pdetab([pdetab(:,3)==l],:);
                        idLoc = find(tmpterm.x==rcols(:,1));
                        cSum = [0 cumsum(rcols(:,2))]+1;
                        rcols = cSum(idLoc):cSum(idLoc)-1;
                        der = tmpterm.D;
                        loc = sum(N:-1:N-der+1) + l + 1; % find this
                        PDE_out.PDE.Crp{loc}.coeff(rrows,rcols) = tmpterm.C;
                    else % boundary value
                        del = poly2double(tmpterm.loc(2)); l = pdetab((tmpterm.x==pdeID),3); % max derivative
                        if ~isfield(tmpterm,'D')
                            tmpterm.D = 0;
                        end
                        der = tmpterm.D;
                        loc = del*(np_all_derivatives-N-1) + sum(N-1:-1:N-der) + l;
                        PDE_out.PDE.Drb{loc}.coeff(rrows,rcols) = tmpterm.C;
                    end
                    k = k+1;
                end
            elseif isfield(tmpterm,'w') % dist term; place in Dzw
                cols = find(tmpterm.w==PDE.w_tab(:,1));
                cols = wpartition(cols):wpartition(cols+1)-1;
                PDE_out.ODE.Bxw(rows,cols) = tmpterm.C;
            else % control term; place in Dzu
                cols = find(tmpterm.u==PDE.u_tab(:,1));
                cols = upartition(cols):upartition(cols+1)-1;
                PDE_out.ODE.Bxu(rows,cols) = tmpterm.C;
            end
        end
    end
end

% now the PDE dynamics
k=1; % reset to track rows in signal v
for i=1:length(x)
    if ~odestate(i) % only go through PDE dynamics now
        for j=1:length(x{i}.term)
            tmpterm = x{i}.term{j};
            vrows = vpartition(k):vpartition(k+1)-1;
            if isfield(tmpterm,'w') % add to v signal Dvw
                vcols = find(tmpterm.w==PDE.w_tab(:,1));
                vcols = wpartition(vcols):wpartition(vcols+1)-1;
                PDE_out.PDE.Dvw(vrows,vcols) = tmpterm.C;
                rows = Xpartition(i-n.nx):Xpartition(i-n.nx+1)-1;
                PDE_out.PDE.Bpv(rows,vcols) = eye(size(vcols,1));
                k=k+1;
            elseif isfield(tmpterm,'u') % add to v signal Dvu
                vcols = find(tmpterm.u==PDE.u_tab(:,1));
                vcols = wpartition(vcols):wpartition(vcols+1)-1;
                PDE_out.PDE.Dvu(vrows,vcols) = tmpterm.C;
                rows = Xpartition(i-n.nx):Xpartition(i-n.nx+1)-1;
                PDE_out.PDE.Bpv(rows,vcols) = eye(size(vcols,1));
                k=k+1;
            elseif isfield(tmpterm,'x')&& ismember(tmpterm.x,odeID) % add to v signal Cv
                vcols = find(tmpterm.x==odeID);
                vcols = wpartition(vcols):wpartition(vcols+1)-1;
                PDE_out.PDE.Cv(vrows,vcols) = tmpterm.C;
                rows = Xpartition(i-n.nx):Xpartition(i-n.nx+1)-1;
                PDE_out.PDE.Bpv(rows,vcols) = eye(size(vcols,1));
                k=k+1;
            else % add to A/Bpb
                if % integral add to A
                else % boundary add to Bpb
                end
            end
        end
    end
end


% interconnection signals
PDE_out.ODE.Cv = ;
PDE_out.ODE.Dvw = ;
PDE_out.ODE.Dvu = ;
PDE_out.ODE.Dvr = ;

% PDE dynamics
PDE_out.PDE.A = ;
PDE_out.PDE.Bpv = ;
PDE_out.PDE.Bpb = ;

% PDE BCs
PDE_out.BC.Ebb = ;
PDE_out.BC.Ebp = ;
PDE_out.BC.Ebv = ;


% finally, use old converter to get the PIE
PIE = convert_PIETOOLS_PDE_terms_legacy(PDE_out);
end