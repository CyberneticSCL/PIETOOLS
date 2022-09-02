function PIE = convert_PIETOOLS_PDE_terms(PDE)
% For now convert the new terms format to old terms format

% extract standard information from the new terms format (domain,
% dimension, vars, etc.)
dom = PDE.dom;
xtab = PDE.x_tab;
odestate = [xtab(:,3)==0];
odeID = xtab(odestate,1);
pdeID = xtab(~odestate,1);
difforder = xtab(:,4);
pdesize = xtab(:,2);
for i=0:max(difforder)
    n.n_pde(i+1) = sum(pdesize(difforder==i));
end
clear xtab pdesize difforder;
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
xpartition = [0 cumsum(PDE.x_tab(odestae,2))]+1;
rpartition = 1;
vpartition = 1;

% calculate nr and nv the interconnection signals
nr = 0; nv = 0;
x = PDE.x; z = PDE.z; y = PDE.y; BC = PDE.BC;
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
k= 1;
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
                    l = PDE.x_tab(find(tmpterm.x==PDE.x_tab(:,1)),3); % max derivative
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    rcols = cumsum( PDE.x_tab([PDE.x_tab(:,3)==l],2) )+1:cumsum();
                    der = tmpterm.D;
                    loc = sum(N:-1:N-der+1) + l + 1; % find this
                    PDE_out.PDE.Crp{loc}.coeff(rrows,rcols) = tmpterm.C;
                else % boundary value
                    PDE_out.PDE.Drb = ;
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
                cols = rpartition(k):rpartition(k+1)-1;
                PDE_out.ODE.Dyr(rows,cols) = tmpterm.C;
                if %integral
                    PDE_out.PDE.Crp = ;
                else % boundary value
                    PDE_out.PDE.Drb = ;
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

% first collect ODE parameters
PDE_out.ODE.A = ;
PDE_out.ODE.Bxw = ;
PDE_out.ODE.Bxu = ;
PDE_out.ODE.Bxr = ;

% y parameters
PDE_out.ODE.Cy = ;
PDE_out.ODE.Dyw = ;
PDE_out.ODE.Dyu = ;
PDE_out.ODE.Dyr = ;

% interconnection signals
PDE_out.ODE.Cv = ;
PDE_out.ODE.Dvw = ;
PDE_out.ODE.Dvu = ;
PDE_out.ODE.Dvr = ;

PDE_out.PDE.Crp = ;
PDE_out.PDE.Drv = ;
PDE_out.PDE.Drb = ;

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