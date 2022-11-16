function PIE = convert_PIETOOLS_PDE_terms(PDE)
% For now convert the new terms format to old terms format
evalin('base','silent_initialize_pde = true;')

% extract standard information from the new terms format (domain,
% dimension, vars, etc.)
dom = PDE.dom;

% first lets separate ode and pde, rearrange pde by differentiability order
% and put them back together
odetab = PDE.x_tab(PDE.x_tab(:,3)==0,:);
pdetab = PDE.x_tab(PDE.x_tab(:,3)~=0,:);
[~,sortidx] = sort(pdetab(:,3),1); % sort pdes by differentiability order
odeID = odetab(:,1); % find the ODE names
pdeID = pdetab(sortidx,1); % find the PDE names in sorted order
[~,reorder] = ismember([odeID;pdeID],PDE.x_tab(:,1)); % find the permutation to reorder PDE.x

% now we reorder PDE.x and reobtain odetab/pdetab/odeID/pdeID in a sorted
% order
PDE.x = PDE.x(reorder); PDE.x_tab = PDE.x_tab(reorder,:);
x = PDE.x; z = PDE.z; y = PDE.y; BC = PDE.BC;

odetab = PDE.x_tab(PDE.x_tab(:,3)==0,:); % find the ODE table
odeID = odetab(:,1);
pdetab = PDE.x_tab(PDE.x_tab(:,3)~=0,:); % find the PDE table
pdeID = pdetab(:,1);
xtab = PDE.x_tab;
odestate = xtab(:,3)==0; % logical array that says whether ith state is ODE or not

% next we find n_pde from the legacy format using differentiability order
difforder = pdetab(:,4);
pdesize = pdetab(:,2);
for i=0:max(difforder)
    n.n_pde(i+1) = sum(pdesize(difforder==i));
end
N = length(n.n_pde)-1; % highest differentiability order

% find frequently used dimensions
np_all_derivatives = sum(1:N+1);
% np_all_derivatives = sum((1:N+1).*n.n_pde);


% next find other dimensions to put in PDE.n of legacy format
n.nx = sum(xtab(:,2).*odestate); % ode size
n.nz = sum(PDE.z_tab(:,2)); % z size
n.ny = sum(PDE.y_tab(:,2)); % y size
n.nw = sum(PDE.w_tab(:,2)); % w size
n.nu = sum(PDE.u_tab(:,2)); % u size


% find partitioning within the signal, if signal lengths are [2,3,4] then
% the partitioning would be [1,3,6,11]
wpartition = [0 cumsum(PDE.w_tab(:,2))']+1;
upartition = [0 cumsum(PDE.u_tab(:,2))']+1;
zpartition = [0 cumsum(PDE.z_tab(:,2))']+1;
ypartition = [0 cumsum(PDE.y_tab(:,2))']+1;
xpartition = [0 cumsum(odetab(:,2))']+1;
Xpartition = [0 cumsum(pdetab(:,2))']+1;
rpartition = 1; % these are yet to be found
vpartition = 1; % these are yet to be found
BCpartition = 1;
for i=1:length(BC)
    BCpartition(end+1) = BCpartition(end)+size(BC{i}.term{1}.C,1);
end

% now we go through each term to find signals the need to go into v or r
% interconnection signal. We will store them in order. For r signals, the
% order is z, y, xdot. For v signals, the order in Xdot, BC
% calculate nr and nv the interconnection signals
nr = 0; nv = 0;

for i=1:length(z) % look for nr terms in z
    for j=1:length(z{i}.term)
        nrlen = size(z{i}.term{j}.C,1);
        if isfield(z{i}.term{j},'x')&& ismember(z{i}.term{j}.x,pdeID) % add to nr only if its a PDE state
            nr=nr+nrlen;
            rpartition(end+1) = rpartition(end)+nrlen;
        end
    end
end
for i=1:length(y) % look for nr terms in y
    for j=1:length(y{i}.term)
        nrlen = size(y{i}.term{j}.C,1);
        if isfield(y{i}.term{j},'x')&& ismember(y{i}.term{j}.x,pdeID) % add to nr only if its a PDE state
            nr=nr+nrlen;
            rpartition(end+1) = rpartition(end)+nrlen;
        end
    end
end
for i=1:length(x)
    if odestate(i) % look only in ode dynamics; look for nr terms
        for j=1:length(x{i}.term)
            nrlen = size(x{i}.term{j}.C,1);
            if isfield(x{i}.term{j},'x')&& ismember(x{i}.term{j}.x,pdeID)
                nr=nr+nrlen;
                rpartition(end+1) = rpartition(end)+nrlen;
            end
        end
    else % it is pde dynamics; look for nv terms
        for j=1:length(x{i}.term)
            nvlen = size(x{i}.term{j}.C,2);
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
        nvlen = size(BC{i}.term{j}.C,2);
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


% place all the parameters in an appropriate structure
PDE_out.n = n;
PDE_out.dom = dom;
PDE_out = initialize_PIETOOLS_PDE_terms_legacy(PDE_out); % this will make everything to be zeros of correct dimensions

% next, place the coefficients from PDE in the correct location of PDE_out
% we will go in the same order as finding r and v signals, that way
% rpartition and vpartition should match the expected signal lengths
% first z, y, xdot. Then Xdot, BC

% Start finding parameters for z
k= 1; % we use counter k for rows in r signal (increase k by one when next rpartition is needed)
for i=1:length(z)
    rows = zpartition(i):zpartition(i+1)-1; % find which rows the coefficients will go to
    for j=1:length(z{i}.term) % now go term by term
        tmpterm = z{i}.term{j};
        if isfield(tmpterm,'x') % either ode or pde coeff; place in Cz or Dzr
            if ismember(tmpterm.x,odeID) % ode state; place in Cz
                cols = find(tmpterm.x==odetab(:,1));
                cols = xpartition(cols):xpartition(cols+1)-1;
                PDE_out.ODE.Cz(rows,cols) = tmpterm.C;
            else % pde state; place in Dzr and simutaneously find Crp/Drb
                rrows = rpartition(k):rpartition(k+1)-1;
                maxder = pdetab((tmpterm.x==pdeID),4); % max derivative of the current PDE term
                rcols = pdetab(pdetab(:,4)==maxder,:); % separate only the pde states with given maxder
                idLoc = find(tmpterm.x==rcols(:,1)); % find where in those pde states the current state is
                cSum = [0 cumsum(rcols(:,2))']+1;
                rcols = cSum(idLoc):cSum(idLoc+1)-1;
                PDE_out.ODE.Dzr(rows,rrows) = eye(length(rows)); % columns will be same as rrows
                if isfield(tmpterm,'I') && ~isempty(tmpterm.I{1}) % integral term, goes to Crp
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    der = tmpterm.D; % derivative of current term
                    loc = sum(N:-1:N-der+1) + maxder + 1;
                    PDE_out.PDE.Crp{loc}.coeff(rrows,rcols) = tmpterm.C;
                else % boundary value term goes to Drb
                    del = (double(tmpterm.loc(1))==dom(2)); % location value
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    der = tmpterm.D;
                    loc = del*(np_all_derivatives-N-1) + sum(N-1:-1:N-der) + maxder;
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
    rows = ypartition(i):ypartition(i+1)-1; % find which rows the coefficients will go to
    for j=1:length(y{i}.term) % now go term by term
        tmpterm = y{i}.term{j};
        if isfield(tmpterm,'x') % either ode or pde coeff; place in Cy or Dyr
            if ismember(tmpterm.x,odeID) % ode state; place in Cy
                cols = find(tmpterm.x==odetab(:,1));
                cols = xpartition(cols):xpartition(cols+1)-1;
                PDE_out.ODE.Cy(rows,cols) = tmpterm.C;
            else % pde state; place in Dyr and simutaneously find Crp/Drb
                rrows = rpartition(k):rpartition(k+1)-1;
                maxder = pdetab((tmpterm.x==pdeID),4); % max derivative of the current PDE term
                rcols = pdetab(pdetab(:,4)==maxder,:); % separate only the pde states with given maxder
                idLoc = find(tmpterm.x==rcols(:,1)); % find where in those pde states the current state is
                cSum = [0 cumsum(rcols(:,2))']+1;
                rcols = cSum(idLoc):cSum(idLoc+1)-1;
                PDE_out.ODE.Dyr(rows,rrows) = eye(length(rows)); % columns will be same as rrows
                if isfield(tmpterm,'I') && ~isempty(tmpterm.I{1}) % integral term, goes to Crp
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    der = tmpterm.D; % derivative of current term
                    loc = sum(N:-1:N-der+1) + maxder + 1;
                    PDE_out.PDE.Crp{loc}.coeff(rrows,rcols) = tmpterm.C;
                else % boundary value term goes to Drb
                    del = (double(tmpterm.loc(1))==dom(2)); % location value
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    der = tmpterm.D;
                    loc = del*(np_all_derivatives-N-1) + sum(N-1:-1:N-der) + maxder;
                    PDE_out.PDE.Drb{loc}.coeff(rrows,rcols) = tmpterm.C;
                end
                k = k+1;
            end
        elseif isfield(tmpterm,'w') % dist term; place in Dyw
            cols = find(tmpterm.w==PDE.w_tab(:,1));
            cols = wpartition(cols):wpartition(cols+1)-1;
            PDE_out.ODE.Dyw(rows,cols) = tmpterm.C;
        else % control term; place in Dyu
            cols = find(tmpterm.u==PDE.u_tab(:,1));
            cols = upartition(cols):upartition(cols+1)-1;
            PDE_out.ODE.Dyu(rows,cols) = tmpterm.C;
        end
    end
end
% find ode parameters, from xdot equations
for i=1:length(x)
    if odestate(i) % only go through ODE dynamics for now
        rows = xpartition(i):xpartition(i+1)-1; % find which rows the coefficients will go to
        for j=1:length(x{i}.term) % now go term by term
            tmpterm = x{i}.term{j};
            if isfield(tmpterm,'x') % either ode or pde coeff; place in A or Bxr
                if ismember(tmpterm.x,odeID) % ode state; place in A
                    cols = find(tmpterm.x==odetab(:,1));
                    cols = xpartition(cols):xpartition(cols+1)-1;
                    PDE_out.ODE.A(rows,cols) = tmpterm.C;
                else % pde state; place in Br and simutaneously find Crp/Drb
                    rrows = rpartition(k):rpartition(k+1)-1;
                    maxder = pdetab((tmpterm.x==pdeID),4); % max derivative of the current PDE term
                    rcols = pdetab(pdetab(:,4)==maxder,:); % separate only the pde states with given maxder
                    idLoc = find(tmpterm.x==rcols(:,1)); % find where in those pde states the current state is
                    cSum = [0 cumsum(rcols(:,2))']+1;
                    rcols = cSum(idLoc):cSum(idLoc+1)-1;
                    PDE_out.ODE.Bxr(rows,rrows) = eye(length(rows)); % columns will be same as rrows
                    if isfield(tmpterm,'I') && ~isempty(tmpterm.I{1}) % integral term, goes to Crp
                        if ~isfield(tmpterm,'D')
                            tmpterm.D = 0;
                        end
                        der = tmpterm.D; % derivative of current term
                        loc = sum(N:-1:N-der+1) + maxder + 1;
                        PDE_out.PDE.Crp{loc}.coeff(rrows,rcols) = tmpterm.C;
                    else % boundary value term goes to Drb
                        del = (double(tmpterm.loc(1))==dom(2)); % location value
                        if ~isfield(tmpterm,'D')
                            tmpterm.D = 0;
                        end
                        der = tmpterm.D;
                        loc = del*(np_all_derivatives-N-1) + sum(N-1:-1:N-der) + maxder;
                        PDE_out.PDE.Drb{loc}.coeff(rrows,rcols) = tmpterm.C;
                    end
                    k = k+1;
                end
            elseif isfield(tmpterm,'w') % dist term; place in Bxw
                cols = find(tmpterm.w==PDE.w_tab(:,1));
                cols = wpartition(cols):wpartition(cols+1)-1;
                PDE_out.ODE.Bxw(rows,cols) = tmpterm.C;
            else % control term; place in Bxu
                cols = find(tmpterm.u==PDE.u_tab(:,1));
                cols = upartition(cols):upartition(cols+1)-1;
                PDE_out.ODE.Bxu(rows,cols) = tmpterm.C;
            end
        end
    end
end

% now the PDE dynamics, Xdot equations
k=1; % reset to track rows in signal v instead of signal r
for i=1:length(x)
    if ~odestate(i) % only go through PDE dynamics now
        for j=1:length(x{i}.term)
            tmpterm = x{i}.term{j};
            rows = Xpartition(i-sum(odestate)):Xpartition(i-sum(odestate)+1)-1; % remove odestate because first i elements should now be only ode dynamics
            if isfield(tmpterm,'x')&& ismember(tmpterm.x,odeID) % add to Bpv and add to v signal Cv
                vrows = vpartition(k):vpartition(k+1)-1;
                PDE_out.PDE.Bpv(rows,vrows) = tmpterm.C; % columns are same as vrows
                vcols = find(tmpterm.x==odetab(:,1));
                vcols = xpartition(vcols):xpartition(vcols+1)-1;
                PDE_out.ODE.Cv(vrows,vcols) = eye(length(vrows));
                k=k+1;
            elseif isfield(tmpterm,'w') % add to v signal Dvw and Bpv
                vrows = vpartition(k):vpartition(k+1)-1;
                PDE_out.PDE.Bpv(rows,vrows) = tmpterm.C; % columns are same as vrows
                vcols = find(tmpterm.w==PDE.w_tab(:,1));
                vcols = wpartition(vcols):wpartition(vcols+1)-1;
                PDE_out.ODE.Dvw(vrows,vcols) = eye(length(vrows));
                k=k+1;
            elseif isfield(tmpterm,'u') % add to v signal Dvu and Bpv
                vrows = vpartition(k):vpartition(k+1)-1;
                PDE_out.PDE.Bpv(rows,vrows) = tmpterm.C; % columns are same as vrows
                vcols = find(tmpterm.u==PDE.u_tab(:,1));
                vcols = upartition(vcols):upartition(vcols+1)-1;
                PDE_out.ODE.Dvu(vrows,vcols) = eye(length(vrows));
                k=k+1;
            else % add to A/Bpb
                lstate = xtab(i,4); % max derivative of LHS state
                rstate = pdetab((tmpterm.x==pdeID),4); % max derivative of RHS term state
                rcols = pdetab(pdetab(:,4)==rstate,:);
                idLoc = find(tmpterm.x==rcols(:,1));
                cSum = [0 cumsum(rcols(:,2))']+1;
                rcols = cSum(idLoc):cSum(idLoc+1)-1;
                rrows = pdetab(pdetab(:,4)==lstate,:);
                idLoc = find(tmpterm.x==rrows(:,1));
                cSum = [0 cumsum(rrows(:,2))']+1;
                rrows = cSum(idLoc):cSum(idLoc+1)-1;
                if ~poly2double(tmpterm.loc(1))% integral add to A
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    der = tmpterm.D; % derivative order of current term
                    if ~isempty(tmpterm.I{1})
                        interval = tmpterm.I{1}'
                        if poly2double(interval(1))%R1 term
                            loc=(N+1)*np_all_derivatives + (lstate)*np_all_derivatives + sum(N:-1:N-der+1) + rstate + 1;
                        end
                        if poly2double(interval(2))%R2 term
                            loc = 2*(N+1)*np_all_derivatives + (lstate)*np_all_derivatives + sum(N:-1:N-der+1) + rstate + 1;
                        end
                    else% R0 term
                        loc = (lstate)*np_all_derivatives + sum(N:-1:N-der+1) + rstate + 1;
                    end
                    PDE_out.PDE.A{loc}.coeff(rrows,rcols) = tmpterm.C;
                else % boundary add to Bpb
                    if ~isfield(tmpterm,'D')
                        tmpterm.D = 0;
                    end
                    der = tmpterm.D;
                    deltaval = (double(tmpterm.loc(1))==dom(2));    % 0 if lower boundary, 1 if upper boundary
                    loc = lstate*(2*np_all_derivatives-2*(N+1)) + deltaval*(np_all_derivatives-N-1) + sum(N-1:-1:N-der) + rstate;
                    PDE_out.PDE.Bpb{loc}.coeff(rrows,rcols) = tmpterm.C;
                end
            end
        end
    end
end


% PDE BCs
for i=1:length(BC)
    for j=1:length(BC{i}.term)
        tmpterm = BC{i}.term{j};
        rows = BCpartition(i):BCpartition(i+1)-1;
        if isfield(tmpterm,'x')&& ismember(tmpterm.x,odeID) % add to Ebv and add to v signal Cv
            vrows = vpartition(k):vpartition(k+1)-1;
            PDE_out.BC.Ebv(rows,vrows) = tmpterm.C; % columns are same as vrows
            vcols = find(tmpterm.x==odetab(:,1));
            vcols = xpartition(vcols):xpartition(vcols+1)-1;
            PDE_out.ODE.Cv(vrows,vcols) = eye(length(vrows));
            k=k+1;
        elseif isfield(tmpterm,'w') % add to v signal Dvw and Ebv
            vrows = vpartition(k):vpartition(k+1)-1;
            PDE_out.BC.Ebv(rows,vrows) = tmpterm.C; % columns are same as vrows
            vcols = find(tmpterm.w==PDE.w_tab(:,1));
            vcols = wpartition(vcols):wpartition(vcols+1)-1;
            PDE_out.ODE.Dvw(vrows,vcols) = eye(length(vrows));
            k=k+1;
        elseif isfield(tmpterm,'u') % add to v signal Dvu and Ebv
            vrows = vpartition(k):vpartition(k+1)-1;
            PDE_out.BC.Ebv(rows,vrows) = tmpterm.C; % columns are same as vrows
            vcols = find(tmpterm.u==PDE.u_tab(:,1));
            vcols = upartition(vcols):upartition(vcols+1)-1;
            PDE_out.ODE.Dvu(vrows,vcols) = eye(length(vrows));
            k=k+1;
        else % add to Ebp/Ebb
            rstate = pdetab((tmpterm.x==pdeID),4); % max derivative of RHS term state
            rcols = pdetab(pdetab(:,4)==rstate,:);
            idLoc = find(tmpterm.x==rcols(:,1));
            cSum = [0 cumsum(rcols(:,2))']+1;
            rcols = cSum(idLoc):cSum(idLoc+1)-1;
            if ~poly2double(tmpterm.loc(1))% integral add to Ebp
                if ~isfield(tmpterm,'D')
                    tmpterm.D = 0;
                end
                der = tmpterm.D; % derivative order of current term
                loc=sum(N:-1:N-der+1) + rstate + 1;
                PDE_out.BC.Ebp{loc}.coeff(rows,rcols) = tmpterm.C;
            else % boundary value, add to Ebb
                if ~isfield(tmpterm,'D')
                    tmpterm.D = 0;
                end
                der = tmpterm.D;
                deltaval = (double(tmpterm.loc(1))==dom(2));    % 0 if lower boundary, 1 if upper boundary
                loc=deltaval*(np_all_derivatives-N-1) + sum(N-1:-1:N-der) + rstate;
                PDE_out.BC.Ebb{loc}.coeff(rows,rcols) = tmpterm.C;
            end
        end
    end
end



% finally, use old converter to get the PIE
PIE = convert_PIETOOLS_PDE_terms_legacy(PDE_out);
PIE = pie_struct(PIE);
end