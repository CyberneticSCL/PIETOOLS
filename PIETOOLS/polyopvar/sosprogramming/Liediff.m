function dV = Liediff(V,PIE)
% DV = LIEDIFF(V,PIE) computes the Lie derivative of a polynomial
% functional V along the PIE defined by PIE
%
% NOTES
% Each variable in V should correspond to a row in the PIE.
% Variable x_{i} in V will be replaced by T(i,:)*x in evaluating the Lie
% derivative of V along the PIE

Top = PIE.T;
RHS = PIE.f;
Vdegs = V.degmat;
fdegs = RHS.degmat;
nvars = numel(V.varname);

if ~isequal(RHS.varname,V.varname)
    error("State variables in PIE and functional must match.")
end

% Decompose the T operator into individual operators acting on each state
% variable (this is really not ideal, but the only option for now...)
Top_cell = cell(nvars,nvars);
if isa(Top,'nopvar') || isa(Top,'ndopvar')
    Top_tmp = ndopvar2dopvar(Top);
else
    Top_tmp = Top;
end
for lindx=1:numel(Top_cell)
    [ridx,cidx] = ind2sub(size(Top_cell),lindx);
    Top_cell{lindx} = dopvar2ndopvar(Top_tmp(ridx,cidx));
end

% Replace the variable x by T*x in each of the monomials defining V. Since
% the monomials are expressed in terms of elements x_{i}, but T(i,:)*x may
% involve all variables again, this is somewhat cumbersome
if nvars>1
    degmat_new = [];
    old_monnum = [];
    rep_mon = [];       % Keep track of how many times the considered monomial appears
    for ii=1:size(Vdegs,1)
        degmat_ii = Vdegs(ii,:);
        dtot_ii = sum(degmat_ii);
        % Build all monomials up to degree at most dtot
        degmat_tmp = eye(nvars);
        for jj=1:dtot_ii
            degmat_tmp = kron(degmat_tmp,ones([nvars,1])) + repmat(eye(nvars),[nvars,1]);
        end
        % Remove duplicate monomials
        [P,degmat_tmp] = uniquerows_integerTable(degmat_tmp);
        rep_mon_ii = sum(P,1)';
        degmat_new = [degmat_new; degmat_new];
        old_monnum = [old_monnum; ii*ones(size(degmat_tmp,1),1)];
        rep_mon = [rep_mon; rep_mon_ii];
    end
    Vdegs = degmat_new;
else
    old_monnum = (1:size(Vdegs,1))';
    rep_mon = ones(size(Vdegs,1),1);
end

% Declare an empty polynomial 
tmp_poly = RHS;
tmp_poly.C.ops = {};
tmp_poly.degmat = zeros(0,nvars);
dV = tmp_poly;
for ii=1:size(Vdegs,1)
    % Loop over each of the monomials in V, replacing x by T*x, and
    % evaluating the derivative d/dt T*x = f(x)
    degs_ii = Vdegs(ii,:);
    Kop_ii = V.C.ops{old_monnum(ii)};
    if rep_mon(ii)>1
        Kop_ii.params = rep_mon(ii)*Kop_ii.params;
    end
    dtot_ii = sum(degs_ii);
    state_idcs_ii = zeros(1,dtot_ii);
    strt_idx = 0;
    for jj=1:numel(degs_ii)
        state_idcs_ii(strt_idx+(1:degs_ii(jj))) = jj;
        strt_idx = strt_idx + V.degmat(ii,jj);
    end
    % Apply the product rule, d/dt(Tx*Tx) = Ax*Tx + Tx*Ax, looping over
    % each factor in the monomial and taking the temporal derivative only
    % along this factor
    dV_ii = tmp_poly;
    for jj=1:dtot_ii
        % Indicate the state variable for which to take the derivative by 0
        state_idcs_jj = state_idcs_ii;
        state_idcs_jj(jj) = 0;
        Kop_jj = combine_terms(Kop_ii,state_idcs_jj);
        % Build the product (T*x)*...*(T*x)*(A*x)*(T*x*...*T*x) for A in 
        % factor j
        Fx_jj = cell(1,dtot_ii);
        for state_num = 1:nvars
            Fx_tmp = tmp_poly;
            Fx_tmp.degmat = [zeros(1,state_num-1),1,zeros(1,nvars-state_num)];
            Fx_tmp.C.ops{1} = Top_cell{ii,state_num};
            is_state_num = state_idcs_jj==state_num;
            Fx_jj{is_state_num} = Fx_tmp;
        end
        % Set factor jj as the right-hand side of the PIE
        dV_jj = tmp_poly;
        for ll=1:size(fdegs,1)
            % Loop over all terms in the PIE
            Fx_ll = RHS;
            Fx_ll.C.ops = Fx_ll.C.ops(state_idcs_ii(jj),ll);
            Fx_ll.degmat = Fx_ll.degmat(ll,:);
            Fx_jj{jj} = Fx_ll;
            % Take the composition of the functional operator with this
            % particular term
            dV_ll = mtimes_functional(Kop_jj,Fx_jj);
            dV_ll = combine_terms(dV_ll);
            dV_jj = dV_jj + dV_ll;
        end
        % Add the derivative for this particular factor to the full
        % derivative
        dV_ii = dV_ii + dV_jj;
    end
    dV = dV + dV_ii;
end

end