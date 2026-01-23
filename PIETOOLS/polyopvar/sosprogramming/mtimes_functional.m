function Cop = mtimes_functional(Kop,Fx,cntr)
% COP = MTIMES_FUNCTIONAL(KOP,FX) computes the composition of a functional 
% operator,
%   KOP*w = sum_{i=1}^{m} int_{a}^{b} ... int_{t_{ord(i,d-1)}}^{b} 
%               K{i}(t1,...,td)*w(t1,...,td) dt_{ord(i,d)} ... dt_{ord(i,1)}
% with a tensor-PI operator
%   (F(x))(s_{1},...,s_{d}) = F{1}(x)(s_{1})*...*F{d}(x)(s_{n})*x(s_{n+1})*...*x(s_{d})
% where Z(x) is a single monomial term

% - cntr:   internal argument to keep track of which call to the function
%           this is. Should not be specified by user!


% NOTE: If Fx has n elements, but Kop integrates over d>n variables, we
% assume the last d-n factors to be identity

Kcell = Kop.params;
idx_mat = Kop.omat;
Kvars = Kop.vars;
Kdom = Kop.dom;

if ~isa(Fx,'cell')
    error("Factors to integrate must be specified as a cell of 'polyopvar' objects.")
elseif numel(Fx)>numel(Kvars)
    error("Number of factors should not exceed the number of variables of integration.")
end
if nargin<=2
    cntr = numel(Fx);
% else
%     if ~numel(is_id)==numel(Fx)
%         error("Number of element of 'is_id' array must match number of elements of 'Fx'.")
%     end
%     is_id = reshape(is_id,[1,numel(Fx)]);
%     %nvars = numel(Fx)-sum(is_id);
end


% Loop over all the terms in the functional defined by Kop
for ii=1:size(idx_mat,1)
    % Extract the order of the dummy variables for the considered term
    idx_ii = idx_mat(ii,:);

    % Find the largest dummy variable for which we have a non-identity
    % factor
    midx = idx_ii(cntr);
    var_max = Kvars(midx);


    % Establish the limits of integration: 
    %       var(idx_ii(m-1)) <= var(m) <= var(idx_ii(m+1))
    if cntr==numel(idx_ii)
        U = Kdom(2);
    else
        U = Kvars(idx_ii(cntr+1));
    end
    if cntr==1
        L = Kdom(1);
    else
        L = Kvars(idx_ii(cntr-1));
    end
    dom_ii = [L,U];

    % Perform the composition of the integral operator with the factor
    % along variable var_max
    Fii = Fx{midx};
    Cii = quad2lin_term(1,Kcell{ii},Fii,dom_ii,var_max);

    % Augment the functional to include the full list of variables
    vars_full = [Kvars(1:midx-1,:); Cii.vars; Kvars(midx+1:end,:)];   % add the variables of Cii to Kop
    nvars_ii = numel(Cii.vars);    
    idx_ii(idx_ii>midx) = idx_ii(idx_ii>midx) + nvars_ii-1; % relabel some of the old variables
    omat_ii = Cii.omat + midx-1;                               % relabel the newly introduced variables
    omat_full = [repmat(idx_ii(1:cntr-1),[size(omat_ii,1),1]), omat_ii, repmat(idx_ii(cntr+1:end),[size(omat_ii,1),1])];
    Cii.vars = vars_full;
    Cii.omat = omat_full;
    Cii.dom = Kdom;

    % Apply the augmented functional to the remaining factors
    if cntr>1
        Cii = mtimes_functional(Cii,Fx,cntr-1);
    end

    

    % Add the kernels to the functional operator Cop
    if ii==1
        % Make sure the kernels are sorted in a consistent order
        omat = Cii.omat;
        [omat_new,idcs] = sortrows_integerTable(omat);
        Cii.omat = omat_new;
        Cii.params = Cii.params(idcs);
        Cop = Cii;
    else
        if ~isequal(Cop.vars,Cii.vars)
            error("Something is going wrong, the variables don't match...")
        end
        omat = [Cop.omat;Cii.omat];
        [P,omat_unique] = uniquerows_integerTable(omat);    % P*omat_unique = omat
        nC_old = size(Cop.omat,1);
        old2new_idcs = P*(1:size(omat_unique,1))';
        Ccell = cell(size(omat_unique,1),1);
        Ccell(old2new_idcs(1:nC_old)) = Cop.params;
        for jj=1:numel(Cii.params)
            idx = old2new_idcs(jj+nC_old);
            if isempty(Ccell{idx})
                Ccell(idx) = Cii.params(jj);
            else
                Ccell{idx} = Ccell{idx} + Cii.params{jj};
            end
        end
        Cop.omat = omat_unique;
        Cop.params = Ccell;
    end
end
% Finally, deal with possible multiplier term,
%       int_{a}^{b} Kop{3}(s)*(Rop1*u)(s)*(Rop2*u)(s) ds
if numel(Kcell)>size(idx_mat,1)
    if numel(Fx)~=2
        var1 = polynomial(Kcell{end}.varname);
        C_mult = quad2lin_term(Kcell{end},Fx{1},Fx{2},dom_ii,var1,Kvars);
    end

    % Make sure the kernels are sorted in a consistent order
    omat = C_mult.omat;
    [omat_new,idcs] = sortrows_integerTable(omat);
    C_mult.omat = omat_new;
    C_mult.params = C_mult.params(idcs);
    for jj=1:numel(Cop.params)
        Cop.params{jj} = Cop.params{jj} + C_mult.params{jj};
    end
end


end