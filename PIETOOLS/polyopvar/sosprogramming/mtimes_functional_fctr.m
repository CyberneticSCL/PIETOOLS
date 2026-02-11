function Cop = mtimes_functional_fctr(Kfun,Ffctr,dom,var1,var2,cntr)
% Cop = MTIMES_FUNCTIONAL_FCTR(KFUN,FFCTR,DOM,VAR1,VAR2) 
% computes the composition of a functional operator,
%   KOP*y = int_{a}^{b} KFUN(s)*y(s) ds
% with a tensor-PI operator
%   (FFCTR(x))(s) = Fop{1}(x)(s) * ... * Fop{d}(x)(s)
% so that
%  KOP*Ffctr(x)(t1,...,td) = 
%   sum_{i=1}^{q} int_{a}^{b} ... int_{t_{ord(i,d-1)}}^{b} 
%               C_{i}(t1,...,td)*w(t1,...,td) dt_{ord(i,d)} ... dt_{ord(i,1)}
% where [a,b] = DOM, s = VAR1 and (t1,...,td) = VAR2;
%
% INPUTS
% - Kfun:   scalar 'polynomial' or 'dpvar' object representing the kernel
%           defining the functional operator;
% - Ffctr:  'polyopvar' representing a distributed polynomial defined by a
%           single monomial (i.e. size(Ffctr.degmat,1)=1) of degree d;
% - dom:    1 x 2 array specifying the domain of integration, [a,b]. Can
%           also be polynomial, i.e. dom = [s1,s2] for pvar objects s1, s2;
% - var1:   1x1 pvar object specifying the variable in Kfun over which to
%           integrate;
% - var2:   d x 1 'polynomial' array of pvar objects, specifying the
%           desired names of the variables (t1,...,td) defining the
%           composite functional operator;
% - cntr:   internal variable to keep track of which call to the function
%           this is. Should not be specified by user!
%
% OUTPUT
% - Cop:    'struct' representing the composition of the functional with 
%           kernel Kfun with the polynomial defined by Fx. Will have fields:
%               params: 1 x q 'polynomial' or 'dpvar' object,
%                       with each element i corresponding to 
%                       the kernels in the ith term of the functional;
%               omat:   q x d array of integers, with row i
%                       specifying the order of the variables in the
%                       integral associated with the ith kernel.
%                       Specifically, if omat(i,:) = [k1,k2,...,kd], 
%                       then term i is defined by the integral
%                           int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b};
%               vars:   d x 1 'polynomial' array specifying the names of
%                       the dummy variables used in definition of the
%                       kernels in Kop.params;
%               dom:    interval [a,b] over which to integrate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes_functional_fctr
%
% Copyright (C) 2026 PIETOOLS Team
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
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 01/25/2026: Initial coding


% Extract the tensor-PI operator which to compose with the integral defined
% by Kfun
if isa(Ffctr.C.ops{1},'nopvar')
    Ffctr.C.ops{1} = Ffctr.C.ops(1);
end

if size(Ffctr.C.ops{1},1)~=1
    % If the function involves multiple terms, we just apply the routine to
    % each term.
    Ffctr_1 = Ffctr;
    Ffctr_1.C.ops{1} = Ffctr.C.ops{1}(1,:);
    Cop = mtimes_functional_fctr(Kfun,Ffctr_1,dom,var1,var2);
    for ii=2:size(Ffctr.C.ops{1},1)
        Ffctr_ii = Ffctr;
        Ffctr_ii.C.ops{1} = Ffctr.C.ops{1}(ii,:);
        Cop_ii = mtimes_functional_fctr(Kfun,Ffctr_ii,dom,var1,var2);

        % Add functional operator associated to term ii
        omat_full = [Cop.omat;Cop_ii.omat];
        params_full = [Cop.params, Cop_ii.params];
        [omat_unique,~,sort_idx] = unique(omat_full,'rows');
        params_unique = params_full*sparse((1:numel(sort_idx)),sort_idx,1);
        Cop.params = params_unique;
        Cop.omat = omat_unique;
    end
    return
end

% We focus on a single PI operator in the tensor product defining Ffctr
if nargin<=5
    cntr = size(Ffctr.C.ops{1},2);
end
Fop = Ffctr.C.ops{1}{1,cntr};

% If Fop.dom = [a,b], but dom = [t1,t2] for t1,t2 in [a,b], we
% decompose the operator Fop into three parts:
%   a PI operator over [a,t1]
%   a PI operator over [t1,t2]
%   a PI operator over [t2,b]
is_split = false;
if ~isequal(Fop.dom(1),dom(1))
    % Extract parameters for integral over [a,t1]
    Fjj_1 = Fop.C{2}; 
    Fjj_1 = coeff2poly(Fjj_1,Fop.dim,Fop.deg,[var1,var2(cntr)]);
    is_split = true;
else
    Fjj_1 = 0;
end
if ~isequal(Fop.dom(2),dom(2))
    Fjj_2 = Fop.C{3}; 
    Fjj_2 = coeff2poly(Fjj_2,Fop.dim,Fop.deg,[var1,var2(cntr)]);
    is_split = true;
else
    Fjj_2 = 0;
end

% For the operator defined on [t1,t2], we can use the "quad2lin" routine to
% compute the composition of the functional with the operator. However, we
% may still need to deal with the other PI operators defining Ffctr.
Fop.dom = dom;
Ffctr.C.ops{1}{cntr} = Fop;
if cntr==1
    Cop = quad2lin_term_old(1,Kfun,Ffctr,dom,var1,var2);
else
    Cop = mtimes_functional_fctr(Kfun,Ffctr,dom,var1,var2,cntr-1);
end


% For the integrals outside of [t1,t2], we remove the operator Fop from
% the product
if ~is_split
    return
else
    Ffctr_rem = Ffctr;
    deg_red = find(cntr<=cumsum(Ffctr.degmat(1,:)),1,'first');
    Ffctr_rem.degmat(:,deg_red) = Ffctr_rem.degmat(:,deg_red)-1;
    Ffctr_rem.C.ops{1} = Ffctr.C.ops{1}(:,[1:cntr-1,cntr+1:end]);    % retain all other operators
    var2_rem = var2([1:cntr-1,cntr+1:end]);             % perform composition only along remaining variables
end

% Add the terms associated with the integral over [a,t1]
if ~all(all(isequal(Fjj_1,0)))
    if isempty(var2_rem)
        Cfun1 = int_simple(Kfun*Fjj_1,var1,dom(1),dom(2));
        Cop.params = [Cop.params, Cfun1];
        Cop.omat = [Cop.omat; -1];     % use -jj to indicate integral [a,t1]
    else
        if cntr==1
            Cop_1 = quad2lin_term_old(1,Kfun*Fjj_1,Ffctr_rem,dom,var1,var2_rem);
        else
            Cop_1 = mtimes_functional_fctr(Kfun*Fjj_1,Ffctr_rem,dom,var1,var2_rem,cntr-1);
        end
        % Add integral terms from Cop_1 to Cop
        Cop.params = [Cop.params, Cop_1.params];
        Cop_1.omat(abs(Cop_1.omat)>=cntr) = Cop_1.omat(abs(Cop_1.omat)>=cntr)+sign(Cop_1.omat(abs(Cop_1.omat)>=cntr));
        omat_1 = [-cntr*ones([size(Cop_1.omat,1),1]),Cop_1.omat];   % use -jj to indicate integral [a,t1]
        Cop.omat = [Cop.omat; omat_1];
    end
end

% Add the terms associated with the integral over [t2,b]
if ~all(all(isequal(Fjj_2,0)))
    if isempty(var2_rem)
        Cfun2 = int_simple(Kfun*Fjj_2,var1,dom(1),dom(2));
        Cop.params = [Cop.params, Cfun2];
        Cop.omat = [Cop.omat; 1i];     % use jj*1i to indicate integral [t2,b]
    else
        if cntr==1
            Cop_2 = quad2lin_term_old(1,Kfun*Fjj_2,Ffctr_rem,dom,var1,var2_rem);
        else
            Cop_2 = mtimes_functional_fctr(Kfun*Fjj_2,Ffctr_rem,dom,var1,var2_rem,cntr-1);
        end
        % Add integral terms from Cop_2 to Cop
        Cop.params = [Cop.params, Cop_2.params];
        Cop_2.omat(abs(Cop_2.omat)>=cntr) = Cop_2.omat(abs(Cop_2.omat)>=cntr)+sign(Cop_2.omat(abs(Cop_2.omat)>=cntr));
        omat_2 = [Cop_2.omat, 1i*cntr*ones([size(Cop_2.omat,1),1])];   % use 1i*jj to indicate integral [t2,b]
        Cop.omat = [Cop.omat; omat_2];
    end
end

end