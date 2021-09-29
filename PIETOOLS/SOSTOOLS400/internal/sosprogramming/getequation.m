function [At,b,Z] = getequation(symexpr,vartable,decvartable,varmat)
%function [At,b,Z] = getequation(symexpr,vartable,decvartable,decvartablename)
%
% GETEQUATION --- Convert a symbolic expression to At, b, and Z
%         used in an SOS program.
%
% Inputs:
% symexpr: This expression may be in dpvar, pvar, or symbolic format. 
% It may be matrix or scalar valued. At present, matrix-valued inputs should 
% only be used with equality constraints
% vartable: List of the independent variables in the sosprogram


% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 4.00.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018, 2021  
%                                      A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5),
%                                      M. Peet (6), D. Jagt (6)
% (1) Department of Engineering Science, University of Oxford, Oxford, U.K.
% (2) Laboratoire de Signaux et Systmes, CentraleSupelec, Gif sur Yvette,
%     91192, France
% (3) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (4) Aerospace and Engineering Mechanics Department, University of
%     Minnesota, Minneapolis, MN 55455-0153, USA.
% (5) Laboratory for Information and Decision Systems, M.I.T.,
%     Massachusetts, MA 02139-4307
% (6) Cybernetic Systems and Controls Laboratory, Arizona State University,
%     Tempe, AZ 85287-6106, USA.
%
% Send bug reports and feedback to: sostools@cds.caltech.edu
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
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change log and developer notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 03/01/02 - SP
% 03/10/02 - SP -- Use Maple
% 07/10/13 - JA&GV -- Use the matlab's symbolic engine
% 09/11/13 - PJS -- Addition for multipoly objects
% 6/30/21  - MP -- dpvar modifications
% 7/12/21 - SS -- fixed typo in line 59, n_vars_dp to n_vars_p (however, this value seems to be unused)


if isa(symexpr,'dpvar')
    vartable = [vartable; varmat];
    n_dvars_prog = numel(decvartable);
    n_vars_prog = numel(vartable);
    cvartable_prog = char(vartable);
    cdvartable_prog = char(decvartable);
    dvarname=symexpr.dvarname;
    varname=symexpr.varname;
    %    cdvarname=char(symexpr.dvarname);
    %    cvarname=char(symexpr.varname);
    degmat=symexpr.degmat;
    C=symexpr.C;
    m= size(symexpr,1);
    n= size(symexpr,2);
    [n_mons,n_vars_p] = size(degmat);
    n_dvars_dp = length(dvarname);
    
    %    if all([m n]==[1 1])
    % scalar case
    %    else
    % matrix case
    C_cell=mat2cell(C,(n_dvars_dp+1)*ones(m,1),n_mons*ones(n,1));
    C_reshape=cell2mat(reshape(C_cell,1,m*n));
    % eliminate empty columns
    col_nums = any(C_reshape,1); %logical vector of non-zero columns
    n_cons=sum(col_nums);
    C_reshape2=C_reshape(:,col_nums);
    
    % delete the monomials corresponding to eliminated constraints
    degmat_big=repmat(degmat,n*m,1);
    degmat_big2=degmat_big(col_nums,:); 
    
    % synchronize varnames in Z to SOS program cvartable (same as in polynomial case)
    
    % Reorder the monomials matrix with variables in the order
    % listed in cvartable and sorted monomials
    %    Z = zeros(n_mons,n_vars_dp);
    %     [~,idx]=ismember(varname,cvartable_prog);
    [~,idx]=ismember(varname,cvartable_prog); % beware the case where varname does not appear in vartable!
    if ~isempty(find(~idx))
        error('the given expression has an independent variable which does not appear in the sosprogram')
    end
%    Z(:,idx) = full(degmat_big2); % Why full?

%        Z = sortNoRepeat([],Z); % shouldn't sorting Z but not At and b defeat
    %    the purpose? not sure what the purpose is, actually.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % separate At and b
%     const_rows = 1:n_dvars_dp+1:size(C_reshape2,1);
    %[C_reshape3]=(unique(C_reshape2,'rows'));
    % b_temp=C_reshape3(:,1);
    % At_temp=C_reshape3(:,2:end)';%???
    %degmat_big3=degmat_big2(col_nums,:); 

    %    b=b_temp;
    b_temp=C_reshape2(1,:);
    At_temp=C_reshape2(2:end,:);%???
    b=b_temp';

    Z = sparse(n_cons,n_vars_prog);
    Z(:,idx) = degmat_big2;

    % synchronize dpvars to SOS program decvartable
    [~,idx]=ismember(dvarname,cdvartable_prog); % also slow??
    if ~isempty(find(~idx))
        error('the given expression has a decision variable which does not appear in the sosprogram')
    end
    %At=sparse(n_dvars_prog,n_cons,nnz(At_temp)); 
    At=sparse([],[],[],n_dvars_prog,n_cons,nnz(At_temp)); %<-- I think this is what we want?
%    At(idx,:)=At_temp;
    At(idx,:)=-At_temp; % main time sink?
    
elseif isa(symexpr,'polynomial')
    % PJS: Handle Polynomial Objects including the Matrix Case
    decvarnum = numel(decvartable);
    vartable = [vartable; varmat];
    cvartable = char(vartable);
    dimp = size(symexpr,2);
    
    % Collect symexpr = g(c)*h(x) where x are the vars in vartable,
    % c are the decision variables, and h(x) is a vector of monoms.
    % in the actual polynomial vars.  Use this to find unique
    % monomials in the polynomial variables.
    [g0,g,h] = collect(symexpr(:),setdiff(symexpr.varname,cvartable));
    g = g(:);
    if ~isequal(g0,0)
        g = [g0;g];
        h = [1; h];
    end
    [nmon,nvar] = size(h.degmat);
    
    % Reorder the monomials matrix with variables in the order
    % listed in cvartable and sorted monomials
    Z = zeros(nmon,nvar);
    [~,idx]=ismember(h.varname,cvartable);
    Z(:,idx) = full(h.degmat);
    Z = sortNoRepeat([],Z);
    
    % Initialize coefficient variables
    [nmon,nvar] = size(Z);
    coeffnts = sparse(nmon*dimp,dimp);
    if decvarnum>0
        coeffnts_decvar = polynomial(sparse(nmon*dimp,dimp));
    end
    
    % Define coefficients
    for i = 1:dimp
        for j = i:dimp
            exprij = symexpr(i,j);
            [g0,g,h] = collect(exprij,setdiff(exprij.varname,cvartable));
            g = g(:);
            if ~isequal(g0,0)
                g = [g0;g];
                h = [1; h];
            end
            coefmatr = double(subs(g,decvartable,zeros(decvarnum,1)));
            monmatr = zeros(length(h),nvar);
            [~,idx]=ismember(h.varname,cvartable);
            monmatr(:,idx) = h.coef'*h.degmat;
            
            for k = 1:length(h);
                s_ijk = coefmatr(k,1);
                s_ijk_decvar = g(k);
                mon_k= monmatr(k,:);
                
                [val,ind_k] = max(sum((Z == kron(ones(nmon,1),mon_k)),2));
                coeffnts((ind_k-1)*dimp+i,j) = s_ijk;
                coeffnts((ind_k-1)*dimp+j,i) = s_ijk;
                
                coeffnts_decvar((ind_k-1)*dimp+i,j) = s_ijk_decvar;
                coeffnts_decvar((ind_k-1)*dimp+j,i) = s_ijk_decvar;
            end
        end
    end
    
    % Constructing At and b matrices
    At = sparse(decvarnum,size(Z,1)*dimp^2);
    if decvarnum>0
        for i=1:nmon
            Mivec = reshape(coeffnts_decvar((i-1)*dimp+1:i*dimp,:),dimp^2,1);
            At(:,(i-1)*dimp^2+1:i*dimp^2) = sparse(-double(jacobian(Mivec,decvartable))');
        end
    end
    b = sparse(coeffnts);
    
else
    % PJS: Original Code to Handle Symbolic Objects
    
    if strcmp(decvartable,'[]');
        decvarnum = 0;
    else
        decvarnum = length(find(decvartable == ','))+1;
    end;
    expr = evalin(symengine,symexpr);
    %vartable = evalin(symengine,vartable);
    if  length(varmat)==2
        vartable = evalin(symengine,vartable);
    else
        vartable = evalin(symengine,[vartable(1:end-1),',',varmat(2:end)]);  %JA&GV 07/06/13
    end
    
    decvartable = evalin(symengine,decvartable);
    charvartable = converttochar([vartable]);
    
    FPexpr = feval(symengine,'expand',expr);
    
    %JA&GV 07/13 commented the below to implement the matrix case
    %FPexpr = feval(symengine,'collect',FPexpr,charvartable);
    
    dimp = size(FPexpr,2);%JA&GV jul/04/2013 sets the dimension of the matrix of polynomials
    
    for i = 1:size(FPexpr,1)
        for j = 1:dimp
            FPexpr(i,j) = feval(symengine,'collect',FPexpr(i,j),charvartable);
        end
    end
    
    if isempty(decvartable)==0%JA&GV 07/13 the code below  addresses the matrix case
        
        
        Zfull = [];
        for i = 1:dimp
            for j = i:dimp
                coefmon = feval(symengine,'poly2list',FPexpr(i,j),charvartable);
                coefmonmatr = subs(coefmon.',decvartable,zeros(1,length(decvartable)));
                Z = double(coefmonmatr(:,2:end));
                Zfull = sortNoRepeat(Zfull,Z);
            end
        end
        Z = Zfull;
        [nmon,nvar] = size(Z);
        coeffnts = sparse(nmon*dimp,dimp);
        %coeffnts = sym(sparse(nmon*dimp,dimp)); %JA edit - sparse doesnt't do anything here?
        coeffnts_decvar = sym(sparse(nmon*dimp,dimp));
        
        
        for i = 1:dimp
            for j = i:dimp
                coefmon = feval(symengine,'poly2list',FPexpr(i,j),charvartable);
                coefmonmatr = subs(coefmon.',symvar(decvartable),zeros(1,length(decvartable)));
                for k = 1:size(coefmonmatr,1)
                    s_ijk = coefmonmatr(k,1);
                    
                    dummy_decvar = reshape(coefmon(k),2,1);
                    s_ijk_decvar = dummy_decvar;
                    
                    mon_k= double(coefmonmatr(k,2:end));
                    [val,ind_k] = max(sum((Zfull == kron(ones(nmon,1),mon_k)),2));
                    coeffnts((ind_k-1)*dimp+i,j) = s_ijk;
                    coeffnts((ind_k-1)*dimp+j,i) = s_ijk;
                    
                    coeffnts_decvar((ind_k-1)*dimp+i,j) = s_ijk_decvar(1);%THIS HAS BEEN CHANGED TO ASSIGN JUST THE FIRST ENTRY OF THE VECTOR GV 17jun2015
                    coeffnts_decvar((ind_k-1)*dimp+j,i) = s_ijk_decvar(1);
                    
                    
                end
                
            end
        end
        
        
    else
        
        
        Zfull = [];
        for i = 1:dimp
            for j = i:dimp
                coefmon = feval(symengine,'poly2list',FPexpr(i,j),charvartable);
                coefmon = coefmon.';
                for k = 1:length(coefmon)
                    dummyvar = reshape(coefmon(k),2,1);
                    Z(k,:) = double(dummyvar(2));
                end
                Zfull = sortNoRepeat(Zfull,Z);
                Z = [];
            end
        end
        Z = sparse(Zfull);
        [nmon,nvar] = size(Z);
        coeffnts = sparse(nmon*dimp,dimp);
        
        for i = 1:dimp
            for j = i:dimp
                coefmon = feval(symengine,'poly2list',FPexpr(i,j),charvartable);
                coefmon = coefmon.';
                for k = 1:length(coefmon)
                    dummyvar = reshape(coefmon(k),2,1);
                    s_ijk = double(dummyvar(1));
                    mon_k= double(dummyvar(2));
                    [val,ind_k] = max(sum((Zfull == kron(ones(nmon,1),mon_k)),2));
                    coeffnts((ind_k-1)*dimp+i,j) = s_ijk;
                    coeffnts((ind_k-1)*dimp+j,i) = s_ijk;
                end
            end
        end
        
    end
    
    % Constructing At and b matrices
    At = sparse(decvarnum,size(Z,1)*dimp^2);
    
    %JA&GV may/2013 uses the jacobian function to derive the expression
    if isempty(decvartable)==0
        for i = 1:size(Z,1)
            Mivec = reshape(coeffnts_decvar((i-1)*dimp+1:i*dimp,:),dimp^2,1);
            At(:,(i-1)*dimp^2+1:i*dimp^2) = sparse(-double(jacobian(Mivec,decvartable))');
        end
    end
    
    %b = sparse(coeffnts);
    b = coeffnts;
end


