function C = mtimes(E,F)
if ~isa(E,'quadpoly')
    E = quadpoly(E);
end
if ~isa(F,'quadpoly') 
    F = quadpoly(F);
end
if (E.isdecVar)&& (F.isdecVar) 
    error('Cannot multiply two decision polynomials');
end

if F.isdecVar % right object is decision var --> multiply transposes instead
    C_transpose = F.'*E.';
    C = C_transpose.';
else  % left object is decision var (or neither is decision var)
    Eldegmat=E.ldegmat;       Fldegmat=F.ldegmat;
    Erdegmat=E.rdegmat;       Frdegmat=F.rdegmat; % we assume there are no repeats in degmat.
    Elvarname = E.lvarname;  Flvarname = F.lvarname;
    Ervarname = E.rvarname;  Frvarname = F.rvarname;
    ntE=size(Erdegmat,1);    ntF=size(Frdegmat,1);
    
    % synchronize the pvarnames
    [pvarnames_all,IE_p,~] = union(Elvarname,Flvarname,'stable');   % returns indices in E
    [~,~,IF_p] = intersect(Flvarname,pvarnames_all,'stable');       % returns indices in F
    np_all=length(pvarnames_all);
    
    % adjust the degmats accordingly
    Edegmat_new=spalloc(ntE,np_all,nnz(Erdegmat));
    Fdegmat_new=spalloc(ntF,np_all,nnz(Frdegmat));
    Edegmat_new(:,IE_p)=Erdegmat;
    Fdegmat_new(:,IF_p)=Frdegmat;
    
    % adjust E and F to use the new varnames and degmats
    F = polynomial(F.coef,Fdegmat_new,pvarnames_all,F.matdim);
    E.degmat=Edegmat_new;
    E.varname = vec(pvarnames_all); % NOTE: if E and F have only 1 varname, need vec to get appropriate dim
    
    % Distinguish three cases:
    if all(size(E)==[1,1]) % scalar dpvar multiplying a pvar
        C = E;
        C.matdim = size(F);
        
        % multiply coefficient matrices
        fcoefs = spalloc(size(F,1),size(F.coefficient,1)*size(F,2),nnz(F.coefficient));
        for i = 1:size(F,2) % maybe implement fancy reshape from pvar mtimes??
            tmp = F.coefficient(:,[(i-1)*size(F,1)+1:i*size(F,1)]);
            fcoefs(:,(i-1)*size(F.coefficient,1)+1:i*size(F.coefficient,1)) = tmp';
        end
%     Alternative implementations
%         fcoefs = [];
%         for i = 1:size(F,2) % maybe implement fancy reshape from pvar mtimes??
%             tmp = F.coefficient(:,[(i-1)*size(F,1)+1:i*size(F,1)]);
%             fcoefs = [fcoefs, tmp'];
%         end
%     Alternative implementation 2
%         fcoefs = Bshape(F); % this is faster only for large dimensions
        C.C = kron(fcoefs,E.C); % Construction of Cnew
        
        % multiply degmats
        %   this assumes both E and F have same pvars
        tmpdegmat = rowwisesum(C.degmat,F.degmat);
        
        % combine the product degmat
        [Adeg, degmat] = combine_degmat(tmpdegmat);
        
        % adjust coefficient matrix to combine rows with same monomials
        C.C = C.C*kron(eye(size(F,2)),Adeg);
        C.degmat = degmat;
        
    elseif all(size(F)==[1,1]) % matrix dpvar times scalar polynomial
        C = E;
        C.matdim = size(E);
        
        % multiply coefficient matrices
        fcoef = F.coefficient;
        ecoef = E.C;
        C.C = kron(ecoef,fcoef'); % Construction of Cnew
        
        % multiply degmats
        %   this assumes both E and F have same pvars
        tmpdegmat = rowwisesum(F.degmat,C.degmat);
        
        % combine the product degmat
        [Adeg, degmat] = combine_degmat(tmpdegmat);
        
        % adjust coefficient matrix to combine rows with same monomials
        C.C = C.C*kron(eye(size(E,2)),Adeg);
        C.degmat = degmat;
        
    else % dpvar matrix (E) times pvar matrix (F)
        if size(E,2)~=size(F,1)
            error('Inner dimensions dont match');
        end
        C = E;
        C.matdim = [E.matdim(1),size(F,2)];
        
        % multiply coefficient matrices
        fcoefs = Bshape(F); % this is faster only for large dimensions
        Cnew = C.C*kron(fcoefs,speye(size(E.degmat,1)));

        % multiply degmats
        %   this assumes both E and F have same pvars
        tmpdegmat = rowwisesum(C.degmat,F.degmat);
        
        % combine the product degmat
        [Adeg, degmat] = combine_degmat(tmpdegmat);
        
        % adjust coefficient matrix to combine rows with same monomials
        C.C = Cnew*kron(eye(size(F,2)),Adeg);
        C.degmat = degmat;
    end
end
end



function [A,d] = combine_degmat(dmat)
% Combine non-unqiue rows of dmat to build d, with A such that A*d=dmat
[d, ~, idxd] = unique(dmat,'rows');
A = sparse(1:length(idxd),idxd,ones(length(idxd),1),size(dmat,1), size(d,1));
% alternative way
% A = zeros(size(dmat,1), size(d,1));
% 
% for j=1:length(idxd)
%     A(j,idxd(j)) = 1;
% end
end
function [ C ] = rowwisesum( A, B ) 
% return the elementwise sum of two matrices, adds matrix A to every
% row of B 
C = kron(B,ones(size(A,1),1)) + repmat(A,size(B,1),1);
end