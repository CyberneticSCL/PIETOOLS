function sos = sosconstr_mp(sos,Type,symexpr)
% SOSCONSTR --- Add a new constraint (equality/inequality) 
%    to an SOS program 
%
% SOSP = sosconstr(SOSP,TYPE,EXPR)
%
% SOSP is the sum of squares program.
% The new constraint is described by TYPE as follows:
%   TYPE = 'eq'   : constraint is an equality, viz., f(x) = 0
%   TYPE = 'ineq' : constraint is an inequality, viz., f(x) >= 0 (an SOS)
% EXPR is the expression in the left hand side of the constraint, i.e., f(x).
%

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 1.00.
%
% Copyright (C) 2002  S. Prajna (1), A. Papachristodoulou (1), P. A. Parrilo (2)
% (1) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (2) Institut fur Automatik - ETH Zurich, CH-8092 Zurich, Switzerland.
%
% Send feedback to: sostools@cds.caltech.edu
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


% 12/24/01 - SP
% 03/01/02 - SP -- New syntax
% 04/06/03 - PJS -- Handle poly objects w/out for-loops
% 07/12/13 - MP -- Allow Matrix-valued equality constraints

if isempty(find(symexpr.coefficient)) % If the xpression is simply zero, return
    return
end

sos.expr.num = sos.expr.num+1;
expr = sos.expr.num;
sos.expr.type{expr} = Type;


if isfield(sos,'symvartable')

  [sos.expr.At{expr},sos.expr.b{expr},sos.expr.Z{expr}] = ...
      getequation(char(symexpr),sos.vartable,sos.decvartable);
    
else

  % Pull out information from the polynomial
  coef = symexpr.coefficient;
  deg = symexpr.degmat;
  var = symexpr.varname;
  %matdim=size(symexpr); % dimensions of the matrix MMP
  
  % Sort variables: decision stacked on non-decision
  [dummy,idxdvar1,idxdvar2] = intersect(var,sos.decvartable);
  ldv = length(idxdvar1);
  [dummy,idxvar1,idxvar2] = intersect(var,sos.vartable);
  var = var([idxdvar1; idxvar1]);
  deg = deg(:,[idxdvar1; idxvar1]);
  
  % Sort terms: Stack terms without decision variables above 
  % terms with decision variables.
  [deg,sortidx] = sortrows(deg);
  coef = coef(sortidx,:);
  if ldv~=0
    [idx,dummy1]=find( deg(:,1:ldv) );
    nvterm = min(idx)-1;
  else
    nvterm = size(deg,1); % last term without a decision variable
  end
  D12 = deg(1:nvterm,ldv+1:end); % terms without decision variables. leaving off the zeros portion (dvar part)
  D22 = deg(nvterm+1:end,ldv+1:end); % terms with dvars, but the non-dvar part
  D21 = deg(nvterm+1:end,1:ldv);   % terms with dvars, the dvar part
  
  % Get the equality constraints    
  if any( sum(D21,2) > 1 )
    error(['The expression is not linear in the decision' ...
	   ' variables']);
  else
    % Extract monomials
    [Ztemp,idx1,idx2] = unique([D12; D22],'rows'); % idx2 maps non-dvar monomials in the non-dvar part to new monomials, then also maps monomials in the dvar part to new monomials.
    lZ = size(Ztemp,1);                            % number of non dvar monomials
    sos.expr.Z{expr} = sparse(lZ,length(sos.vartable));   
    sos.expr.Z{expr}(:,idxvar2) = Ztemp;           % assigns set of monomials in this expression - will yeild number of equality constraints
    
    % Form the equality constraints
    % Ac^T xdec = b
    ndvars=length(sos.decvartable);
    ncoef = size(coef,2); %MMP
    Ac = sparse(ndvars,lZ*ncoef);       % NOTE: Ac^T 1 row for each non-dec monomials, 1 column for each decvar. If the coefficients are going to be matrix-valued, however, we need lots more rows. 1 for each col of coeff.
    b= sparse(lZ*ncoef,1);
    
    % If symexpr has not been combined down, need to check that
    % idx2(1:nvterm) are all unique....
    %  MMP - but this hasn't been done!

if nvterm~=0
        [ipositC1 jpositC1 C1coef]=find(coef(1:nvterm,:));  % iposit is term number, jposit is matrix element number
        if size(jpositC1,2)>1     %Fixes Matlab stupidity
            jpositC1=jpositC1';ipositC1=ipositC1';C1coef=C1coef';
        end
        positb=(idx2(ipositC1)-1)*ncoef+jpositC1;  %idx indicates the non-decision variable monomial(position in Ztemp) in the term i. 
        b(positb)=C1coef; %b(idx2(1:nvterm)) = coef(1:nvterm);
end        

    if ldv~=0
        [ipositC jpositC vcoef]=find(coef(nvterm+1:end,:));  % iposit is term number, jposit is matrix element numberc
        if size(jpositC,2)>1     %Fixes Matlab stupidity
            jpositC=jpositC';ipositC=ipositC';vcoef=vcoef';
        end
        %to where do the coefficients get mapped?
        [dv,dummy]=find(D21');      %gives the variable number associated to each term in degmat. this gives the i position data in Ac
        if size(dv,2)>1     %Fixes Matlab stupidity
            dv=dv';
        end
        ipositAc=idxdvar2(dv(ipositC));        %idxdvar2 indicates the decision variable(position in varname) associated with term i 
        jpositAc=(idx2(nvterm+ipositC)-1)*ncoef+jpositC;  %idx indicates the non-decision variable monomial(position in Ztemp) in the term i. 
        tempidx = sub2ind(size(Ac) , ipositAc  , jpositAc );
        Ac(tempidx) = vcoef;
        

      % for each entry here, I want to map kkk->ncoef*(kkk-1)+1:ncoef*kkk
        

      % For each element in var2, we create a vector of entries in Ac
      % position dv
    end

    % We should probably eliminate all repeated or trivial constraints
    [tempAb]=unique([Ac' b],'rows');
    Ac=tempAb(:,1:(end-1))';
    b=tempAb(:,end);
    sos.expr.At{expr} = -Ac; 
    sos.expr.b{expr} = b;
  end

% $$$   % Extract degmat and sort the terms according to degmat     
% $$$   [dummy,idxvar1,idxvar2] = intersect(symexpr.varname,sos.vartable);
% $$$   [Ztemp,idx1,idx2] = unique(sparse(symexpr.degmat(:,idxvar1)),'rows');
% $$$   sos.expr.Z{expr} = sparse(size(Ztemp,1),length(sos.vartable));
% $$$   sos.expr.Z{expr}(:,idxvar2) = Ztemp;
% $$$   
% $$$   % Now get the coefficients
% $$$   Ac = sparse(length(sos.decvartable),size(sos.expr.Z{expr},1));
% $$$   b = sparse(size(sos.expr.Z{expr},1),1);
% $$$   
% $$$   [dummy,idxdecvar1,idxdecvar2] = ... 
% $$$       intersect(symexpr.varname,sos.decvartable);
% $$$   for i = 1:size(Ztemp,1)
% $$$     
% $$$     for j = find(idx2==i)'         % These are the terms corresponding to the i-th row of Ztemp
% $$$       k = find(symexpr.degmat(j,idxdecvar1)==1);
% $$$       if isempty(k)
% $$$ 	b(i) = symexpr.coefficient(j);
% $$$       elseif length(k)>1
% $$$ 	error('The expression is not linear in the decision variables');
% $$$       else
% $$$ 	Ac(idxdecvar2(k),i) = symexpr.coefficient(j);
% $$$       end;
% $$$     end;
% $$$   end;
% $$$   
% $$$   sos.expr.At{expr} = -Ac; 
% $$$   sos.expr.b{expr} = b;
% $$$ end  

end



