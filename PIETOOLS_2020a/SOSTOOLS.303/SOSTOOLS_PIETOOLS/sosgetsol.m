function p = sosgetsol(sos,V,digit)
% SOSGETSOL --- Get the solution from a solved SOS program 
%
% SOL = sosgetsol(SOSP,VAR,DIGIT) 
%
% SOL is the solution from a solved sum of squares program SOSP,
% obtained through substituting all the decision variables
% in VAR by the numerical values which are the solutions to
% the corresponding semidefinite program. 
%
% The third argument DIGIT (optional) will determine the 
% accuracy of SOL in terms of the number of digits. Default 
% value is 5.
%

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 2.00.
%
% Copyright (C)2002, 2004  S. Prajna (1), A. Papachristodoulou (1), P. Seiler (2),
%                          P. A. Parrilo (3)
% (1) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (2) Mechanical and Industrial Engineering Department - University of Illinois 
%     Urbana, IL 61801, USA
% (3) Institut fur Automatik - ETH Zurich, CH-8092 Zurich, Switzerland.
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
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%

% 12/27/01 - SP
% 02/21/02 - SP -- Symbolic polynomial
% 03/01/02 - SP -- New syntax
% 03/15/02 - SP -- Fast
% 06/09/13 - MP -- Matrix-valued polynomials
% 06/25/20 - Sachin -- fixed bug related constant polynomials

if nargin == 2
    digit = 5;   % Default
end;

if isfield(sos,'symvartable')
    
    p = mysymsubs(V,sos.symdecvartable,sos.solinfo.RRx(1:length(sos.symdecvartable)),digit);

else
    
    [dummy,idxdecvar1,idxdecvar2] = intersect(V.varname,sos.decvartable);
    idxvar = setdiff(1:length(V.varname),idxdecvar1);
    coeffs = [V.degmat(:,idxdecvar1) 1-sum(V.degmat(:,idxdecvar1),2)]*[sos.solinfo.RRx(idxdecvar2);1];
                                                %MMP 6/9/2013 updated to
                                                %allow for matrix-valued
                                                %polynomials and to fix
                                                %problem with terms with no
                                                %decision variables
    coeffs = V.coefficient.*repmat(coeffs,1,size(V.coefficient,2));         % 01/31/02
    varname = V.varname(idxvar);
    if isempty(idxvar)
        degmat = [];
    else
        degmat = V.degmat(:,idxvar);
    end;
 %   p = set(V,'varname',varname,'degmat',degmat,'coefficient',coeffs);
%     degmat 
%     coeffs
%     size(V)
    if isempty(degmat)
      coeffs=sum(coeffs,1); % modified by sachin - 6/25/2020 original version sum(coeffs)
      p=polynomial(coeffs);
      p=reshape(p,size(V));
  else
    p=polynomial(coeffs,degmat,varname,size(V));
   end
    p=combine(p);
end;
