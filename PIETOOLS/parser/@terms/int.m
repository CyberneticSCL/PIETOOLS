function intTerms = int(objA, var, lim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms-class method that integrates objA
% Input: 
% objA - terms class objects
% var - polynomial variable wrt which integration is performed
% lim - limits of integration, 1x2 pvar or double
% Output:
% intTerms - terms class object \int_{lim{1}}^{lim{2}} objA d(var)
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
if any(isequal(lim,var))
    error('Limits of integration and integration variable are the same. Change the variable naming');
end
if ~poly2double(lim(2))&&~poly2double(lim(1))
    error('Integral should have at least one limit of integration as 0 or 1');
end
if poly2double(lim(2))&&(double(lim(2))~=1)
    error('Upper limit of integration can only be a "pvar" variable or 1');
end
if poly2double(lim(1))&&(double(lim(1))~=0)
    error('Lower limit of integration can only be a "pvar" variable or 0');
end

try
    intterms = double(objA.operator.R.R2)|double(objA.operator.R.R1); %R1 and R2 must be zero
    if any(intterms(:))
        error('Double integrals are currently not supported')
    end
catch
    error('Double integrals are currently not supported')
end

% isdot_A = []; isout_A=[]; 
tempstatevec = objA.statevec;
isdot_A = isdot(tempstatevec); isout_A=isout(tempstatevec); 
% for i=1:length(tempstatevec)
%     isdot_A = [isdot_A; tempstatevec(i).diff_order(1)*ones(tempstatevec(i).veclength,1)];
%     isout_A = [isout_A; strcmp(tempstatevec(i).type,'out')*ones(tempstatevec(i).veclength,1)];
% end
if any((isdot_A|isout_A))
    error("Integration of vectors with outputs or time-derivative of state is not allowed");
end


opvar T; T.var2 = var;
if poly2double(lim(2))
T.R.R2 = objA.operator.R.R0;
end
if poly2double(lim(1))
T.R.R1 = objA.operator.R.R0;    
end

intTerms = terms(T,tempstatevec);
end