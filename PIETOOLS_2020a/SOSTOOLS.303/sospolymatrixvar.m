function [sos,P] = sospolymatrixvar(sos,ZSym,n,matrixstr)
% SOSMATRIXVAR --- Declare a polynomial matrix variable P in the sos program
% SOS of
%
% [SOSP,P] = sosmatrixvar(SOSP,ZSym,n,matrixstr)
%
% SOSP is the sum of squares program.
% P is the new polynomial matrix.
% n is the dimension of the matrix P: n(1) x n(2)
% ZSym is the vector of monomials contained in VAR. Decision
% variables corresponding to those monomials will be assigned
% automatically by SOSPOLYVAR.
% matrixstr is a char string with the option 'symmetric' when required

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 3.03.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018  A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5)
% (1) Department of Engineering Science, University of Oxford, Oxford, U.K.
% (2) Laboratoire de Signaux et Systmes, CentraleSupelec, Gif sur Yvette,
%     91192, France
% (3) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (4) Aerospace and Engineering Mechanics Department, University of
%     Minnesota, Minneapolis, MN 55455-0153, USA.
% (5) Laboratory for Information and Decision Systems, M.I.T.,
%     Massachusetts, MA 02139-4307
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


% JA&GV - 6/13/2013



% Original Code
if nargin == 4
    if matrixstr=='symmetric'
        if n(1)==n(2)
            if isfield(sos,'symvartable')
                P = sym(zeros(n(1),n(2)));
            else
                % Code for multipoly: PJS 9/9/2013
                P = polynomial(zeros(n(1),n(2)));
            end
            for i = 1:n(1)
                for j = i:n(1)
                    [sos,var] = sospolyvar(sos,ZSym);
                    P(i,j) = var;
                    P(j,i) = var;
                    clear var;
                end
            end
        else
            disp(['''symmetric''' ' option used, matrix must be square.']);
            P = [];
            return
        end
    else
        disp(['Matrix structure ' matrixstr ' is not defined.' ]);
        P = [];
        return
        
    end
else
    if isfield(sos,'symvartable')
        P = sym(zeros(n(1),n(2)));
    else
        % Code for multipoly: PJS 9/9/2013
        P = polynomial(zeros(n(1),n(2)));
    end
    for i = 1:n(1)
        for j = 1:n(2)
            [sos,var] = sospolyvar(sos,ZSym);
            P(i,j) = var;
            clear var;
        end
    end
    
end