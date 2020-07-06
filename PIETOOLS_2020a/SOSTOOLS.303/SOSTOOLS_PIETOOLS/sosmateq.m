function sos = sosmateq(sos,symexpr)
% SOSEQ --- Add a new equality constraint f(x) = 0
%    to an SOS program 
%
% SOSP = soseq(SOSP,EXPR)
%
% SOSP is the sum of squares program.
% EXPR is the expression on the left hand side of the constraint, i.e., f(x).
%
% EXPR can be a column vector. In this case, several equality
% constraints will be added simultaneously to the sum of
% squares program.
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


    sos = sosconstr_mp(sos,'eq',symexpr);
