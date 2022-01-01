function [d] = degbalance(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [d] = degbalance(P,opts) provides degrees to be selected for a poslpivar Q such that
% components of P and Q have similar degrees
% 
% INPUT
%   P: PI opvar class object
%   opts: options describing the type of poslpivar (pure, diag)
% 
% OUTPUT 
%   d: vector of degrees that is the input to sos_posopvar
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - degbalance
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 7_26_2019

T = getdeg(P);
Q_max_s = T(2,3);
R0_max_s = T(4,3);
R1_max_s = T(5,3);
R1_max_t = T(5,4);

if nargin ==1
    dup = 0;
    d1 = ceil(R0_max_s/2);
    d2t = ceil(0.5*R1_max_t);
    d2s = ceil(0.5*max([R1_max_s-d1, R1_max_s-Q_max_s,0]));
    d3t = ceil(0.5*R1_max_s);
    d3s = ceil(0.5*max([R1_max_t-d1, R1_max_t-Q_max_s,0]));
    d2_total = max([Q_max_s-1,d2s,d2t]); 
    d3_total = max([Q_max_s-1,d3s,d3t]);
    d = {d1,[d2s d2t d2_total],[d3s d3t d3_total]};
elseif nargin==2
    if opts.pure==1&& opts.full==0
        d1 = 0;
        d2t = R1_max_t;
        d2s = max(ceil((R1_max_s-1-d2t)/2),Q_max_s-1-d2t);
        d = [d1, max([d2t,0]), max([d2s,0])];
    elseif opts.pure==1
        d1 = 0;
        d2tot = Q_max_s -1;
        d2s = R1_max_s -1 -d2tot;
        d2t = d2tot-d2s;
        d = [d1, max([d2t,0]), max([d2s,0])];
    elseif opts.full==0
        if R1_max_t>R1_max_s
            d1 = R0_max_s/2;
            d3t = max(R1_max_s, Q_max_s-1);
            d3s = (R1_max_t-R1_max_s)/2-0.5;
            d3 = min(R1_max_t,(R1_max_t+R1_max_s-1)/2);
            d = {d1, [0 0 0], [d3s d3t d3]};
        else
            d1 = R0_max_s/2;
            d2t = max(R1_max_t, Q_max_s-1);
            d2s = (R1_max_s-R1_max_t)/2-0.5;
            d2 = min(R1_max_s,(R1_max_t+R1_max_s-1)/2);
            d = {d1, [d2s d2t d2],[0 0 0]};
        end
    else
        d1 = R0_max_s/2;
        d2t = max(R1_max_t,Q_max_s-1);
        d2s = 0;
        d3t = max(R1_max_s,Q_max_s-1);
        d3s = 0;
        d = {d1,[d2s d2t],[d3s d3t]};
    end
end
end