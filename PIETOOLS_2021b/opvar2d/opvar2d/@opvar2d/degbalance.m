function [d] = degbalance(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [d] = degbalance(P,opts) provides degrees to be selected for a poslpivar Q such that
% components of P and Q have similar degrees
% 
% INPUT
%   P: PI opvar2d class object
%   opts: options describing the type of poslpivar (pure, diag)
% 
% OUTPUT 
%   d: vector of degrees that is the input to sos_posopvar
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - degbalance
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 02_27_2021
%   ^ Based heavily on "@opvar"-degbalance code by SS ^

T = getdeg(P);

R0x_max_s1 = T(2,5);
Rxxo_max_s1 = T(6,5);
Rxxa_max_s1 = T(7,5);
Rxxa_max_t1 = T(7,6);
Rxy_max_s1 = T(9,5);


R0y_max_s2 = T(3,7);
Ryyo_max_s2 = T(15,7);
Ryya_max_s2 = T(15,7);
Ryya_max_t2 = T(16,8);
Rxy_max_s2 = T(9,7);


R02_max_s1 = T(4,5);
R02_max_s2 = T(4,7);

Rx2o_max_s1 = T(10,5);      Rx2o_max_t1 = T(10,6);      Rx2o_max_s2 = T(10,7);
Rx2a_max_s1 = T(11,5);      Rx2a_max_t1 = T(11,6);      Rx2a_max_s2 = T(11,7);
Rx2b_max_s1 = T(12,5);      Rx2b_max_t1 = T(12,6);      Rx2b_max_s2 = T(12,7);

Ry2o_max_s1 = T(18,5);      Ry2o_max_s2 = T(18,7);      Ry2o_max_t2 = T(18,8);
Ry2a_max_s1 = T(19,5);      Ry2a_max_s2 = T(19,7);      Ry2a_max_t2 = T(19,8);
Ry2b_max_s1 = T(20,5);      Ry2b_max_s2 = T(20,7);      Ry2b_max_t2 = T(20,8);

R22oo_max_s1 = T(28,5);     R22oo_max_s2 = T(28,7);
R22ao_max_s1 = T(29,5);     R22ao_max_t1 = T(29,6);     R22ao_max_s2 = T(29,7);
R22bo_max_s1 = T(30,5);     R22bo_max_t1 = T(30,6);     R22bo_max_s2 = T(30,7);
R22oa_max_s1 = T(31,5);     R22oa_max_s2 = T(31,7);     R22oa_max_t2 = T(31,8);
R22ob_max_s1 = T(32,5);     R22ob_max_s2 = T(32,7);     R22ob_max_t2 = T(32,8);
R22aa_max_s1 = T(33,5);     R22aa_max_t1 = T(33,6);     R22aa_max_s2 = T(33,7);     R22aa_max_t2 = T(33,8);
R22ba_max_s1 = T(34,5);     R22ba_max_t1 = T(34,6);     R22ba_max_s2 = T(34,7);     R22ba_max_t2 = T(34,8);
R22ab_max_s1 = T(35,5);     R22ab_max_t1 = T(35,6);     R22ab_max_s2 = T(35,7);     R22ab_max_t2 = T(35,8);
R22bb_max_s1 = T(36,5);     R22bb_max_t1 = T(36,6);     R22bb_max_s2 = T(36,7);     R22bb_max_t2 = T(36,8);



if nargin ==1
    dup = 0;
    
    dx1 = ceil(Rxxo_max_s1/2);
    dx2t = ceil(0.5*Rxxa_max_t1);
    dx2s = ceil(0.5*max([Rxxa_max_s1-dx1, Rxxa_max_s1-R0x_max_s1, Rxxa_max_s1-Rxy_max_s1,0]));
    dx3t = ceil(0.5*Rxxa_max_s1);
    dx3s = ceil(0.5*max([Rxxa_max_t1-dx1, Rxxa_max_t1-R0x_max_s1, Rxxa_max_t1-Rxy_max_s1,0]));
    dx2_tot = max([R0x_max_s1-1,Rxy_max_s1-1,dx2s,dx2t]); 
    dx3_tot = max([R0x_max_s1-1,Rxy_max_s1-1,dx3s,dx3t]);
    dx = {dx1,[dx2s dx2t dx2_tot],[dx3s dx3t dx3_tot]};
    
    
    dy1 = ceil(Ryyo_max_s2/2);
    dy2t = ceil(0.5*Ryya_max_t2);
    dy2s = ceil(0.5*max([Ryya_max_s2-dy1, Ryya_max_s2-R0y_max_s2, Ryya_max_s2-Rxy_max_s2,0]));
    dy3t = ceil(0.5*Ryya_max_s2);
    dy3s = ceil(0.5*max([Ryya_max_t2-dy1, Ryya_max_t2-R0y_max_s2, Ryya_max_t2-Rxy_max_s2,0]));
    dy2_tot = max([R0y_max_s2-1,Rxy_max_s2-1,dy2s,dy2t]); 
    dy3_tot = max([R0y_max_s2-1,Rxy_max_s2-1,dy3s,dy3t]);
    dy = {dy1,[dy2s dy2t dy2_tot],[dy3s dy3t dy3_tot]};
    
    
    d2oos1 = ceil(R22oo_max_s1/2);
    d2oos2 = ceil(R22oo_max_s2/2);
    d2oos1s2 = max([Rx2o_max_s1-1,Ry2o_max_s2-1,d2oos1,d2oos2]); 
    
    d2aot1 = ceil(0.5*max([Rx2o_max_t1 R22ao_max_t1]));
    d2aos1 = ceil(0.5*max([R22ao_max_s1-d2oos1, R22ao_max_s1-R02_max_s1, R22ao_max_s1-Ry2o_max_s1,0]));
    d2aos2 = ceil(0.5*max([R22ao_max_s2-d2oos2, R22ao_max_s2-Rx2o_max_s2,0]));
    d2bot1 = ceil(0.5*max([Rx2o_max_s1 R22ao_max_s1]));
    d2bos1 = ceil(0.5*max([R22ao_max_t1-d2oos1, R22ao_max_t1-R02_max_s1, R22ao_max_t1-Ry2o_max_s1,0]));
    d2bos2 = ceil(0.5*max([R22bo_max_s2-d2oos2, R22bo_max_s2-Rx2o_max_s2,0]));
    d2aos1t1 = max([R02_max_s1-1, Ry2o_max_s1-1, Ry2a_max_s1-1, Ry2b_max_s1-1, d2aos1, d2aot1]);
    d2bos1t1 = max([R02_max_s1-1, Ry2o_max_s1-1, Ry2a_max_s1-1, Ry2b_max_s1-1, d2bos1, d2bot1]);
    d2ao_tot = max([d2aos1, d2aos2, d2aos1t1, d2oos1s2]);
    d2bo_tot = max([d2bos1, d2bos2, d2bos1t1, d2oos1s2]);
    
    d2oat2 = ceil(0.5*max([Ry2o_max_t2 R22oa_max_t2]));
    d2oas2 = ceil(0.5*max([R22oa_max_s2-d2oos2, R22oa_max_s2-R02_max_s2, R22oa_max_s2-Rx2o_max_s2,0]));
    d2oas1 = ceil(0.5*max([R22oa_max_s1-d2oos1, R22oa_max_s1-Ry2o_max_s1,0]));
    d2obt2 = ceil(0.5*max([Ry2o_max_s2 R22oa_max_s2]));
    d2obs2 = ceil(0.5*max([R22oa_max_t2-d2oos2, R22oa_max_t2-R02_max_s2, R22oa_max_t2-Rx2o_max_s2,0]));
    d2obs1 = ceil(0.5*max([R22bo_max_s1-d2oos1, R22ob_max_s1-Ry2o_max_s1,0]));
    d2oas2t2 = max([R02_max_s2-1, Rx2o_max_s2-1, Rx2a_max_s2-1, Rx2b_max_s2-1, d2oas2, d2oat2]);
    d2obs2t2 = max([R02_max_s2-1, Rx2o_max_s2-1, Rx2a_max_s2-1, Rx2b_max_s2-1, d2obs2, d2obt2]);
    d2oa_tot = max([d2oas2, d2oas1, d2oas2t2, d2oos1s2]);
    d2ob_tot = max([d2obs2, d2obs1, d2obs2t2, d2oos1s2]);
    
    d2aat1 = ceil(0.5*max([R22aa_max_t1 R22ao_max_t1]));
    d2aas1 = ceil(0.5*max([R22aa_max_s1-d2oos1, R22aa_max_s1-R02_max_s1,0]));
    d2aas1t1 = max([R02_max_s1-1, Ry2o_max_s1-1, Ry2a_max_s1-1, Ry2b_max_s1-1, d2aas1, d2aat1]);
    d2aat2 = ceil(0.5*max([R22aa_max_t2 R22oa_max_t2]));
    d2aas2 = ceil(0.5*max([R22aa_max_s2-d2oos2, R22aa_max_s2-R02_max_s2,0]));
    d2aas2t2 = max([R02_max_s2-1, Rx2o_max_s2-1, Rx2a_max_s2-1, Rx2b_max_s2-1, d2aas2, d2aat2]);
    d2aat1t2 = max([d2aat1, d2aat2]);
    d2aas1s2 = max([d2aas1, d2aas2]);
    d2aa_tot = max([d2aat1t2, d2aas1s2, d2aas1t1, d2aas2t2]);
    
    d2bat1 = ceil(0.5*max([R22ba_max_s1 R22bo_max_s1]));
    d2bas1 = ceil(0.5*max([R22ba_max_t1-d2oos1, R22ba_max_t1-R02_max_s1,0]));
    d2bas1t1 = max([R02_max_s1-1, Ry2o_max_s1-1, Ry2a_max_s1-1, Ry2b_max_s1-1, d2bas1, d2bat1]);
    d2bat2 = ceil(0.5*max([R22ba_max_t2 R22oa_max_t2]));
    d2bas2 = ceil(0.5*max([R22ba_max_s2-d2oos2, R22ba_max_s2-R02_max_s2,0]));
    d2bas2t2 = max([R02_max_s2-1, Rx2o_max_s2-1, Rx2a_max_s2-1, Rx2b_max_s2-1, d2bas2, d2bat2]);
    d2bat1t2 = max([d2bat1, d2bat2]);
    d2bas1s2 = max([d2bas1, d2bas2]);
    d2ba_tot = max([d2bat1t2, d2bas1s2, d2bas1t1, d2bas2t2]);
    
    d2abt1 = ceil(0.5*max([R22ab_max_t1 R22ao_max_t1]));
    d2abs1 = ceil(0.5*max([R22ab_max_s1-d2oos1, R22ab_max_s1-R02_max_s1,0]));
    d2abs1t1 = max([R02_max_s1-1, Ry2o_max_s1-1, Ry2a_max_s1-1, Ry2b_max_s1-1, d2abs1, d2abt1]);
    d2abt2 = ceil(0.5*max([R22ab_max_s2 R22ob_max_s2]));
    d2abs2 = ceil(0.5*max([R22ab_max_t2-d2oos2, R22ab_max_t2-R02_max_s2,0]));
    d2abs2t2 = max([R02_max_s2-1, Rx2o_max_s2-1, Rx2a_max_s2-1, Rx2b_max_s2-1, d2abs2, d2abt2]);
    d2abt1t2 = max([d2abt1, d2abt2]);
    d2abs1s2 = max([d2abs1, d2abs2]);
    d2ab_tot = max([d2abt1t2, d2abs1s2, d2abs1t1, d2abs2t2]);
    
    d2bbt1 = ceil(0.5*max([R22bb_max_s1 R22bo_max_s1]));
    d2bbs1 = ceil(0.5*max([R22bb_max_t1-d2oos1, R22bb_max_t1-R02_max_s1,0]));
    d2bbs1t1 = max([R02_max_s1-1, Ry2o_max_s1-1, Ry2a_max_s1-1, Ry2b_max_s1-1, d2bbs1, d2bbt1]);
    d2bbt2 = ceil(0.5*max([R22bb_max_s2 R22ob_max_s2]));
    d2bbs2 = ceil(0.5*max([R22bb_max_t2-d2oos2, R22bb_max_t2-R02_max_s2, 0]));
    d2bbs2t2 = max([R02_max_s2-1, Rx2o_max_s2-1, Rx2a_max_s2-1, Rx2b_max_s2-1, d2bbs2, d2bbt2]);
    d2bbt1t2 = max([d2bbt1, d2bbt2]);
    d2bbs1s2 = max([d2bbs1, d2bbs2]);
    d2bb_tot = max([d2bbt1t2, d2bbs1s2, d2bbs1t1, d2bbs2t2]);
    
    d2 = {[d2oos1;d2oos2;d2oos1s2],[d2oas1;d2oas2;d2oat2;d2oas2t2;d2oa_tot],[d2obs1;d2obs2;d2obt2;d2obs2t2;d2ob_tot];
        [d2aos1;d2aot1;d2aos1t1;d2aos2;d2ao_tot],[d2aas1,d2aas2,d2aas1s2;d2aat1,d2aat2,d2aat1t2;d2aas1t1,d2aas2t2,d2aa_tot],[d2abs1,d2abs2,d2abs1s2;d2abt1,d2abt2,d2abt1t2;d2abs1t1,d2abs2t2,d2ab_tot];
        [d2bos1;d2bot1;d2bos1t1;d2bos2;d2bo_tot],[d2bas1,d2bas2,d2bas1s2;d2bat1,d2bat2,d2bat1t2;d2bas1t1,d2bas2t2,d2ba_tot],[d2bbs1,d2bbs2,d2bbs1s2;d2bbt1,d2bbt2,d2bbt1t2;d2bbs1t1,d2bbs2t2,d2bb_tot]};
    
    d.dx = dx;
    d.dy = dy;
    d.d2 = d2;
    
elseif nargin==2
    if opts.pure==1&& opts.full==0
        dx1 = 0;
        dx2t = Rxxa_max_t1;
        dx2s = max(ceil((Rxxa_max_s1-1-dx2t)/2),R0x_max_s1-1-dx2t);
        dx = [dx1, max([dx2t,0]), max([dx2s,0])];
    elseif opts.pure==1
        dx1 = 0;
        d2tot = R0x_max_s1 -1;
        dx2s = Rxxa_max_s1 -1 -d2tot;
        dx2t = d2tot-dx2s;
        dx = [dx1, max([dx2t,0]), max([dx2s,0])];
    elseif opts.full==0
        if Rxxa_max_t1>Rxxa_max_s1
            dx1 = Rxxo_max_s1/2;
            dx3t = max(Rxxa_max_s1, R0x_max_s1-1);
            dx3s = (Rxxa_max_t1-Rxxa_max_s1)/2-0.5;
            d3 = min(Rxxa_max_t1,(Rxxa_max_t1+Rxxa_max_s1-1)/2);
            dx = {dx1, [0 0 0], [dx3s dx3t d3]};
        else
            dx1 = Rxxo_max_s1/2;
            dx2t = max(Rxxa_max_t1, R0x_max_s1-1);
            dx2s = (Rxxa_max_s1-Rxxa_max_t1)/2-0.5;
            d2 = min(Rxxa_max_s1,(Rxxa_max_t1+Rxxa_max_s1-1)/2);
            dx = {dx1, [dx2s dx2t d2],[0 0 0]};
        end
    else
        dx1 = Rxxo_max_s1/2;
        dx2t = max(Rxxa_max_t1,R0x_max_s1-1);
        dx2s = 0;
        dx3t = max(Rxxa_max_s1,R0x_max_s1-1);
        dx3s = 0;
        dx = {dx1,[dx2s dx2t],[dx3s dx3t]};
    end
end
end