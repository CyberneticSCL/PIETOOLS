function [Xop,eps,deg_fctr_final] = mrdivide(P2op,P1op,deg_fctr,tol,deg_fctr_max)
% [Xop,eps,deg_fctr_final] = mrdivide(P2op,P1op,deg_fctr,tol,deg_fctr_max)
% computes an opvar2d object Xop=P2op/P1op, such that Xop*P1op = P2op.
%
% INPUTS:
% - P2op:         [m0,n0; mx,nx; my,ny; m2,n2] opvar2d object.     
% - P1op:         [p0,n0; px,nx; py,ny; p2,n2] opvar2d object.
% - deg_fctr:     4x2 'double' array of integers, specifying a factor with
%                 which to increase the degrees of the monomials appearing
%                 in P1op to set the degrees of the monomials appearing in
%                 Xop. Here, we represent P1op and P2op using monomials
%                   Zx(x),   Ztt(tt);
%                   Zy(y),   Znu(nu);
%                   Zxy(xy), Zttnu(tt,nu)
%                 and we represent Xop using monomials
%                   iZx(x),   iZtt(tt);
%                   iZy(y),   iZnu(nu);
%                   iZxy(xy), iZttnu(tt,nu)
%                 Then the degrees of
%                   iZx in x,   iZtt in tt;
%                   iZy in y,   iZnu in nu;
%                   iZxy in x,  iZttnu in tt;
%                   iZxy in y,  iZttnu in nu
%                 are given by the maximal degrees of
%                   Zx in x and Ztt in tt;
%                   Zy in y and Znu in nu;
%                   Zxy in x and Zttnu in tt;
%                   Zxy in y and Zttnu in nu
%                 scaled with factors
%                   deg_fctr(1,1),  deg_fctr(1,2);
%                   deg_fctr(2,1),  deg_fctr(2,2);
%                   deg_fctr(3,1),  deg_fctr(3,2);
%                   deg_fctr(4,1),  deg_fctr(4,2);
% - tol:          positive scalar specifying a desired accuracy of the
%                 inverse. An estimate of the accuracy is determined as a
%                 maximal absolute value of the parameters defining the
%                 difference Xop*P1op-P2op. If the current guess of Xop
%                 is not sufficiently accurate, the deg_fctr will be
%                 increased, and a new estimate will be computed.
%                   Defaults to 1e-8.
% - deg_fctr_max: 4x2 'double' array of integers, specifying the maximal
%                 factor with which each monomial degree can be
%                 multiplied in searching for an accurate inverse. If no
%                 operator Xop has been found such that the error in
%                 Xop*P1op-P2op is below the tolerance with deg_fctr =
%                 deg_fctr_max, the current best guess will be returned.
%
% OUTPUTS:
% - Xop:            [m0,p0; mx,px; my,py; m2,p2] opvar2d object, solving
%                   the system Xop*P1op = P2op. This operator is computed
%                   by representing P1op and P2op by standard monomial
%                   functions/operators, scaled by matrices, solving a
%                   matrix system X*A = B to establish the matrices
%                   defining Xop. NOTE that polynomial functions in general
%                   do not have polynomial inverses, so the solution may
%                   only be approximate in most cases.
% - eps:            Estimate of the error in the solution, computed as the
%                   maximal absolute value of the elements of X*A-B.
% - deg_fctr_final: Final value of deg_fctr corresponding to the returned
%                   solution Xop.
%
%
% NOTES:
% The function assumes the parameters in the operator Xop to be polynomial,
% just like the parameters in P2op and P1op. The degree of the parameters
% in Xop is unknown, and increasingly large degrees will be used to
% approximate an accurate solution Xop. As a result, the function may be
% computationally demanding, and the output may not be accurate, depending
% on the combination of P2op and P1op.   
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
%
% Initial coding DJ - 07/11/2022
% 08/03/2022, DJ: Call exact inverse function if P1op is separable.


% % % Extract the inputs
if nargin<2
    error(['At least two arguments P1 and P2 are necessary to compute P2/P1']);
elseif nargin==2
    tol = 1e-8;
    deg_fctr = 2*ones(4,2);
    deg_fctr_max = [10,10; 10,10; 5,5; 5,5];
elseif nargin==3
    tol = 1e-8;
    deg_fctr_max = [10,10; 10,10; 5,5; 5,5];
elseif nargin==4
    if all(size(tol)==[4,2])
        deg_fctr_max = tol;
        tol = 1e-8;
    else
        deg_fctr_max = [10,10; 10,10; 5,5; 5,5];
    end
elseif nargin>5
    error('At most 5 inputs are accepted')
end
% Keep track of how many iterations have already been performed.
global iter_num_mrdivide_opvar2d


% % % Check and process the inputs
if ~isa(P1op,'opvar2d')
    error('Input P1 in P2/P1 must be an opvar2d class object.')
end
P1dim = P1op.dim;
var1 = P1op.var1;     var2 = P1op.var2;
if isa(P2op,'double') && all(size(P2op)==1)
    P2dim = [P1dim(:,2), P1dim(:,2)];
    P2op = P2op * opvar2d(eye(size(P1op,2)),P2dim,P1op.I,var1,var2);
elseif ~isa(P2op,'opvar2d')
    if size(P2op,2)~=size(P1op,2)
        error('The in dimensions of operators P1 and P2 in P2/P1 must match.')
    end
    P2dim = [P1dim(:,2),P1dim(:,2)];
    try P2op = opvar2d(P2op,P2dim,P1op.I,var1,var2);    % Not quite....
    catch
        error('Please specify P2 as an opvar2d object when calling P2/P1.')
    end
else
    P2dim = P2op.dim;
    if any(P1dim(:,2)~=P2dim(:,2))
        error('The input dimensions of operators P1 and P2 in P2/P1 must match.')
    elseif any(~isequal(P1op.I,P2op.I))
        error('Cannot compute P2/P1: Domains of the operators P1 and P2 don''t match.')
    elseif any(~isequal(P1op.var1,P2op.var1)) || any(~isequal(P1op.var2,P2op.var2))
        error('Cannot compute P2/P1: Spatial variables of the operators P1 and P2 don''t match.')
    end
end

if tol<=0 || ~isreal(tol)
    error('The provided tolerance must be a real scalar value, greater than zero.')
end
if numel(deg_fctr)==1
    deg_fctr = deg_fctr * ones(4,2);
elseif any(size(deg_fctr)~=[4,2])
    error('The degree increase factors should be specified as a 4x2 array.')
end
if any(deg_fctr<=0)
    error('The degree increase factors should be strictly positive.')
end
if numel(deg_fctr_max)==1
    deg_fctr_max = deg_fctr_max * ones(4,2);
elseif any(size(deg_fctr_max)~=[4,2])
    error('The maximal degree increase factors should be specified as a 4x2 array.')
end
if any(deg_fctr_max<=0)
    error('The maximal degree increase factors should be strictly positive.')
end
if ~isa(P1op,'opvar2d')
    error('Input argument must be an opvar2d class object.')
end


P1dim = P1op.dim;
xx = var1(1);   tt = var2(1);
yy = var1(2);   nu = var2(2);


% % % First, we establish what monomials appear in the operator Pop

% % Test for separability of the kernels
sep_Rxx = any(P1dim(2,:)==0) || all(all(isequal(P1op.Rxx{2},P1op.Rxx{3})));
sep_Rx2 = P1dim(2,1)==0 || P1dim(4,2)==0 || all(all(isequal(P1op.Rx2{2},P1op.Rx2{3})));
sep_R2x = P1dim(4,1)==0 || P1dim(2,2)==0 || all(all(isequal(P1op.R2x{2},P1op.R2x{3})));

sep_Ryy = any(P1dim(3,:)==0) || all(all(isequal(P1op.Ryy{2},P1op.Ryy{3}))); 
sep_Ry2 = P1dim(3,1)==0 || P1dim(4,2)==0 || all(all(isequal(P1op.Ry2{2},P1op.Ry2{3})));
sep_R2y = P1dim(4,1)==0 || P1dim(3,2)==0 || all(all(isequal(P1op.R2y{2},P1op.R2y{3})));

sep_R22_ao = any(P1dim(4,:)==0) || all(all(isequal(P1op.R22{2,1},P1op.R22{3,1})));
sep_R22_oa = any(P1dim(4,:)==0) || all(all(isequal(P1op.R22{1,2},P1op.R22{1,3})));
sep_R22_aa_x = any(P1dim(4,:)==0) || (all(all(isequal(P1op.R22{2,2},P1op.R22{3,2}))) && all(all(isequal(P1op.R22{2,3},P1op.R22{3,3}))));
sep_R22_aa_y = any(P1dim(4,:)==0) || (all(all(isequal(P1op.R22{2,2},P1op.R22{2,3}))) && all(all(isequal(P1op.R22{3,2},P1op.R22{3,3}))));


% % If the operator is separable and maps only 1D functions, see if we can
% % use the algebraic inverse.
if sep_Rxx && sep_Ryy && all(P1dim(end,:)==0)

    isconstant_Rxx_o = isa(P1op.Rxx{1},'double') || isdouble(P1op.Rxx{1});
    isconstant_Ryy_o = isa(P1op.Ryy{1},'double') || isdouble(P1op.Ryy{1});

    if isconstant_Rxx_o && isconstant_Ryy_o
        Rxx0 = double(P1op.Rxx{1});
        sgm_x = svd(Rxx0);
        Ryy0 = double(P1op.Ryy{1});
        sgm_y = svd(Ryy0);
        if (isempty(Rxx0) || sgm_x(end)/sgm_x(1) > tol && length(sgm_x) >= size(Rxx0,2)) && ...
                (isempty(Ryy0) || sgm_y(end)/sgm_y(1) > tol && length(sgm_y) >= size(Ryy0,2))
            % Compute Xop using exact inverse of P1op
            [P1inv] = inv_opvar2d_separable(P1op',tol); % P1inv * P1op' = I --> P1op * P1inv' = I
            Xop = P2op * P1inv';        % Xop*P1op = P2op --> Xop = P2op * P1inv
            
            % Return a (heuristic) estimate of the accuracy of the result.
            ERRop = Xop * P1op - P2op;
            eps = 0;
            Rparams = {'R00','R0x','R0y','R02';
                       'Rx0','Rxx','Rxy','Rx2';
                       'Ry0','Ryx','Ryy','Ry2';
                       'R20','R2x','R2y','R22'};
            for kk=1:numel(Rparams)
                PR = ERRop.(Rparams{kk});
                if isa(PR,'double')
                    eps = max([eps,max(max(abs(PR)))]);
                elseif isa(PR,'polynomial')
                    eps = max([eps,max(max(abs(PR.coeff)))]);
                elseif isa(PR,'cell')
                    for ll=1:numel(PR)
                        PRR = PR{ll};
                        if isa(PRR,'double')
                            eps = max([eps,max(max(abs(PRR)))]);
                        elseif isa(PRR,'polynomial')
                            eps = max([eps,max(max(abs(PRR.coeff)))]);
                        else
                            error(['Parameter ',Rparams{kk},'{',num2str(ll),'} is not appropriately defined.']);
                        end
                    end
                else
                    error(['Parameter ',Rparams{kk},' is not appropriately defined.']);
                end
            end
            deg_fctr_final = deg_fctr;
            return
        end
    end
end


% % Establish degree matrices for the monomial bases Z defining Pop
[~,maxdegs_P1,~] = getdeg(P1op);
Zx_deg_P1 = max([maxdegs_P1.Rx0(2,1);maxdegs_P1.Rxx{1}(2,1);maxdegs_P1.Rxx{2}(2,1);maxdegs_P1.Rxx{3}(2,1);
                 maxdegs_P1.Rxy(2,1);maxdegs_P1.Rx2{1}(2,1);maxdegs_P1.Rx2{2}(2,1);maxdegs_P1.Rx2{3}(2,1)]);
Zy_deg_P1 = max([maxdegs_P1.Ry0(1,2);maxdegs_P1.Ryy{1}(1,2);maxdegs_P1.Ryy{2}(1,2);maxdegs_P1.Ryy{3}(1,2);
                 maxdegs_P1.Ryx(1,2);maxdegs_P1.Ry2{1}(1,2);maxdegs_P1.Ry2{2}(1,2);maxdegs_P1.Ry2{3}(1,2)]);
Zxy_deg_x_P1 = max([maxdegs_P1.R20(2,1);     maxdegs_P1.R2x{1}(2,1);  maxdegs_P1.R2x{2}(2,1);maxdegs_P1.R2x{3}(2,1);
                    maxdegs_P1.R2y{1}(2,1);  maxdegs_P1.R2y{2}(2,1);  maxdegs_P1.R2y{3}(2,1);
                    maxdegs_P1.R22{1,1}(2,1);maxdegs_P1.R22{1,2}(2,1);maxdegs_P1.R22{1,3}(2,1);
                    maxdegs_P1.R22{2,1}(2,1);maxdegs_P1.R22{2,2}(2,1);maxdegs_P1.R22{2,3}(2,1)
                    maxdegs_P1.R22{3,1}(2,1);maxdegs_P1.R22{3,2}(2,1);maxdegs_P1.R22{3,3}(2,1)]);
Zxy_deg_y_P1 = max([maxdegs_P1.R20(1,2);     maxdegs_P1.R2x{1}(1,2);  maxdegs_P1.R2x{2}(1,2);maxdegs_P1.R2x{3}(1,2);
                    maxdegs_P1.R2y{1}(1,2);  maxdegs_P1.R2y{2}(1,2);  maxdegs_P1.R2y{3}(1,2);
                    maxdegs_P1.R22{1,1}(1,2);maxdegs_P1.R22{1,2}(1,2);maxdegs_P1.R22{1,3}(1,2);
                    maxdegs_P1.R22{2,1}(1,2);maxdegs_P1.R22{2,2}(1,2);maxdegs_P1.R22{2,3}(1,2)
                    maxdegs_P1.R22{3,1}(1,2);maxdegs_P1.R22{3,2}(1,2);maxdegs_P1.R22{3,3}(1,2)]);

Ztt_deg_P1 = max([maxdegs_P1.R0x(2,1);maxdegs_P1.Rxx{2}(3,1);maxdegs_P1.Rxx{3}(3,1);
                  maxdegs_P1.Ryx(2,1);maxdegs_P1.R2x{2}(3,1);maxdegs_P1.R2x{3}(3,1)]);
Znu_deg_P1 = max([maxdegs_P1.R0y(1,2);maxdegs_P1.Ryy{2}(1,3);maxdegs_P1.Ryy{3}(1,3);
                  maxdegs_P1.Rxy(1,2);maxdegs_P1.R2y{2}(1,3);maxdegs_P1.R2y{3}(1,3)]);
Zttnu_deg_tt_P1 = max([maxdegs_P1.R02(2,1);     maxdegs_P1.Rx2{2}(3,1);  maxdegs_P1.Rx2{3}(3,1);
                       maxdegs_P1.Ry2{2}(2,1);  maxdegs_P1.R2y{3}(2,1);
                       maxdegs_P1.R22{2,1}(3,1);maxdegs_P1.R22{2,2}(3,1);maxdegs_P1.R22{2,3}(3,1)
                       maxdegs_P1.R22{3,1}(3,1);maxdegs_P1.R22{3,2}(3,1);maxdegs_P1.R22{3,3}(3,1)]);
Zttnu_deg_nu_P1 = max([maxdegs_P1.R02(1,2);     maxdegs_P1.Rx2{2}(1,2);  maxdegs_P1.Rx2{3}(1,2);
                       maxdegs_P1.Ry2{2}(1,3);  maxdegs_P1.Ry2{3}(1,3);
                       maxdegs_P1.R22{1,2}(1,3);maxdegs_P1.R22{1,3}(1,3);
                       maxdegs_P1.R22{2,2}(1,3);maxdegs_P1.R22{2,3}(1,3)
                       maxdegs_P1.R22{3,2}(1,3);maxdegs_P1.R22{3,3}(1,3)]);
            
% Establish maximal degrees for operator P2
[~,maxdegs_P2,~] = getdeg(P2op);
Zx_deg_P2 = max([maxdegs_P2.Rx0(2,1);maxdegs_P2.Rxx{1}(2,1);maxdegs_P2.Rxx{2}(2,1);maxdegs_P2.Rxx{3}(2,1);
                 maxdegs_P2.Rxy(2,1);maxdegs_P2.Rx2{1}(2,1);maxdegs_P2.Rx2{2}(2,1);maxdegs_P2.Rx2{3}(2,1)]);
Zy_deg_P2 = max([maxdegs_P2.Ry0(1,2);maxdegs_P2.Ryy{1}(1,2);maxdegs_P2.Ryy{2}(1,2);maxdegs_P2.Ryy{3}(1,2);
                 maxdegs_P2.Ryx(1,2);maxdegs_P2.Ry2{1}(1,2);maxdegs_P2.Ry2{2}(1,2);maxdegs_P2.Ry2{3}(1,2)]);
Zxy_deg_x_P2 = max([maxdegs_P2.R20(2,1);     maxdegs_P2.R2x{1}(2,1);  maxdegs_P2.R2x{2}(2,1);maxdegs_P2.R2x{3}(2,1);
                    maxdegs_P2.R2y{1}(2,1);  maxdegs_P2.R2y{2}(2,1);  maxdegs_P2.R2y{3}(2,1);
                    maxdegs_P2.R22{1,1}(2,1);maxdegs_P2.R22{1,2}(2,1);maxdegs_P2.R22{1,3}(2,1);
                    maxdegs_P2.R22{2,1}(2,1);maxdegs_P2.R22{2,2}(2,1);maxdegs_P2.R22{2,3}(2,1)
                    maxdegs_P2.R22{3,1}(2,1);maxdegs_P2.R22{3,2}(2,1);maxdegs_P2.R22{3,3}(2,1)]);
Zxy_deg_y_P2 = max([maxdegs_P2.R20(1,2);     maxdegs_P2.R2x{1}(1,2);  maxdegs_P2.R2x{2}(1,2);maxdegs_P2.R2x{3}(1,2);
                    maxdegs_P2.R2y{1}(1,2);  maxdegs_P2.R2y{2}(1,2);  maxdegs_P2.R2y{3}(1,2);
                    maxdegs_P2.R22{1,1}(1,2);maxdegs_P2.R22{1,2}(1,2);maxdegs_P2.R22{1,3}(1,2);
                    maxdegs_P2.R22{2,1}(1,2);maxdegs_P2.R22{2,2}(1,2);maxdegs_P2.R22{2,3}(1,2)
                    maxdegs_P2.R22{3,1}(1,2);maxdegs_P2.R22{3,2}(1,2);maxdegs_P2.R22{3,3}(1,2)]);

Ztt_deg_P2 = max([maxdegs_P2.R0x(2,1);maxdegs_P2.Rxx{2}(3,1);maxdegs_P2.Rxx{3}(3,1);
                  maxdegs_P2.Ryx(2,1);maxdegs_P2.R2x{2}(3,1);maxdegs_P2.R2x{3}(3,1)]);
Znu_deg_P2 = max([maxdegs_P2.R0y(1,2);maxdegs_P2.Ryy{2}(1,3);maxdegs_P2.Ryy{3}(1,3);
                  maxdegs_P2.Rxy(1,2);maxdegs_P2.R2y{2}(1,3);maxdegs_P2.R2y{3}(1,3)]);
Zttnu_deg_tt_P2 = max([maxdegs_P2.R02(2,1);     maxdegs_P2.Rx2{2}(3,1);  maxdegs_P2.Rx2{3}(3,1);
                       maxdegs_P2.Ry2{2}(2,1);  maxdegs_P2.R2y{3}(2,1);
                       maxdegs_P2.R22{2,1}(3,1);maxdegs_P2.R22{2,2}(3,1);maxdegs_P2.R22{2,3}(3,1)
                       maxdegs_P2.R22{3,1}(3,1);maxdegs_P2.R22{3,2}(3,1);maxdegs_P2.R22{3,3}(3,1)]);
Zttnu_deg_nu_P2 = max([maxdegs_P2.R02(1,2);     maxdegs_P2.Rx2{2}(1,2);  maxdegs_P2.Rx2{3}(1,2);
                       maxdegs_P2.Ry2{2}(1,3);  maxdegs_P2.Ry2{3}(1,3);
                       maxdegs_P2.R22{1,2}(1,3);maxdegs_P2.R22{1,3}(1,3);
                       maxdegs_P2.R22{2,2}(1,3);maxdegs_P2.R22{2,3}(1,3)
                       maxdegs_P2.R22{3,2}(1,3);maxdegs_P2.R22{3,3}(1,3)]);
            
                
% % % Now we start working on Xop, assuming no separability
Xdim = [P2dim(:,1), P1dim(:,1)];

% sep_iRxx = sep_Rxx;     sep_iRx2 = sep_R2x;     sep_iR2x = sep_Rx2;
% sep_iRyy = sep_Ryy;     sep_iRy2 = sep_R2y;     sep_iR2y = sep_Ry2;
% sep_iR22_ao = sep_R22_ao;           sep_iR22_oa = sep_R22_oa;
% sep_iR22_aa_x = sep_R22_aa_x;       sep_iR22_aa_y = sep_R22_aa_y;

sep_iRxx = false;     sep_iRx2 = false;     sep_iR2x = false;
sep_iRyy = false;     sep_iRy2 = false;     sep_iR2y = false;
sep_iR22_ao = false;        sep_iR22_oa = false;
sep_iR22_aa_x = false;      sep_iR22_aa_y = false;


% % % We decompose the operator Xop as
% % Xop = iIZL(x,y)^T * iHmat * iIZRop
% %
% %  = [iZ0                           ]^T *
% %    [    iIZx(x)                   ]
% %    [            iIZy(y)           ]
% %    [                    iIZxy(x,y)]
% % 
% % [ iH00 0      iH0x   iH0x   0      iH0y   iH0y   0       0       0       0       0       iH02    iH02    iH02    iH02   ] * 
% % [ iHx0 iHxx_o iHxx_a iHxx_b 0      iHxy   iHxy   0       0       0       iHx2_o  iHx2_o  iHx2_a  iHx2_b  iHx2_a  iHx2_b ]
% % [ iHy0 0      iHyx   iHyx   iHyy_o iHyy_a iHyy_b 0       iHy2_o  iHy2_o  0       0       iHy2_a  iHy2_a  iHy2_b  iHy2_b ]
% % [ iH20 iH2x_o oH2x_a iH2x_b iH2y_o iH2y_a iH2y_b iH22_oo iH22_ao iH22_bo iH22_oa iH22_ob iH22_aa iH22_ba iH22_ab iH22_bb]
% %
% %    [iZ0                                                                                         ] 
% %    [    I                                                                                       ] 
% %    [    int_{a}^{x}dtt iIZtt(tt)                                                                ] 
% %    [    int_{x}^{b}dtt iIZtt(tt)                                                                ] 
% %    [                             I                                                              ] 
% %    [                             int_{c}^{y}dnu iIZnu(nu)                                       ] 
% %    [                             int_{y}^{d}dnu iIZnu(nu)                                       ] 
% %    [                                                      I                                     ] 
% %    [                                                      int_{a}^{x}dtt iIZttnu_tt(tt)         ] 
% %    [                                                      int_{x}^{b}dtt iIZttnu_tt(tt)         ] 
% %    [                                                      int_{c}^{y}dnu iIZttnu_nu(nu)         ] 
% %    [                                                      int_{y}^{d}dnu iIZttnu_nu(nu)         ]
% %    [                                                      int_{a,c}^{x,y}dtt dnu iIZttnu(tt,nu) ] 
% %    [                                                      int_{x,c}^{b,y}dtt dnu iIZttnu(tt,nu) ] 
% %    [                                                      int_{a,y}^{x,d}dtt dnu iIZttnu(tt,nu) ] 
% %    [                                                      int_{x,y}^{b,d}dtt dnu iIZttnu(tt,nu) ] 
% %
% % Where e.g. iIZx(x) = I_{imx} \otimes iZx(x),      iIZy(y) = I_{imx} \otimes iZy(y)
% %            iIZtt(tt) = I_{inx} \otimes iZtt(tt),  iIZnu(nu) = I_{iny} \otimes iZnu(nu)
% % for dimension iPdim = [im0, in0; imx, inx; imy, iny; im2, in2];
% % We first establish the degree matrices of the different monomials in
% %this decomposition

% Increase the monomial degrees by a factor deg_fctr to establish monomial 
% degrees for the operator Xop.
deg_boost = ones(4,2);
%deg_fctr = ones(4,2);

iZx_deg = round(deg_fctr(1,1)*(max(Zx_deg_P2-Zx_deg_P1,Zx_deg_P1) + deg_boost(1,1)));
iZy_deg = round(deg_fctr(2,1)*(max(Zy_deg_P2-Zy_deg_P1,Zy_deg_P1) + deg_boost(2,1)));
iZxy_deg_x = round(deg_fctr(3,1)*(max(Zxy_deg_x_P2-Zxy_deg_x_P1,Zxy_deg_x_P1) + deg_boost(3,1)));
iZxy_deg_y = round(deg_fctr(4,1)*(max(Zxy_deg_y_P2-Zxy_deg_y_P1,Zxy_deg_y_P1) + deg_boost(4,1)));

iZtt_deg = round(deg_fctr(1,2)*(max(Ztt_deg_P2-Ztt_deg_P1,Ztt_deg_P1) + deg_boost(1,2)));
iZnu_deg = round(deg_fctr(2,2)*(max(Znu_deg_P2-Znu_deg_P1,Znu_deg_P1) + deg_boost(2,2)));
iZttnu_deg_tt = round(deg_fctr(3,2)*(max(Zttnu_deg_tt_P2-Zttnu_deg_tt_P1,Zttnu_deg_tt_P1) + deg_boost(3,2)));
iZttnu_deg_nu = round(deg_fctr(4,2)*(max(Zttnu_deg_nu_P2-Zttnu_deg_nu_P1,Zttnu_deg_nu_P1) + deg_boost(4,2)));

iZx_deg = max([iZx_deg,iZtt_deg]);
iZtt_deg = max([iZx_deg,iZtt_deg]);
iZy_deg = max([iZy_deg,iZnu_deg]);
iZnu_deg = max([iZy_deg,iZnu_deg]);
iZxy_deg_x = max([iZxy_deg_x,iZttnu_deg_tt]);
iZttnu_deg_tt = max([iZxy_deg_x,iZttnu_deg_tt]);
iZxy_deg_y = max([iZxy_deg_y,iZttnu_deg_nu]);
iZttnu_deg_nu = max([iZxy_deg_y,iZttnu_deg_nu]);

% If the inverse has no parameters mapping to e.g. L2[x], include only the
% constant monomial.
if Xdim(2,1)==0
    iZx_deg = 1;
    deg_fctr(1,1) = 1;
end
if Xdim(2,2)==0
    iZtt_deg = 1;
    deg_fctr(1,2) = 1;
end
if Xdim(3,1)==0
    iZy_deg = 1;
    deg_fctr(2,1) = 1;
end
if Xdim(3,2)==0
    iZnu_deg = 1;
    deg_fctr(2,2) = 1;
end
if Xdim(4,1)==0
    iZxy_deg_x = 1;
    iZxy_deg_y = 1;
    deg_fctr([3;4],1) = [1;1];
end
if Xdim(4,2)==0
    iZttnu_deg_tt = 1;
    iZttnu_deg_nu = 1;
    deg_fctr([3;4],2) = [1;1];
end

% iZx_deg = max([iZx_deg,iZtt_deg]);
% iZtt_deg = max([iZx_deg,iZtt_deg]);
% iZy_deg = max([iZy_deg,iZnu_deg]);
% iZnu_deg = max([iZy_deg,iZnu_deg]);
% iZxy_deg_x = max([iZxy_deg_x,iZttnu_deg_tt]);
% iZttnu_deg_tt = max([iZxy_deg_x,iZttnu_deg_tt]);
% iZxy_deg_y = max([iZxy_deg_y,iZttnu_deg_nu]);
% iZttnu_deg_nu = max([iZxy_deg_y,iZttnu_deg_nu]);

iZtt_dmat = monomials(1,0:iZtt_deg);
iZnu_dmat = monomials(1,0:iZnu_deg);
iZttnu_tt_dmat = monomials(1,0:iZttnu_deg_tt);
iZttnu_nu_dmat = monomials(1,0:iZttnu_deg_nu);
iZttnu_dmat = [repmat(iZttnu_tt_dmat,[length(iZttnu_nu_dmat),1]),kron(iZttnu_nu_dmat,ones([length(iZttnu_tt_dmat),1]))];

iZx_dmat = monomials(1,0:iZx_deg);
iZy_dmat = monomials(1,0:iZy_deg);
iZxy_x_dmat = monomials(1,0:iZxy_deg_x);
iZxy_y_dmat = monomials(1,0:iZxy_deg_y);
iZxy_dmat = [repmat(iZxy_x_dmat,[length(iZxy_y_dmat),1]),kron(iZxy_y_dmat,ones([length(iZxy_x_dmat),1]))];

inZx = size(iZx_dmat,1);
inZy = size(iZy_dmat,1);
inZxy = size(iZxy_dmat,1);
inZL_arr = [1; inZx; inZy; inZxy];

inZtt = size(iZtt_dmat,1);
inZnu = size(iZnu_dmat,1);
inZttnu = size(iZttnu_dmat,1);
inZttnu_tt = size(iZttnu_tt_dmat,1);
inZttnu_nu = size(iZttnu_nu_dmat,1);
inZR_arr = [1, inZtt, inZnu, inZttnu];

% Construct the acutal polynomials
iZtt = polynomial(speye(inZtt),iZtt_dmat,tt.varname,[inZtt,1]);
iZnu = polynomial(speye(inZnu),iZnu_dmat,nu.varname,[inZnu,1]);
iZttnu = polynomial(speye(inZttnu),iZttnu_dmat,[tt.varname;nu.varname],[inZttnu,1]);
iZttnu_tt = polynomial(speye(inZttnu_tt),iZttnu_tt_dmat,tt.varname,[inZttnu_tt,1]);
iZttnu_nu = polynomial(speye(inZttnu_nu),iZttnu_nu_dmat,nu.varname,[inZttnu_nu,1]);


% % % Next, build the operator iIZRop defining the inverse 
% % % Pop_inv = iIZL(x,y)^T * iHmat * iIZRop
% Rows mapping to R^m
iIZRop_0 = opvar2d();
iIZRop_0.R00 = eye(Xdim(1,2));

% Rows mapping to L2[x]
iIZRop_x = opvar2d();
iIZRop_x.Rxx{1} = [eye(Xdim(2,2)); zeros(2*Xdim(2,2)*inZtt,Xdim(2,2))];
iIZRop_x.Rxx{2} = [zeros(Xdim(2,2)); kron(eye(Xdim(2,2)),iZtt); zeros(Xdim(2,2)*inZtt,Xdim(2,2))];
iIZRop_x.Rxx{3} = [zeros(Xdim(2,2)*(1+inZtt),Xdim(2,2)); kron(eye(Xdim(2,2)),iZtt)];

% Rows mapping to L2[y]
iIZRop_y = opvar2d();
iIZRop_y.Ryy{1} = [eye(Xdim(3,2)); zeros(2*Xdim(3,2)*inZnu,Xdim(3,2))];
iIZRop_y.Ryy{2} = [zeros(Xdim(3,2)); kron(eye(Xdim(3,2)),iZnu); zeros(Xdim(3,2)*inZnu,Xdim(3,2))];
iIZRop_y.Ryy{3} = [zeros(Xdim(3,2)*(1+inZnu),Xdim(3,2)); kron(eye(Xdim(3,2)),iZnu)];

% Rows mapping to L2[x,y]
iIZRop_2 = opvar2d();
iIZRop_2.R22{1,1} = [eye(Xdim(4,2)); zeros(Xdim(4,2)*(2*inZttnu_tt + 2*inZttnu_nu + 4*inZttnu),Xdim(4,2))];
iIZRop_2.R22{2,1} = [zeros(Xdim(4,2)); kron(eye(Xdim(4,2)),iZttnu_tt); zeros(Xdim(4,2)*(inZttnu_tt + 2*inZttnu_nu + 4*inZttnu),Xdim(4,2))];
iIZRop_2.R22{3,1} = [zeros(Xdim(4,2)*(1 + inZttnu_tt),Xdim(4,2)); kron(eye(Xdim(4,2)),iZttnu_tt); zeros(Xdim(4,2)*(2*inZttnu_nu + 4*inZttnu),Xdim(4,2))];
iIZRop_2.R22{1,2} = [zeros(Xdim(4,2)*(1 + 2*inZttnu_tt),Xdim(4,2)); kron(eye(Xdim(4,2)),iZttnu_nu); zeros(Xdim(4,2)*(inZttnu_nu + 4*inZttnu),Xdim(4,2))];
iIZRop_2.R22{1,3} = [zeros(Xdim(4,2)*(1 + 2*inZttnu_tt + inZttnu_nu),Xdim(4,2)); kron(eye(Xdim(4,2)),iZttnu_nu); zeros(Xdim(4,2)*4*inZttnu,Xdim(4,2))];
iIZRop_2.R22{2,2} = [zeros(Xdim(4,2)*(1 + 2*inZttnu_tt + 2*inZttnu_nu),Xdim(4,2)); kron(eye(Xdim(4,2)),iZttnu); zeros(Xdim(4,2)*3*inZttnu,Xdim(4,2))];
iIZRop_2.R22{3,2} = [zeros(Xdim(4,2)*(1 + 2*inZttnu_tt + 2*inZttnu_nu + inZttnu),Xdim(4,2)); kron(eye(Xdim(4,2)),iZttnu); zeros(Xdim(4,2)*2*inZttnu,Xdim(4,2))];
iIZRop_2.R22{2,3} = [zeros(Xdim(4,2)*(1 + 2*inZttnu_tt + 2*inZttnu_nu + 2*inZttnu),Xdim(4,2)); kron(eye(Xdim(4,2)),iZttnu); zeros(Xdim(4,2)*inZttnu,Xdim(4,2))];
iIZRop_2.R22{3,3} = [zeros(Xdim(4,2)*(1 + 2*inZttnu_tt + 2*inZttnu_nu + 3*inZttnu),Xdim(4,2)); kron(eye(Xdim(4,2)),iZttnu)];

% Combine
iIZRop = blkdiag(iIZRop_0,iIZRop_x,iIZRop_y,iIZRop_2);
iIZRop.I = P1op.I;
%iIZRop.var1 = P1op.var1; iIZRop.var2 = P1op.var2;

% % Decomposing the inverse as Xop = iIZL(x,y)^T * iH_mat * iIZRop,
% % we get
% %     Xop*P1op = iIZR(x,y)^T * iHmat * (iIZRop * P1op)
% % We compute this product ZPop = iIZRop * P1op
ZPop = poly_opvar2d(iIZRop * P1op);
ZPop = struct(ZPop);    % Convert opvar2d to struct to be able to adjust varnames
ZPdim = ZPop.dim;

% % Then decompose the operator ZPop in a similar way as Xop:
% %     ZPop = IZL(x,y) * Hmat * IZRop 
% % For this, we first establish appropriate monomials to define IZL(x,y)
% % and IZRop. For each set of monomials, also establish "pindcs" that
% % indicates where each monomial that appears in each parameter in ZPop
% % appears in the combined monomial vector

% % To perform the decomposition, in the integral parameters, we replace 
% % the variables x,y with the dummy variables tt,nu
ZPop.R0x.varname(ismember(ZPop.R0x.varname,xx.varname)) = tt.varname;
ZPop.Ryx.varname(ismember(ZPop.Ryx.varname,xx.varname)) = tt.varname;
ZPop.R02.varname(ismember(ZPop.R02.varname,xx.varname)) = tt.varname;

ZPop.Ry2{1}.varname(ismember(ZPop.Ry2{1}.varname,xx.varname)) = tt.varname;
ZPop.Ry2{2}.varname(ismember(ZPop.Ry2{2}.varname,xx.varname)) = tt.varname;
ZPop.Ry2{3}.varname(ismember(ZPop.Ry2{3}.varname,xx.varname)) = tt.varname;

ZPop.R0y.varname(ismember(ZPop.R0y.varname,yy.varname)) = nu.varname;
ZPop.Rxy.varname(ismember(ZPop.Rxy.varname,yy.varname)) = nu.varname;
ZPop.R02.varname(ismember(ZPop.R02.varname,yy.varname)) = nu.varname;

ZPop.Rx2{1}.varname(ismember(ZPop.Rx2{1}.varname,yy.varname)) = nu.varname;
ZPop.Rx2{2}.varname(ismember(ZPop.Rx2{2}.varname,yy.varname)) = nu.varname;
ZPop.Rx2{3}.varname(ismember(ZPop.Rx2{3}.varname,yy.varname)) = nu.varname;

% Constant monomial Z0 = 1
%Z0_dmat = zeros(1,0);
Z0L_pindcs = {ones(size(ZPop.R00.degmat,1),1),...
                ones(size(ZPop.R0x.degmat,1),1),...
                  ones(size(ZPop.R0y.degmat,1),1),...
                    ones(size(ZPop.R02.degmat,1),1)};
Z0R_pindcs = {ones(size(ZPop.R00.degmat,1),1);
                ones(size(ZPop.Rx0.degmat,1),1);
                  ones(size(ZPop.Ry0.degmat,1),1);
                    ones(size(ZPop.R20.degmat,1),1)};
                  
% Zx(x)
Zx_params = {'Rx0','Rxx','Rxy','Rx2'};
[Zx_dmat, Zx_pindcs] = collect_monomials(ZPop,Zx_params,xx.varname);
nZx = size(Zx_dmat,1);

% Zy(y)
Zy_params = {'Ry0','Ryx','Ryy','Ry2'};
[Zy_dmat, Zy_pindcs] = collect_monomials(ZPop,Zy_params,yy.varname);
nZy = size(Zy_dmat,1);

% Zxy(x,y)
Zxy_params = {'R20','R2x','R2y','R22'};
[Zxy_dmat, Zxy_pindcs] = collect_monomials(ZPop,Zxy_params,[xx.varname;yy.varname]);
nZxy = size(Zxy_dmat,1);

% Ztt(tt)
Ztt_params = {'R0x';'Rxx';'Ryx';'R2x'};
[Ztt_dmat, Ztt_pindcs] = collect_monomials(ZPop,Ztt_params,tt.varname);
nZtt = size(Ztt_dmat,1);

% Znu(nu)
Znu_params = {'R0y';'Rxy';'Ryy';'R2y'};
[Znu_dmat, Znu_pindcs] = collect_monomials(ZPop,Znu_params,nu.varname);
nZnu = size(Znu_dmat,1);

% Zttnu(tt,nu)
Zttnu_params = {'R02';'Rx2';'Ry2';'R22'};
[Zttnu_dmat, Zttnu_pindcs] = collect_monomials(ZPop,Zttnu_params,[tt.varname;nu.varname]);
nZttnu = size(Zttnu_dmat,1);

% Zttnu_tt_indcs = Zttnu_dmat(:,2)==0;
% Zttnu_tt_dmat = Zttnu_dmat(Zttnu_tt_indcs,1);
% Zttnu_nu_indcs = Zttnu_dmat(:,1)==0;
% Zttnu_nu_dmat = Zttnu_dmat(Zttnu_nu_indcs,2);

% Collect values for later use
nZL_arr = [1; nZx; nZy; nZxy];
nZR_arr = [1, nZtt, nZnu, nZttnu];

ZL_indcs = [Z0L_pindcs; Zx_pindcs; Zy_pindcs; Zxy_pindcs];
ZR_indcs = [Z0R_pindcs, Ztt_pindcs, Znu_pindcs, Zttnu_pindcs];


% % % Now, build the matrix Hmat such that
% % %       ZPop = IZL(x,y) * Hmat * IZRop 
Rparams = {'R00','R0x','R0y','R02';
           'Rx0','Rxx','Rxy','Rx2';
           'Ry0','Ryx','Ryy','Ry2';
           'R20','R2x','R2y','R22'};

Hcell = cell(size(Rparams));
for k=1:numel(Rparams)
    [ii,jj] = ind2sub(size(Rparams),k);
    
    % Parameter of interest
    PR = ZPop.(Rparams{k});
    % Row and column dimension of matrix-valued parameter
    nr = ZPdim(ii,1);            nc = ZPdim(jj,2);
    % Size of left (e.g. Zx) and right (e.g. Ztt) monomial vectors associated to this parameter
    nZL = nZL_arr(ii);          nZR = nZR_arr(jj);
    % Indices of monomials in left and right monomial vectors associated to this parameter
    ZL_indx = ZL_indcs{ii,jj};  ZR_indx = ZR_indcs{ii,jj}; 
       
    if isempty(PR)
        Hcell{k} = zeros(nr*nZL,nc*nZR);
        
    elseif isa(PR,'double') %<-- shouldn't really be possible
        
        Cvals = reshape(PR,[],1);           % Reshape coefficients into a vector
        ridcs = repmat((1:nr)',[nc,1]);     % Row index in matrix-values object to which each coefficient corresponds
        cidcs = kron((1:nc)',ones(nr,1));   % Column index in matrix-values object to which each coefficient corresponds
        zidcs = ones(nr,nc);                % Index of monomial to which coefficient corresponds
        
        zLidcs = ZL_indx(zidcs);            % Index in unique monomials on left-hand side to which each coefficient corresponds
        Hridcs = zLidcs + (ridcs-1)*nZL;    % Row index for each coefficient in matrix H..
        zRidcs = ZR_indx(zidcs);            % Index in unique monomials on right-hand side to which each coefficient corresponds
        Hcidcs = zRidcs + (cidcs-1)*nZR;    % Column index for each coefficient in matrix H..
        
        Hcell{k} = sparse(Hridcs,Hcidcs,Cvals,nr*nZL,nc*nZR);
          
    elseif isa(PR,'polynomial')
        ndegs = size(PR.degmat,1);
        % Extract coefficient values and associated row and column indices
        % within the matrix-valued parameter associated to this coefficient
        C_old = PR.coeff;
        ridx_old = repmat((1:nr),[ndegs,nc]);
        cidx_old = repmat(kron((1:nc),ones(1,nr)),[ndegs,1]);

        Cvals = reshape(C_old,[],1);            % Reshape coefficients into a vector
        ridcs = reshape(ridx_old,[],1);         % Row index in matrix-valued object to which each coefficient corresponds
        cidcs = reshape(cidx_old,[],1);         % Column index in matrix-valued object to which each coefficient corresponds
        zidcs = repmat((1:ndegs)',[nr*nc,1]);   % Index of monomial to which coefficient corresponds

        zLidcs = ZL_indx(zidcs);            % Index in unique monomials on left-hand side to which each coefficient corresponds
        Hridcs = zLidcs + (ridcs-1)*nZL;    % Row index for each coefficient in matrix H..
        zRidcs = ZR_indx(zidcs);            % Index in unique monomials on right-hand side to which each coefficient corresponds
        Hcidcs = zRidcs + (cidcs-1)*nZR;    % Column index for each coefficient in matrix H..
        
        % Construct the matrix
        Hcell{k} = sparse(Hridcs,Hcidcs,Cvals,nr*nZL,nc*nZR);
        
    elseif isa(PR,'cell')
        PP = PR;
        ZL_indcs_k = ZL_indx;   ZR_indcs_k = ZR_indx;
        Hcell_k = cell(size(PP));
        for l=1:numel(PP)
            PR = PP{l};
            ZL_indx = ZL_indcs_k{l};     ZR_indx = ZR_indcs_k{l};
        
            if isempty(PR)
                Hcell_k{l} = zeros(nr*nZL, nc*nZR);
            elseif isa(PR,'double')        
                Cvals = reshape(PR,[],1);           % Reshape coefficients into a vector
                ridcs = repmat((1:nr)',[nc,1]);     % Row index in matrix-valued object to which each coefficient corresponds
                cidcs = kron((1:nc)',ones(nr,1));   % Column index in matrix-valued object to which each coefficient corresponds
                zidcs = ones(nr,nc);                % Index of monomial to which coefficient corresponds

                zLidcs = ZL_indx(zidcs);            % Index in unique monomials on left-hand side to which each coefficient corresponds
                Hridcs = zLidcs + (ridcs-1)*nZL;    % Row index for each coefficient in matrix H..
                zRidcs = ZR_indx(zidcs);            % Index in unique monomials on right-hand side to which each coefficient corresponds
                Hcidcs = zRidcs + (cidcs-1)*nZR;    % Column index for each coefficient in matrix H..

                Hcell_k{l} = sparse(Hridcs,Hcidcs,Cvals,nr*nZL,nc*nZR);
                
            elseif isa(PR,'polynomial')        
                ndegs = size(PR.degmat,1);
                C_old = PR.coeff;
                ridx_old = repmat((1:nr),[ndegs,nc]);
                cidx_old = repmat(kron((1:nc),ones(1,nr)),[ndegs,1]);

                Cvals = reshape(C_old,[],1);            % Reshape coefficients into a vector
                ridcs = reshape(ridx_old,[],1);         % Row index in matrix-valued object to which each coefficient corresponds
                cidcs = reshape(cidx_old,[],1);         % Column index in matrix-valued object to which each coefficient corresponds
                zidcs = repmat((1:ndegs)',[nr*nc,1]);   % Index of monomial to which coefficient corresponds

                zLidcs = ZL_indx(zidcs);            % Index in unique monomials on left-hand side to which each coefficient corresponds
                Hridcs = zLidcs + (ridcs-1)*nZL;    % Row index for each coefficient in new matrix H..
                zRidcs = ZR_indx(zidcs);            % Index in unique monomials on right-hand side to which each coefficient corresponds
                Hcidcs = zRidcs + (cidcs-1)*nZR;    % Column index for each coefficient in new matrix H..

                Hcell_k{l} = sparse(Hridcs,Hcidcs,Cvals,nr*nZL,nc*nZR);    
                
            else
                error('The opvar2d parameters must be specified as ''double'' or ''polynomial'' class objects');
            end
            Hcell{k} = Hcell_k;
        end
    else
        error('The opvar2d parameters must be specified as ''double'' or ''polynomial'' class objects');
    end    
end

% % Combine the different blocks into a single matrix
nIZL = ZPdim(:,1).*nZL_arr(:);
nIZR = ZPdim(:,2)'.*nZR_arr(:)';

Hmat = [Hcell{1,1}, zeros(nIZL(1),nIZR(2)), Hcell{1,2},    Hcell{1,2},    zeros(nIZL(1),nIZR(3)), Hcell{1,3},    Hcell{1,3},    ...
                zeros(nIZL(1),5*nIZR(4)), repmat(Hcell{1,4},[1,4]);
         Hcell{2,1}, [Hcell{2,2}{:}], zeros(nIZL(2),nIZR(3)), Hcell{2,3},    Hcell{2,3},    ...
                zeros(nIZL(2),3*nIZR(4)), Hcell{2,4}{1}, Hcell{2,4}{1}, Hcell{2,4}{2}, Hcell{2,4}{3}, Hcell{2,4}{2}, Hcell{2,4}{3};
         Hcell{3,1}, zeros(nIZL(3),nIZR(2)), Hcell{3,2},    Hcell{3,2},    [Hcell{3,3}{:}], ...
                zeros(nIZL(3),nIZR(4)),   Hcell{3,4}{1}, Hcell{3,4}{1}, zeros(nIZL(3),2*nIZR(4)), Hcell{3,4}{2}, Hcell{3,4}{2}, Hcell{3,4}{3}, Hcell{3,4}{3};
         Hcell{4,1}, Hcell{4,2}{1},      Hcell{4,2}{2}, Hcell{4,2}{3}, Hcell{4,3}{1},      Hcell{4,3}{2}, Hcell{4,3}{3}, ...
                [Hcell{4,4}{[1,2,3,4,7,5,6,8,9]}]];

% Further get rid of redundant monomials 
% (we really shouldn't be including them in the first place...)
dum = 0;        log_arr = true(size(Hmat,2),1);
dimval_arr = [1,2*ones(1,3),3*ones(1,3),4*ones(1,9)];
ZR_retain_indcs = cell(1,16);
nZR_retain = zeros(1,16);
for Ridx = 1:16
    nZ_k = nZR_arr(dimval_arr(Ridx));
    nc = P1dim(dimval_arr(Ridx),2);
    if nc==0
        continue
    end
    log_arr_k = true(nZ_k,1);
    for zidx = 1:nZ_k
        Zcidcs = dum + (0:nc-1)*nZ_k + zidx;        % Column indcs in Hmat associated to monomial specified by zidx
%        log_arr_k(zidx) = any(any(Hmat(:,Zcidcs))); % Does this monomial actually contribute?
    end
    
    ZR_retain_indcs{Ridx} = log_arr_k;
    nZR_retain(Ridx) = sum(log_arr_k);
    log_arr(dum+1:Zcidcs(end)) = repmat(log_arr_k,[nc,1]);
    dum = Zcidcs(end);
end
nIZR_arr = P2dim([1;2*ones(3,1);3*ones(3,1);4*ones(9,1)],2).*nZR_retain';
        
Hmat = Hmat(:,log_arr);

% % % At this point, we can represent
% % %   Xop*P1op = iIZR(x,y)^T * iHmat * (iIZRop * P1op)
% % %            = iIZR(x,y)^T * iHmat * IZL(x,y) * Hmat * IZRop
% % % We consider the central polynomial matrix, IZL(x,y), and represent it
% % % using a single vector of monomials ZZ(x,y), so that
% % IZL(x,y) 
% %  := [ I                                                                  ]
% %     [   I_{3} \otimes IZx(x)                                             ]
% %     [                        I_{3} \otimes IZy(y)                        ]
% %     [                                             I_{9}\otimes IZxy(x,y) ]
% %
% % = [ IZZ0(x,y)                                                                          ]^T * 
% %   [           I_{3} \otimes IZZx(x,y)                                                  ]     
% %   [                                   I_{3} \otimes IZZy(x,y)                          ]       
% %   [                                                           I_{9} \otimes IZZxy(x,y) ]     
% %                 [ C0                  ]
% %                 [    Cx_o             ]
% %                 [    Cx_a             ]
% %                 [    Cx_b             ]
% %                 [         Cy_o        ]
% %                 [         Cy_a        ]
% %                 [         Cy_b        ] 
% %                 [              Cxy_oo ]
% %                 [              Cxy_ao ]
% %                 [              Cxy_bo ]
% %                 [              Cxy_oa ]
% %                 [              Cxy_ob ]
% %                 [              Cxy_aa ]   
% %                 [              Cxy_ba ]   
% %                 [              Cxy_ab ]
% %                 [              Cxy_bb ]
% %
% % where IZZ0(x,y) = I_m0 \otimes (I_{nZ0} \otimes ZZ(x,y)), 
% %       IZZx(x,y) = I_mx \otimes (I_{nZx} \otimes ZZ(x,y)),
% %       IZZy(x,y) = I_my \otimes (I_{nZy} \otimes ZZ(x,y)),
% %       IZZxy(x,y) = I_m2 \otimes (I_{nZxy} \otimes ZZ(x,y)),
% % for row dimensions [m0; mx; my; m2] of Pop, and sizes
% % [nZ0; nZx; nZy; nZxy] of the monomial bases Z0, Zx, Zy and Zxy. 

% % First, construct the degree matrix for ZZ

% Combine all monomials of each product into a single vector
[ZZ_dmat,~,indx2] = uniquerows_integerTable([[0,0];[Zx_dmat,zeros(nZx,1)];[zeros(nZy,1),Zy_dmat];Zxy_dmat]);
nZZ = size(ZZ_dmat,1);

% Establish indices in the unique set that correspond to each of the
% non-unique monomials
ZZ_pindcs = cell(length(nZL_arr),1);
nnZL_arr = [0;cumsum(nZL_arr)];
for k=1:length(ZZ_pindcs)
    ZZ_pindcs{k} = indx2(nnZL_arr(k)+1 : nnZL_arr(k+1));
end

% % Next, build the C matrix
nrC = nZZ*sum(ZPdim(1:4,1));
ridcs_0 = reshape(ZZ_pindcs{1} + nZZ*(0:ZPdim(1,1)-1),[],1);
ridcs_x = nZZ*ZPdim(1,1) + reshape(ZZ_pindcs{2} + nZZ*(0:ZPdim(2,1)-1),[],1);
ridcs_y = nZZ*sum(ZPdim(1:2,1)) + reshape(ZZ_pindcs{3} + nZZ*(0:ZPdim(3,1)-1),[],1);
ridcs_xy = nZZ*sum(ZPdim(1:3,1)) + reshape(ZZ_pindcs{4} + nZZ*(0:ZPdim(4,1)-1),[],1);
ridcs = [ridcs_0; ridcs_x; ridcs_y; ridcs_xy];

ncC = sum(ZPdim(:,1).*nZL_arr);
cidcs = (1:ncC)';

Cmat = sparse(ridcs,cidcs,1,nrC,ncC);
 
% % % With that, we can represent
% % % IZL(x,y) = IIZZ(x,y)^T * Cmat, and thus
% %    Xop*P1op = iIZL(x,y)^T * iHmat * (iIZRop*P1op)      
% %             = iIZL(x,y)^T * iHmat * IZL(x,y) * Hmat * IZRop  
% %             = iIZL(x,y)^T * iHmat * IIZZ(x,y) * Cmat * Hmat * IZRop
% %             = iIZL(x,y)^T * iHmat * IIZZ(x,y) * CHmat * IZRop
% % where
CHmat = Cmat*Hmat;
% % We want to move IIZZ(x,y) through iHmat, aiming to represent
% %
% % iIZL(x,y)^T * iHmat * IIZZ(x,y) = IZZZL(x,y) * iGmat,
% %
% % In particular, we define
% %
% % IZZZL(x,y)^T = 
% %    [ZZZ0(x,y)^T                                         ]                                                               ] 
% %    [            IZZZx(x,y)^T                            ]                                                           ] 
% %    [                         IZZZy(x,y)^T               ]
% %    [                                      IZZZxy(x,y)^T ] 
% % 
% % where e.g. IZZZ0(x,y) = I_{inx} \otimes ZZZx(x,y) with
% % ZZZx(x,y) = unique(ZZ(x,y) \otimes Zx(x)).
% %
% % First, we construct the new monomial bases IZZZL, by simply taking the
% % Kronecker product of the monomials ZZ(x,y) with the monomials iZx(x),
% % iZy(y) and iZxy(x,y), and retaining only the unique monomials:

%ZZZ0_dmat = ZZ_dmat;
nZZZ0 = nZZ;
ZZZ0_pindcs = (1:nZZ)';

ZZZx_dmat = kron(ZZ_dmat,ones(inZx,1)) + [repmat(iZx_dmat,[nZZ,1]),zeros([nZZ*inZx,1])];
[ZZZx_dmat,~,ZZZx_pindcs] = uniquerows_integerTable(ZZZx_dmat);     % ZZZ_dmat(pindcs,:) = ZZZ_dmat_non_unique
nZZZx = size(ZZZx_dmat,1);
ZZZx_xindcs = (ZZZx_dmat(:,2)==0);   % Rows in ZZZx_dmat associated to monomials only in x
%ZZZx = polynomial(speye(nZZZx),ZZZx_dmat,[xx.varname;yy.varname],[nZZZx,1]);

ZZZy_dmat = kron(ZZ_dmat,ones(inZy,1)) + [zeros([nZZ*inZy,1]),repmat(iZy_dmat,[nZZ,1])];
[ZZZy_dmat,~,ZZZy_pindcs] = uniquerows_integerTable(ZZZy_dmat);     % ZZZ_dmat(pindcs,:) = ZZZ_dmat_non_unique
nZZZy = size(ZZZy_dmat,1);
ZZZy_yindcs = (ZZZy_dmat(:,1)==0);
%ZZZy = polynomial(speye(nZZZy),ZZZy_dmat,[xx.varname;yy.varname],[nZZZy,1]);

ZZZxy_dmat = kron(ZZ_dmat,ones(inZxy,1)) + repmat(iZxy_dmat,[nZZ,1]);
[ZZZxy_dmat,~,ZZZxy_pindcs] = uniquerows_integerTable(ZZZxy_dmat);     % ZZZ_dmat(pindcs,:) = ZZZ_dmat_non_unique
nZZZxy = size(ZZZxy_dmat,1);
%ZZZxy = polynomial(speye(nZZZxy),ZZZxy_dmat,[xx.varname;yy.varname],[nZZZxy,1]);

nZZZL_arr = [nZZZ0; nZZZx; nZZZy; nZZZxy];
nIZZZL_arr = nZZZL_arr.*Xdim(:,1);
ZZZL_indcs = {ZZZ0_pindcs; ZZZx_pindcs; ZZZy_pindcs; ZZZxy_pindcs};

% % Given these monomial vectors, we look for a matrix Fmat such that
% %  IZZZL(x,y)^T * Fmat * IZRop = P2op
% % Then, if we can represent 
% %     iIZL(x,y)^T * iHmat * IIZZ(x,y) = IZZZL(x,y)^T * iGmat,
% % we would have to find the values of iGmat such that 
% %     iGmat * CHmat = Fmat


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% First for components mapping to R^n
nr = P2dim(1,1);        nc = P2dim(1,2);            % Number of rows and columns of parameter P2.R00
nZL = nZZZL_arr(1);
nIZL = nIZZZL_arr(1);   nIZR = sum(nIZR_arr);

Fmat_00 = double(P2op.R00);
Fmat_0x_o = zeros(nr,nIZR_arr(2));
Fmat_0x_a = coeff_quad_decomp(P2op.R0x,{},xx.varname,zeros(1,0),Ztt_dmat(ZR_retain_indcs{3},:));
Fmat_0x_b = coeff_quad_decomp(P2op.R0x,{},xx.varname,zeros(1,0),Ztt_dmat(ZR_retain_indcs{4},:));
Fmat_0y_o = zeros(nr,nIZR_arr(5));
Fmat_0y_a = coeff_quad_decomp(P2op.R0y,{},yy.varname,zeros(1,0),Znu_dmat(ZR_retain_indcs{6},:));
Fmat_0y_b = coeff_quad_decomp(P2op.R0y,{},yy.varname,zeros(1,0),Znu_dmat(ZR_retain_indcs{7},:));
Fmat_02_o = zeros(nr,sum(nIZR_arr(8:12)));
Fmat_02_aa = coeff_quad_decomp(P2op.R02,{},[xx.varname;yy.varname],zeros(1,0),Zttnu_dmat(ZR_retain_indcs{13},:));
Fmat_02_ba = coeff_quad_decomp(P2op.R02,{},[xx.varname;yy.varname],zeros(1,0),Zttnu_dmat(ZR_retain_indcs{14},:));
Fmat_02_ab = coeff_quad_decomp(P2op.R02,{},[xx.varname;yy.varname],zeros(1,0),Zttnu_dmat(ZR_retain_indcs{15},:));
Fmat_02_bb = coeff_quad_decomp(P2op.R02,{},[xx.varname;yy.varname],zeros(1,0),Zttnu_dmat(ZR_retain_indcs{16},:));

Fmat_0 = [Fmat_00, ...
          Fmat_0x_o, Fmat_0x_a, Fmat_0x_b, ...
          Fmat_0y_o, Fmat_0y_a, Fmat_0y_b, ...
          Fmat_02_o, Fmat_02_aa, Fmat_02_ba, Fmat_02_ab, Fmat_02_bb];
FFmat_0 = spalloc(nIZL,nIZR,nnz(Fmat_0));
FFmat_0(1:nZL:nIZL,:) = Fmat_0;

% Next for components mapping to L2[x]
% In theory, ZR_retain indices should exclude e.g. variation in nu from
% iZttnu where appropriate
% However, ZZZx_dmat is a degmat in both x and y...
nIZL = nIZZZL_arr(2);
ZZZx_x_dmat = ZZZx_dmat(:,1);       % Extract x-degrees from ZZZx_dmat
ZZZx_x_dmat(~ZZZx_xindcs,:) = -1;   % Assign rows that correspond to nonzero y contributions an unreasonable degree

Fmat_x0 = coeff_quad_decomp(P2op.Rx0,xx.varname,{},ZZZx_x_dmat,zeros(nZR_retain(1),0));
Fmat_xx_o = coeff_quad_decomp(P2op.Rxx{1},xx.varname,{},ZZZx_x_dmat,zeros(nZR_retain(2),0));
Fmat_xx_a = coeff_quad_decomp(P2op.Rxx{2},xx.varname,tt.varname,ZZZx_x_dmat,Ztt_dmat(ZR_retain_indcs{3},:));
Fmat_xx_b = coeff_quad_decomp(P2op.Rxx{3},xx.varname,tt.varname,ZZZx_x_dmat,Ztt_dmat(ZR_retain_indcs{4},:));
Fmat_xy_o = zeros(nIZL,nIZR_arr(5));
Fmat_xy_a = coeff_quad_decomp(P2op.Rxy,xx.varname,yy.varname,ZZZx_x_dmat,Znu_dmat(ZR_retain_indcs{6},:));
Fmat_xy_b = coeff_quad_decomp(P2op.Rxy,xx.varname,yy.varname,ZZZx_x_dmat,Znu_dmat(ZR_retain_indcs{7},:));
Fmat_x2_oo = zeros(nIZL,nIZR_arr(8));
Fmat_x2_ao = zeros(nIZL,nIZR_arr(9));
Fmat_x2_bo = zeros(nIZL,nIZR_arr(10));
Fmat_x2_oa = coeff_quad_decomp(P2op.Rx2{1},xx.varname,yy.varname,ZZZx_x_dmat,Zttnu_dmat(ZR_retain_indcs{11},2));
Fmat_x2_ob = coeff_quad_decomp(P2op.Rx2{1},xx.varname,yy.varname,ZZZx_x_dmat,Zttnu_dmat(ZR_retain_indcs{12},2));
Fmat_x2_aa = coeff_quad_decomp(P2op.Rx2{2},xx.varname,[tt.varname;yy.varname],ZZZx_x_dmat,Zttnu_dmat(ZR_retain_indcs{13},:));
Fmat_x2_ba = coeff_quad_decomp(P2op.Rx2{3},xx.varname,[tt.varname;yy.varname],ZZZx_x_dmat,Zttnu_dmat(ZR_retain_indcs{14},:));
Fmat_x2_ab = coeff_quad_decomp(P2op.Rx2{2},xx.varname,[tt.varname;yy.varname],ZZZx_x_dmat,Zttnu_dmat(ZR_retain_indcs{15},:));
Fmat_x2_bb = coeff_quad_decomp(P2op.Rx2{3},xx.varname,[tt.varname;yy.varname],ZZZx_x_dmat,Zttnu_dmat(ZR_retain_indcs{16},:));

FFmat_x = [Fmat_x0, ...
           Fmat_xx_o, Fmat_xx_a, Fmat_xx_b, ...
           Fmat_xy_o, Fmat_xy_a, Fmat_xy_b, ...
           Fmat_x2_oo, Fmat_x2_ao, Fmat_x2_bo, Fmat_x2_oa, Fmat_x2_ob, Fmat_x2_aa, Fmat_x2_ba, Fmat_x2_ab, Fmat_x2_bb];

% Next for components mapping to L2[y]
nIZL = nIZZZL_arr(3);
ZZZy_y_dmat = ZZZy_dmat(:,2);       % Extract y-degrees from ZZZy_dmat
ZZZy_y_dmat(~ZZZy_yindcs,:) = -1;   % Assign rows that correspond to nonzero x contributions an unreasonable degree

Fmat_y0 = coeff_quad_decomp(P2op.Ry0,yy.varname,{},ZZZy_y_dmat,zeros(nZR_retain(1),0));
Fmat_yx_o = zeros(nIZL,nIZR_arr(2));
Fmat_yx_a = coeff_quad_decomp(P2op.Ryx,yy.varname,xx.varname,ZZZy_y_dmat,Ztt_dmat(ZR_retain_indcs{3},:));
Fmat_yx_b = coeff_quad_decomp(P2op.Ryx,yy.varname,xx.varname,ZZZy_y_dmat,Ztt_dmat(ZR_retain_indcs{4},:));
Fmat_yy_o = coeff_quad_decomp(P2op.Ryy{1},yy.varname,{},ZZZy_y_dmat,zeros(nZR_retain(5),0));
Fmat_yy_a = coeff_quad_decomp(P2op.Ryy{2},yy.varname,nu.varname,ZZZy_y_dmat,Znu_dmat(ZR_retain_indcs{6},:));
Fmat_yy_b = coeff_quad_decomp(P2op.Ryy{3},yy.varname,nu.varname,ZZZy_y_dmat,Znu_dmat(ZR_retain_indcs{7},:));
Fmat_y2_oo = zeros(nIZL,nIZR_arr(8));
Fmat_y2_ao = coeff_quad_decomp(P2op.Ry2{1},yy.varname,xx.varname,ZZZy_y_dmat,Zttnu_dmat(ZR_retain_indcs{9},1));
Fmat_y2_bo = coeff_quad_decomp(P2op.Ry2{1},yy.varname,xx.varname,ZZZy_y_dmat,Zttnu_dmat(ZR_retain_indcs{10},1));
Fmat_y2_oa = zeros(nIZL,nIZR_arr(11));
Fmat_y2_ob = zeros(nIZL,nIZR_arr(12));
Fmat_y2_aa = coeff_quad_decomp(P2op.Ry2{2},yy.varname,[xx.varname;nu.varname],ZZZy_y_dmat,Zttnu_dmat(ZR_retain_indcs{13},:));
Fmat_y2_ba = coeff_quad_decomp(P2op.Ry2{2},yy.varname,[xx.varname;nu.varname],ZZZy_y_dmat,Zttnu_dmat(ZR_retain_indcs{14},:));
Fmat_y2_ab = coeff_quad_decomp(P2op.Ry2{3},yy.varname,[xx.varname;nu.varname],ZZZy_y_dmat,Zttnu_dmat(ZR_retain_indcs{15},:));
Fmat_y2_bb = coeff_quad_decomp(P2op.Ry2{3},yy.varname,[xx.varname;nu.varname],ZZZy_y_dmat,Zttnu_dmat(ZR_retain_indcs{16},:));

FFmat_y = [Fmat_y0, ...
           Fmat_yx_o, Fmat_yx_a, Fmat_yx_b, ...
           Fmat_yy_o, Fmat_yy_a, Fmat_yy_b, ...
           Fmat_y2_oo, Fmat_y2_ao, Fmat_y2_bo, Fmat_y2_oa, Fmat_y2_ob, Fmat_y2_aa, Fmat_y2_ba, Fmat_y2_ab, Fmat_y2_bb];

% Finally for components mapping to L2[x,y]
Fmat_20 = coeff_quad_decomp(P2op.R20,[xx.varname;yy.varname],{},ZZZxy_dmat,zeros(1,0));
Fmat_2x_o = coeff_quad_decomp(P2op.R2x{1},[xx.varname;yy.varname],{},ZZZxy_dmat,zeros(nZR_retain(2),0));
Fmat_2x_a = coeff_quad_decomp(P2op.R2x{2},[xx.varname;yy.varname],tt.varname,ZZZxy_dmat,Ztt_dmat(ZR_retain_indcs{3},:));
Fmat_2x_b = coeff_quad_decomp(P2op.R2x{3},[xx.varname;yy.varname],tt.varname,ZZZxy_dmat,Ztt_dmat(ZR_retain_indcs{4},:));
Fmat_2y_o = coeff_quad_decomp(P2op.R2y{1},[xx.varname;yy.varname],{},ZZZxy_dmat,zeros(nZR_retain(5),0));
Fmat_2y_a = coeff_quad_decomp(P2op.R2y{2},[xx.varname;yy.varname],nu.varname,ZZZxy_dmat,Znu_dmat(ZR_retain_indcs{6},:));
Fmat_2y_b = coeff_quad_decomp(P2op.R2y{3},[xx.varname;yy.varname],nu.varname,ZZZxy_dmat,Znu_dmat(ZR_retain_indcs{7},:));
Fmat_22_oo = coeff_quad_decomp(P2op.R22{1,1},[xx.varname;yy.varname],{},ZZZxy_dmat,zeros(nZR_retain(8),0));
Fmat_22_ao = coeff_quad_decomp(P2op.R22{2,1},[xx.varname;yy.varname],tt.varname,ZZZxy_dmat,Zttnu_dmat(ZR_retain_indcs{9},1));
Fmat_22_bo = coeff_quad_decomp(P2op.R22{3,1},[xx.varname;yy.varname],tt.varname,ZZZxy_dmat,Zttnu_dmat(ZR_retain_indcs{10},1));
Fmat_22_oa = coeff_quad_decomp(P2op.R22{1,2},[xx.varname;yy.varname],nu.varname,ZZZxy_dmat,Zttnu_dmat(ZR_retain_indcs{11},2));
Fmat_22_ob = coeff_quad_decomp(P2op.R22{1,3},[xx.varname;yy.varname],nu.varname,ZZZxy_dmat,Zttnu_dmat(ZR_retain_indcs{12},2));
Fmat_22_aa = coeff_quad_decomp(P2op.R22{2,2},[xx.varname;yy.varname],[tt.varname;nu.varname],ZZZxy_dmat,Zttnu_dmat(ZR_retain_indcs{13},:));
Fmat_22_ba = coeff_quad_decomp(P2op.R22{3,2},[xx.varname;yy.varname],[tt.varname;nu.varname],ZZZxy_dmat,Zttnu_dmat(ZR_retain_indcs{14},:));
Fmat_22_ab = coeff_quad_decomp(P2op.R22{2,3},[xx.varname;yy.varname],[tt.varname;nu.varname],ZZZxy_dmat,Zttnu_dmat(ZR_retain_indcs{15},:));
Fmat_22_bb = coeff_quad_decomp(P2op.R22{3,3},[xx.varname;yy.varname],[tt.varname;nu.varname],ZZZxy_dmat,Zttnu_dmat(ZR_retain_indcs{16},:));

FFmat_2 = [Fmat_20, ...
           Fmat_2x_o, Fmat_2x_a, Fmat_2x_b, ...
           Fmat_2y_o, Fmat_2y_a, Fmat_2y_b, ...
           Fmat_22_oo, Fmat_22_ao, Fmat_22_bo, Fmat_22_oa, Fmat_22_ob, Fmat_22_aa, Fmat_22_ba, Fmat_22_ab, Fmat_22_bb];


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % VALUES OF FMAT FOR INVERSE
% % First for ZZZ0 * Imat0 * ZRop
% zLI_idx = 1;            zRI_idx = 1;            % Index associated to constant monomial 1
% nr = P1dim(1,2);         nc = P1dim(1,2);         % Number of rows and columns of identity parameter P.R00
% nZL = nZZZL_arr(1);      nZR = nZR_retain(1);   % Number of LHS and RHS monomials
% Ivals = ones(nr,1);
% Irindcs = (0:nr-1)*nZL + zLI_idx;
% Icindcs = (0:nc-1)*nZR + zRI_idx;
% Imat_0 = sparse(Irindcs,Icindcs,Ivals,nr*nZL,nc*nZR);
% 
% % Next for ZZZx * Imatx * ZRop
% zLI_idx = find(all(ZZZx_dmat==0,2),1);          % Index associated to constant monomial 1 in ZZZx (should be 1)
% zRI_idx = find(all(Ztt_dmat(ZR_retain_indcs{2},:)==0,2),1); % Index associated to constant monomial 1 in Ztt (should be 1)     
% nr = P1dim(2,2);         nc = P1dim(2,2);         % Number of rows and columns of identity parameter P.Rxx{1}
% nZL = nZZZL_arr(2);      nZR = nZR_retain(2);   % Number of LHS and RHS monomials
% Ivals = ones(nr,1);
% Irindcs = (0:nr-1)*nZL + zLI_idx;
% Icindcs = (0:nc-1)*nZR + zRI_idx;
% Imat_x = sparse(Irindcs,Icindcs,Ivals,nr*nZL,nc*nZR);
% Imat_x = [Imat_x, zeros(size(Imat_x,1),nc*sum(nZR_retain(3:4)))];  % Set contribution of P.Rxx{2} and P.Rxx{3} to zero
% 
% % Next for ZZZy * Imaty * ZRop
% zLI_idx = find(all(ZZZy_dmat==0,2),1);          % Index associated to constant monomial 1 in ZZZy (should be 1)
% zRI_idx = find(all(Znu_dmat(ZR_retain_indcs{5},:)==0,2),1); % Index associated to constant monomial 1 in Znu (should be 1)  
% nr = P1dim(3,2);         nc = P1dim(3,2);         % Number of rows and columns of identity parameter P.Ryy{1}
% nZL = nZZZL_arr(3);      nZR = nZR_retain(5);   % Number of LHS and RHS monomials
% Ivals = ones(nr,1);
% Irindcs = (0:nr-1)*nZL + zLI_idx;
% Icindcs = (0:nc-1)*nZR + zRI_idx;
% Imat_y = sparse(Irindcs,Icindcs,Ivals,nr*nZL,nc*nZR);
% Imat_y = [Imat_y, zeros(size(Imat_y,1),nc*sum(nZR_retain(6:7)))];  % Set contribution of P.Ryy{2} and P.Ryy{3} to zero
% 
% % Finally for ZZZxy * Imaty * ZRop
% zLI_idx = find(all(ZZZxy_dmat==[0,0],2),1);         % Index associated to constant monomial 1 in ZZZxy (should be 1)
% zRI_idx = find(all(Zttnu_dmat(ZR_retain_indcs{8},:)==[0,0],2),1);   % Index associated to constant monomial 1 in Zttnu (should be 1)  
% nr = P1dim(4,2);         nc = P1dim(4,2);             % Number of rows and columns of identity parameter P.R22{1}
% nZL = nZZZL_arr(4);      nZR = nZR_retain(8);       % Number of LHS and RHS monomials
% Ivals = ones(nr,1);
% Irindcs = (0:nr-1)*nZL + zLI_idx;
% Icindcs = (0:nc-1)*nZR + zRI_idx;
% Imat_xy = sparse(Irindcs,Icindcs,Ivals,nr*nZL,nc*nZR);
% Imat_xy = [Imat_xy, zeros(size(Imat_xy,1),nc*sum(nZR_retain(9:16)))];  % Set contribution of P.R22{i,j} for i~=1 or j~=1 to zero

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Then, we want to find iGmat such that iGmat * CHmat = Imat
% % We do not have a very "nice" expression for iGmat, but we do know it
% % has a particular structure:
% %
% % iGmat = 
% % [ iG00 0      iG0x   iG0x   0      iG0y   iG0y   0       0       0       0       0       iG02    iG02    iG02    iG02   ] * 
% % [ iGx0 iGxx_o iGxx_a iGxx_b 0      iGxy   iGxy   0       0       0       iGx2_o  iGx2_o  iGx2_a  iGx2_b  iGx2_a  iGx2_b ]
% % [ iGy0 0      iGyx   iGyx   iGyy_o iGyy_a iGyy_b 0       iGy2_o  iGy2_o  0       0       iGy2_a  iGy2_a  iGy2_b  iGy2_b ]
% % [ iG20 iG2x_o iG2x_a iG2x_b iG2y_o iG2y_a iG2y_b iG22_oo iG22_ao iG22_bo iG22_oa iG22_ob iG22_aa iG22_ba iG22_ab iG22_bb]
% %
% % To account for this structure, we decompose the matrix as
% %
% % iGmat
% % = [ iG00 iG0x iG0y iG02
% %   [                      iGx0 iGxx_o iGxx_a iGxx_b iGxy iGx2_o iGx2_a iGx2_b   
% %   [                                                                           iGy0 iGyx iGyy_o iGyy_a iGyy_b iGy2_o iGy2_a iGy2_b   
% %   [
% % *
% %   [ I0 
% %   [ 0   0    Ix Ix
% %   [ 0   0    0  0   0    Iy Iy
% %   [ 0   0    0  0   0    0  0   0     0      0      0      0      Ixy Ixy Ixy  Ixy ]
% %   [ I0  0
% %   [ 0   Ix_0
% %   [ 0   0    Ix
% %   [ 0   0    0  Ix
% %   [ 0   0    0  0   0    Iy Iy
% %   [ 0   0    0  0   0    0  0   0     0      0      Ixy_nu Ixy_nu                 ]
% %   [ 0   0    0  0   0    0  0   0     0      0      0      0      Ixy 0   Ixy 0   ]
% %   [ 0   0    0  0   0    0  0   0     0      0      0      0      0   Ixy 0   Ixy ]
% %   [ I0  0
% %   [ 0   0    Ix Ix
% %   [ 0   0    0  0   Iy_0
% %   [ 0   0    0  0   0    Iy
% %   [ 0   0    0  0   0    0  Iy
% %   [ 0   0    0  0   0    0  0   0     Ixy_tt Ixy_tt 0      0                      ]
% %   [ 0   0    0  0   0    0  0   0     0      0      0      0      Ixy Ixy 0   0   ]
% %   [ 0   0    0  0   0    0  0   0     0      0      0      0      0   0   Ixy Ixy ]
% %   [ I0
% %   [ 0   Ix_0
% %   [ 0   0    Ix
% %   [ 0   0    0  Ix
% %   [ 0   0    0  0   Iy
% %   [ 0   0    0  0   0    Iy
% %   [ 0   0    0  0   0    0  Iy
% %   [ 0   0    0  0   0    0  0   Ixy_0
% %   [ 0   0    0  0   0    0  0   0     Ixy_tt
% %   [ 0   0    0  0   0    0  0   0     0      Ixy_tt
% %   [ 0   0    0  0   0    0  0   0     0      0      Ixy_nu
% %   [ 0   0    0  0   0    0  0   0     0      0      0      Ixy_nu
% %   [ 0   0    0  0   0    0  0   0     0      0      0      0      Ixy
% %   [ 0   0    0  0   0    0  0   0     0      0      0      0      0  Ixy
% %   [ 0   0    0  0   0    0  0   0     0      0      0      0      0  0  Ixy
% %   [ 0   0    0  0   0    0  0   0     0      0      0      0      0  0  0  Ixy ]
% %
% % = [ iG_blocks_0                                      ] * [ iG_pmat_0  ]
% %   [             iG_blocks_x                          ]   [ iG_pmat_x  ]
% %   [                         iG_blocks_y              ]   [ iG_pmat_y  ]
% %   [                                     iG_blocks_xy ]   [ iG_pmat_xy ]
% %
% % = iG_blocks * iG_pmat
% %
% % We construct the matrices iG_pmat_, by first building the individual
% % identity matrices
I0 = speye(Xdim(1,2)*nZZ);             
Ix = speye(Xdim(2,2)*inZR_arr(2)*nZZ);
Iy = speye(Xdim(3,2)*inZR_arr(3)*nZZ);
Ixy = speye(Xdim(4,2)*inZR_arr(4)*nZZ);

% Ix_0
Ix_0 = speye(Xdim(2,2)*nZZ);
Ox_0 = zeros(size(Ix,1),size(Ix_0,2));

% Iy_0
Iy_0 = speye(Xdim(3,2)*nZZ);
Oy_0 = zeros(size(Iy,1),size(Iy_0,2));

% Ixy_0
Ixy_0 = speye(Xdim(4,2)*nZZ);
Oxy_02 = zeros(size(Ixy,1),size(Ixy_0,2));

% Ixy_nu
Ixy_nu = speye(Xdim(4,2)*inZttnu_nu*nZZ);
Oxy_nu_02 = zeros(size(Ixy,1),size(Ixy_nu,2));

% Ixy_tt
Ixy_tt = speye(Xdim(4,2)*inZttnu_tt*nZZ);
Oxy_tt_02 = zeros(size(Ixy,1),size(Ixy_tt,2));

Oxy_0_y2 = zeros(size(Ixy_tt,1),size(Ixy_0,2));
Oxy_0_x2 = zeros(size(Ixy_nu,1),size(Ixy_0,2));
Oxy_nu_y2 = zeros(size(Ixy_tt,1),size(Ixy_nu,2));
Oxy_tt_x2 = zeros(size(Ixy_nu,1),size(Ixy_tt,2));

% % Then, combine into the matirces iG_pmat_
iG_pmat_0 = blkdiag(I0,[Ox_0,Ix,Ix],[Oy_0,Iy,Iy],[[Oxy_02,Oxy_tt_02,Oxy_tt_02,Oxy_nu_02,Oxy_nu_02],repmat(Ixy,[1,4])]);

% Account for potential separability of kernels
if sep_iRxx
    Ixx_block = blkdiag(Ix_0,[Ix,Ix]);
else
    Ixx_block = blkdiag(Ix_0,Ix,Ix);
end
if sep_iRx2
    Ix2_a_block = [Ixy,Ixy,Ixy,Ixy];
else
    Ix2_a_block = [blkdiag(Ixy,Ixy),blkdiag(Ixy,Ixy)];
end
iG_pmat_x = blkdiag(I0,Ixx_block,[Oy_0,Iy,Iy],...
                        blkdiag([Oxy_0_x2,Oxy_tt_x2,Oxy_tt_x2,Ixy_nu,Ixy_nu],Ix2_a_block));

% Account for potential separability of kernels
if sep_iRyy
    Iyy_block = blkdiag(Iy_0,[Iy,Iy]);
else
    Iyy_block = blkdiag(Iy_0,Iy,Iy);
end
if sep_iRy2
    Iy2_a_block = [Ixy,Ixy,Ixy,Ixy];
else
    Iy2_a_block = blkdiag([Ixy,Ixy],[Ixy,Ixy]);
end
iG_pmat_y = blkdiag(I0,[Ox_0,Ix,Ix],Iyy_block,...
                        blkdiag([Oxy_0_y2,Ixy_tt,Ixy_tt,Oxy_nu_y2,Oxy_nu_y2],Iy2_a_block));

% Account for potential separability of kernels
if sep_iR2x
    I2x_block = blkdiag(Ix_0,[Ix,Ix]);
else
    I2x_block = blkdiag(Ix_0,Ix,Ix);
end
if sep_iR2y
    I2y_block = blkdiag(Iy_0,[Iy,Iy]);
else
    I2y_block = blkdiag(Iy_0,Iy,Iy);
end
if sep_iR22_ao
    I22_ao_block = [Ixy_tt,Ixy_tt];
else
    I22_ao_block = blkdiag(Ixy_tt,Ixy_tt);
end
if sep_iR22_oa
    I22_oa_block = [Ixy_nu,Ixy_nu];
else
    I22_oa_block = blkdiag(Ixy_nu,Ixy_nu);
end
if sep_iR22_aa_x && sep_iR22_aa_y
    I22_aa_block = [Ixy,Ixy,Ixy,Ixy];
elseif sep_iR22_aa_x
    I22_aa_block = blkdiag([Ixy,Ixy],[Ixy,Ixy]);
elseif sep_iR22_aa_y
    I22_aa_block = [blkdiag(Ixy,Ixy),blkdiag(Ixy,Ixy)];
else
    I22_aa_block = blkdiag(Ixy,Ixy,Ixy,Ixy);
end
iG_pmat_xy = blkdiag(I0,I2x_block,I2y_block,...
                        Ixy_0,I22_ao_block,I22_oa_block,I22_aa_block);

% % % Then we have
% % Pop_inv*Pop = iIZL(x,y)^T * iHmat * (iIZRop * Pop)
% %             = iIZL(x,y)^T * iHmat * IZZ(x,y) * CHmat * IZRop,
% %             = IZZZL(x,y)^T * iGmat * CHmat * IZRop,
% %             = IZZZL(x,y)^T * iG_blocks * iG_pmat * CHmat * IZRop
% %             = IZZZL(x,y)^T * iG_blocks * PCHmat * IZRop,
% % where
% % PCHmat = iG_pmat * CHmat;
% % We want iG_blocks to be such that iG_blocks * PCHmat = Imat.
% % Here, we note that 
% %
% % iG_blocks
% % = [ iG_blocks_0                                      ]
% %   [             iG_blocks_x                          ]
% %   [                         iG_blocks_y              ]
% %   [                                     iG_blocks_xy ]
% %
% % Exploiting this block-diagonal structure, the first block iG_pmat_0
% % will operate only on iG_pmat_0 * CHmat, defining the first ... rows
% % of the desired identity matrix. In particular, we decompose Imat as
% %     Imat = [ IImat_0;  IImat_x;  IImat_y;  IImat_xy ];
% % where
% rdim_Imat = [size(Imat_0,1);size(Imat_x,1);size(Imat_y,1);size(Imat_xy,1)];
% cdim_Imat = [size(Imat_0,2),size(Imat_x,2),size(Imat_y,2),size(Imat_xy,2)];
% FFmat_0 = [Imat_0,zeros(rdim_Imat(1),sum(cdim_Imat(2:4)))];
% FFmat_x = [zeros(rdim_Imat(2),cdim_Imat(1)), Imat_x, zeros(rdim_Imat(2),sum(cdim_Imat(3:4)))];
% FFmat_y = [zeros(rdim_Imat(3),sum(cdim_Imat(1:2))), Imat_y, zeros(rdim_Imat(3),cdim_Imat(4))];
% FFmat_2 = [zeros(rdim_Imat(4),sum(cdim_Imat(1:3))), Imat_xy];
% % Then we want 
% %     iG_blocks_.. * PCHmat_.. = IImat_..;
% % where
PCHmat_0 = iG_pmat_0 * CHmat;
PCHmat_x = iG_pmat_x * CHmat;
PCHmat_y = iG_pmat_y * CHmat;
PCHmat_xy = iG_pmat_xy * CHmat;

% % Here we can get rid of any columns that are guaranteed and desired to
% % be zero in the product iGmat_blocks_.. * PCHmat_..
cretain_indcs_0 = any(PCHmat_0,1) | any(FFmat_0,1);
cretain_indcs_x = any(PCHmat_x,1) | any(FFmat_x,1);
cretain_indcs_y = any(PCHmat_y,1) | any(FFmat_y,1);
cretain_indcs_xy = any(PCHmat_xy,1) | any(FFmat_2,1);

PCHmat_0 = PCHmat_0(:,cretain_indcs_0);     FFmat_0 = FFmat_0(:,cretain_indcs_0);
PCHmat_x = PCHmat_x(:,cretain_indcs_x);     FFmat_x = FFmat_x(:,cretain_indcs_x);
PCHmat_y = PCHmat_y(:,cretain_indcs_y);     FFmat_y = FFmat_y(:,cretain_indcs_y);
PCHmat_xy = PCHmat_xy(:,cretain_indcs_xy);  FFmat_2 = FFmat_2(:,cretain_indcs_xy);


% % This leaves the big challenge: establishing parameters iHmat, such that
% % iG_blocks * PCHmat = Imat, where iGmat = iG_blocks * iG_pmat is such
% % that iIZL(x,y)^T * iHmat * IIZZ(x,y) = IZZZL(x,y)^T * iGmat...
% % We note here that, if e.g. Rxx{1} is scalar, we can represent
% %
% % iZx(x)^T * iHxx^o * (I_{inZtt} \otimes ZZ(x,y))^T
% % = (ZZ(x,y)^T \otimes iZx(x)^T) * 
% %   [ I_{nZZ} \otimes iHxx^o(:,1), ..., I_{nZZ} \otimes iHxx^o(:,inZtt) ]
% %
% % = (ZZ(x,y)^T \otimes iZx(x)^T) * 
% %   (I_{inZx*nZZ} \otimes iHxx^o(:)^T) * 
% %     [ [e_1       0  ...  0      ]     [ 0        0  ...  0      ] ] }                  }  
% %     [ [ :        :   .   :      ]     [ :        .   :   :      ] ] }} inZtt*inZx rows }  
% %     [ [ 0        0  ...  0      ]     [e_1       0  ...  0      ] ] }                  }  
% %     [ [e_2       0  ...  0      ] ... [ 0        0  ...  0      ] ]                    }  
% %     [ [ :        :   .   :      ]     [ :        :   .   :      ] ]                    } } inZx*inZtt*inZx rows  
% %     [ [e_{inZx}  0  ...  0      ]     [ 0        0  ...  0      ] ]                    }  
% %     [ [ :        :   .   :      ]     [ :        :   .   :      ] ]                    }  
% %     [ [ 0        0  ...  0      ]     [e_{inZx}  0  ...  0      ] ]                    }  
% %     [                                                             ]
% %     [ [ 0       e_1 ...  0      ]     [ 0        0  ...  0      ] ]
% %     [ [ :        :   .   :      ] ... [ :        :   .   :      ] ]
% %     [ [ 0        0  ...  0      ]     [ 0       e_1 ...  0      ] ]
% %     [ [ :        :   .   :      ]     [ :        :   .   :      ] ]
% %     [                                                             ]
% %     [            :                 .             :                ]
% %     [                                                             ]
% %     [ [ :        :   .   :      ]     [ :        :   .   :      ] ]
% %     [ [ 0        0  ... e_{inZx}]     [ 0        0  ...  0      ] ]
% %     [ [ :        :   .   :      ] ... [ :        :   .   :      ] ]
% %     [ [ 0        0  ...  0      ]     [ 0        0  ... e_{inZx}] ]
% %
% % where e_j is the jth standard basis function in \R^{inZx}
% % Each block has nZZ columns, and there are inZtt blocks
% % We note that this decomposition is valid only for the non-unique
% % combination of monomials (ZZ(x,y)^T \otimes iZx(x)^T). Retaining only
% % the unique monomials, 

% % Let's start with 0 components iG_blocks_0
%inIZL_arr = iPdim(:,1)'.*inZL_arr;
inIZR_arr = Xdim(:,2)'.*inZR_arr;  % Number of columns in iH associated to R00, R0x, R0y and R02
ncols_Xop_0 = inIZR_arr;

inr = Xdim(1,1);           % Number of rows to which R00 etc. map
nZL = inZL_arr(1);          % Number of monomials in LHS monomial of inverse
nZZZL = nZZZL_arr(1);       % Number of monomials in ZZZ_0(x,y)
ncols_Xop_tot = sum(ncols_Xop_0);   % Total number of columns in iH associated to R00, R0x, R0y and R02
zLidcs = ZZZL_indcs{1};     % Position of monomials from ZZ(x,y)\otimes Z0 in the unique vector of monomials ZZZ_0(x,y)

nvals_iGvec = nZL * ncols_Xop_tot;          % Number of values in Gvec_0 = [iH00(1:nZL,1)',...,iH00(1:nZL,end)',iH0x(1:nZL,1)',...,iH02(1:nZL,end)']
ncE = nZZ * ncols_Xop_tot;                  % Number of columns in Emat
Ecidcs = kron((1:ncE)',ones([nZL,1]));      % Column indices in Emat associated to each coefficient
iGvec_idcs = reshape((repmat((1:nZL)',[nZZ,1]) + nZL*(0:ncols_Xop_tot-1)),[],1);   % Element in Gvec associated to each coefficient

% For practical purposes, we decompose Emat into
% nZZZL blocks of nvals_Gvec rows
EPCHmat_0 = [];
EIImat_0 = [];
for k = 1:nZZZL
    zidcs_k = find(zLidcs==k);
    rc_idcs_k = reshape(zidcs_k + length(zLidcs)*(0:ncols_Xop_tot-1),[],1);
    Emat_0_k = sparse(iGvec_idcs(rc_idcs_k), Ecidcs(rc_idcs_k),1,nvals_iGvec,ncE);
    
    EPCHmat_0_k = Emat_0_k * PCHmat_0;
    EPCHmat_0 = [EPCHmat_0, EPCHmat_0_k];
    EIImat_0 = [EIImat_0, FFmat_0(k+nZZZL*(0:inr-1)',:)];
end

% Then, perform the matrix inverse...
iGvec_block_0 = EIImat_0 / EPCHmat_0;
%eps_0 = max(max(abs(iGvec_block_0 * EPCHmat_0 - EIImat_0))/min(P1dim(1,:))); % Estimate of error in the result
eps_0 = max(max(abs(iGvec_block_0 * EPCHmat_0 - EIImat_0)));
iGvec_cell_0 = mat2cell(iGvec_block_0,ones(inr,1));
for k = 1:inr
    iGvec_cell_0{k} = reshape(iGvec_cell_0{k},nZL,[]);
end
Xop_params_0 = cell2mat(iGvec_cell_0);


% % Next for the rows mapping to L_2[x]
% Number of columns in iHp associated to Rx0, Rxx^o, Rxx^a, Rxx^b, Rxy, 
%       Rx2^o, Rx2^a and Rx2^b
ncols_Xop_x = [inIZR_arr(1), Xdim(2,2), inIZR_arr(2), inIZR_arr(2), inIZR_arr(3), ...
                                inZttnu_nu*Xdim(4,2), inIZR_arr(4), inIZR_arr(4)];
if sep_iRxx
    ncols_Xop_x(4) = 0;
end
if sep_iRx2
    ncols_Xop_x(8) = 0;
end

inr = Xdim(2,1);           % Number of rows to which Rx0 etc. map
nZL = inZL_arr(2);          % Number of monomials in LHS monomial of inverse
nZZZL = nZZZL_arr(2);       % Number of monomials in ZZZ_0(x,y)
ncols_Xop_tot = sum(ncols_Xop_x);   % Total number of columns in iH associated to Rx0, Rxx, Rxy and Rx2
zLidcs = ZZZL_indcs{2};     % Position of monomials from ZZ(x,y)\otimes Z0 in the unique vector of monomials ZZZ_x(x,y)

nvals_iGvec = nZL * ncols_Xop_tot;  % Number of values in Gvec_0 = [iHx0(1:nZL,1)',...,iHx0(1:nZL,end)',iHxx^o(1:nZL,1)',...,iHx2^b(1:nZL,end)']
ncE = nZZ * ncols_Xop_tot;          % Number of columns in Emat

Ecidcs = kron((1:ncE)',ones([nZL,1]));              % Column indices in Emat associated to each coefficient
iGvec_idcs = reshape((repmat((1:nZL)',[nZZ,1]) + nZL*(0:ncols_Xop_tot-1)),[],1);   % Element in Gvec associated to each coefficient

% For practical purposes, we decompose Emat into
% nZZZL blocks of nvals_Gvec rows
EPCHmat_x = [];
EIImat_x = [];
for k = 1:nZZZL
    zidcs_k = find(zLidcs==k);
    rc_idcs_k = reshape(zidcs_k + length(zLidcs)*(0:ncols_Xop_tot-1),[],1);
    Emat_x_k = sparse(iGvec_idcs(rc_idcs_k), Ecidcs(rc_idcs_k),1,nvals_iGvec,ncE);
    
    EPCHmat_x_k = Emat_x_k * PCHmat_x;
    EPCHmat_x = [EPCHmat_x, EPCHmat_x_k];
    EIImat_x = [EIImat_x, FFmat_x(k+nZZZL*(0:inr-1)',:)];
end

% Then, perform the matrix inverse...
iGvec_block_x = EIImat_x / EPCHmat_x;
%eps_x = max(max(abs(iGvec_block_x * EPCHmat_x - EIImat_x))/min(P1dim(2,:))); % Estimate of error in the result
eps_x = max(max(abs(iGvec_block_x * EPCHmat_x - EIImat_x)));
iGvec_cell_x = mat2cell(iGvec_block_x,ones(inr,1));
for k = 1:inr
    iGvec_cell_x{k} = reshape(iGvec_cell_x{k},nZL,[]);
end
Xop_params_x = cell2mat(iGvec_cell_x);


% % Next for the rows mapping to L_2[y]
% Number of columns in iHp associated to Ry0, Ryx, Ryy^o, Ryy^a, Ryy^b, 
%       Ry2^o, Ry2^a and Ry2^b
ncols_Xop_y = [inIZR_arr(1), inIZR_arr(2), Xdim(3,2), inIZR_arr(3), inIZR_arr(3), ...
                                inZttnu_tt*Xdim(4,2), inIZR_arr(4), inIZR_arr(4)];
if sep_iRyy
    ncols_Xop_y(5) = 0;
end
if sep_iRy2
    ncols_Xop_y(8) = 0;
end
inr = Xdim(3,1);           % Number of rows to which Rx0 etc. map
nZL = inZL_arr(3);          % Number of monomials in LHS monomial of inverse
nZZZL = nZZZL_arr(3);       % Number of monomials in ZZZ_y(x,y)
ncols_Xop_tot = sum(ncols_Xop_y);   % Total number of columns in iH associated to Ry0, Ryx, Ryy and Ry2
zLidcs = ZZZL_indcs{3};     % Position of monomials from ZZ(x,y)\otimes Zy in the unique vector of monomials ZZZ_y(x,y)

nvals_iGvec = nZL * ncols_Xop_tot;          % Number of values in Gvec_0 = [iHy0(1:nZL,1)',...,iHy0(1:nZL,end)',iHyx(1:nZL,1)',...,iHy2^b(1:nZL,end)']
ncE = nZZ * ncols_Xop_tot;                  % Number of columns in Emat
Ecidcs = kron((1:ncE)',ones([nZL,1]));      % Column indices in Emat associated to each coefficient
iGvec_idcs = reshape((repmat((1:nZL)',[nZZ,1]) + nZL*(0:ncols_Xop_tot-1)),[],1);   % Element in Gvec associated to each coefficient

% For practical purposes, we decompose Emat into
% nZZZL blocks of nvals_Gvec rows
EPCHmat_y = [];
EIImat_y = [];
for k = 1:nZZZL
    zidcs_k = find(zLidcs==k);
    rc_idcs_k = reshape(zidcs_k + length(zLidcs)*(0:ncols_Xop_tot-1),[],1);
    Emat_y_k = sparse(iGvec_idcs(rc_idcs_k), Ecidcs(rc_idcs_k),1,nvals_iGvec,ncE);
    
    EPCHmat_y_k = Emat_y_k * PCHmat_y;
    EPCHmat_y = [EPCHmat_y, EPCHmat_y_k];
    EIImat_y = [EIImat_y, FFmat_y(k+nZZZL*(0:inr-1)',:)];
end

% Then, perform the matrix inverse...
iGvec_block_y = EIImat_y / EPCHmat_y;
%eps_y = max(max(abs(iGvec_block_y * EPCHmat_y - EIImat_y))/min(P1dim(3,:))); % Estimate of error in the result
eps_y = max(max(abs(iGvec_block_y * EPCHmat_y - EIImat_y)));
iGvec_cell_y = mat2cell(iGvec_block_y,ones(inr,1));
for k = 1:inr
    iGvec_cell_y{k} = reshape(iGvec_cell_y{k},nZL,[]);
end
Xop_params_y = cell2mat(iGvec_cell_y);


% % Finally, for the rows mapping to L_2[x,y]
% Number of columns in iHp associated to R20, R2x^o, R2x^a, R2x^b, Ryy^o, Ryy^a, Ryy^b, 
%       R22^oo, R22^a0, R22^bo, R22^oa, R22^ob, R22^aa, R22^ba, R22^ab and R22^bb
ncols_Xop_xy = [inIZR_arr(1), Xdim(2,2), inIZR_arr(2), inIZR_arr(2), Xdim(3,2), inIZR_arr(3), inIZR_arr(3), ...
                                Xdim(4,2), inZttnu_tt*Xdim(4,2), inZttnu_tt*Xdim(4,2),...
                                inZttnu_nu*Xdim(4,2), inZttnu_nu*Xdim(4,2), inIZR_arr(4)*ones(1,4)];
if sep_iR2x
    ncols_Xop_xy(4) = 0;
end
if sep_iR2y
    ncols_Xop_xy(7) = 0;
end
if sep_iR22_ao
    ncols_Xop_xy(10) = 0;
end
if sep_iR22_oa
    ncols_Xop_xy(12) = 0;
end
if sep_iR22_aa_x
    ncols_Xop_xy([14,16]) = [0,0];
end
if sep_iR22_aa_y
    ncols_Xop_xy([15,16]) = [0,0];
end
inr = Xdim(4,1);           % Number of rows to which Rx0 etc. map
nZL = inZL_arr(4);          % Number of monomials in LHS monomial of inverse
nZZZL = nZZZL_arr(4);       % Number of monomials in ZZZ_xy(x,y)
ncols_Xop_tot = sum(ncols_Xop_xy);   % Total number of columns in iH associated to Ry0, Ryx, Ryy and Ry2
zLidcs = ZZZL_indcs{4};     % Position of monomials from ZZ(x,y)\otimes Zxy in the unique vector of monomials ZZZ_y(x,y)

nvals_iGvec = nZL * ncols_Xop_tot;          % Number of values in Gvec_0 = [iHy0(1:nZL,1)',...,iHy0(1:nZL,end)',iHyx(1:nZL,1)',...,iHy2^b(1:nZL,end)']
ncE = nZZ * ncols_Xop_tot;                  % Number of columns in Emat
Ecidcs = kron((1:ncE)',ones([nZL,1]));      % Column indices in Emat associated to each coefficient
iGvec_idcs = reshape((repmat((1:nZL)',[nZZ,1]) + nZL*(0:ncols_Xop_tot-1)),[],1);   % Element in Gvec associated to each coefficient

% For practical purposes, we decompose Emat into
% nZZZL blocks of nvals_Gvec rows
EPCHmat_xy = [];
EIImat_xy = [];
for k = 1:nZZZL
    zidcs_k = find(zLidcs==k);
    rc_idcs_k = reshape(zidcs_k + length(zLidcs)*(0:ncols_Xop_tot-1),[],1);
    Emat_xy_k = sparse(iGvec_idcs(rc_idcs_k), Ecidcs(rc_idcs_k),1,nvals_iGvec,ncE);
    
    %test = Emat_xy((k-1)*nvals_iGvec + (1:nvals_iGvec)',:) - Emat_xy_k
    EPCHmat_xy_k = Emat_xy_k * PCHmat_xy;
    EPCHmat_xy = [EPCHmat_xy, EPCHmat_xy_k];
    EIImat_xy = [EIImat_xy, FFmat_2(k+nZZZL*(0:inr-1)',:)];
end

% Then, perform the matrix inverse...
iGvec_block_xy = EIImat_xy / EPCHmat_xy;
%eps_xy = max(max(abs(iGvec_block_xy * EPCHmat_xy - EIImat_xy))/min(P1dim(4,:))); % Estimate of error in the result
eps_xy = max(max(abs(iGvec_block_xy * EPCHmat_xy - EIImat_xy)));
iGvec_cell_xy = mat2cell(iGvec_block_xy,ones(inr,1));
for k = 1:inr
    iGvec_cell_xy{k} = reshape(iGvec_cell_xy{k},nZL,[]);
end
Xop_params_xy = cell2mat(iGvec_cell_xy);


% % % If the obtained parameters do not produce a sufficiently accurate
% % % inverse, try again with more monomials
eps_arr = [eps_0;eps_x;eps_y;eps_xy];
eps = max(eps_arr);
if any(eps>tol)
    if isempty(iter_num_mrdivide_opvar2d)
        iter_num_mrdivide_opvar2d = 1;
    else
        iter_num_mrdivide_opvar2d = iter_num_mrdivide_opvar2d + 1;
    end
    if any(any(deg_fctr==deg_fctr_max)) || iter_num_mrdivide_opvar2d>=10
        fprintf(['\n An accurate inverse (up to specified tolerance) could not be obtained with a degree increase factor of \n'])
        display(deg_fctr)
        fprintf([' Returning the current best guess of the inverse.\n '])
        deg_fctr_final = deg_fctr;
    else
        deg_fctr_new = min(round(deg_fctr.*2),deg_fctr_max);
        [Xop,eps,deg_fctr_final] = mrdivide(P2op,P1op,deg_fctr_new,tol,deg_fctr_max);
    return
    end
else
    deg_fctr_final = deg_fctr;
end


% % % Finally, we can extract the different parameter coefficients from
% % % iH_params, and construct the inverse operator
Xop = opvar2d([], Xdim, P1op.I, var1, var2);

% Determine which columns in iH_params correspond to which parameter
nncols_Xop_0 = cumsum([0,ncols_Xop_0]);
nncols_Xop_x = cumsum([0,ncols_Xop_x]);
nncols_Xop_y = cumsum([0,ncols_Xop_y]);
nncols_Xop_xy = cumsum([0,ncols_Xop_xy]);

% Construct the LHS and RHS monomials of Pop_inv
iIZT0 = polynomial(speye(Xdim(1,1)));
iIZTx = kron(eye(Xdim(2,1)), polynomial(speye(inZx),iZx_dmat,xx.varname,[1,inZx]));
iIZTy = kron(eye(Xdim(3,1)), polynomial(speye(inZy),iZy_dmat,yy.varname,[1,inZy]));
iIZTxy = kron(eye(Xdim(4,1)), polynomial(speye(inZxy),iZxy_dmat,[xx.varname;yy.varname],[1,inZxy]));

iIZ0 = polynomial(speye(Xdim(1,2)));
iIZtt_0 = kron(eye(Xdim(2,2)), polynomial(1,0,tt.varname,[1,1]));
iIZtt = kron(eye(Xdim(2,2)), polynomial(speye(inZtt),iZtt_dmat,tt.varname,[inZtt,1]));
iIZnu_0 = kron(eye(Xdim(3,2)), polynomial(1,0,nu.varname,[1,1]));
iIZnu = kron(eye(Xdim(3,2)), polynomial(speye(inZnu),iZnu_dmat,nu.varname,[inZnu,1]));
iIZttnu_0 = kron(eye(Xdim(4,2)), polynomial(1,[0,0],[tt.varname;nu.varname],[1,1]));
iIZttnu_nu = kron(eye(Xdim(4,2)), polynomial(speye(inZttnu_nu),iZttnu_nu_dmat,nu.varname,[inZttnu_nu,1]));
iIZttnu_tt = kron(eye(Xdim(4,2)), polynomial(speye(inZttnu_tt),iZttnu_tt_dmat,tt.varname,[inZttnu_tt,1]));
iIZttnu = kron(eye(Xdim(4,2)), polynomial(speye(inZttnu),iZttnu_dmat,[tt.varname;nu.varname],[inZttnu,1]));

iIZL_cell = {iIZT0; iIZTx; iIZTy; iIZTxy};
iIZR_cell = {iIZ0, {iIZtt_0,iIZtt}, {iIZnu_0,iIZnu},...
              {iIZttnu_0,iIZttnu_tt,iIZttnu_nu,iIZttnu}};

% % Multiply the parameter matrices with the appropriate monomial vectors
% R0
if ~isempty(Xop_params_0)
Xop.R00 = iIZL_cell{1} * Xop_params_0(:, nncols_Xop_0(1)+1:nncols_Xop_0(2)) * iIZR_cell{1};
Xop.R0x = iIZL_cell{1} * Xop_params_0(:, nncols_Xop_0(2)+1:nncols_Xop_0(3)) * subs(iIZR_cell{2}{2},tt,xx);
Xop.R0y = iIZL_cell{1} * Xop_params_0(:, nncols_Xop_0(3)+1:nncols_Xop_0(4)) * subs(iIZR_cell{3}{2},nu,yy);
Xop.R02 = iIZL_cell{1} * Xop_params_0(:, nncols_Xop_0(4)+1:nncols_Xop_0(5)) * subs(iIZR_cell{4}{4},[tt;nu],[xx;yy]);
end

% Rx
if ~isempty(Xop_params_x)
Xop.Rx0    = iIZL_cell{2} * Xop_params_x(:, nncols_Xop_x(1)+1:nncols_Xop_x(2)) * iIZR_cell{1};
Xop.Rxx{1} = iIZL_cell{2} * Xop_params_x(:, nncols_Xop_x(2)+1:nncols_Xop_x(3)) * iIZR_cell{2}{1};
Xop.Rxx{2} = iIZL_cell{2} * Xop_params_x(:, nncols_Xop_x(3)+1:nncols_Xop_x(4)) * iIZR_cell{2}{2};
if sep_iRxx
    Xop.Rxx{3} = Xop.Rxx{2};
else
    Xop.Rxx{3} = iIZL_cell{2} * Xop_params_x(:, nncols_Xop_x(4)+1:nncols_Xop_x(5)) * iIZR_cell{2}{2};
end
Xop.Rxy =    iIZL_cell{2} * Xop_params_x(:, nncols_Xop_x(5)+1:nncols_Xop_x(6)) * subs(iIZR_cell{3}{2},nu,yy);
Xop.Rx2{1} = iIZL_cell{2} * Xop_params_x(:, nncols_Xop_x(6)+1:nncols_Xop_x(7)) * subs(iIZR_cell{4}{3},nu,yy);
Xop.Rx2{2} = iIZL_cell{2} * Xop_params_x(:, nncols_Xop_x(7)+1:nncols_Xop_x(8)) * subs(iIZR_cell{4}{4},nu,yy);
if sep_iRx2
    Xop.Rx2{3} = Xop.Rx2{2};
else
Xop.Rx2{3} = iIZL_cell{2} * Xop_params_x(:, nncols_Xop_x(8)+1:nncols_Xop_x(9)) * subs(iIZR_cell{4}{4},nu,yy);
end
end

% Ry
if ~isempty(Xop_params_y)
Xop.Ry0 =    iIZL_cell{3} * Xop_params_y(:, nncols_Xop_y(1)+1:nncols_Xop_y(2)) * iIZR_cell{1};
Xop.Ryx =    iIZL_cell{3} * Xop_params_y(:, nncols_Xop_y(2)+1:nncols_Xop_y(3)) * subs(iIZR_cell{2}{2},tt,xx);
Xop.Ryy{1} = iIZL_cell{3} * Xop_params_y(:, nncols_Xop_y(3)+1:nncols_Xop_y(4)) * iIZR_cell{3}{1};
Xop.Ryy{2} = iIZL_cell{3} * Xop_params_y(:, nncols_Xop_y(4)+1:nncols_Xop_y(5)) * iIZR_cell{3}{2};
if sep_iRyy
    Xop.Ryy{3} = Xop.Ryy{2};
else
Xop.Ryy{3} = iIZL_cell{3} * Xop_params_y(:, nncols_Xop_y(5)+1:nncols_Xop_y(6)) * iIZR_cell{3}{2};
end
Xop.Ry2{1} = iIZL_cell{3} * Xop_params_y(:, nncols_Xop_y(6)+1:nncols_Xop_y(7)) * subs(iIZR_cell{4}{2},tt,xx);
Xop.Ry2{2} = iIZL_cell{3} * Xop_params_y(:, nncols_Xop_y(7)+1:nncols_Xop_y(8)) * subs(iIZR_cell{4}{4},tt,xx);
if sep_iRy2
    Xop.Ry2{3} = Xop.Ry2{2};
else
Xop.Ry2{3} = iIZL_cell{3} * Xop_params_y(:, nncols_Xop_y(8)+1:nncols_Xop_y(9)) * subs(iIZR_cell{4}{4},tt,xx);
end
end

% R2
if ~isempty(Xop_params_xy)
Xop.R20 =    iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(1)+1:nncols_Xop_xy(2)) * iIZR_cell{1};
Xop.R2x{1} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(2)+1:nncols_Xop_xy(3)) * iIZR_cell{2}{1};
Xop.R2x{2} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(3)+1:nncols_Xop_xy(4)) * iIZR_cell{2}{2};
if sep_iR2x
    Xop.R2x{3} = Xop.R2x{2};
else
Xop.R2x{3} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(4)+1:nncols_Xop_xy(5)) * iIZR_cell{2}{2};
end
Xop.R2y{1} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(5)+1:nncols_Xop_xy(6)) * iIZR_cell{3}{1};
Xop.R2y{2} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(6)+1:nncols_Xop_xy(7)) * iIZR_cell{3}{2};
if sep_iR2y
    Xop.R2y{3} = Xop.R2y{2};
else
Xop.R2y{3} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(7)+1:nncols_Xop_xy(8)) * iIZR_cell{3}{2};
end
Xop.R22{1,1} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(8)+1:nncols_Xop_xy(9)) * iIZR_cell{4}{1};
Xop.R22{2,1} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(9)+1:nncols_Xop_xy(10)) * iIZR_cell{4}{2};
if sep_iR22_ao
    Xop.R22{3,1} = Xop.R22{2,1};
else
Xop.R22{3,1} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(10)+1:nncols_Xop_xy(11)) * iIZR_cell{4}{2};
end
Xop.R22{1,2} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(11)+1:nncols_Xop_xy(12)) * iIZR_cell{4}{3};
if sep_iR22_oa
    Xop.R22{1,3} = Xop.R22{1,2};
else
Xop.R22{1,3} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(12)+1:nncols_Xop_xy(13)) * iIZR_cell{4}{3};
end
Xop.R22{2,2} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(13)+1:nncols_Xop_xy(14)) * iIZR_cell{4}{4};
if sep_iR22_aa_x
    Xop.R22{3,2} = Xop.R22{2,2};
else
Xop.R22{3,2} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(14)+1:nncols_Xop_xy(15)) * iIZR_cell{4}{4};
end
if sep_iR22_aa_y
    Xop.R22{2,3} = Xop.R22{2,2};
else
Xop.R22{2,3} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(15)+1:nncols_Xop_xy(16)) * iIZR_cell{4}{4};
end
if sep_iR22_aa_x
    Xop.R22{3,3} = Xop.R22{2,3};
elseif sep_iR22_aa_y
    Xop.R22{3,3} = Xop.R22{3,2};
else
Xop.R22{3,3} = iIZL_cell{4} * Xop_params_xy(:, nncols_Xop_xy(16)+1:nncols_Xop_xy(17)) * iIZR_cell{4}{4};
end
end

% Convert dummy variables to primary variables in parameters above diagonal
% Pop_inv.R0x.varname(ismember(Pop_inv.R0x.varname,tt.varname)) = xx.varname;
% Pop_inv.Ryx.varname(ismember(Pop_inv.Ryx.varname,tt.varname)) = xx.varname;
% Pop_inv.R02.varname(ismember(Pop_inv.R02.varname,tt.varname)) = xx.varname;
% 
% Pop_inv.Ry2{1}.varname(ismember(Pop_inv.Ry2{1}.varname,tt.varname)) = xx.varname;
% Pop_inv.Ry2{2}.varname(ismember(Pop_inv.Ry2{2}.varname,tt.varname)) = xx.varname;
% Pop_inv.Ry2{3}.varname(ismember(Pop_inv.Ry2{3}.varname,tt.varname)) = xx.varname;
% 
% Pop_inv.R0y.varname(ismember(Pop_inv.R0y.varname,nu.varname)) = yy.varname;
% Pop_inv.Rxy.varname(ismember(Pop_inv.Rxy.varname,nu.varname)) = yy.varname;
% Pop_inv.R02.varname(ismember(Pop_inv.R02.varname,nu.varname)) = yy.varname;
% 
% Pop_inv.Rx2{1}.varname(ismember(Pop_inv.Rx2{1}.varname,nu.varname)) = yy.varname;
% Pop_inv.Rx2{2}.varname(ismember(Pop_inv.Rx2{2}.varname,nu.varname)) = yy.varname;
% Pop_inv.Rx2{3}.varname(ismember(Pop_inv.Rx2{3}.varname,nu.varname)) = yy.varname;

% Get rid of terms below tolerance
clean_tol = 1e-15;
Xop = clean_opvar(Xop,clean_tol);

end



function [Z_dmat, p_indcs] = collect_monomials(Pop,params,var)

if ischar(var)
    var = {var};
elseif ispvar(var)
    var = var.varname;
end
nvars = length(var);

% Combine all the degree matrices in the desired variable into a single
% matrix
Z_dmat = zeros(0,nvars);
nZ_arr = [];
n_params = ones(numel(params),2);
for k=1:numel(params)
    PR = Pop.(params{k});
    
    if isempty(PR) || isa(PR,'double')
%         nZ_arr = [nZ_arr; 0];
%     elseif isa(PR,'double')
        Z_dmat = [Z_dmat; zeros(1,nvars)];
        nZ_arr = [nZ_arr; 1];
    elseif isa(PR,'polynomial')
        [isvar,Loc] = ismember(PR.varname,var);
        %[cindx,Loc] = ismember(PR.varname,var);
        if ~any(isvar)
            Z_dmat = [Z_dmat; zeros(size(PR.degmat,1),nvars)];
            nZ_arr = [nZ_arr; size(PR.degmat,1)];
        else
            Zk_dmat = zeros(size(PR.degmat,1),nvars);
            Zk_dmat(:,Loc(isvar)) = PR.degmat(:,isvar);
            Z_dmat = [Z_dmat; Zk_dmat];
            nZ_arr = [nZ_arr; size(PR.degmat,1)];
        end
    elseif isa(PR,'cell')
        PP = PR;
        n_params(k,:) = size(PP);
        for l = 1:numel(PP)
            PR = PP{l};
            if isempty(PR) || isa(PR,'double')
%                 nZ_arr = [nZ_arr; 0];
%             elseif isa(PR,'double')
                Z_dmat = [Z_dmat; zeros(1,nvars)];
                nZ_arr = [nZ_arr; 1];
            elseif isa(PR,'polynomial')
                [isvar,Loc] = ismember(PR.varname,var);
                if ~any(isvar)
                    Z_dmat = [Z_dmat; zeros(size(PR.degmat,1),nvars)];
                    nZ_arr = [nZ_arr; size(PR.degmat,1)];
                else
                    Zk_dmat = zeros(size(PR.degmat,1),nvars);
                    Zk_dmat(:,Loc(isvar)) = PR.degmat(:,isvar);
                    Z_dmat = [Z_dmat; Zk_dmat];
                    nZ_arr = [nZ_arr; size(PR.degmat,1)];
                end
            else
                error('Parameters in opvar2d object should be specified as ''double'' or ''polynomial'' class objects');
            end
        end
    else
        error('Parameters in opvar2d object should be specified as ''double'' or ''polynomial'' class objects');
    end
end

% Establish a unique set of degrees
[Z_dmat,~,indx2] = uniquerows_integerTable(Z_dmat);     % Z_dmat(indx2,:) = Z_dmat_non_unique

% % For each parameter k, also establish indcs p_indcs{k} in the unique
% % degmat that correspond to the original degmat of this parameter
nnZx = [0;cumsum(nZ_arr)];

p_indcs = cell(size(params));
dum = 0;
for k=1:numel(params)
    if any(n_params(k,:)>1)
        p_indcs_k = cell(n_params(k,:));
        for l=1:numel(p_indcs_k)
            dum = dum+1;
            p_indcs_k{l} = indx2(nnZx(dum)+1 : nnZx(dum+1));
        end
    else
        dum = dum+1;
        p_indcs_k = indx2(nnZx(dum)+1 : nnZx(dum+1));
    end
    p_indcs{k} = p_indcs_k;
end 
    


end



function [A_unique,indx1,indx2] = uniquerows_integerTable(A)
% Establish unique rows in the integer array A.
% Produces the same result as 
%   [A_unique,indx1,indx2] = unique(A,'rows');
% but uses sortrows_integerTable.

[Asort, indx1] = sortrows_integerTable(A);  % Asort = A(indx1,:);
[~,indx2] = sort(indx1);                    % A = Asort(indx2,:);
     

retain_indcs = [any(1); any(Asort(1:end-1,:) - Asort(2:end,:),2)];
A_unique = Asort(retain_indcs,:);
indx1 = indx1(retain_indcs);            % A_unique = A(indx1,:);

retain_indcs2 = cumsum(retain_indcs);   % Asort = A_unique(retain_indcs2,:);
indx2 = retain_indcs2(indx2);           % A = A_unique(indx2,:);

end



function [Cmat] = coeff_quad_decomp(P,varL,varR,ZL_degs,ZR_degs)

% ZL_poly = polynomial(eye(nZL),ZL_degs,varL,[nZL,1])
% ZR_poly = polynomial(eye(nZR),ZR_degs,varR,[nZR,1])
% P = kron(eye(nr),ZL_poly)'*Cmat*kron(eye(nc),ZR_poly)

% if isempty(ZL_degs)
%     ZL_degs = zeros(1,length(varL));
% end
% if isempty(ZR_degs)
%     ZR_degs = zeros(1,length(varR));
% end

if isempty(P)
    Cmat = zeros(size(P,1)*size(ZL_degs,1),size(P,2)*size(ZR_degs,1));
    return
end
if size(ZR_degs,1)==0
    Cmat = zeros(size(P,1)*size(ZL_degs,1),0);
    return
end


if isa(P,'double')
    P = polynomial(P);
end

% Extract the inputs
degmat = P.degmat;
coeff = P.coeff;
varname = P.varname;
[nr,nc] = size(P);

% If certain do not appear, add them with constant contribution (deg=0)
if length(varname) < length(varL)+length(varR)
    degmat = [degmat,zeros(size(degmat,1),(length(varL)+length(varR)) - length(varname))];
    varname = [varname; varL(~ismember(varL,varname)); varR(~ismember(varR,varname))];
end

%if (size(ZL_degs,2)+size(ZR_degs,2))~=size(degmat,2)
%    error('The polynomial includes more or fewer variables than that appear in the proposed monomial basis.')
%end

if ispvar(varL)
    varL = varL.varname;
end
if ispvar(varR)
    varR = varR.varname;
end

ndegs = size(degmat,1);
nZL = size(ZL_degs,1);      nZR = size(ZR_degs,1);

% For each coefficient, determine which monomial in ZL_degs it is
% associated to.
if nZL==0
    ZL_indcs = ones(ndegs,nr*nc);
    nZL = 1;
    logval_L = true;
else
    [~,varL_indcs] = ismember(varL,varname);
    degmat_L = degmat(:,varL_indcs);
    [logval_L,ZL_indcs] = ismember(degmat_L,ZL_degs,'rows');
    ZL_indcs = repmat(ZL_indcs,[1,nr*nc]);
end
% For each coefficient, determine which monomial in ZR_degs it is
% associated to.
if nZR==0
    ZR_indcs = ones(ndegs,nr*nc);
    nZR = 1;
    logval_R = true;
else
    [~,varR_indcs] = ismember(varR,varname);
    degmat_R = degmat(:,varR_indcs);
    [logval_R,ZR_indcs] = ismember(degmat_R,ZR_degs,'rows');
    ZR_indcs = repmat(ZR_indcs,[1,nr*nc]);
end

if any(~logval_L) || any(~logval_R)
    error('The polynomial includes at least one monomial that does not appear in the proposed monomial basis.')
end

% For each coefficient, establish a row and column number in the
% matrix-valued object.
row_indcs = repmat(1:nr,[ndegs,nc]);
col_indcs = repmat(kron((1:nc),ones(1,nr)),[ndegs,1]);

% Determine row and column indices for each coefficient in the new
% coefficient matrix
Cr_indcs = (row_indcs-1)*nZL + ZL_indcs;
Cc_indcs = (col_indcs-1)*nZR + ZR_indcs;

% Build the new matrix C
Cmat = sparse(Cr_indcs(:),Cc_indcs(:),coeff(:),nr*nZL,nc*nZR);

end