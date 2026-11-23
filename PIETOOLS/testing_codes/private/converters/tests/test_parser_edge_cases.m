function results = test_parser_edge_cases(groups)
% TEST_PARSER_EDGE_CASES Edge cases of the PDE parser and of the PIE
% structures: integrals, products, sums and differences, substitution,
% derivatives, delays, outputs, BCs, constraints and inputs, concatenation
% of PIEs, and comparison of PIEs.
%
% Each 'eq' case declares the same system in two independent ways and
% checks that the converted PIEs are equal; for nonlinear PIEs, for which
% == need not identify x*x_s with x_s*x, the right-hand sides f are also
% compared by evaluation at random polynomial states. Each 'err' case
% declares an unsupported system, and passes if the parser or converter
% throws an error of its own (not an internal MATLAB error such as a
% dimension mismatch). Each 'ok' case must convert.
% The 1D systems live on [a,b] = [-1,2], with BCs x_s(a)=0 and x(b)=0, so
% that neither the domain nor the boundary values are standard.
%
% Run as: results = test_parser_edge_cases;  or for some groups only,
%         results = test_parser_edge_cases({'I','S'});
%
% INPUT
% - groups: (optional) cellstr, the groups of cases to run, from 'I'
%           (integrals), 'M' (products), 'P' (sums and differences), 'S'
%           (substitution), 'D' (derivatives), 'E' (delays), 'B' (outputs,
%           BCs, constraints, inputs), 'C' (concatenation of PIEs, PIEs
%           built from operators), 'N' (comparison of nonlinear PIEs); all
%           by default;
%
% OUTPUT
% - results: 'struct' array with fields name, passed and message, one
%            element per case; for an 'err' case, message is the error.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_parser_edge_cases
%
% Copyright (C) 2026 PIETOOLS Team
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding DJ - 11/23/2026

if nargin<1
    groups = {'I','M','P','S','D','E','B','C','N'};
end
pvar s t s1 s2 th
a = -1;     b = 2;
R = cell(0,1);

if ismember('I',groups)
fprintf('\n===== int =====\n');
R{end+1} = eq_case('I1 [a,s] then [a,b]','y',@(V) int(int(V.x,s,[a,s]),s,[a,b]),@(V) int((b-s)*V.x,s,[a,b]),a,b);
R{end+1} = eq_case('I2 [s,b] then [a,b]','y',@(V) int(int(V.x,s,[s,b]),s,[a,b]),@(V) int((s-a)*V.x,s,[a,b]),a,b);
R{end+1} = eq_case('I3 [a,s] then [s,b]','y',@(V) int(int(V.x,s,[a,s]),s,[s,b]),@(V) (b-s)*int(V.x,s,[a,s])+int((b-s)*V.x,s,[s,b]),a,b);
R{end+1} = eq_case('I3b [s,b] then [a,s]','y',@(V) int(int(V.x,s,[s,b]),s,[a,s]),@(V) int((s-a)*V.x,s,[a,s])+(s-a)*int(V.x,s,[s,b]),a,b);
R{end+1} = eq_case('I4 triple [a,s],[a,s],[a,b]','y',@(V) int(int(int(V.x,s,[a,s]),s,[a,s]),s,[a,b]),@(V) int((b-s)^2/2*V.x,s,[a,b]),a,b);
R{end+1} = eq_case('I4b triple [a,s],[s,b],[a,b]','y',@(V) int(int(int(V.x,s,[a,s]),s,[s,b]),s,[a,b]),@(V) int(((b-s)^2/2+(b-s)*(s-a))*V.x,s,[a,b]),a,b);
R{end+1} = eq_case('I5 vector state, matrix kernel','y2',@(V) int(int([1 s;0 2]*V.x2,s,[a,s]),s,[a,b]),@(V) int((b-s)*[1 s;0 2]*V.x2,s,[a,b]),a,b);
R{end+1} = eq_case('I6 vector state, 1x2 kernel','y21',@(V) int(int([1 s]*V.x2,s,[a,s]),s,[a,b]),@(V) int((b-s)*[1 s]*V.x2,s,[a,b]),a,b);
R{end+1} = eq_case('I7 sum of terms inside','y',@(V) int(int(V.x+s*V.x,s,[a,s]),s,[a,b]),@(V) int((b-s)*(1+s)*V.x,s,[a,b]),a,b);
R{end+1} = eq_case('I8 full + partial inside','y',@(V) int(int(V.x,s,[a,b])+int(V.x,s,[a,s]),s,[a,b]),@(V) int((2*b-a-s)*V.x,s,[a,b]),a,b);
R{end+1} = eq_case('I9 int of x_s over [a,b]','X',@(V) int(diff(V.x,s),s,[a,b]),@(V) subs(V.x,s,b)-subs(V.x,s,a),a,b);
R{end+1} = eq_case('I10 int of x_s over [a,s]','y',@(V) int(diff(V.x,s),s,[a,s]),@(V) V.x-subs(V.x,s,a),a,b);
R{end+1} = eq_case('I11 int of x(a) over [a,b]','y',@(V) int(subs(V.x,s,a),s,[a,b]),@(V) (b-a)*subs(V.x,s,a),a,b);
R{end+1} = eq_case('I11b int of x(a) over [a,s]','y',@(V) int(subs(V.x,s,a),s,[a,s]),@(V) (s-a)*subs(V.x,s,a),a,b);
R{end+1} = eq_case('I12 int of ODE state over s','X',@(V) int(V.X,s,[a,b]),@(V) (b-a)*V.X,a,b);
R{end+1} = eq_case('I17 reversed [b,a]','X',@(V) int(V.x,s,[b,a]),@(V) -1*int(V.x,s,[a,b]),a,b);
R{end+1} = eq_case('I18 reversed [s,a], kernel s','y',@(V) int(s*V.x,s,[s,a]),@(V) -1*int(s*V.x,s,[a,s]),a,b);
R{end+1} = eq_case('I18b inner [b,s], outer [a,s]','y',@(V) int(int(V.x,s,[b,s]),s,[a,s]),@(V) -1*(int((s-a)*V.x,s,[a,s])+(s-a)*int(V.x,s,[s,b])),a,b);
R{end+1} = eq_case('I19 input: [a,s] then [a,b], ODE eq','Xw',@(V) int(int(V.w,s,[a,s]),s,[a,b]),@(V) int((b-s)*V.w,s,[a,b]),a,b);
R{end+1} = eq_case('I20 input: reversed [b,s]','yw',@(V) int(V.w,s,[b,s]),@(V) -1*int(V.w,s,[s,b]),a,b);
R{end+1} = eq_case('I21 input: [a,s] then [s,b]','yw',@(V) int(int(V.w,s,[a,s]),s,[s,b]),@(V) (b-s)*int(V.w,s,[a,s])+int((b-s)*V.w,s,[s,b]),a,b);
R{end+1} = eq_case('I22 input: [a,s] then [a,b], PDE eq','yw',@(V) int(int(V.w,s,[a,s]),s,[a,b]),@(V) int((b-th)*subs(V.w,s,th),th,[a,b]),a,b);
R{end+1} = eq_case('I23 input: int(s*w) in PDE eq = s*int(w)','yw',@(V) int(s*V.w,s,[a,b]),@(V) s*int(V.w,s,[a,b]),a,b);
R{end+1} = eq_case('I24 input: int(s*w) in ODE eq','Xw',@(V) int(s*V.w,s,[a,b]),@(V) int(int(V.w,s,[s,b]),s,[a,b])+a*int(V.w,s,[a,b]),a,b);
R{end+1} = eq_case('I26 input: int(s*w,[a,s])','yw',@(V) int(s*V.w,s,[a,s]),@(V) s*int(V.w,s,[a,s])-int(int(V.w,s,[a,s]),s,[a,s]),a,b);
R{end+1} = eq_case('I13 2D int over [s1;s2] = nested','2DX',@(V) int(V.x,[s1;s2],[0,1;0,1]),@(V) int(int(V.x,s1,[0,1]),s2,[0,1]),0,1);
R{end+1} = eq_case('I14 2D partial in s1, full in s2','2Dy',@(V) int(V.x,[s1;s2],[0,s1;0,1]),@(V) int(int(V.x,s1,[0,s1]),s2,[0,1]),0,1);
R{end+1} = eq_case('I15 2D order [s2;s1]','2DX',@(V) int(V.x,[s2;s1],[0,1;0,1]),@(V) int(V.x,[s1;s2],[0,1;0,1]),0,1);
R{end+1} = eq_case('I16 2D partial in s1 and s2, either order','2Dyy',@(V) int(int(V.x,s1,[0,s1]),s2,[0,s2]),@(V) int(int(V.x,s2,[0,s2]),s1,[0,s1]),0,1);
R{end+1} = eq_case('I16b 2D [0,s1] then [s1,1]','2Dyy',@(V) int(int(V.x,s1,[0,s1]),s1,[s1,1]),@(V) (1-s1)*int(V.x,s1,[0,s1])+int((1-s1)*V.x,s1,[s1,1]),0,1);
R{end+1} = eq_case('I16c 2D nested, then over s2','2DX',@(V) int(int(int(V.x,s1,[0,s1]),s1,[0,1]),s2,[0,1]),@(V) int((1-s1)*V.x,[s1;s2],[0,1;0,1]),0,1);
R{end+1} = eq_case('I16d 2D split in s1, then over s2 in one call','2Dy',@(V) int(int(V.x,s1,[0,s1]),[s1;s2],[s1,1;0,1]),@(V) (1-s1)*int(int(V.x,s1,[0,s1]),s2,[0,1])+int(int((1-s1)*V.x,s1,[s1,1]),s2,[0,1]),0,1);
R{end+1} = err_case('I30 limits inside the domain, [a,0]','X',@(V) int(V.x,s,[a,0]),a,b);
R{end+1} = err_case('I31 integration over t','X',@(V) int(V.x,t,[0,1]),a,b);
R{end+1} = err_case('I32 integration of a nonlinear term','X',@(V) int(V.x*V.x,s,[a,b]),a,b);
R{end+1} = err_case('I33 partial integral in an ODE equation','X',@(V) int(V.x,s,[a,s]),a,b);
R{end+1} = err_case('I34 nested, outer limits [a,0]','X',@(V) int(int(V.x,s,[a,s]),s,[a,0]),a,b);
R{end+1} = err_case('I35 nested, outer limit not matching','y',@(V) int(int(V.x,s,[a,s]),s,[0,b]),a,b);
R{end+1} = err_case('I25 partial integral of input in ODE equation','Xw',@(V) int(V.w,s,[a,s]),a,b);
end

if ismember('M',groups)
fprintf('\n===== mtimes =====\n');
R{end+1} = eq_case('M1 x*2 = 2*x','y',@(V) V.x*2,@(V) 2*V.x,a,b);
R{end+1} = eq_case('M1b x*s = s*x','y',@(V) V.x*s,@(V) s*V.x,a,b);
R{end+1} = eq_case('M3 2*I*(x2.*x2) = 2*(x2.*x2)','y2',@(V) [2 0;0 2]*(V.x2.*V.x2),@(V) 2*(V.x2.*V.x2),a,b);
R{end+1} = eq_case('M5 (x*x)*s = s*(x*x)','y',@(V) (V.x*V.x)*s,@(V) s*(V.x*V.x),a,b);
R{end+1} = eq_case('M6 s*x*x = s*(x*x)','y',@(V) s*V.x*V.x,@(V) s*(V.x*V.x),a,b);
R{end+1} = eq_case('M9 x*x*x = x^3','y',@(V) V.x*V.x*V.x,@(V) V.x^3,a,b);
R{end+1} = ok_case('M7 x*int(x) (nonlinear, integral factor)','y',@(V) V.x*int(V.x,s,[a,b]),a,b);
R{end+1} = eq_case('M4 diag(2,3)*(x2.*x2)','y2',@(V) [2 0;0 3]*(V.x2.*V.x2),@(V) ([2 0;0 3]*V.x2).*V.x2,a,b);
R{end+1} = eq_case('M8 [1;s]*(x*x)','y12',@(V) [1;s]*(V.x*V.x),@(V) ([1;s]*V.x).*([1;1]*V.x),a,b);
R{end+1} = eq_case('M11 [1 2;3 4]*(x2.*x2)','y2',@(V) [1 2;3 4]*(V.x2.*V.x2),...
    @(V) ([1;3]*([1 0]*V.x2)).*([1;1]*([1 0]*V.x2))+([2;4]*([0 1]*V.x2)).*([1;1]*([0 1]*V.x2)),a,b);
R{end+1} = eq_case('M12 [1 s]*(x2.*x2)','y21',@(V) [1 s]*(V.x2.*V.x2),@(V) ([1 0]*V.x2).*([1 0]*V.x2)+([0 s]*V.x2).*([0 1]*V.x2),a,b);
R{end+1} = eq_case('M13 [1;s]*(x*x*x)','y12',@(V) [1;s]*(V.x*V.x*V.x),@(V) ([1;s]*V.x).*([1;1]*V.x).*([1;1]*V.x),a,b);
R{end+1} = err_case('M2 [1 2]*x for scalar x','y',@(V) [1 2]*V.x,a,b);
R{end+1} = err_case('M10 x2*[1 2] (matrix on the right)','y2',@(V) V.x2*[1 2],a,b);
end

if ismember('P',groups)
fprintf('\n===== plus / minus =====\n');
R{end+1} = eq_eqn_case('P1 (y_t==y_ss)-s*x',@(V) (diff(V.y,t)==diff(V.y,s,2))-s*V.x,@(V) diff(V.y,t)==diff(V.y,s,2)-s*V.x,a,b);
R{end+1} = eq_eqn_case('P2 s*x-(y_t==y_ss) = y_t==y_ss+s*x',@(V) s*V.x-(diff(V.y,t)==diff(V.y,s,2)),@(V) diff(V.y,t)==diff(V.y,s,2)+s*V.x,a,b);
R{end+1} = eq_eqn_case('P3 s*x+(y_t==y_ss) = y_t==y_ss-s*x',@(V) s*V.x+(diff(V.y,t)==diff(V.y,s,2)),@(V) diff(V.y,t)==diff(V.y,s,2)-s*V.x,a,b);
R{end+1} = eq_eqn_case('P4 (y_t==y_ss+x)-x = y_t==y_ss',@(V) (diff(V.y,t)==diff(V.y,s,2)+V.x)-V.x,@(V) diff(V.y,t)==diff(V.y,s,2),a,b);
R{end+1} = eq_eqn_case('P5 -(y_t==y_ss) = the same equation',@(V) -(diff(V.y,t)==diff(V.y,s,2)),@(V) diff(V.y,t)==diff(V.y,s,2),a,b);
R{end+1} = eq_eqn_case('P6 y_ss==y_t (time derivative on the right)',@(V) diff(V.y,s,2)==diff(V.y,t),@(V) diff(V.y,t)==diff(V.y,s,2),a,b);
R{end+1} = eq_eqn_case('P7 0+(y_t==y_ss)+0',@(V) 0+(diff(V.y,t)==diff(V.y,s,2))+0,@(V) diff(V.y,t)==diff(V.y,s,2),a,b);
R{end+1} = err_eqn_case('P8 equation + equation',@(V) (diff(V.y,t)==diff(V.y,s,2))+(diff(V.y,t)==V.y),a,b);
R{end+1} = err_eqn_case('P9 x + 1',@(V) diff(V.y,t)==diff(V.y,s,2)+V.x+1,a,b);
R{end+1} = err_eqn_case('P10 two time derivatives in one equation',@(V) diff(V.y,t)+diff(V.x,t)==diff(V.y,s,2),a,b);
R{end+1} = err_eqn_case('P11 x2 (2 rows) + (y_t==y_ss) (1 row)',@(V) V.x2+(diff(V.y,t)==diff(V.y,s,2)),a,b);
end

if ismember('S',groups)
fprintf('\n===== subs =====\n');
R{end+1} = eq_case('S2 subs(s*x,s,b) = b*x(b)','X',@(V) subs(s*V.x,s,b),@(V) b*subs(V.x,s,b),a,b);
R{end+1} = eq_case('S5 subs(int(x,[a,s]),s,b) = int(x,[a,b])','X',@(V) subs(int(V.x,s,[a,s]),s,b),@(V) int(V.x,s,[a,b]),a,b);
R{end+1} = eq_case('S6 subs(int(x,[s,b]),s,a) = int(x,[a,b])','X',@(V) subs(int(V.x,s,[s,b]),s,a),@(V) int(V.x,s,[a,b]),a,b);
R{end+1} = eq_case('S7 subs(int(s*x,[a,s]),s,b)','X',@(V) subs(int(s*V.x,s,[a,s]),s,b),@(V) int(s*V.x,s,[a,b]),a,b);
R{end+1} = eq_case('S8 subs(s*int(x,[a,s]),s,b)','X',@(V) subs(s*int(V.x,s,[a,s]),s,b),@(V) b*int(V.x,s,[a,b]),a,b);
R{end+1} = eq_case('S9 subs(int(x,[a,s]),s,a) = 0','X',@(V) subs(int(V.x,s,[a,s]),s,a),@(V) 0,a,b);
R{end+1} = eq_case('S14 subs(int(int(x,[a,s]),[a,s]),s,b)','X',@(V) subs(int(int(V.x,s,[a,s]),s,[a,s]),s,b),@(V) int((b-s)*V.x,s,[a,b]),a,b);
R{end+1} = eq_case('S11 delay of a partial integral','y',@(V) subs(int(V.x,s,[a,s]),t,t-1),@(V) int(subs(V.x,t,t-1),s,[a,s]),a,b);
R{end+1} = eq_case('S12 delay of a full integral','X',@(V) subs(int(V.x,s,[a,b]),t,t-1),@(V) int(subs(V.x,t,t-1),s,[a,b]),a,b);
R{end+1} = eq_case('S15 delay of a partial integral with kernel','y',@(V) int(s*subs(V.x,t,t-1),s,[a,s]),@(V) subs(int(s*V.x,s,[a,s]),t,t-1),a,b);
R{end+1} = err_case('S1 subs(x,s,0.5) (interior point)','y',@(V) subs(V.x,s,0.5),a,b);
R{end+1} = err_case('S4 subs(x,t,t+1) (time advance)','y',@(V) subs(V.x,t,t+1),a,b);
R{end+1} = err_case('S9b subs(int(x,[a,s]),s,0.5)','X',@(V) subs(int(V.x,s,[a,s]),s,0.5),a,b);
R{end+1} = err_case('S10 subs(x,s,s2)','y',@(V) subs(V.x,s,s2),a,b);
end

if ismember('D',groups)
fprintf('\n===== diff =====\n');
R{end+1} = eq_case('D1 d/ds int(x,[a,s]) = x','y',@(V) diff(int(V.x,s,[a,s]),s),@(V) V.x,a,b);
R{end+1} = eq_case('D2 d/ds int(x,[s,b]) = -x','y',@(V) diff(int(V.x,s,[s,b]),s),@(V) -1*V.x,a,b);
R{end+1} = eq_case('D3 d/ds int(x,[a,b]) = 0','y',@(V) diff(int(V.x,s,[a,b]),s),@(V) 0,a,b);
R{end+1} = eq_case('D4 d/ds s*int(x,[a,s])','y',@(V) diff(s*int(V.x,s,[a,s]),s),@(V) int(V.x,s,[a,s])+s*V.x,a,b);
R{end+1} = eq_case('D5 d/ds int(s*x,[a,s]) = s*x','y',@(V) diff(int(s*V.x,s,[a,s]),s),@(V) s*V.x,a,b);
R{end+1} = eq_case('D6 d/ds (s*x)','y',@(V) diff(s*V.x,s),@(V) V.x+s*diff(V.x,s),a,b);
R{end+1} = eq_case('D7 d/ds x(a) = 0','y',@(V) diff(subs(V.x,s,a),s),@(V) 0,a,b);
R{end+1} = eq_case('D8 d/ds int(int(x,[a,s]),[a,s])','y',@(V) diff(int(int(V.x,s,[a,s]),s,[a,s]),s),@(V) int(V.x,s,[a,s]),a,b);
R{end+1} = eq_case('D11 d2/ds2 int(x,[a,s]) = x_s','y',@(V) diff(int(V.x,s,[a,s]),s,2),@(V) diff(V.x,s),a,b);
R{end+1} = eq_case('D12 d/ds (x*x) = 2 x x_s','y',@(V) diff(V.x*V.x,s),@(V) 2*(V.x*diff(V.x,s)),a,b);
R{end+1} = eq_case('D13 d/ds int(x,[b,s]) = x','y',@(V) diff(int(V.x,s,[b,s]),s),@(V) V.x,a,b);
R{end+1} = eq_case('D14 d/ds int((s-th)x(th),[a,s])','y',@(V) diff(int((s-th)*subs(V.x,s,th),th,[a,s]),s),@(V) int(V.x,s,[a,s]),a,b);
R{end+1} = err_case('D9 d/ds w (input)','yw',@(V) diff(V.w,s),a,b);
R{end+1} = err_case('D10 x_t on the right of another equation','y',@(V) diff(V.x,t),a,b);
end

if ismember('E',groups)
fprintf('\n===== delays =====\n');
R{end+1} = ok_case('E1 delayed x_s(b)','X',@(V) subs(subs(diff(V.x,s),s,b),t,t-1),a,b);
R{end+1} = ok_case('E2 delayed input in an integral','Xw',@(V) int(subs(V.w,t,t-1),s,[a,b]),a,b);
R{end+1} = ok_case('E3 two different delays','X',@(V) subs(subs(V.x,s,a),t,t-1)+subs(subs(V.x,s,a),t,t-2),a,b);
R{end+1} = eq_case('E4 the same delay twice','X',@(V) subs(subs(V.x,s,a),t,t-1)+subs(subs(V.x,s,a),t,t-1),@(V) 2*subs(subs(V.x,s,a),t,t-1),a,b);
R{end+1} = ok_case('E7 delayed x_s in an integral','X',@(V) int(subs(diff(V.x,s),t,t-1),s,[a,b]),a,b);
R{end+1} = ok_case('E8 delayed ODE state','X',@(V) subs(V.X,t,t-1),a,b);
R{end+1} = ok_case('E9 delayed partial integral of an input','yw',@(V) int(subs(V.w,t,t-1),s,[a,s]),a,b);
R{end+1} = case_E5(a,b);
end

if ismember('B',groups)
fprintf('\n===== outputs, BCs, constraints, inputs =====\n');
R{end+1} = case_B1(a,b);
R{end+1} = case_B2(a,b);
R{end+1} = case_B3(a,b);
R{end+1} = case_K1(a,b);
R{end+1} = err_case('Q1 coefficient in a non-spatial variable p','y',@(V) polynomial({'p'})*V.x,a,b);
R{end+1} = eq_case('V1 vector input, [1 s] kernel, [a,s]','yw2',@(V) int([1 s]*V.w2,s,[a,s]),@(V) int([1 0]*V.w2,s,[a,s])+int([0 s]*V.w2,s,[a,s]),a,b);
R{end+1} = eq_case('V2 vector input, nested, opposite sides','yw2',@(V) int(int([1 s]*V.w2,s,[a,s]),s,[s,b]),@(V) (b-s)*int([1 s]*V.w2,s,[a,s])+int((b-s)*[1 s]*V.w2,s,[s,b]),a,b);
end

if ismember('C',groups)
fprintf('\n===== concatenation / piess =====\n');
R{end+1} = case_C1();
R{end+1} = case_C3();
end

if ismember('N',groups)
fprintf('\n===== npie_struct eq =====\n');
R{end+1} = case_N1();
R{end+1} = case_N2();
R{end+1} = case_N3();
R{end+1} = case_N4(a,b);
end

results = [R{:}];
n_failed = sum(~[results.passed]);
fprintf('\n%d of %d cases passed.\n',numel(results)-n_failed,numel(results));
if nargout==0
    clear results
end
end



%% Helpers
function r = report(name,passed,message)
% R = REPORT(NAME,PASSED,MESSAGE) prints the outcome of a case and returns
% it as a result.
% INPUT
% - name:       char, the name of the case;
% - passed:     logical, whether the case passed;
% - message:    char, a detail of the outcome (may be empty);
% OUTPUT
% - r:          'struct' with fields name, passed and message.
r = struct('name',name,'passed',logical(passed),'message',message);
if passed && isempty(message)
    fprintf('passed  %s\n',name);
elseif passed
    fprintf('passed  %s: %s\n',name,message);
else
    fprintf('FAILED  %s: %s\n',name,message);
end
end

function str = where(ME)
% STR = WHERE(ME) returns ' [file:line]' of the deepest frame of the error
% ME that is not in this test file.
% INPUT
% - ME:     'MException' object;
% OUTPUT
% - str:    char, empty if no such frame exists.
str = '';
this_file = mfilename('fullpath');
for k=1:numel(ME.stack)
    if ~strcmp(erase(ME.stack(k).file,'.m'),this_file)
        [~,f] = fileparts(ME.stack(k).file);
        str = sprintf(' [%s:%d]',f,ME.stack(k).line);
        return
    end
end
end

function [V,lhs,base] = build(kind,a,b)
% [V,LHS,BASE] = BUILD(KIND,A,B) declares the components of a test system.
% INPUT
% - kind:   char: 'y' (states x, y), 'y2' (2x1 states x2, y), 'y21' (2x1
%           state x2, scalar y), 'y12' (scalar x, 2x1 state y), 'yw'
%           (input w, state y), 'yw2' (2x1 input w2, state y), 'X' (state
%           x, ODE state X), 'Xw' (input w, ODE state X), '2Dy' (x(s1,s2),
%           y(s1)), '2DX' (x(s1,s2), X), '2Dyy' (x(s1,s2), y(s1,s2));
% - a, b:   'double' scalars, the domain [a,b] of the 1D systems;
% OUTPUT
% - V:      'struct' with the declared components ('pde_struct' objects);
% - lhs:    function handle, lhs(r) is the equation of y (or X) with the
%           term r added to its right-hand side;
% - base:   'pde_struct' with the other equations and the BCs ([] if none).
pvar s t s1 s2
clear('stateNameGenerator');
V = struct();
heat = @(x) [diff(x,t)==diff(x,s,2); subs(diff(x,s),s,a)==0; subs(x,s,b)==0];
switch kind
    case {'y','y2','y21','y12','yw','yw2'}
        % The component in the term, then the state y
        if any(strcmp(kind,{'y','y12'})),  V.x = pde_var(s,[a,b]);         base = heat(V.x);   end
        if any(strcmp(kind,{'y2','y21'})), V.x2 = pde_var(2,s,[a,b]);      base = heat(V.x2);  end
        if strcmp(kind,'yw'),   V.w = pde_var('in',1,s,[a,b]);  base = [];          end
        if strcmp(kind,'yw2'),  V.w2 = pde_var('in',2,s,[a,b]); base = [];          end
        ny = 1 + any(strcmp(kind,{'y2','y12'}));
        V.y = pde_var(ny,s,[a,b]);
        BCy = [subs(V.y,s,a)==0; subs(V.y,s,b)==0];
        if isempty(base), base = BCy; else, base = [base; BCy]; end
        lhs = @(r) diff(V.y,t)==diff(V.y,s,2)+r;
    case 'X'
        V.x = pde_var(s,[a,b]);     V.X = pde_var(1);
        base = heat(V.x);
        lhs = @(r) diff(V.X,t)==-V.X+r;
    case 'Xw'
        V.w = pde_var('in',1,s,[a,b]);  V.X = pde_var(1);
        base = [];
        lhs = @(r) diff(V.X,t)==-V.X+r;
    case {'2Dy','2DX','2Dyy'}
        % A 2D heat equation on [0,1]^2, and the state y or X
        V.x = pde_var([s1;s2],[0,1;0,1]);
        base = [diff(V.x,t)==diff(V.x,s1,2)+diff(V.x,s2,2); subs(diff(V.x,s1),s1,0)==0; subs(V.x,s1,1)==0;
                subs(diff(V.x,s2),s2,0)==0; subs(V.x,s2,1)==0];
        if strcmp(kind,'2Dy')
            V.y = pde_var(s1,[0,1]);
            base = [base; subs(V.y,s1,0)==0; subs(V.y,s1,1)==0];
            lhs = @(r) diff(V.y,t)==diff(V.y,s1,2)+r;
        elseif strcmp(kind,'2DX')
            V.X = pde_var(1);
            lhs = @(r) diff(V.X,t)==-V.X+r;
        else
            V.y = pde_var([s1;s2],[0,1;0,1]);
            base = [base; subs(V.y,s1,0)==0; subs(V.y,s1,1)==0; subs(V.y,s2,0)==0; subs(V.y,s2,1)==0];
            lhs = @(r) diff(V.y,t)==diff(V.y,s1,2)+diff(V.y,s2,2)+r;
        end
end
end

function P = add_base(eqn,base)
% P = ADD_BASE(EQN,BASE) appends the equations BASE to EQN.
% INPUT
% - eqn:    'pde_struct', an equation;
% - base:   'pde_struct' or [], the other equations;
% OUTPUT
% - P:      'pde_struct', the equations together.
if isempty(base), P = eqn; else, P = [eqn; base]; end
end

function PIE = cv(P)
% PIE = CV(P) converts the PDE P to a PIE without printing.
% INPUT
% - P:      'pde_struct', the PDE;
% OUTPUT
% - PIE:    'pie_struct' or 'npie_struct', the PIE.
evalc('PIE = convert_PIETOOLS_PDE(P,''silent'');');
end

function [e,err] = nl_eval_eq(N1,N2)
% [E,ERR] = NL_EVAL_EQ(N1,N2) compares the right-hand sides f of two
% nonlinear PIEs with the same state, by evaluation at random polynomial
% values (quadratic in the first variable) of the fundamental state.
% INPUT
% - N1, N2: 'npie_struct' objects;
% OUTPUT
% - e:      logical, true if the values agree to 1e-8;
% - err:    'double', the largest difference of the coefficients.
names = N1.f.varname;
v = polynomial(N1.space.vars(1));
err = 0;
for trial=1:3
    xvals = cell(numel(names),1);
    for i=1:numel(names)
        c = randn(3,1);
        xvals{i} = (c(1)+c(2)*v+c(3)*v^2)*ones(N1.f.dims.in(i),1);
    end
    Y = polynomial(apply_poly(N1.f,xvals) - apply_poly(N2.f,xvals));
    if ~isempty(Y.coefficient)
        err = max(err,full(max(abs(Y.coefficient(:)))));
    end
end
e = err<1e-8;
end

function y = apply_cop(P,xf)
% Y = APPLY_COP(P,XF) applies the 'copvar' operator P to polynomial
% functions, block by block.
% INPUT
% - P:      'copvar' object with m x n blocks;
% - xf:     n x 1 cell of 'polynomial' functions, one per block column;
% OUTPUT
% - y:      m x 1 cell of 'polynomial' functions, P*xf by block row.
y = cell(size(P.C,1),1);
for i=1:size(P.C,1)
    y{i} = polynomial(0);
    for j=1:size(P.C,2)
        if ~isempty(P.C{i,j})
            y{i} = y{i} + apply_sopvar(P.C{i,j},xf{j});
        end
    end
end
end

function m = max_coeff(p)
% M = MAX_COEFF(P) returns the largest absolute coefficient of P.
% INPUT
% - p:      'polynomial' or 'double';
% OUTPUT
% - m:      'double', 0 for an empty polynomial.
p = polynomial(p);
if isempty(p.coefficient), m = 0; else, m = full(max(abs(p.coefficient(:)))); end
end

function r = eq_case(name,kind,r1,r2,a,b)
% R = EQ_CASE(NAME,KIND,R1,R2,A,B) checks that the terms r1(V) and r2(V),
% each placed in the equation of y (or X) of the system KIND, give the
% same PIE.
% INPUT
% - name:   char, the name of the case;
% - kind:   char, the system (see build);
% - r1, r2: function handles of the components V, the term and its
%           reference;
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result (fields name, passed, message).
try
    [V,lhs,base] = build(kind,a,b);
    P1 = add_base(lhs(r1(V)),base);
    P2 = add_base(lhs(r2(V)),base);
catch ME
    r = report(name,false,['declaring: ',ME.message,where(ME)]);     return
end
try, PIE1 = cv(P1);
catch ME, r = report(name,false,['converting: ',ME.message,where(ME)]);   return
end
try, PIE2 = cv(P2);
catch ME, r = report(name,false,['converting the reference: ',ME.message,where(ME)]);   return
end
warn_state = warning('off','npie_struct:eqUndecided');
e = PIE1==PIE2;
warning(warn_state);
if ~e && isa(PIE1,'npie_struct') && isa(PIE2,'npie_struct')
    % Polynomials not shown equal by ==: compare f by evaluation
    [e,err] = nl_eval_eq(PIE1,PIE2);
    r = report(name,e,sprintf('f compared by evaluation, max difference %.1e',err));
    return
end
r = report(name,e,'');
end

function r = err_case(name,kind,r1,a,b)
% R = ERR_CASE(NAME,KIND,R1,A,B) checks that declaring or converting the
% term r1(V) in the system KIND throws an error of PIETOOLS, not an
% internal MATLAB error.
% INPUT
% - name:   char, the name of the case;
% - kind:   char, the system (see build);
% - r1:     function handle of the components V, the term;
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result, with the error message.
try
    [V,lhs,base] = build(kind,a,b);
    cv(add_base(lhs(r1(V)),base));
    r = report(name,false,'converted, but an error was expected');
catch ME
    r = report(name,~startsWith(ME.identifier,'MATLAB:'),['error: ',ME.message,where(ME)]);
end
end

function r = ok_case(name,kind,r1,a,b)
% R = OK_CASE(NAME,KIND,R1,A,B) checks that the term r1(V) in the system
% KIND converts.
% INPUT
% - name:   char, the name of the case;
% - kind:   char, the system (see build);
% - r1:     function handle of the components V, the term;
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result.
try
    [V,lhs,base] = build(kind,a,b);
    cv(add_base(lhs(r1(V)),base));
    r = report(name,true,'');
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = eq_eqn_case(name,f1,f2,a,b)
% R = EQ_EQN_CASE(NAME,F1,F2,A,B) checks that the equations f1(V) and
% f2(V) of y, in the system 'y', give the same PIE.
% INPUT
% - name:   char, the name of the case;
% - f1, f2: function handles of the components V, the equation and its
%           reference;
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result.
try
    [V,~,base] = build('y',a,b);
    r = report(name,cv([f1(V); base])==cv([f2(V); base]),'');
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = err_eqn_case(name,f1,a,b)
% R = ERR_EQN_CASE(NAME,F1,A,B) checks that the equation f1(V) in the
% system 'y', with a 2x1 state x2 added, throws an error of PIETOOLS.
% INPUT
% - name:   char, the name of the case;
% - f1:     function handle of the components V, the equation;
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result, with the error message.
pvar s
try
    [V,~,base] = build('y',a,b);
    V.x2 = pde_var(2,s,[a,b]);
    cv([f1(V); base]);
    r = report(name,false,'converted, but an error was expected');
catch ME
    r = report(name,~startsWith(ME.identifier,'MATLAB:'),['error: ',ME.message,where(ME)]);
end
end



%% Cases that do not fit the helpers
function r = case_E5(a,b)
% R = CASE_E5(A,B) checks that a delay in a boundary condition and in an
% output converts.
% INPUT
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result.
pvar s t
name = 'E5 delay in a BC and in an output';
try
    clear('stateNameGenerator');
    x = pde_var(s,[a,b]); z = pde_var('out',1);
    cv([diff(x,t)==diff(x,s,2); subs(diff(x,s),s,a)==subs(subs(x,s,b),t,t-1); subs(x,s,b)==0; z==subs(subs(x,s,b),t,t-1)]);
    r = report(name,true,'');
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = case_B1(a,b)
% R = CASE_B1(A,B) checks an output with a nested integral against the
% output with the single integral of the same kernel.
% INPUT
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result.
pvar s t
name = 'B1 output with a nested integral';
try
    clear('stateNameGenerator');
    x = pde_var(s,[a,b]); z = pde_var('out',1);
    base = [diff(x,t)==diff(x,s,2); subs(diff(x,s),s,a)==0; subs(x,s,b)==0];
    r = report(name,cv([base; z==int(int(x,s,[a,s]),s,[a,b])])==cv([base; z==int((b-s)*x,s,[a,b])]),'');
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = case_B2(a,b)
% R = CASE_B2(A,B) checks a nonlocal BC with a nested integral against the
% BC with the single integral of the same kernel.
% INPUT
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result.
pvar s t
name = 'B2 nonlocal BC with a nested integral';
try
    clear('stateNameGenerator');
    x = pde_var(s,[a,b]);
    P1 = [diff(x,t)==diff(x,s,2); subs(diff(x,s),s,a)==0; subs(x,s,b)==int(int(x,s,[a,s]),s,[a,b])];
    P2 = [diff(x,t)==diff(x,s,2); subs(diff(x,s),s,a)==0; subs(x,s,b)==int((b-s)*x,s,[a,b])];
    r = report(name,cv(P1)==cv(P2),'');
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = case_B3(a,b)
% R = CASE_B3(A,B) checks a BC with a reversed integral against the
% negated integral.
% INPUT
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result.
pvar s t
name = 'B3 BC with a reversed integral';
try
    clear('stateNameGenerator');
    x = pde_var(s,[a,b]);
    P1 = [diff(x,t)==diff(x,s,2); subs(diff(x,s),s,a)==0; subs(x,s,b)==int(x,s,[b,a])];
    P2 = [diff(x,t)==diff(x,s,2); subs(diff(x,s),s,a)==0; subs(x,s,b)==-1*int(x,s,[a,b])];
    r = report(name,cv(P1)==cv(P2),'');
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = case_K1(a,b)
% R = CASE_K1(A,B) checks the constraint 0 = x2 - int_a^s x1 on the
% domain: for PDE states x1 (with x1_s(a)=0, x1(b)=0) and x2, with
% fundamental state xf = (x1_ss, x2_ss), T*xf must give (x1,x2), and K*xf
% must vanish if x2 = int_a^s x1, and not if x2 differs by p(s)(s-a)(s-b).
% INPUT
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result.
pvar s t
name = 'K1 constraint with a partial integral, K checked on states';
try
    clear('stateNameGenerator');
    x1 = pde_var(s,[a,b]); x2 = pde_var(s,[a,b]);
    PIE = cv([diff(x1,t)==diff(x1,s,2); diff(x2,t)==diff(x2,s,2)+x1; subs(diff(x1,s),s,a)==0; subs(x1,s,b)==0;
              subs(x2,s,a)==0; subs(x2,s,b)==0; x2-int(x1,s,[a,s])==0]);
    idx = PIE.pde2pie.x.idx(:,1);
    err_T = 0;  err_K = 0;  min_K = inf;
    for deg = 3:5
        % x1 of degree deg with x1_s(a)=0, x1(b)=0 and int_a^b x1 = 0, so
        % that x2 = int_a^s x1 satisfies x2(a)=x2(b)=0
        k = 0:deg;
        M = [k.*a.^max(k-1,0); b.^k; (b.^(k+1)-a.^(k+1))./(k+1)];
        c = null(M)*randn(deg-2,1);
        X1 = polynomial(0);
        for j=0:deg, X1 = X1 + c(j+1)*s^j; end
        X2 = int(X1,s,a,s);
        for violate = [false true]
            X2v = X2 + violate*(1+randn*s)*(s-a)*(s-b);
            xf = {diff(X1,s,2); diff(X2v,s,2)};     xf = xf(idx);
            x_pde = {X1; X2v};                      x_pde = x_pde(idx);
            Tx = apply_cop(PIE.T,xf);
            err_T = max(err_T,max(cellfun(@(p,q) max_coeff(p-q),Tx,x_pde)));
            Kx = max(cellfun(@max_coeff,apply_cop(PIE.K,xf)));
            if violate, min_K = min(min_K,Kx); else, err_K = max(err_K,Kx); end
        end
    end
    r = report(name,err_T<1e-8 && err_K<1e-8 && min_K>1e-4,...
        sprintf('max|T*xf-x| %.1e, max|K*xf| (satisfied) %.1e, min max|K*xf| (violated) %.1e',err_T,err_K,min_K));
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = case_C1()
% R = CASE_C1() checks horizontal and vertical concatenation, and the PIE
% built from the operators, for a PIE whose variables are not in
% alphabetical order (state on zz, input on aa).
% OUTPUT
% - r:      'struct', the result.
pvar t zz aa
name = 'C1 unsorted variable names: horzcat, vertcat, piess';
try
    clear('stateNameGenerator');
    x = pde_var(zz,[0,1]); w = pde_var('in',1,aa,[0,2]); z = pde_var('out',1);
    PIE = cv([diff(x,t)==diff(x,zz,2)+int(w,aa,[0,2]); subs(x,zz,0)==0; subs(x,zz,1)==0; z==int(x,zz,[0,1])]);
    Ph = initialize([PIE,PIE]);     Pv = initialize([PIE;PIE]);
    P5 = piess(PIE.T,PIE.A,PIE.Bw,PIE.Cz,PIE.Dzw);
    ok = isequal(Ph.space.vars,PIE.space.vars) && isequal(Pv.space.vars,PIE.space.vars) && P5==PIE;
    r = report(name,ok,'');
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = case_C3()
% R = CASE_C3() checks that concatenating PIEs in different variables
% throws an error of PIETOOLS.
% OUTPUT
% - r:      'struct', the result, with the error message.
pvar t s s1
name = 'C3 [PIE(s),PIE(s1)]';
try
    clear('stateNameGenerator');
    x = pde_var(s,[0,1]);   x1 = pde_var(s1,[0,1]);
    P1 = cv([diff(x,t)==diff(x,s,2); subs(x,s,0)==0; subs(x,s,1)==0]);
    P2 = cv([diff(x1,t)==diff(x1,s1,2); subs(x1,s1,0)==0; subs(x1,s1,1)==0]);
    [P1,P2]; %#ok<VUNUS>
    r = report(name,false,'concatenated, but an error was expected');
catch ME
    r = report(name,~startsWith(ME.identifier,'MATLAB:'),['error: ',ME.message,where(ME)]);
end
end

function r = case_N1()
% R = CASE_N1() checks that a nonlinear PIE equals itself when converted
% twice, and differs from the PIE with another nonlinear coefficient.
% OUTPUT
% - r:      'struct', the result.
pvar s t
name = 'N1 nonlinear PIE: == itself, ~= other coefficient';
try
    clear('stateNameGenerator');
    x = pde_var(s,[0,1]); w = pde_var('in',1);
    PDE1 = [diff(x,t)==diff(x,s,2)+x*x+w; subs(x,s,0)==0; subs(x,s,1)==0];
    PDE2 = [diff(x,t)==diff(x,s,2)+2*(x*x)+w; subs(x,s,0)==0; subs(x,s,1)==0];
    warn_state = warning('off','npie_struct:eqUndecided');
    ok = cv(PDE1)==cv(PDE1) && ~(cv(PDE1)==cv(PDE2));
    warning(warn_state);
    r = report(name,ok,'');
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = case_N2()
% R = CASE_N2() checks that a linear PDE with a boundary input, converted
% to an 'npie_struct' directly, equals npie_struct of its 'pie_struct'.
% OUTPUT
% - r:      'struct', the result.
pvar s t
name = 'N2 boundary input: npie_struct(P)==N';
try
    clear('stateNameGenerator');
    x = pde_var(s,[0,1]); w = pde_var('in',1);
    PDE = [diff(x,t)==diff(x,s,2); subs(x,s,0)==w; subs(x,s,1)==0];
    N = convert_PIETOOLS_PDE(PDE,{'silent','nonlinear'});
    r = report(name,npie_struct(cv(PDE))==N,'');
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = case_N3()
% R = CASE_N3() checks, as case N2, the two 'npie_struct' objects of a 2D
% PDE.
% OUTPUT
% - r:      'struct', the result.
pvar t s1 s2
name = 'N3 2D: npie_struct(P)==N';
try
    clear('stateNameGenerator');
    x = pde_var([s1;s2],[0,1;0,1]);
    PDE = [diff(x,t)==diff(x,s1,2)+diff(x,s2,2); subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
    N = convert_PIETOOLS_PDE(PDE,{'silent','nonlinear'});
    r = report(name,npie_struct(cv(PDE))==N,'');
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end

function r = case_N4(a,b)
% R = CASE_N4(A,B) checks that == of the PIEs of x*x_s and x_s*x, which
% are equal but not recognized as such, returns false with the warning
% 'npie_struct:eqUndecided', and that their f agree by evaluation.
% INPUT
% - a, b:   'double' scalars, the domain;
% OUTPUT
% - r:      'struct', the result.
name = 'N4 x*x_s vs x_s*x: == warns, f equal by evaluation';
pvar s
try
    [V,lhs,base] = build('y',a,b);
    N1 = cv(add_base(lhs(V.x*diff(V.x,s)),base));
    N2 = cv(add_base(lhs(diff(V.x,s)*V.x),base));
    lastwarn('');
    warn_state = warning('off','npie_struct:eqUndecided');
    e = N1==N2;
    warning(warn_state);
    [~,warn_id] = lastwarn;
    [e_val,err] = nl_eval_eq(N1,N2);
    r = report(name,(e || strcmp(warn_id,'npie_struct:eqUndecided')) && e_val,...
        sprintf('== %d, warning "%s", f difference %.1e',e,warn_id,err));
catch ME
    r = report(name,false,[ME.message,where(ME)]);
end
end
