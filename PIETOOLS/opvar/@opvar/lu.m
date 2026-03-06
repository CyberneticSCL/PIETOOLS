function [L,U,info] = lu(A,opts)
% Given a 4-PI operator P = [A, Bop; 
%                            Cop, Dop]
% This function returns L and U such that P = L*U
% More specifically,
% P = [A,    0]  *  [I, A^{-1}B]
%     [C,   LD]     [0,      UD]
% where LD x(s) = R0(s)*x(s) + int_0^s K1(s,t)x(t)dt 
% and   UD x(s) =     I*x(s) + int_s^1 K2(s,t)x(t)dt.

end