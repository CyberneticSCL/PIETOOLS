function C = polyopvar_times_v2(A,B)                         
% C = POLYOPVAR_TIMES_V2(A,B) returns the 'polyopvar' object C
% representing the elementwise product of two distributed polynomials,
% by calling `times`.
%
% It is used for p1*g in LocalStability, which is the product of
% functional distributed polynomials described by Lemma 5 in the
% Automatica paper.
%
% INPUTS
% - A:     m x n 'polyopvar' object;
% - B:     m x n 'polyopvar' object;
%
% OUTPUTS
% - C:     'polyopvar' object representing A.*B;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - polyopvar_times_v2
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CR, 09/28/2026: LocalStability copy of polyopvar/times.m.
% CR, 09/28/2026: Preserve zero entries in integration-order matrices while
%                   reordering the positive entries required by Lemma 5.
% DJ, 10/15/2026: Updated to use new polyopvar functions;

% DJ, 10/14/2026: Update to new polyopvar class

C = times(A,B);                                                             % DJ, 10/14/2026

end
