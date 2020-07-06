function solver = solver
% This function sets the solver called by sossolve
%
% Choose whether SeDuMi or SDPT3 is used.

solver = 'SeDuMi';
%solver = 'SDPT3';

