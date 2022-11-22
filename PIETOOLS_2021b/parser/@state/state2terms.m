function out = state2terms(obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that converts objA to terms format
% Input: 
% objA - state class objects
% Output:
% out - A*[objA], where A is an identity PI operator

s.type = '.'; s.subs = 'veclength';
opvar T; T.R.R0 = eye(sum(subsref(obj,s)));
termsOut = obj;
out = terms(T,termsOut);
end