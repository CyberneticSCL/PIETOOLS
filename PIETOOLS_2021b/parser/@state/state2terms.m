function out = state2terms(obj)
s.type = '.'; s.subs = 'veclength';
opvar T; T.R.R0 = eye(subsref(obj,s));
termsOut = obj;
out = terms(T,termsOut);
end