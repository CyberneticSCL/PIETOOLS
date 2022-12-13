function strOut = stateNameGenerator()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% internal method to generate unique ID for state class objects

persistent n

if isempty(n)
    n=0;
end
n=n+1;

strOut = n;
end