function strOut = stateNameGenerator()
persistent n

if isempty(n)
    n=0;
end
n=n+1;

strOut = n;
end