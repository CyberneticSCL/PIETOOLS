function e = end(F,k,n)
m = F.dim(1);
p = F.dim(2);

if n == 1
    e = m*p;
else
    if k == 1
        e = m;
    else
        e = p;
    end
end
end
