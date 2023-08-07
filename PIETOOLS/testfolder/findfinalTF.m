function varargout = findfinalTF(Cq1,invTA,Bq2,si)
N = length(si); 
nz = size(Cq1,1)/2; nw = size(Bq2,2)/2;
R0 = invTA{1}; R1 = invTA{2}; R2 = invTA{3};
var1 = pvar('s');

ds = si(2)-si(1);
diagint=0; lowint=0; upint = 0;
for i=1:N
    diagint = diagint+double(subs(Cq1,var1,si(i)))*R0{i}*double(subs(Bq2,var1,si(i)))*ds;
    for j=i:N
        lowint = lowint+double(subs(Cq1,var1,si(i)))*R1{i,j}*double(subs(Bq2,var1,si(j)))*ds^2;
    end
    for k=1:i
        upint = upint+double(subs(Cq1,var1,si(i)))*R2{i,k}*double(subs(Bq2,var1,si(k)))*ds^2;
    end
end

tmp = diagint+lowint+upint; % TF is size n_grid*2*nz*2*nw;

varargout{1} = tmp(1:nz,1:nw);
varargout{2} = tmp(1:nz,nw+1:2*nw);
end