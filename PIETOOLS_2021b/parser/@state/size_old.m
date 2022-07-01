function varargout = size(objA,dim)
veclen = length(objA);
if nargin==1
    varargout{1} = [veclen 1];
elseif nargin==2
    if ~(isa(dim,'double') && isreal(dim) &&  all(size(dim)==[1 1]) &&...
            all(ceil(dim)==floor(dim)) )
        error('Dimension must be a positive integer scalar');
    else
        if dim==1            
            out = veclen;
        elseif dim==2
            out = 1;
        else
            error('Dimension must be a positive integer scalar 1 or 2');
        end
    end
    varargout{1} = out;  
end 
end