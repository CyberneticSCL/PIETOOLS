function disp(obj, type)
if nargin==2 && strcmp(type,"latex")
    latex(obj,"show");
else
    fprintf('%s\n', char(obj));
end
end

function s = char(obj)
    s = format(obj);
end