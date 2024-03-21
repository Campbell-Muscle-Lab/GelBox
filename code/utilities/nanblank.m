function output = nanblank(values)
[m,n] = size(values)
for i = 1:m
    mask(i,1) = contains(values(i,:),'NaN');
end
if nnz(mask)
    output = string(values);
    output(mask) = " ";
    output = char(output);
else
    output = values;
end
end