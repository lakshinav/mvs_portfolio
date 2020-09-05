function z=lags(x,l)
m = length(x);
X = cell(1, l);
[X{:}] = ndgrid(x);
X = X(end : -1 : 1);
z = cat(l+1, X{:});
z = reshape(z, [m^l, l]);
end