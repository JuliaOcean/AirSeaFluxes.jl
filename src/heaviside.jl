function heaviside(x)
# HEAVISIDE
# y = HEAVISIDE(x) returns 0 when x<0 and 1 when x>=0
y=zeros(size(x));
y[x.>=0].=1;
return y
end

