function holtslag_psim(xi)
if xi>0
    y=1+5*xi;
else
    y=(1-16*xi)^(-1/4);
end
return y
end
