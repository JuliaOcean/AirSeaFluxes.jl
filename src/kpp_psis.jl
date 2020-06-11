function kpp_psis(xi)
if xi>=0
  y=1+5*xi;
elseif (xi<0) & (xi>=-1.)
  y=(1-16*xi)^(-1/2);
else
  y=(-28.86-98.96*xi)^(-1/3);
end
return y
end
