function kpp_psim(xi)

if xi>=0
  y=1+5*xi;
elseif (xi<0) & (xi>=-0.2)
  y=(1-16*xi)^(-1/4);
else
  y=(1.26-8.38*xi)^(-1/3);
end
return y
end
