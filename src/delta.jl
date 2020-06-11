function delta(n)
# DELTA calculates the sun's declination angle
# y = DELTA(jd) gets the Julian day (jd) and calculate declination of the sun

d=(asin(-sin(23.45/180*pi).*cos(360/365.25*(n+10)/180*pi)));

return d

end
