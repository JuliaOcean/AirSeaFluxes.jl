function angle_of_incidence(lat,lon,jd,time)
# angle_of_incidence calculates the cosine of the sun incident angle
#  y = ANGLE_OF_INCIDENCE(lat,lon,jd,time) calculates the cosine of the sun incident angle for 
#   a given location (lat,lon), julian day (jd) and time in hours and fraction of hours.

# inclination angle
delta=(asin(-sin(23.45/180*pi).*cos(360/365.25*(jd+10)/180*pi)));
# lat to radians
lamda=lat/180*pi;
# time to radians
omega=time/24*2*pi+lon/360*2*pi+pi;

y=cos(delta)*cos(lamda)*cos(omega)+sin(delta)*sin(lamda);

return y

end
