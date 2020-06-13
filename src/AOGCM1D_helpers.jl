
"""
    angle_of_incidence(lat,lon,jd,time)

Calculate the cosine of the sun incident angle for a given location (lat,lon),
julian day (jd) and time in hours and fraction of hours.
"""
function angle_of_incidence(lat, lon, jd, time)
    # inclination angle
    delta = (asin(-sin(23.45 / 180 * pi) .* cos(360 / 365.25 * (jd + 10) / 180 * pi)))
    # lat to radians
    lamda = lat / 180 * pi
    # time to radians
    omega = time / 24 * 2 * pi + lon / 360 * 2 * pi + pi

    return cos(delta) * cos(lamda) * cos(omega) + sin(delta) * sin(lamda)
end

"""
    delta(jd)

Calculates the Sun's declination angle as a function of Julian day (jd)
"""
delta(n) = (asin(-sin(23.45 / 180 * pi) .* cos(360 / 365.25 * (n + 10) / 180 * pi)))


"""
    heaviside(x)

Heaviside function : 0 when x<0 and 1 when x>=0
"""
function heaviside(x)
    y = zeros(size(x))
    y[x.>=0] .= 1
    return y
end
