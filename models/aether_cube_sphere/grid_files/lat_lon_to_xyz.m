function[xyz] = lat_lon_to_xyz(lat, lon)

xyz(1) = cos(lat) * cos(lon);
xyz(2) = cos(lat) * sin(lon);
xyz(3) = sin(lat);
