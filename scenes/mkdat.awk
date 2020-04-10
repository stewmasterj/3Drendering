#!/usr/bin/awk -f
BEGIN {
 nx = 4;
 ny = nx;
 nz = nx;
 for (i=-nx; i<=nx; ++i) {
  for (j=-ny; j<=ny; ++j) {
   for (k=-nz; k<=nz; ++k) {
    dist = sqrt((i)**2 +(j)**2 +(k)**2);
    d = (dist<0.5*nx ? 0.866 : 0.5 ); # point radius
    t = (dist<0.5*nx ? 0.0 : 200.0 ); # transperency [0-255]
    r = (dist<0.5*nx ? 255 : 0 ); # red [0-255]
    g = 0; #(dist<5 ? 255 : 0 );  # green [0-255]
    b = (dist<0.5*nx ? 0 : 255 ); # blue [0-255]
    print i, j, k, d, t, r, g, b
   }
  }
 }
}

