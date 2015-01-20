#ifndef TRANSFORM_HEADER
#define TRANSFORM_HEADER

void wgs2gcj(double wgsLat, double wgsLng, double *gcjLat, double *gcjLng);
void gcj2wgs(double gcjLat, double gcjLng, double *wgsLat, double *wgsLnt);
void gcj2wgs_exact(double gcjLat, double gcjLng, double *wgsLat, double *wgsLnt);
double distance(double latA, double lngA, double latB, double lngB);

void bd_encrypt(double gg_lat, double gg_lon, double *bd_lat, double *bd_lon);
void bd_decrypt(double bd_lat, double bd_lon, double *gg_lat, double *gg_lon);

#endif
