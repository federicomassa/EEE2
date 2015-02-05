double absval2(double x) {
  if (x >= 0) return x;
  else return (-x);
}

int myround(double x) {
  int sign;
  if (x >= 0) sign = 1;
  else sign = -1;
  double result;
  if (absval2(x) - double(int(absval(x))) < 0.5) result = double(floor(absval2(x)));
  else result = double(ceil(absval2(x)));
  return result*double(sign);
}

//questo vale per EEE perché la distanza totale tra le strip non corrisponde
//alla lunghezza della camera
int xcell(double x,double L, const int imax) {
  double dist_strip = L - L/double(imax);  
double x2 = x + dist_strip/2;
  
  int n = myround(x2/dist_strip*double(imax-1));
  if (n < 0) n = 0;
  if (n > imax -1) n = imax -1;
  return n;
    }

int ycell(double y,double L, const int imax) {  
double y2 = y + L/2;
  
  int n = myround(y2/L*double(imax-1));
  if (n < 0) n = 0;
  if (n > imax -1) n = imax -1;
  return n;
    }


double xstrip(double x) {
  double strd = 82.0/24.0;
  double start = -(82-strd)/2;
  int c = xcell(x,82,24);
  
  return (start+double(c)*strd);
}
