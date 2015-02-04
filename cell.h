int cell(double x,double L, const int imax) {
  double x2 = x + L/2;
  
  int n = round(x2/L*double(imax-1));
  if (n < 0) n = 0;
  if (n > imax -1) n = imax -1;
  return n;
    }
