cppFunction('NumericVector qalloC(NumericVector p, double min_phi,
                                  double max_phi, double a, bool lower_tail,
                                   bool log_p) {
   int n = p.size();
   NumericVector out(n);

   if(log_p){
     p = exp(p);
   }

   if(!lower_tail){
     p = 1 - p;
   }

  if(a == 1){
    for(int i = 0; i < n; ++i) {
    if((p[i] > 1) | (p[i] < 0)){
        out[i] = NAN;
      } else {
        out[i] = exp(p[i]*(log(max_phi)-log(min_phi))+log(min_phi));
      }
    }
  } else {

      for(int i = 0; i < n; ++i) {
    if((p[i] > 1) | (p[i] < 0)){
        out[i] = NAN;
      } else {
        out[i] = pow((p[i]*(pow(max_phi, (1-a)) - pow(min_phi, (1-a))) +
        pow(min_phi, (1-a))), 1/(1-a));
      }
    }
  }
  return out;
}')
