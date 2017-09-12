data {
  int G;
  int P;
  int T;
  vector[T] volumes_cor_xf;
  matrix[G,P] genes;
}

parameters {
  vector[G] weights;
  real<lower=0> sigma;
}

model {
  matrix[G,P] weighted_vols = diag_pre_multiply(weights, genes);
  vector[T] corr_xf_difs;
  real l1 = 0;
  real l2 = 0;
  
  {
    int count = 1; 
    
    for(i in 1:(P - 1)){
      for(j in (i + 1):P){
	vector[G] x = col(weighted_vols, i);
	vector[G] y = col(weighted_vols, j);
	real EXY = 0;
	real EX = 0;
	real EXX = 0;
	real EY = 0;
	real EYY = 0;
	real corr;
	
	for(ind in 1:G){
	  EX = EX + x[ind]/G;
	  EY = EY + y[ind]/G;
	  EXX = EXX + pow(x[ind], 2)/G;
	  EYY = EYY + pow(y[ind], 2)/G;
	  EXY = EXY + x[ind] * y[ind] /G;

	  l1 = l1 + fabs(weights[ind]);
	  l2 = l2 + pow(weights[ind], 2);
	}

	corr = (EXY - EX * EY) / sqrt(EXX - pow(EX,2)) / sqrt(EYY - pow(EY,2));
	corr_xf_difs[count] =
	  volumes_cor_xf[count] - .5 * log( (1 + corr) / (1 - corr) );

	if(count > T) reject("T is counting too high");
       	count = count + 1;
      }
    }
  }

  target += normal_lpdf(sum(fabs(corr_xf_difs)) | 0, sigma);
  target += exponential_lpdf(1 * l1 + 5 * l2 | 1/pow(sigma,2) );
  target += normal_lpdf( sigma | 0, 3);
}
  
