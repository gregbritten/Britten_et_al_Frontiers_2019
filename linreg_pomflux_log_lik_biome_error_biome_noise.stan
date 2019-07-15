data{
  int N_R;       //number of biome
	int np[N_R+1]; //linear indices for biomes 
	int N_p;       //total number of stations
	int N;         //total number of observations
	int ni[N_p+1]; //linear indices for station data
	vector[N] y;   //y variable
	vector[N] x;   //x variable
}
parameters{
	real<lower=-2,upper=8> logFeu_bar[N_R]; //number of station means for intercept equal to number of biomes	
	real<lower=-2,upper=8> b_bar[N_R];      //number of station means for slope equal to number of biomes 
	real<lower=-2,upper=8> logFeu[N_p];     //p intercepts for individual stations
	real<lower=0,upper=3> b[N_p];	          //p slopes for individual stations
	real<lower=1e-15,upper=5> logFeu_sd[N_R];     //standard deviation for individual intercepts
	real<lower=1e-15,upper=3> b_sd[N_R];           //standard deviation for individual slopes
	real<lower=1e-15,upper=10> sigma[N_R];                  //standard deviation of measurement noise for individual observations
}
model{
	for(j in 1:N_R){
		for(i in (np[j]+1):np[j+1]){                                  //loop over stations
			logFeu[i]            ~ normal(logFeu_bar[j],logFeu_sd[j]);  //indv int drawn from dist with biome-specific mean
			b[i]                 ~ normal(b_bar[j],     b_sd[j]);       //indv slope from dist with biome-specific mean
			y[(ni[i]+1):ni[i+1]] ~ normal(logFeu[i] - b[i]*x[(ni[i]+1):ni[i+1]], sigma[j]); //likelihood of data, station by station
		}
	}
}
//generated quantities{
//	vector[N] log_p_y_theta;
//	int j;
//	j = 1;
//	for(i in 1:p){
//		for(k in 1:n[i]){
//			log_p_y_theta[j] = normal_lpdf(y[j] | logJ0[i] - b[i]*x[j], sigma);
//			j = j+1;
//		}
//	}
//}


