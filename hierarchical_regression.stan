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
	real<lower=-10,upper=10> logFeu_bar[N_R]; //number of station means for intercept equal to number of biomes	
	real<lower=-10,upper=10> b_bar[N_R];      //number of station means for slope equal to number of biomes 
	real<lower=-10,upper=10> logFeu[N_p];     //p intercepts for individual stations
	real<lower=0,upper=3> b[N_p];	          //p slopes for individual stations
	real<lower=1e-15,upper=10> logFeu_sd;     //standard deviation for individual intercepts
	real<lower=1e-15,upper=3> b_sd;           //standard deviation for individual slopes
	real<lower=1e-15> sigma;                  //standard deviation of measurement noise for individual observations
}
model{
	for(j in 1:N_R){
		for(i in (np[j]+1):np[j+1]){                                                                    //loop over stations
			logFeu[i]            ~ normal(logFeu_bar[j],logFeu_sd);                      //indv int drawn from dist with biome-specific mean
			b[i]                 ~ normal(b_bar[j],     b_sd);                           //indv slope from dist with biome-specific mean
			y[(ni[i]+1):ni[i+1]] ~ normal(logFeu[i] - b[i]*x[(ni[i]+1):ni[i+1]], sigma); //likelihood of data, station by station
		}
	}
}
