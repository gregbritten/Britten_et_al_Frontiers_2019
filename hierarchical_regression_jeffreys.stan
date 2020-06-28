data{
    int N_R;   
	int np[N_R+1]; 
	int N_p;       
	int N;         
	int ni[N_p+1]; 
	vector[N] y;   
	vector[N] x;   
}
parameters{
	real<lower=-10,upper=10> a_bar[N_R]; 
	real<lower=-10,upper=10> b_bar[N_R];      
	real<lower=-10,upper=10> a[N_p];     
	real<lower=0,upper=3> b[N_p];	          
	real a_sd;     
	real b_sd;           
	real<lower=1e-15> sigma;                  
}
transformed parameters{
	real a_sd_jeffreys = exp(a_sd);
	real b_sd_jeffreys = exp(b_sd);
}
model{
	for(j in 1:N_R){
		for(i in (np[j]+1):np[j+1]){                                                     //loop over stations
			a[i]                 ~ normal(a_bar[j], a_sd_jeffreys);                      //indv int drawn from dist with biome-specific mean
			b[i]                 ~ normal(b_bar[j], b_sd_jeffreys);                      //indv slope from dist with biome-specific mean
			y[(ni[i]+1):ni[i+1]] ~ normal(a[i] - b[i]*x[(ni[i]+1):ni[i+1]], sigma); //likelihood of data, station by station
		}
	}
}
