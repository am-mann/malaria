// applicable to weekly cases data by symptoms onset date
functions {
	real[] SEIR(real t, 
							real[] y, 
							real[] theta, 
							real[] x_r,
							int[] x_i
	) {
		int TotalPop = x_i[1];            // total population
		real L_h = x_r[1];               //  
		real beta_h = x_r[2];               //  
		real mu_h = x_r[3];          // 
		real sigma_h = x_r[4];          // 
		real gamma_h = x_r[5];          // 
		real beta_v = x_r[6];          // 
		real mu_v = x_r[7];          // 
		real PIE = x_r[8];
		
		real temp[12];
		temp[1:12] = x_r[16+1:16+12];

		real y0_vars[8];

		row_vector[8] dydt;       //  an array containing the derivative of all the compartment arrange in the order S, E1, E2, I1, I2, R, and Inci 
		
		real Nh;
		real b;
		real lambda_h;
		real lambda_v;
		
		real Sh = y[1]; 
		real Eh = y[2];
		real Ih = y[3]; 
		real Rh = y[4]; 
		real Sv =  y[5];
		real Ev = y[6];
		real Iv =  y[7];
		real CumInci = y[8];
		
		real dSh_dt;
		real dEh_dt;
		real dIh_dt;
		real dRh_dt;
		real dSv_dt;
		real dEv_dt;
		real dIv_dt;
		real dInci_dt;
		
		// real IhR_0 = theta[1];
		// real IhU_0 = theta[2];   
		real i0 = theta[1];
		real q = theta[2];
		real Iv0 = theta[3];
		real w_h = theta[4];
		real delta = theta[5];
		real q_birth = theta[6];
		real q_kappa = theta[7];
		real q_mu_v = theta[8];
		//print("kappa_v", kappa_v);
		//real tc = theta[10];

		y0_vars = x_r[9:16];
		int temp_index;
    temp_index = bin_search(floor(t), 0, 11)+1;
    if(temp_index > 11)
    temp_index =11;
    //print(t);
    //print(temp_index);
    //print("temp_index", temp_index);
    real T = (t-floor(t))*temp[temp_index+1]+(1-(t-floor(t)))*temp[temp_index];
    //print("temp", temp[temp_index]);
    //print("temperature", T, " and temp_index", temp_index);
    //turn bin_search into a floor function and then make a smoothed out T
    
		Nh = Sh + Eh + Ih + Rh;
		//b = q* 25 * (25 - t_min_b) * pow(t_max_b - 25, 0.5);
		//b =  30*q* T* (T - 13.8) * pow(40 - T, 0.5); // data for Aedes aegypti from Mordecai, 2019
		b =  30*q*0.000202* T* (T - 13.8) * pow(40 - T, 0.5);
		real ro;
		// ro=0.0025*T^2 - 0.094*T+1.0257; Ngwarakana
		ro=0.0011*T^2 - 0.049*T+0.609; //Mordecai
		//ro=A*T^2 + B*T+C;
		//ASK: How do you force ro to be positive
    //mu_v = -30*log(ro);
    mu_v = q_mu_v*30*ro;
		//b=q*(1 + t_min*cos(2*PIE*t_max*(t + tc) ) );
		lambda_h = b*beta_h*Iv/Nh;
		lambda_v = b*beta_v*(Ih)/Nh;
		
		real EFP = 0.00816* T* (T-14.7) * pow(34.4 - T, 0.5); 
		//real EFP =30*(A*T^2+B*T+C); 
    real B_T = EFP/mu_v; 
    // real MDR =30*0.0000952* T* (T-15.4) * pow(34.7 - T, 0.5); // Anopheles gambiae from Mordecai 2019
    real MDR =0.0000783* T* (T-11.6) * pow(39.1 - T, 0.5); // Anopheles gambiae from Mordecai 2019
    //real p_EA = -4.77+0.453*T+0.009*T^2;// Mordecai 2013
    real p_EA = -0.00765*(T-14.5)*(T-34.1);// Anopheles gambae from Mordecai 2019
    real birth = q_birth*B_T*p_EA*MDR*30;
    //print("q_birth: ",  q_birth);
    //print("q_mu_v: ", q_mu_v);
    //print("EFP: ", EFP);
    //print("T: ", T);
		//print("MDR:  ", MDR);
		//print("p_EA:  ", p_EA);
		//print("mu_v:  ", mu_v);
		//print("b:  ", b);
		//real mu_v;
		//death = T_ = get_temp(t,T);
    //%p=-0.000828*T_^2 + 0.0367*T_+0.522; (mordecai)
    //p=-tc*T^2 + 0.0367*T+0.522;
    //mu_v = 30*(1-p);
    real kappa_v = 30*0.000106*q_kappa* T* (T - 14.1) * pow(34.2 - T, 0.5); // data from modercai
    if (t < 10){
      //print("birth:  ", birth);
      //print("q_mu_v: ", q_mu_v);
      //print("mu_v: ", mu_v);
      //print("mu_v: ", 0.0011*(T^2) - (0.049*T)+0.609);
      //print("T: ", T^2);
		//  print("Sh: ", round(Sh), "Eh: ", round(Eh), "Ih: ", round(Ih), "Rh: ", round(Rh), "Nh: ", round(Nh), "Sv: ", round(Sv), "Ev: ", round(Ev), "Iv: ", round(Iv));
		//  print("birth: ", round(birth), "biting rate: ", round(b), "lambda_v*Sv: ", round(lambda_v*Sv), " mu_v: ", mu_v, "kappa_v: ", round(kappa_v));
		//  print("time: ",t);
    }
		dSh_dt = L_h - lambda_h*Sh + w_h*Rh - mu_h*Sh;    
		dEh_dt = lambda_h*Sh - sigma_h*Eh - mu_h*Eh;
		dIh_dt = sigma_h*Eh - gamma_h*Ih - (mu_h + delta)*Ih;
		dRh_dt = gamma_h*Ih - w_h*Rh - mu_h*Rh;
		dSv_dt = birth - lambda_v*Sv - mu_v*Sv;
		dEv_dt = lambda_v*Sv - (mu_v+kappa_v)*Ev;
		dIv_dt = kappa_v*Ev - mu_v*Iv;
		dInci_dt = sigma_h*Eh;
		return {dSh_dt, dEh_dt, dIh_dt, dRh_dt, dSv_dt, dEv_dt, dIv_dt, dInci_dt};
		
	}
 int bin_search(real x, int min_val, int max_val){
    // This assumes that min_val >= 0 is the minimum integer in range, 
    //  max_val > min_val,
    // and that x has already been rounded. 
    //  It should find the integer equivalent to x.
    int range = (max_val - min_val+1)/2; // We add 1 to make sure that truncation doesn't exclude a number
    int mid_pt = min_val + range;
    //print("mid_pt", mid_pt);
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1)/2; 
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range; 
        }
    }
    return out;
  }
}



data {
	// Structure
	int TotalPop;
	int n_months;       // total number of weeks
	real t0;          //starting time
	real ts[n_months];  // time bins
	int doprint;
	int inference;
	real L_h;
	real beta_h;
	real mu_h;
	real sigma_h;
	real gamma_h;
	real beta_v;
	real mu_v;
	real PIE;
	real temperature[n_months];
	
	
	// Data to fit
	// number of days with reported incidence
	int Reported_cases[n_months];
	
		// Priors
	real p_i0[2];
	//real p_t_min[2];
	real p_q[2];
	real p_Iv0[2];
	//real p_t_max[2];
	real p_wh[2];
	real p_delta[2];
	real p_q_birth[2];

	real p_q_kappa[2];

	//real p_Lv[2];
	//real p_tc[2];
	real p_phi;
	real p_q_mu_v[2];

	// Fixed parameters
	  real y0_vars[8];  
}
	

transformed data {
  real x_r[16+12];
  int x_i[1];

  x_i[1] = TotalPop;
  
  x_r[1] =  L_h;
  x_r[2] =  beta_h;
  x_r[3] =  mu_h;
  x_r[4] =  sigma_h;
  x_r[5] =  gamma_h;
  x_r[6] =  beta_v;
  x_r[7] =  mu_v;
  x_r[8] =  PIE;
  
  x_r[9:16] = y0_vars;
  x_r[16+1:16+12] = temperature;
}

parameters{
	// real<lower=1, upper=10000> IhR_0;   
	// real<lower=1, upper=100000> IhU_0;   
	real<lower=1, upper=20000> i0;  
	//real<lower=-0.01, upper=0.01> A; 
	//real<lower=-0.05, upper=0.5> B; 
	//real<lower=-1, upper=1> C; 
	real<lower=0, upper=10> q_birth;
	real<lower=0, upper=20> q; 
	real<lower=0, upper=20> q_kappa; //0.000111 in Mordecai
	real<lower=1, upper=50000> Iv0; 
	real<lower=0, upper=2> w_h; 
	real<lower=0, upper=2> delta; 
	real<lower=0, upper=10> q_mu_v;
	//real<lower=-1, upper=1> tc; 
	real<lower=0> phi; // 
}

transformed parameters {
	// transformed parameters
	real y[n_months, 8]; // raw ODE output
	real y0_temp[8];  // initial state
	real y0[8];  //
	
	// change of format for integrate_ode_rk45
	real theta[8];     // vector of parameters
	
	
	// ode outputs
	// cumulative incidence
	vector[n_months] Cum_inci; // overall case incidence by day
	
	// incidence
	vector[n_months] Inci; // overall case incidence by day
	
		
	// set up the initial conditions:
	y0_temp[1] = y0_vars[1] - (y0_vars[2]*i0 + y0_vars[3]*i0);
	y0_temp[2] = y0_vars[2]*i0;
	y0_temp[3] = y0_vars[3]*i0;
	y0_temp[4] = y0_vars[4];
	y0_temp[5] = 5000;  // Sv0 is assumed
	y0_temp[6] = Iv0/2; // Iv0 is fitted
	y0_temp[7] = Iv0/2;  //y0_vars[6];
	y0_temp[8] = y0_vars[8];
  //print("Sh0: ", y0_temp[1], "  Ih0: ", y0_temp[2], " Eh0: ", y0_temp[3], " Rh0: ", y0_temp[4]);
  //print("Sv0: ", y0_temp[5], "  Iv0: ", y0_temp[6], " Ev0: ", y0_temp[7]);
	
	y0 = to_array_1d(y0_temp);
	
	
	
	// change of format for integrate_ode_rk45
	theta[1:8] =  { i0, q, Iv0, w_h, delta, q_birth, q_kappa, q_mu_v};

	// // // run ODE solver
	//y = integrate_ode_rk45( SEIR, // ODE function
	//												y0, // initial states
	//												t0, // t0
	//												ts, // evaluation dates (ts)
	//												theta, // parameters
	//												x_r, // real data
	//												x_i, // integer data
	//												1.0E-6, 1.0E-6, 1.0E3); // tolerances and maximum steps
	 y = integrate_ode_rk45(SEIR, y0, t0, ts, theta, x_r, x_i); // tolerances and maximum steps
	
	// extract and format ODE results (1.0E-9 correction to avoid negative values due to unprecise estimates of zeros as tolerance is 1.0E-10)
	
	
	// Extracting incidence from the output of the ode solver	
	Cum_inci = to_vector(y[,8]) ;  
	
	Inci[1] =  Cum_inci[1]; 	  Inci[2:n_months] =  Cum_inci[2:n_months] - Cum_inci[1:(n_months-1)];
	
}

model {
	
	// priors
	i0 ~ lognormal(p_i0[1], p_i0[2]);
	//t_min ~ normal(p_t_min[1], p_t_min[2]);
	q ~ normal(p_q[1], p_q[2]);
	q_kappa ~ normal(p_q_kappa[1], p_q_kappa[2]);
	q_birth ~ normal(p_q_birth[1], p_q_birth[2]);
	Iv0 ~ lognormal(p_Iv0[1], p_Iv0[2]);
	//t_max ~ normal(p_t_max[1], p_t_max[2]);
	w_h ~ lognormal(p_wh[1], p_wh[2]);
	delta ~ lognormal(p_delta[1], p_delta[2]);
	//L_v ~ lognormal(p_Lv[1], p_Lv[2]);
	//tc ~ normal(p_tc[1], p_tc[2]);
	phi ~ exponential(p_phi);

	// debug
	if(doprint == 1) {
		//print("gamma: ", gamma);
		print("y[5,]: ",y[5,]);
		
	}

		if (inference!=0){
			Reported_cases ~ neg_binomial_2(Inci, 1/phi);
			
		}
}

generated quantities{
	real pred_cases[n_months];
	pred_cases = neg_binomial_2_rng(Inci, 1/phi);

}

