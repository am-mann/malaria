// applicable to weekly cases data by symptoms onset date
functions {
  array[] real SEIR(real t,
                    array[] real y,
                    array[] real theta,
                    array[] real x_r,
                    array[] int x_i) {
    int TotalPop = x_i[1];
    real L_h = x_r[1];
    real beta_h = x_r[2];
    real mu_h = x_r[3];
    real sigma_h = x_r[4];
    real gamma_h = x_r[5];
    real beta_v = x_r[6];
    real mu_v = x_r[7];
    real PIE = x_r[8];

    array[12] real temp;
    temp[1:12] = x_r[16:27];

    array[7] real y0_vars;

    row_vector[7] dydt;

    real Nh;
    real b;
    real lambda_h;
    real lambda_v;

    real Sh = y[1];
    real Eh = y[2];
    real Ih = y[3];
    real Rh = y[4];
    real Sv = y[5];
    real Iv = y[6];
    real CumInci = y[7];

    real dSh_dt;
    real dEh_dt;
    real dIh_dt;
    real dRh_dt;
    real dSv_dt;
    real dIv_dt;
    real dInci_dt;

    real i0 = theta[1];
    real q = theta[2];
    real Iv0 = theta[3];
    real w_h = theta[4];
    real delta = theta[5];
    real birth = theta[6];
    real q_mu = theta[7];

    y0_vars = x_r[9:15];
    int temp_index = bin_search(floor(t), 0, 11) + 1;
    if (temp_index > 11) temp_index = 11;

    real T = (t - floor(t)) * temp[temp_index + 1] + (1 - (t - floor(t))) * temp[temp_index];

    Nh = Sh + Eh + Ih + Rh;
    b = 30 * q * T * (T - 13.8) * pow(40 - T, 0.5);
    real ro = q_mu * (0.0011*T^2 - 0.049*T + 0.609);
    mu_v = 30 * ro;
    lambda_h = b * beta_h * Iv / Nh;
    lambda_v = b * beta_v * Ih / Nh;

    dSh_dt = L_h - lambda_h * Sh + w_h * Rh - mu_h * Sh;
    dEh_dt = lambda_h * Sh - sigma_h * Eh - mu_h * Eh;
    dIh_dt = sigma_h * Eh - gamma_h * Ih - (mu_h + delta) * Ih;
    dRh_dt = gamma_h * Ih - w_h * Rh - mu_h * Rh;
    dSv_dt = birth - lambda_v * Sv - mu_v * Sv;
    dIv_dt = lambda_v * Sv - mu_v * Iv;
    dInci_dt = sigma_h * Eh;
    return {dSh_dt, dEh_dt, dIh_dt, dRh_dt, dSv_dt, dIv_dt, dInci_dt};
  }

  int bin_search(real x, int min_val, int max_val) {
    int range = (max_val - min_val + 1) / 2;
    int mid_pt = min_val + range;
    int out;
    while (range > 0) {
      if (x == mid_pt) {
        out = mid_pt;
        range = 0;
      } else {
        range = (range + 1) / 2;
        mid_pt = x > mid_pt ? mid_pt + range : mid_pt - range;
      }
    }
    return out;
  }
}


data {
  // Structure
  int TotalPop;
  int n_months;       // total number of weeks
  real t0;            // starting time
  array[n_months] real ts;  // time bins
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
  array[n_months] real temperature;

  // Data to fit
  array[n_months] int Reported_cases;

  // Priors
  array[2] real p_i0;
  array[2] real p_q;
  array[2] real p_Iv0;
  array[2] real p_wh;
  array[2] real p_delta;
  array[2] real p_birth;
  array[2] real p_q_mu;
  real p_phi;

  // Fixed parameters
  array[7] real y0_vars;
}

transformed data {
  array[27] real x_r;
  array[1] int x_i;

  x_i[1] = TotalPop;

  x_r[1] = L_h;
  x_r[2] = beta_h;
  x_r[3] = mu_h;
  x_r[4] = sigma_h;
  x_r[5] = gamma_h;
  x_r[6] = beta_v;
  x_r[7] = mu_v;
  x_r[8] = PIE;

  x_r[9:15] = y0_vars;
  x_r[16:27] = temperature;
}


parameters{
	// real<lower=1, upper=10000> IhR_0;   
	// real<lower=1, upper=100000> IhU_0;   
	real<lower=1, upper=10000> i0;  
	real<lower=0, upper=10> q_mu; 
	//real<lower=0, upper=0.02> q_MDR;
	//real<lower=5, upper=20> t_min_MDR; 
	//real<lower=10, upper=35> t_max_MDR; 
	real<lower=0, upper=1> q; 
	real<lower=1, upper=15000> birth; 
	//real<lower=0, upper=0.01> q_kappa; //0.000111 in Mordecai
	real<lower=1, upper=50000> Iv0; 
	real<lower=0, upper=2> w_h; 
	real<lower=0, upper=2> delta; 
	//real<lower=0.01, upper=5> q_death; 
	//real<lower=1000, upper=30000> birth;
	//real<lower=-1, upper=1> tc; 
	real<lower=0> phi; // 
}

transformed parameters {
  array[n_months, 7] real y; // raw ODE output
  array[7] real y0_temp;  // initial state
  array[7] real y0;        //

  // change of format for integrate_ode_rk45
  array[7] real theta;     // vector of parameters

  // ode outputs
  vector[n_months] Cum_inci; // overall case incidence by day
  vector[n_months] Inci;     // overall case incidence by day

  // set up the initial conditions:
  y0_temp[1] = y0_vars[1] - (y0_vars[2] * i0 + y0_vars[3] * i0);
  y0_temp[2] = y0_vars[2] * i0;
  y0_temp[3] = y0_vars[3] * i0;
  y0_temp[4] = y0_vars[4];
  y0_temp[5] = 5000;        // Sv0 is assumed
  y0_temp[6] = Iv0 / 2;     // Iv0 is fitted
  y0_temp[7] = y0_vars[7];

  y0 = to_array_1d(y0_temp);

  // change of format for integrate_ode_rk45
  theta = { i0, q, Iv0, w_h, delta, birth, q_mu};

  y = integrate_ode_rk45(SEIR, y0, t0, ts, theta, x_r, x_i);

  // Extracting incidence from the output of the ode solver
  Cum_inci = to_vector(y[, 7]);
  Inci[1] = Cum_inci[1];
  Inci[2:n_months] = Cum_inci[2:n_months] - Cum_inci[1:(n_months - 1)];
}

model {
	
	// priors
	i0 ~ lognormal(p_i0[1], p_i0[2]);
	//A ~ normal(p_A[1], p_A[2]);
	//B ~ normal(p_B[1], p_B[2]);
	//C ~ normal(p_C[1], p_C[2]);
    q_mu  ~ normal(p_q_mu[1], p_q_mu[2]);
	//q_kappa ~ normal(p_q_kappa[1], p_q_kappa[2]);
	//q_MDR ~ normal(p_q_MDR[1], p_q_MDR[2]);
	//t_min_MDR ~ normal(p_t_min_MDR[1], p_t_min_MDR[2]);
	//t_max_MDR ~ normal(p_t_max_MDR[1], p_t_max_MDR[2]);
	Iv0 ~ lognormal(p_Iv0[1], p_Iv0[2]);
	//t_max ~ normal(p_t_max[1], p_t_max[2]);
	w_h ~ lognormal(p_wh[1], p_wh[2]);
	delta ~ lognormal(p_delta[1], p_delta[2]);
	birth ~ lognormal(p_birth[1], p_birth[2]);
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
	array[n_months] real pred_cases;
	pred_cases = neg_binomial_2_rng(Inci, 1/phi);

}

