// Me playing with logistic model as practice. CCM 6/19/2012


DATA_SECTION
  init_int num_years
  init_vector catches(1,num_years)
  init_int num_obs
  init_vector log_pop_obs(1,num_obs)
  init_vector years_obs(1,num_obs)
  init_number sd_obs
  init_int check
  // !! cout << check << endl;
  // !! cout << num_years << endl;
  int j
  
INITIALIZATION_SECTION
  logr -2.302585
  logK 8.517193

PARAMETER_SECTION
  init_bounded_number logr(-4.60517,-0.6931472)
  init_bounded_number logK(6.907755,9.615805)
  number r
  number K 
  vector N(1,num_years)
  number fpen
  number temp 			// used in intermediate calc
  objective_function_value NLL
  sdreport_number sd_K

PROCEDURE_SECTION
  r=mfexp(logr);
  K=mfexp(logK);
  j=1;
  NLL=0;
  N[1]=K;
  sd_K=K;
  fpen=0.0;
  for (int i=2; i<=num_years; i++)
  {
  N(i)=N(i-1)+N(i-1)*r*(1-N(i-1)/K)-catches(i-1);
  N(i)=posfun(N(i), 10, fpen);
  NLL+=100*fpen;
  }
 for(int j=1; j<=num_obs;j++)
    {
    temp=square(log_pop_obs(j)-log(N(years_obs(j))))/(2*square(sd_obs));
    NLL+=temp;
    }


REPORT_SECTION
//   cout << "NLL is " << NLL << endl;
//   cout << "r is " << r << endl;
//   cout << "K is " << K << endl;
//   cout << "fpen is " << fpen << endl;

