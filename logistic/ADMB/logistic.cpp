#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <logistic.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  num_years.allocate("num_years");
  catches.allocate(1,num_years,"catches");
  num_obs.allocate("num_obs");
  log_pop_obs.allocate(1,num_obs,"log_pop_obs");
  years_obs.allocate(1,num_obs,"years_obs");
  sd_obs.allocate("sd_obs");
  check.allocate("check");
}

void model_parameters::initializationfunction(void)
{
  logr.set_initial_value(-2.302585);
  logK.set_initial_value(8.517193);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  logr.allocate(-4.60517,-0.6931472,"logr");
  logK.allocate(6.907755,9.615805,"logK");
  r.allocate("r");
  #ifndef NO_AD_INITIALIZE
  r.initialize();
  #endif
  K.allocate("K");
  #ifndef NO_AD_INITIALIZE
  K.initialize();
  #endif
  N.allocate(1,num_years,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  temp.allocate("temp");
  #ifndef NO_AD_INITIALIZE
  temp.initialize();
  #endif
  NLL.allocate("NLL");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  sd_K.allocate("sd_K");
}

void model_parameters::userfunction(void)
{
  NLL =0.0;
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
