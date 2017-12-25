#ifndef MATHUTIL_H
#define MATHUTIL_H
#include <cmath>
#include <random>
const double pi = 3.14159265358979323846;

double truncate_num(double num, int places)
{
double result=0.0;
int sig = num>0? 1:-1;
unsigned int temp=(num*pow(10,places))*sig;
 
result=(((double)temp)/pow(10,places)*sig);
return result;

}
void get_cumulative_Poisson_distribution( double *F, int l,  double lambda)
{
  try
  {
    double fact=1;
    double elambda = exp(-1*lambda) ;
    F[0]= elambda;
    
    for ( int i = 1 ; i <= l ; ++i )
    {
      double t1=pow(lambda,i);
      fact=fact*i;
      F[i]=F[i-1]+(elambda*(t1/fact));
      //F[i]= truncate_num(F[i],2);
      
    }
    //std::cout<<F[l]<<std::endl;
    //int x=0;
    //std::cin>>x;
   }
catch (const char* msg) {
     std::cerr << msg << std::endl;
   }

}

void compute_prob_poisson_distribution(unsigned int nb_kmers,unsigned int read_coverage, unsigned int *kmer_coverage, double *prob_vals, unsigned int *flags)
{
 try
  { 
      
    for(unsigned int i=0;i<nb_kmers;i++)
    {
       bool flag=true;
       if(read_coverage > 1024)
         {
           read_coverage=1024;
           flags[i]=0;
           flag=false;
         }
       unsigned int kmer_coverage_value = kmer_coverage[i];
       if(kmer_coverage_value >read_coverage)
        {
            kmer_coverage_value = read_coverage;
            flags[i]=0;
            flag=false;
        }
       double *F = (double *)malloc( sizeof( double ) * (kmer_coverage_value+1)) ;
       get_cumulative_Poisson_distribution(F,kmer_coverage_value, read_coverage);
       prob_vals[i]=F[kmer_coverage_value];
       if(flag)
          {flags[i]=1;}
       //prob_vals[i]=truncate_num(F[kmer_coverage_value],2);
    }

   }
  catch (const char* msg) {
     std::cerr << msg << std::endl;
   }     
   
}
/*
void compute_prob_poisson_distribution(unsigned int nb_kmers,unsigned int read_coverage, unsigned int *kmer_coverage, double *prob_vals)
{
    
    for(unsigned int i=0;i<nb_kmers;i++)
    {
       unsigned int kmer_coverage_value = kmer_coverage[i];
       double *F = (double *)malloc( sizeof( double ) * (kmer_coverage_value+1)) ;
       F[0]=exp(-1*read_coverage) ;
       for(int j=1;j<=kmer_coverage_value;++j)
       {
         F[j]=F[j-1]+ pow(M_E, kmer_coverage_value * log(read_coverage) - read_coverage - lgamma(kmer_coverage_value + 1.0));
       }
       prob_vals[i]=F[kmer_coverage_value];
     
    }

        
   
}
*/
void compute_prob_poisson_distribution(unsigned int read_coverage, unsigned int kmer_coverage, double &prob_val)
{
   
     prob_val= pow(M_E, kmer_coverage * log(read_coverage) - read_coverage - lgamma(kmer_coverage + 1.0));
   
}
/*
void compute_prob_poisson_distribution(unsigned int nb_kmers,unsigned int read_coverage, unsigned int *kmer_coverage, double *prob_vals)
{
    
    for(unsigned int i=0;i<nb_kmers;i++)
        prob_vals[i] = pow(M_E, kmer_coverage[i] * log(read_coverage) - read_coverage - lgamma(kmer_coverage[i] + 1.0));
   
}
*/


#endif