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

void compute_prob_poisson_distribution(unsigned int nb_kmers,int read_coverage, int *kmer_coverage, double *prob_vals, unsigned int *flags,std::string ns)
{
 try
  { 
    bool zeroflag=false;
    bool nonzeroflag=false;
    bool adjusted=false;
    int start_non_zero_index=0;
    int end_non_zero_index=0;  
    int start_zero_index=0;
    int end_zero_index=0; 
    if(read_coverage > 128)// 1024)
    {
           read_coverage= 128;//1024;
           adjusted=true;
    }
    double *median_prob = (double *)malloc( sizeof( double ) * (read_coverage+1)) ;
    get_cumulative_Poisson_distribution(median_prob,read_coverage, read_coverage);   
    
    for(unsigned int i=0;i<nb_kmers;i++)
    {
       bool flag=true;
       if(adjusted)
        {
           //read_coverage= 128;//1024;
           flags[i]=0;
           flag=false;
        }
       int kmer_coverage_value = kmer_coverage[i];
       if(kmer_coverage_value >read_coverage)
       {
            //kmer_coverage_value = read_coverage;
            flags[i]=0;
            flag=false;
            prob_vals[i]=median_prob[read_coverage];
          /*if(std::isinf(prob_vals[i]))
            {
            std::cout<<prob_vals[i]<<std::endl;
            int y=0;
            std::cin>>y;}*/            
            if(!nonzeroflag)
                start_non_zero_index=i;
            else
                end_non_zero_index=i;
            nonzeroflag=true;
            continue;
       }
       if(kmer_coverage_value==-1)
       { 
            flags[i]=0;
            flag=false;
            prob_vals[i]=1;//previously was 0
            if(!nonzeroflag)
                start_non_zero_index=i;
            else
                end_non_zero_index=i;
            nonzeroflag=true;
            continue;

       }
       if(kmer_coverage_value==0)
       {
            flags[i]=0;
            flag=false;
           
           if(!zeroflag)
             start_zero_index=i;
           else
             end_zero_index=i;
           prob_vals[i]=0;
             zeroflag=true;
           continue;

        }
       
       double *F = (double *)malloc( sizeof( double ) * (kmer_coverage_value+1)) ;
       get_cumulative_Poisson_distribution(F,kmer_coverage_value, read_coverage);
       prob_vals[i]=F[kmer_coverage_value];
       if(!nonzeroflag)
          start_non_zero_index=i;
       else
          end_non_zero_index=i;
       nonzeroflag=true;
       if(flag)
          {flags[i]=1;}
       //prob_vals[i]=truncate_num(F[kmer_coverage_value],2);
       

    }
   double correct_prob=0.0; 
  
   if(zeroflag)
   {
    /*
    for(unsigned int i=0;i<nb_kmers;i++)
    {
      std::cout<<prob_vals[i]<<",";

    }
   std::cout<<std::endl;
   std::cout<<std::endl;
   std::cout<<"start zero"<<start_zero_index<<std::endl;
   std::cout<<"end zero"<<end_zero_index<<std::endl;
   std::cout<<"start non zero"<<start_non_zero_index<<std::endl;
   std::cout<<"end non zero"<<end_non_zero_index<<std::endl;*/
   
   for(unsigned int i=end_zero_index;i>start_zero_index;i--)
   {
      if(i >end_non_zero_index)
      {continue;}
      if(i<start_non_zero_index)
      {continue;} 
      if (prob_vals[i]!=0)
          {
               correct_prob=prob_vals[i];
               continue;
           }
      else
          {
               correct_prob=prob_vals[i+1];
          }

      if(ns=="keepValN")
      prob_vals[i]=correct_prob;
      else if(ns=="keepN")
      prob_vals[i]=1;


   }

  /*
   std::cout<<std::endl;
   std::cout<<std::endl;
   for(unsigned int i=0;i<nb_kmers;i++)
    {
      std::cout<<prob_vals[i]<<",";

    }
     std::cout<<std::endl;
     int x=0;
     std::cin>>x;*/
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
void get_cumulative_binomial_distribution( std::vector<double> & F, int l, double p )
{
    // p is the probability of getting 1.
    int i ;
    double coef = 1 ;
    double exp = pow( 1 - p, l ) ;
    F[0] = pow( 1 - p, l ) ;
    
    for ( i = 1 ; i <= l ; ++i )
    {
        coef = coef / i * ( l - i + 1 ) ;
        exp =  exp / ( 1 - p ) * p ;
        F[i] = F[i - 1] + coef * exp ;
    }
    
}


#endif
