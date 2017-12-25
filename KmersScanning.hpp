#ifndef KMERSCANNING_H
#define KMERSCANNING_H
#include "KmerUtil.hpp"
#include "MathUtil.hpp"
#include "Utility.hpp"
#include "ReadsParsing.hpp"
#include <string>
#include <vector>
#include <algorithm>
#include <typeinfo>
#include <fstream>


void help_me()
{

       std::cout<<"--- LightTrimmer can not process your dataset !!! "<<std::endl;

       std::cout<<"--- maximum supported read length for this version = "<<MAX_READ_LENGTH<<std::endl;

       std::cout<<"--- try different values for k [kmer size] & g [gap size] or different dataset"<<std::endl;

       std::cout<<std::endl;
       exit(1);


}

int compare_counts (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}


void pass_one_reads_trimming(std::string read,kmersUtil &kmers_util,double *prob_vals, unsigned int *count_vals, unsigned int *flags, unsigned int &val)
{
    int i;
    unsigned int count=0;
    unsigned int nb_kmers= read.length()-kmer_size + 1;
    unsigned int median=0;
    unsigned int *counts = (unsigned int *)malloc( sizeof( unsigned int ) * (nb_kmers)) ;
    unsigned int *counts_tmp = (unsigned int *)malloc( sizeof( unsigned int ) * (nb_kmers)) ;
        
    kmercode_length kmercode=0,kmercode_can=0;
    for(i=0; i<kmer_size; ++i)
    {
      kmercode=kmercode*4+kmers_util.nt2int(read[i]);
    }
      kmercode_can=kmers_util.get_canonical_kmer_code(kmercode);
      count = kmers_util.get_from_table(kmercode_can);
      count_vals[0]=count;
      //if(count !=0)
      counts[0]=count; 
      counts_tmp[0]=count;
    for(i=1; i<read.length()-kmer_size+1; ++i)
    {
        kmercode=(kmercode*4+kmers_util.nt2int(read[i+(kmer_size-1)])& kmer_mask);
        kmercode_can=kmers_util.get_canonical_kmer_code(kmercode);
        count = kmers_util.get_from_table(kmercode_can);
    //  if(count !=0)
        counts[i]=count;
        counts_tmp[i]=count;
        count_vals[i]=count;
    }

    qsort( counts,nb_kmers, sizeof( unsigned int ), compare_counts) ;
    
    if(((nb_kmers+1)%2)==0)
     {
         int median_index= ((nb_kmers+1)/2);
         median = counts[median_index-1];
     }
     else
     {
        int median_index= ((nb_kmers+1)/2);
            median = round((counts[median_index-1]+counts[median_index])/2);
     }
          
    val=median;
    compute_prob_poisson_distribution(nb_kmers,median,counts_tmp,prob_vals,flags);

}//end_method
void init(const std::vector<std::string> read_files, std::string counting_table_name,int nb_threads,bool verbose)
{

try
{
//************************** Step 0: preprocessing the kmers counting table **************************************************
    std::cout<<"--- Start Reading your kmers counting file. "<<std::endl;
    FILE *countingtable = NULL ;
    kmersUtil kmers_util;
    int nb_kmers=0;
    countingtable = fopen(counting_table_name.c_str(), "r" ) ;
    if ( countingtable == NULL )
    {
	std::cerr<<"--- could not open file "<<counting_table_name<<std::endl;
	exit( 1 ) ;
    }
    char buffer[100] ;
    int i=0;
    while ( fscanf( countingtable, "%s", buffer ) != EOF )
    {
        unsigned int count = atoi( &buffer[1] ) ;
	fscanf( countingtable, "%s", buffer ) ;
        kmercode_length kmercode=0,kmercode_can=0;
        for(i=0;buffer[i];i++)
           {
              kmercode=kmercode*4+kmers_util.nt2int(buffer[i]);
           }
        if(i!=kmer_size)
           {
             	std::cerr<<"--- counting module uses different kmer size, k = "<<i<<std::endl;
	        exit( 1 ) ;
           }
        kmercode_can=kmers_util.get_canonical_kmer_code(kmercode);
        kmers_util.add2table(kmercode_can,count); 
        nb_kmers++;

    }
     
   fclose( countingtable ) ;

   if(verbose)
    {
              std::cout<<"--- total number of counting kmers = "<<nb_kmers<< std::endl;       

    }

//************************** Step 1: preprocessing the kmers counting table **************************************************
    
    reads_parsing reads(read_files);
    uint64_t total_bases=0;uint64_t total_reads=0;
    std::string read_seq="";int read_length=0;
    uint64_t average_len=0;
    double gapped_kmer=0.0;
    std::ofstream kmers_prob_file("kmers_prob.txt",std::ios::out);
    std::ofstream kmers_count_file("kmers_count.txt",std::ios::out);
    std::ofstream kmers_correct_file("kmers_correct.txt",std::ios::out);
    std::ofstream kmers_all_file("kmers_info_all.txt",std::ios::out);
    std::cout<<"--- Start Reading your dataset. "<<std::endl;
    std::cout<<std::endl;
    start_time(1);
    if(nb_threads==1)
    {
          if (kmers_prob_file.is_open()&& kmers_count_file.is_open()&&kmers_correct_file.is_open()&&kmers_all_file.is_open())
           {
              kmers_all_file<<"R"<<"  "<<"K"<<"  "<<"C"<<"  "<<"M"<<"  "<<"P"<<"  "<<"Co"<<"  "<<"Ca"<<std::endl;
              while(reads.get_next_sequence(read_seq,read_length))
               {        

                  total_reads++;
                  total_bases+=read_seq.length();
                  std::vector<std::string> reads_without_ns;
                  kmers_util.get_reads_spiliting_around_Ns(read_seq,read_length,reads_without_ns);
                 
                  for(int i=0;i<reads_without_ns.size();i++)
                  {
                      double *prob_vals = (double *)malloc( sizeof( double ) * (reads_without_ns[i].length()-kmer_size+1)) ;
                      unsigned int *count_vals = (unsigned int *)malloc( sizeof( unsigned int) * (reads_without_ns[i].length()-kmer_size+1)) ;
                      unsigned int *flags = (unsigned int *)malloc( sizeof( unsigned int) * (reads_without_ns[i].length()-kmer_size+1)) ;
                      unsigned int val=0;
                      pass_one_reads_trimming(reads_without_ns[i], kmers_util,prob_vals,count_vals,flags,val); 
                      
                      for(int j=0;j<reads_without_ns[i].length()-kmer_size+1;j++)
                         {

                            std::string kmer_seq=reads_without_ns[i].substr(j,kmer_size);
                            bool flag=true;
                            for (int l=0;l<kmer_seq.length();l++)
                                {
                                   if(islower(kmer_seq[l]))
                                      flag=false;
                                      break;
                                  
                                }
                            kmers_prob_file <<prob_vals[j]<<",";
                            kmers_count_file <<count_vals[j]<<",";

                           if(flag)
                            {
                               kmers_correct_file <<"("<<j<<","<<"1"<<")"<<",";
                               kmers_all_file<<total_reads<<"  "<<j<<"  "<<count_vals[j]<<"  "<<val<<"  "<<prob_vals[j]<<"  "<<"1 "<<"  "<<flags[j]<<std::endl;
                            }
                           else
                            {
                               kmers_correct_file <<"("<<j<<","<<"0"<<")"<<",";
                               kmers_all_file<<total_reads<<"  "<<j<<"  "<<count_vals[j]<<"  "<<prob_vals[j]<<"  "<<"0 "<<"  "<<flags[j]<<std::endl;

                            }
                            
                       	
                         }
                         
                  }
                  kmers_prob_file <<std::endl;
                  kmers_count_file<<std::endl;
                  kmers_correct_file<<std::endl;
                  
                  
                 
               }             

             kmers_prob_file.close();
             kmers_count_file.close();
             kmers_correct_file.close();
             kmers_all_file.close();
          }
          else
          {
             std::cerr << "--- can't open kmers probability file: kmers_prob.txt"<< std::endl;
             std::cerr << "--- can't open kmers count file: kmers_count.txt"<< std::endl;
             std::cerr << "--- can't open kmers correct file: kmers_correct.txt"<< std::endl;
             std::cerr << "--- can't open kmers correct file: kmers_info_all.txt"<< std::endl;
             exit(1);
          }

    
    }

    end_time(1);
    std::vector<double> untrust(kmer_size+1);
    std::vector<int> threshold(kmer_size+1);
    if(total_reads >0)
    {
              gapped_kmer= static_cast<double>(1.0/static_cast<double>(gap_size));
              average_len=total_bases/total_reads;
              if((total_bases%total_reads)>(total_reads/2))
              average_len++;
    }

   if(verbose)
    {
              std::cout<<"--- average read length = "<<average_len<< std::endl;
              std::cout<<"--- total number of reads ="<<total_reads<<std::endl;       

    }
}
   
catch (const char* msg) {
     std::cerr << msg << std::endl;
   }
}
#endif
