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
#include <map>

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

// Dust scoring scheme as given by:
// Morgulis A. "A fast and symmetric DUST implementation to Mask
// Low-Complexity DNA Sequences". J Comp Bio.
double calculateDustScore(const std::string& seq)
{
    std::map<std::string, int> scoreMap;
    
    // Cannot calculate dust scores on very short reads
    if(seq.size() < 3)
        return 0.0f;

    // Slide a 3-mer window over the sequence and insert the sequences into the map
    for(size_t i = 0; i < seq.size() - 3; ++i)
    {
        std::string triMer = seq.substr(i, 3);
        scoreMap[triMer]++;
    }

    // Calculate the score by summing the square of every element in the map
    double sum = 0;
    std::map<std::string, int>::iterator iter = scoreMap.begin();
    for(; iter != scoreMap.end(); ++iter)
    {
        int tc = iter->second;
        double score = (double)(tc * (tc - 1)) / 2.0f;
        sum += score;
    }
    return sum / (seq.size() - 2);
}

// Returns the window over seq with the highest dust score
double maxDustWindow(const std::string& seq, size_t windowSize=64, size_t minWindow=64)
{
    double maxScore = 0.0f;
    for(size_t i = 0; i < seq.size(); i += 1)
    {
        size_t r = seq.size() - i;
        size_t w = r < windowSize ? r : windowSize;
        if(w >= minWindow) // don't calculate score for small windows
        {
            double s = calculateDustScore(seq.substr(i, w));
            if(s > maxScore)
                maxScore = s;
        }
    }
    return maxScore;
}


void pass_one_reads_trimming(std::string read,kmersUtil &kmers_util,double *prob_vals, int *count_vals, unsigned int *flags, int &val, unsigned int &val2,double dust, std::string ns)
{
    
    int i;
    int count=0;
    unsigned int nb_kmers=read.length()-kmer_size + 1;
    unsigned int nb_kmers1=0;
             int median=0;
    int *counts = (int *)malloc( sizeof( int ) * (nb_kmers)) ;
    int *counts_tmp = (int *)malloc( sizeof( int ) * (nb_kmers)) ;
    bool flag =false;
    int encode_char=0;    
    kmercode_length kmercode=0,kmercode_can=0;
    for(i=0; i<kmer_size; ++i)
    {
      encode_char=kmers_util.nt2int(read[i]);
      if(encode_char==-1)
        {
           encode_char=0;
           flag=true;
        }
     

      kmercode=kmercode*4+encode_char;

    }

   if(flag) //kmer has N's
     {
          count=0;
          flag=false;
      }
   else
     {
         std::string kmer_seq=read.substr(0,kmer_size);
         double globalDustScore = maxDustWindow(kmer_seq,kmer_size,kmer_size);
         if(globalDustScore > static_cast<double>(dust))
            count=-1;
         else 
            {
                kmercode_can=kmers_util.get_canonical_kmer_code(kmercode);
                count = kmers_util.get_from_table(kmercode_can);
            }
   
     }
    count_vals[0]=count;
    if(count >1)
      {
         counts[nb_kmers1]=count;
         nb_kmers1++;
      } 
    counts_tmp[0]=count;
    for(i=1; i<read.length()-kmer_size+1; ++i)
    {
      encode_char=kmers_util.nt2int(read[i+(kmer_size-1)]);
      if(encode_char==-1)
        {  
           encode_char=0;
           flag=true;
        }

      kmercode=(kmercode*4+encode_char& kmer_mask);
     
      if(flag) //kmer has N's
      {
          count=0;
          flag=false;
      }
     else
     {
         std::string kmer_seq=read.substr(i,kmer_size);
         double globalDustScore = maxDustWindow(kmer_seq,kmer_size,kmer_size);
         if(globalDustScore > static_cast<double>(dust))
            count=-1;
         else
            { 
                kmercode_can=kmers_util.get_canonical_kmer_code(kmercode);
                count = kmers_util.get_from_table(kmercode_can);
            }
    
     }

     if(count >1)
      {
            counts[nb_kmers1]=count;
            nb_kmers1++;

      }
        counts_tmp[i]=count;
        count_vals[i]=count;
    }
   
  if(nb_kmers1!=0)
  { 
    qsort(counts,nb_kmers1, sizeof( unsigned int ), compare_counts) ;
    
    if(((nb_kmers1+1)%2)==0)
     {
         int median_index= ((nb_kmers1+1)/2);
         median = counts[median_index-1];
     }
     else
     {
        int median_index= ((nb_kmers1+1)/2);
            median = round((counts[median_index-1]+counts[median_index])/2);
     }

    val=median;
    val2=nb_kmers1; 
    compute_prob_poisson_distribution(nb_kmers,median,counts_tmp,prob_vals,flags,ns);
  }
}//end_method
bool poly_as_ts_reads_trimming(std::string read,int poly,std::string &trimmed_read,kmersUtil &kmers_util)
{
   int poly_a=0;
   int poly_t=0;
   int encode_char=0;
   int index=0;
   bool flag=false;
   for(int i=0;i<read.length();i++)
   {
        encode_char=kmers_util.nt2int(read[i]);
        if(encode_char==0)
         { 
           if(i==0)
           poly_a++;
           else
           if(kmers_util.nt2int(read[i-1])==0)
           poly_a++;
           else
             {
                index=i;
                break;

             }
         }
        else
        if(encode_char==3)
          {
            if(i==0)
            poly_t++;
            else
            if(kmers_util.nt2int(read[i-1])==3)
            poly_t++;
            else
             {
                index=i;
                break;
             }
          }
        else
            {
               index=i;
               break;
            }
       
   }//for
 
   if((poly_a==read.length())||(poly_t==read.length()))
   {
           trimmed_read="";
           flag=true;
           return flag;
   }
 
   if(poly_a>=poly)
   {
       trimmed_read=read.substr(index);
       flag=true;    
   }
   else
   if(poly_t>=poly)
   {
       trimmed_read=read.substr(index);
       flag=true;

   }

 if(flag==false)
 {
   for(int i=read.length()-1;i>=0;i--)
   {
        encode_char=kmers_util.nt2int(read[i]);
        if(encode_char==0)
          { 
            if(i==read.length()-1)
               poly_a++; 
            else
            if(kmers_util.nt2int(read[i+1])==0) 
             poly_a++;
            else
             {
               index=i;
               break;


             }
          }
        else
        if(encode_char==3)
          {
            if(i==read.length()-1)
               poly_t++;
             else
             if(kmers_util.nt2int(read[i+1])==3)
             poly_t++;
             else
             {
               index=i;
               break;

             }
          }
        else
            {
              index=i;
              break;
            }
   }//for

   if(poly_a>=poly)
   {
          trimmed_read=read.substr(0,index+1);
          flag=true;
   }
   else
   if(poly_t>=poly)
   {
          trimmed_read=read.substr(0,index+1);
          flag=true;
   }
   
   }//if
  else
   {
       poly_a=0;
       poly_t=0;
       index=0;
       for(int i=trimmed_read.length()-1;i>=0;i--)
       {
         encode_char=kmers_util.nt2int(trimmed_read[i]);
         if(encode_char==0)
          {
            if(i==trimmed_read.length()-1)
               poly_a++; 
            else
            if(kmers_util.nt2int(trimmed_read[i+1])==0)
             poly_a++;
            else
             {
               index=i;
               break;


             }
          }
        else
        if(encode_char==3)
          {
            if(i==trimmed_read.length()-1)
               poly_t++;
             else
            if(kmers_util.nt2int(trimmed_read[i+1])==3)
             poly_t++;
             else
             {
               index=i;
               break;

             }
          }
        else
            {
              index=i;
              break;
            }

       }//for
     
   if((poly_a==trimmed_read.length())||(poly_t==trimmed_read.length()))
   {
           trimmed_read="";
           //flag=true;
           return flag;
   }
   if(poly_a>=poly)
   {
          trimmed_read=trimmed_read.substr(0,index+1);
          //flag=true;
   }
   else
   if(poly_t>=poly)
   {
          trimmed_read=trimmed_read.substr(0,index+1);
          //flag=true;
   }

   }//else

   if(flag==false)
      trimmed_read=read;

   return flag;
}

void init(const std::vector<std::string> read_files, std::string counting_table_name,double dust, int poly,std::string ns,int window, bool verbose)
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
    uint64_t total_bases=0;uint64_t total_reads=0;uint64_t total_trimmed_reads=0;uint64_t total_lowcomplex_reads=0;
    std::string read_seq="";int read_length=0;
    uint64_t average_len=0;
    std::ofstream kmers_prob_file("kmers_prob.txt",std::ios::out);
    std::ofstream kmers_prob_file_C("kmers_prob_C.txt",std::ios::out);
    std::ofstream kmers_prob_file_I("kmers_prob_I.txt",std::ios::out);
    std::ofstream kmers_count_file("kmers_count.txt",std::ios::out);
    std::ofstream kmers_correct_file("kmers_correct.txt",std::ios::out);
    std::ofstream kmers_all_file("kmers_info_all.txt",std::ios::out);
    std::cout<<"--- Start Reading your dataset. "<<std::endl;
    std::cout<<std::endl;
    start_time(1);
    if (kmers_prob_file.is_open()&& kmers_count_file.is_open()&&kmers_correct_file.is_open()&&kmers_all_file.is_open()&&kmers_prob_file_C.is_open()&&kmers_prob_file_I.is_open())
           {
              kmers_all_file<<"R"<<","<<"K"<<","<<"C"<<","<<"M"<<","<<"P"<<","<<"Co"<<","<<"Ca"<<","<<"Cbp"<<std::endl;
              while(reads.get_next_sequence(read_seq,read_length))
               {        
                  
                  total_reads++;
                  total_bases+=read_seq.length();
                  bool flag=false;
                  std::string trimmed_read="";
             if(poly!=-1)
               {
                  //*********** (Polys A|T) *******************************************

                  flag=poly_as_ts_reads_trimming(read_seq,poly,trimmed_read,kmers_util);
                  if(flag)
                  { total_trimmed_reads++; }
                  //read_length=trimmed_read.length();

                }
              else
               {
                   trimmed_read=read_seq;

               }
                  read_length=trimmed_read.length();
                  if(read_length<kmer_size)
                   { continue;}
              
              if(dust!=-1)
               {
                  //********** (filter low complexity reads)***************************
                   if(window==0)
                      window=read_length;

                   double globalDustScore = maxDustWindow(trimmed_read,window,window);
                   if(globalDustScore >static_cast<double>( dust))
                   {
                       total_lowcomplex_reads++;
                       continue;

                   }
                }

                  double *prob_vals = (double *)malloc( sizeof( double ) * (read_length-kmer_size+1)) ;
                  int *count_vals = (int *)malloc( sizeof(int) * (read_length-kmer_size+1)) ;
                  unsigned int *flags = (unsigned int *)malloc( sizeof( unsigned int) * (read_length-kmer_size+1)) ;
                  int val=0;
                  unsigned int nb_kmers_used=0;
                  pass_one_reads_trimming(trimmed_read, kmers_util,prob_vals,count_vals,flags,val,nb_kmers_used,10.0,ns); 
                  for(int j=0;j<read_length-kmer_size+1;j++)
                    {

                            std::string kmer_seq=trimmed_read.substr(j,kmer_size);
                            int kc=1;
                            int bp=1;
                            if(islower(trimmed_read[j]))
                                bp=0;
                            if(trimmed_read[j]=='n'||trimmed_read[j]=='N')
                                bp=-1;
                            for (int l=0;l<kmer_seq.length();l++)
                                {

                                    if(islower(kmer_seq[l]))
                                        {
                                           kc=0;
                                           break;
                                        }                                     
                                    
                                  
                                }
                             for (int l=0;l<kmer_seq.length();l++)
                                {

                                   if(kmer_seq[l]=='n'||kmer_seq[l]=='N')
                                     {
                                        kc=-1;
                                        break;

                                     }
                                }   
                            kmers_prob_file <<prob_vals[j]<<",";
                            kmers_count_file <<count_vals[j]<<",";
                            kmers_correct_file <<"("<<j<<","<<kc<<")"<<",";
                            kmers_all_file<<total_reads<<","<<j<<","<<count_vals[j]<<","<<val<<","<<prob_vals[j]<<","<<kc<<","<<flags[j]<<","<<bp<<std::endl;
                            if(kc==1)
                            {kmers_prob_file_C<<prob_vals[j]<<",";kmers_prob_file_I<<"-1"<<",";}
                            if(kc==0)
                            {kmers_prob_file_I<<prob_vals[j]<<",";kmers_prob_file_C<<"-1"<<",";}
                            
                       	
                         }
                         
                  
                  kmers_prob_file <<std::endl;
                  kmers_count_file<<std::endl;
                  kmers_correct_file<<std::endl;
                  kmers_prob_file_C<<std::endl;
                  kmers_prob_file_I<<std::endl;
                  
                 
               }             

             kmers_prob_file.close();
             kmers_count_file.close();
             kmers_correct_file.close();
             kmers_all_file.close();
             kmers_prob_file_C.close();
             kmers_prob_file_I.close();
          }
          else
          {
             std::cerr << "--- can't open kmers probability file: kmers_prob.txt"<< std::endl;
             std::cerr << "--- can't open kmers count file: kmers_count.txt"<< std::endl;
             std::cerr << "--- can't open kmers correct file: kmers_correct.txt"<< std::endl;
             std::cerr << "--- can't open kmers correct file: kmers_info_all.txt"<< std::endl;
             exit(1);
          }

    end_time(1);
    std::vector<double> untrust(kmer_size+1);
    std::vector<int> threshold(kmer_size+1);
    if(total_reads >0)
    {
           
              average_len=total_bases/total_reads;
              if((total_bases%total_reads)>(total_reads/2))
              average_len++;
    }

   if(verbose)
    {
              std::cout<<"--- average read length = "<<average_len<< std::endl;
              std::cout<<"--- total number of reads ="<<total_reads<<std::endl;
              std::cout<<"--- total number of poly(A|T)'s trimmed reads = "<<total_trimmed_reads<<std::endl;       
              std::cout<<"--- total number of low complexity trimmed reads= "<<total_lowcomplex_reads<<std::endl;
    }
}
   
catch (const char* msg) {
     std::cerr << msg << std::endl;
   }
}
#endif
