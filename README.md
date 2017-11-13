# LightTrimmer
Light Weight Trimming Algorithm for RNA-Seq reads

#### System requirements 
64-bit machine with g++ version 4.7 or higher, [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads),and [zlib](http://en.wikipedia.org/wiki/Zlib) libraries.

#### Installation 
1. Clone the [GitHub repo](https://github.com/SaraEl-Metwally/LightTrimmer), e.g. with `git clone https://github.com/SaraEl-Metwally/LightTrimmer.git`
2. Run `make` in the repo directory for **k <= 31**  or `make k=kmersize` for **k > 31**, e.g. `make k=49`.

#### Quick usage guide
``` 
./LightTrimmer -k [kmer size] -g [gap size] -c [kmers counting file] -t [threads] -o [output prefix] [input files] --verbose 

``` 

``` 
* [-k] kmer size                [default: 31]
* [-g] gap size                 [default: 1]
* [-c] kmers counting file      [default: JellyFish format]
* [-t] number of threads        [default: 1]
* [-o] output prefix file name  [default: LightTrimmer]
``` 
#### Notes
- The maximum read length for this version is ``` 1024 bp```.
- The maximum supported read files for this version is ```100``` files.

#### Read files 
LightTrimmer accepts multiple input files of the sequencing reads given in ***fasta/fastq*** format. Also, LightTrimmer can read directly the input files compressed with gzip ***fasta.gz/fastq.gz***.

#### Outputs
The output of LightTrimmer till now is the set of kmers probabilities from Poisson distribution, in the file:

```kmers_prob.txt``` 
