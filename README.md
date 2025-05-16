# MSMuTect
Indel, Allele and Mutation caller, specifically designed to call mutations in microsatellite regions using a pair of sequencing files (normal and tumor samples)

# Installation
### Binary
There is a prebuilt x86_64 linux binary available in releases.  
Note: The binary is slightly slower than the 'Local' option.   
Download the binary from the following link:
[fill in link]
### Local
If on a different platform, or to achieve maximum performance, do the following:  
git clone https://github.com/MaruvkaLab/MSMuTect_4  
cd MSMuTect_4
pip3 install -r requirements.txt  
bash build.sh  #optional; improves performance
When running, use MSMuTect_4/msmutect.sh everywhere the documentation says 'msmutect'

# A Note About Versions
The versions are very confused due to a discontinuous history of development. Practically speaking
you should use either version 4.0 (which is fastest) or version (4.1) which is ~2x slower but 
is much more accurate for longer motif repeats. Both have executable binaries available on the github

# Usage
## Locus File
We strongly recommend using one of the precompiled locus files available at ##. They can be subsampled
as desired as long as the order is maintained

If want to use your own locus file, there are a couple of steps you must do:
1. First, the loci file must be sorted properly. This should work on all unix systems:   
sort -t $'\t' -k1,1 -k5n,5 -k4nr,4 -V [original loci file] > [new loci file]    
2. If using version>=4.1, you must also correct the loci file so that the motif is written as it
appears on the forward strand. You can do this with ./tests/scratch/post_process_loci_file.py 
(change the paths as needed)

## Running the Software

msmutect [flags]
  
To fully analyze all loci:  
msmutect -T [tumor_bam.bam] -N [normalbam.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -m -A -H  

To find mutations in the most efficient runtime possible:  
msmutect -T [tumor_bam.bam] -N [normalbam.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -m  

To call indels and alleles for an individual file:  
msmutect -S [sequence_file.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -A -H  

To call indels but not alleles for an individual file:  
msmutect -S [sequence_file.bam] -l [loci_file.phobos] -O [output_prefix] -c [number of cores to use] -H  

To call mutations using results from previously called indels (can only use a single core):
msmutect --from_file -l [loci_file.phobos]  -N [normal.hist.tsv] -T [tumor.hist.tsv]  -O [output_prefix] -m -A 

To see all flags, such as changing parameters, outputting vcf files, etc., run 'msmutect --help'

MSMuTect will create temporary files when running, with names like .tmp_10242_1721809243.1243694_25529.  
It deletes them at the end. If, for some reason, msmutect is interrupted, these files will not be deleted. They can be safely removed 

### Understanding the 'Call' Column
M = Mutation  
NM = Not Mutation
AN = No Alleles. Either the tumor sample or the normal sample lacks alleles  
RR = Reversion to Reference. The normal sample held an alternative allele, and the tumor had a mutation of the reference allele  
GV = Germline Variation. There are too many SNPs in the vicinity of the locus to confidently say that the indels were of the MS motif
FFT = Failed Fisher Test. Passed other tests to be called a mutation, but failed the Fisher's exact test of significance      
INS = Insufficient Support. The normal sample has multiple alleles, but one of them has insufficient support, indicating a noisy locus  
TMA = Too Many Alleles. Normal sample has too many alleles, and hence the locus is too noisy to call    

MSMuTect can generate a vcf file with the results in addition to the regular output (a tsv file).   
The vcf file will include every locus that had an alternate allele in either the tumor or the normal sample.   
However, every locus that is not called as a mutation will be marked as filtered in the filter column of the vcf file, with the resulting call (ex. NM, INS, AN, etc.) as the filter
# Changes, Pull Requests, and Compiling the Binary Executable
Users should typically install MSMuTect as described in the "Installation" section. However, if you 
wish to make changes, this is also practicable. The code is straightforward and the entry point is in src/Entry/main.py.  
When making changes, you don't need to recompile the code every time: simply don't run build.sh and msmutect.sh will use the python files directly.  
If you add any interesting features, we would be grateful if you could share them with us by emailing us or opening a pull request.  
If you did run build.sh and want to revert the .pyx files to python files, run reverse_rename.sh.   
If for some reason you would like to compile a binary executable, this is a little tricky and not recommended. The most important thing 
to understand is that there are two levels of compilation:  
1. Compiling every individual module. This is purely to improve performance. It is not neccesary to run msmutect.sh
2. Compiling the binary executable so that it can be run by itself.  

The steps to create the binary executable are as follows:   
1. Run build.sh to create compile the individual python modules
2. Go to pyinstaller_build
3. Update the paths in build.sh (this is the build.sh in pyinstaller_build, not the one in the top level directory)
4. Run build.sh. The resulting executable will be put in pyinstaller_build/dist



# Publication and Citation
For orginal paper, see 
YE  Maruvka, Mouw K,  et al, Analysis of somatic microsatellite indels identifies driver events in human tumors
This version is known as version 4.1

# Authors
Avraham Kahan, Dr. Yosef Maruvka, Gaia Frant, and the Maruvka Lab at Technion  
For questions, suggestions, or concerns, open an issue on github or email k.avraham@technion.ac.il or yosi.maruvka@bfe.technion.ac.il

