## Instance 
This lesson aims to download a handful of datasets onto a cloud instance, count
kmers in them, and compare kmer spectra before and after some bioinformatic 
filters.

Launch an m3-xlarge instance on Amazon EC2 using the Ubuntu 14.04 server image.

Hints on starting up can be found here:
https://github.com/datacarpentry/cloud-genomics/blob/gh-pages/lessons/1.logging-onto-cloud.md

And a sketch of a lesson, if you'd like a programming exercise / kmer interpretation
puzzle, can be found at http://angus.readthedocs.org/en/2016/automation.html

## Installing the tools

We will need to set up khmer :

```bash
sudo apt-get update
sudo apt-get install -y python-pip python-dev  
sudo easy_install -U setuptools
sudo pip install khmer
```

And kmerspectrumanalyzer 
```bash
sudo apt-get install -y git python-matplotlib python-scipy jellyfish
git clone http://github.com/wltrimbl/kmerspectrumanalyzer
```

and SRAtools
```
wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.7/sratoolkit.2.5.7-ubuntu64.tar.gz
tar xvf sratoolkit.2.5.7-ubuntu64.tar.gz
export PATH=$PATH:$HOME/sratoolkit.2.5.7-ubuntu64/bin
```

Now clone  
```bash
cd && git clone http://github.com/dib-lab/khmer
```
and we need some parts of these tools in our PATH
```bash
export PATH=$PATH:$HOME/khmer/sandbox
export PATH=$PATH:$HOME/kmerspectrumanalyzer/src
```

One more piece -- trimmomatic.
```bash
sudo apt-get install -y default-jre unzip
cd && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.35.zip
unzip Trimmomatic-0.35.zip
trimmomatic=$HOME/Trimmomatic-0.35/trimmomatic-0.35.jar
```

## Check toolkit for completeness
Ok.  Let us test the kit.
```bash
fastq-dump --help && echo OK
error-correct-pass2.py --help && echo OK 
countkmer21.sh  && echo OK
java -jar $trimmomatic --help && echo OK
touch /mnt/littlebunnyfoofoo && echo OK
```

If we can't create a file in /mnt we aren't going to get very far, so 
we change permissions on /mnt
```bash
sudo chown ubuntu /mnt
```

## Downloading sequence data
Now the tools are in place, it is time some some sequence data and do something with it.

I'd like to draw your attention to two datasets that we can get from SRA.

```bash
cd /mnt
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR519/SRR519926/SRR519926.sra #   86M
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR036/SRR036919/SRR036919.sra
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR447/SRR447649/SRR447649.sra #   168M
fastq-dump --split-3 SRR519926.sra 
fastq-dump --split-spot SRR519926.sra
fastq-dump --split-3 SRR036919.sra 
# This step creates SRR447649_1.fastq and SRR447649_2.fastq 
fastq-dump --split-3 SRR447649.sra 
# This puts the same data into a single file SRR447649.fastq 
fastq-dump --split-spot SRR447649.sra 
```
This downloads two smallish datasets from SRA in SRA's format, and uncompresses them into one or two files, depending on whether the sequencing run produced paired reads or not.

We get output like
```bash
Rejected 1373828 READS because of filtering out non-biological READS
Read 1373828 spots for SRR447649.sra
Written 1373828 spots for SRR447649.sra
```

## Counting kmers

I'll show you two ways to count the kmers in these datasets: first, using khmer.
First we create a kmer hash from the two
```bash
load-into-counting.py -x 1e9 -k 21 SRR447649.4G.kh SRR447649_1.fastq SRR447649_2.fastq 
abundance-dist.py -s SRR447649.4G.kh  SRR447649.fastq SRR447649.4G.21
```

```bash
Saving k-mer countgraph to SRR447649.4G.kh
Loading kmers from sequences in ['SRR447649_1.fastq', 'SRR447649_2.fastq']
making countgraph
consuming input SRR447649_1.fastq
consuming input SRR447649_2.fastq
Total number of unique k-mers: 43156383
saving SRR447649.4G.kh
```

Now we have one kmer spectrum, let us look at it.
```bash
plotkmerspectrum.py SRR447649.21 -g 3
```
This creates a pdf file in /mnt.  

And, since it just takes a minute, let us count our other two datasets
```bash
fastq-dump --split-spot SRR519926
load-into-counting.py -x 1e9 -k 21 SRR519926.4G.kh SRR519926.fastq 
abundance-dist.py -s SRR519926.4G.kh  SRR519926.fastq SRR519926.4G.21
load-into-counting.py -x 1e9 -k 21 SRR036919.4G.kh SRR036919.fastq 
abundance-dist.py -s SRR036919.4G.kh  SRR036919.fastq SRR036919.21
```

The data carpentry 
(cloud genomics class)[https://github.com/JasonJWilliamsNY/cloud-genomics/blob/master/lessons/3.single-analysis.md] has a recipe for Q-value trimming using Trimmomatic.

```bash 
mkdir /mnt/SRR519926_trimmed
cd /mnt/SRR519926_trimmed
java -jar $trimmomatic  PE -phred33 -trimlog trimlog.txt ../SRR519926_1.fastq ../SRR519926_2.fastq p1.fq u1.fq p2.fq u2.fq LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:50 2>&1 | tee cmd.txt
```
This gives us four files for paired and unpaired reads post-trimming p1.fq u1.fq p2.fq u2.fq

Now let us combine the four output files and count kmers
```bash
cat p1.fq u1.fq p2.fq u2.fq | countkmer21.sh > SRR519926_trimmed.21
```

```bash
plotkmerspectrum.py SRR519926_trimmed.21 ../SRR519926.21 -g 5
plotkmerspectrum.py SRR519926_trimmed.21 ../SRR519926.21 -g 3
plotkmerspectrum.py SRR519926_trimmed.21 ../SRR519926.21 -g 6
plotkmerspectrum.py SRR519926_trimmed.21 ../SRR519926.21 -g 20
plotkmerspectrum.py SRR519926_trimmed.21 ../SRR519926.21 -g 1
```
These commands do two things; they create pdf graphs comparing the two spectra, 
and a one-line statistical summaries in the file kmers.log.

Now if we compare the kmer spectrum before and after, we find several things:
* The trimming reduced our total depth.  
* The trimming dramatically reduced the number and fraction of singleton observations.
* The dataset has modest (2%) levels of adapter contamination that were not addressed by the above recipe, and large fractions (%) of unique, presumptively erroneous kmers.

Let us try a different scrubbing recipe.  First, we need a file with the right contaminating adapters

```bash
cat > ~/adapters.fa <<EOF
>TruSeqUniversalAdapter-P5_R
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>P7_indexing74_R
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAACTCATCTCGTATGCCGTCTTCTGCTTG
EOF
```

And now we use a different functionality in trimmomatic to remove the adapters
```bash
mkdir /mnt/SRR519926_adapterscrub
cd /mnt/SRR519926_adapterscrub
java -jar $trimmomatic PE -phred33 -trimlog trimlog.txt ../SRR519926_1.fastq ../SRR519926_2.fastq p1.fq u1.fq p2.fq u2.fq ILLUMINACLIP:$HOME/adapters.fa:2:30:10 2>&1 | tee cmd.txt
cat p1.fq u1.fq p2.fq u2.fq | countkmer21.sh > SRR447649_scrubbed.21
```

Comparing this kmer spectrum, 
* Adapter scrubbing barely reduced the total depth on the genome.
* Adapter scrubbing only slightly reduced the number and fraction of singleton observations.



