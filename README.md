snpedia-23andme
===============

Find health information from 23andMe's raw data.

The command
  ./snpedia-23andme.py [23andMe-genome-file]
reads 23andMe-genome-file, downloads the SNP information from that file from
SNPedia and saves it in the zip file snpedia-archive.zip, to put less strain
on SNPedia in case of multiple runs. Then outputs a list of all SNP genotypes
sorted according to SNPedia's magnitude scale.
