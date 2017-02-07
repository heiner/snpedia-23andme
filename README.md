snpedia-23andme
===============

**Please note:** This code is old and has probably ceased to work. For alternatives,
check out e.g., [r/slatestarcodex](https://www.reddit.com/r/slatestarcodex/comments/5jnnvz/recommended_3rd_party_analysis_of_23andme_data/)
 on Reddit.

Find health information from 23andMe's raw data.

The command
```sh
  ./snpedia-23andme.py [23andMe-genome-file]
```
reads 23andMe-genome-file, downloads the SNP information from that file from
SNPedia and saves it in the zip file snpedia-archive.zip, to put less strain
on SNPedia in case of multiple runs. Then outputs a list of all SNP genotypes
sorted according to SNPedia's magnitude scale.
