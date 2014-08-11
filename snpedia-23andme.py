#!/usr/bin/env python

"""snpedia-23andme.py [23andMe-genome-file]

Find health information from 23andMe's raw data.

Reads 23andMe-genome-file, downloads the SNP information from that file from
SNPedia and saves it in the zip file snpedia-archive.zip, to put less strain
on SNPedia in case of multiple runs. Then print a list of all SNP genotypes
sorted according to SNPedia's magnitude scale.
"""

import sys
import csv
import re
import zipfile

import pprint
pp = pprint.PrettyPrinter(indent=4).pprint

from wikitools import api, category, page, wiki

def read_snpedia(string):
    """Parse and return the information from an SNPedia HTML that we are interested in."""
    match = re.search(r"\[\[Orientation::(\w+)\]\]", string)
    if not match:
        return {}
    orientation = match.group(1)

    match = re.search(r"\[\[Max Magnitude::(\d+)\]\]", string)
    if not match:
        return {}
    max_magnitude = float(match.group(1))

    match = re.search(r"\[\[Magnitude\|Mag\]\] .+? \|- \s*", string,
                      re.MULTILINE | re.VERBOSE | re.DOTALL)
    string = string[match.end(0):]

    genotypes = {}
    while True:
        match = re.match(r"""
\| \s* \[\[(?P<variant>\w+\([ATCG];?[ATCG]\)) \| (?P<genotype>[^\]]*)\]\] \s*
\| (?:\s*\w+=\"[^\"]*\")* \s* \| \s* (?P<magnitude>[\d.]*) \s*
\| (?P<comment>.*?)
\|-\s*""", string, re.MULTILINE | re.VERBOSE | re.DOTALL)
        if not match:
            break
        string = string[match.end(0):]

        genotype = match.group("genotype")
        genotype = re.sub("\W", "", genotype)
        if orientation == "minus":
            genotype = switch_orientation(genotype)

        magnitude = match.group("magnitude")
        magnitude = float(magnitude) if magnitude else None

        genotypes[genotype] = dict(magnitude=magnitude,
                                   variant=match.group("variant"),
                                   comment=match.group("comment").strip())

    return dict(max_magnitude=max_magnitude,
                original_orientation=orientation,
                genotypes=genotypes)

def switch_orientation(genotype):
    """Switch orientation of the genotype string"""
    substitutions = { "A": "T",
                      "T": "A",
                      "C": "G",
                      "G": "C" }
    return "".join(map(substitutions.get, genotype))



if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    with open(sys.argv[1], 'r') as genomefile, \
            zipfile.ZipFile("snpedia-archive.zip", mode="a",
                            compression=zipfile.ZIP_DEFLATED) as archive:
        site = wiki.Wiki("http://bots.snpedia.com/api.php")
        snpsfile = csv.DictReader(genomefile, delimiter="\t",
                              fieldnames=["rsid", "chromosome",
                                          "position", "genotype"])

        if "snpedia_rsids" in archive.namelist():
            snpedia_rsids = [line.rstrip() for line in archive.read("snpedia_rsids")]
        else:
            sys.stdout.write("Get list of SNPs on SNPedia ... ")
            sys.stdout.flush()
            snps = category.Category(site, "Is_a_snp")
            snpedia_rsids = [article.title.lower()
                             for article in snps.getAllMembersGen(namespaces=[0])]
            snpedia_rsids.sort()
            archive.writestr("snpedia_rsids", "\n".join(snpedia_rsids))
            sys.stdout.write("done\n")
        ##pp(snpedia_snp)
        # rsid = "rs11614913"
        # pp(read_snpedia(archive.read(rsid)))

        namelist = set(archive.namelist())
        matches = []

        for snp in snpsfile:
            rsid = snp["rsid"]
            if rsid[0] == "#":
                # skip comments
                continue
            #pp(snp)

            if rsid in namelist:
                #pass
                sys.stdout.write("SNP " + rsid + " present in zip file\n")
                sys.stdout.flush()
            elif rsid in snpedia_rsids:
                sys.stdout.write("Query SNP " + rsid + " ... ")
                sys.stdout.flush()
                pagehandle = page.Page(site, rsid)
                archive.writestr(rsid, pagehandle.getWikiText(expandtemplates=True))
                sys.stdout.write("done.\n")
            else:
                continue
                #sys.stdout.write("SNP " + rsid + " not on SNPedia.\n")

            snp_info = read_snpedia(archive.read(rsid))
            if not "max_magnitude" in snp_info:
                continue

            pp(snp_info)
            genotype = snp["genotype"]
            if genotype in snp_info["genotypes"]:
                matches.append(snp_info["genotypes"][genotype])

        print "Sorted matched genotypes:"
        matches.sort(key=lambda g: g["magnitude"], reverse=True)
        pp(matches)
