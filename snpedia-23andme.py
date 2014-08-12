#!/usr/bin/env python

"""./snpedia-23andme.py [23andMe-genome-file]

Find health information from 23andMe's raw data.

Reads 23andMe-genome-file, downloads the SNP information from that file from
SNPedia and saves it in the zip file snpedia-archive.zip, to put less strain
on SNPedia in case of multiple runs. Then print a list of all SNP genotypes
sorted according to SNPedia's magnitude scale.
"""

import sys
import csv
import json
import re
import zipfile

import pprint
pp = pprint.PrettyPrinter(indent=4).pprint

sys.path.append("wikitools") # Hack for using wikitools as a git submodule
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

def test_and_store_genotype(snpinfo, rsid, genotype, matches=[]): # mutable matches list!
    if snpinfo[rsid] and genotype in snpinfo[rsid]["genotypes"]:
        snpmatch = dict(snpinfo[rsid])
        snpmatch["rsid"] = rsid
        snpmatch["genotype"] = genotype
        snpmatch["magnitude"] = snpmatch["genotypes"][genotype]["magnitude"]
        snpmatch["link"] = "http://www.snpedia.com/index.php/" + rsid.capitalize()
        matches.append(snpmatch)
    return matches

def puts(s):
    sys.stdout.write(s)
    sys.stdout.flush()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    with open(sys.argv[1], 'r') as genomefile, \
            zipfile.ZipFile("snpedia-archive.zip", mode="a",
                            compression=zipfile.ZIP_DEFLATED) as ziparchive:
        site = wiki.Wiki("http://bots.snpedia.com/api.php")
        snpsfile = csv.DictReader(genomefile, delimiter="\t",
                              fieldnames=["rsid", "chromosome",
                                          "position", "genotype"])

        if "snpedia_rsids" in ziparchive.namelist():
            snpedia_rsids = {line.rstrip() for line in ziparchive.read("snpedia_rsids")}
        else:
            puts("Get list of SNPs on SNPedia ... ")
            snps = category.Category(site, "Is_a_snp")
            snpedia_rsids = {article.title.lower()
                             for article in snps.getAllMembersGen(namespaces=[0])}
            ziparchive.writestr("snpedia_rsids", "\n".join(sorted(snpedia_rsids)))
            puts("done\n")

        try:
            with open("snpedia-archive.json", "r") as snpinfofile:
                snpinfo = json.load(snpinfofile)
        except (IOError, ValueError):
            snpinfo = {}

        namelist = set(ziparchive.namelist())

        counter = 0
        puts("  ")
        for snp in snpsfile:
            rsid = snp["rsid"]
            if rsid[0] == "#":
                # skip comments
                continue

            counter += 1
            if counter % 10000 == 0:
                puts(" %6i\n  " % counter)
            elif counter % 500 == 0:
                puts(" ")
            elif counter % 100 == 0:
                puts(".")

            genotype = snp["genotype"]

            # Check json file first
            if rsid in snpinfo:
                #print "SNP", rsid, "present in json file"
                matches = test_and_store_genotype(snpinfo, rsid, genotype)
                continue

            # If not found there, check zip file or snpedia.com
            if rsid in namelist:
                pass
                #puts("SNP " + rsid + " present in zip file\n")
            elif rsid in snpedia_rsids:
                puts("Query SNPedia for " + rsid + " ... ")
                pagehandle = page.Page(site, rsid)
                ziparchive.writestr(rsid, pagehandle.getWikiText(expandtemplates=True))
                puts("done.\n")
            else:
                #puts("SNP " + rsid + " not on SNPedia.\n")
                continue

            # If we made it here, the zip file has some information
            snpinfo[rsid] = read_snpedia(ziparchive.read(rsid))
            if not snpinfo[rsid]:
                continue

            matches = test_and_store_genotype(snpinfo, rsid, genotype)

        puts(" %6i\n" % counter)

        matches.sort(key=lambda g: g["magnitude"], reverse=True)

        with open(sys.argv[1] + ".json", "w") as resultfile:
            json.dump(matches, resultfile, indent=2, separators=(',', ': '))

        # with open(sys.argv[1] + ".csv", "w") as resultfile:
        #     writer = csv.writer(resultfile)
        #     writer.writerow(["magnitude", "comment", "genotype", "rsid", "link"])
        #     for match in matches:
        #         genotype_info = match["genotypes"][match["genotype"]]
        #         writer.writerow([
        #                 match["magnitude"],
        #                 genotype_info["comment"],
        #                 match["genotype"],
        #                 match["rsid"],
        #                 "http://www.snpedia.com/index.php/" + match["rsid"].capitalize()])
        # print "Wrote output to " + sys.argv[1] + ".csv"

        with open("snpedia-archive.json", "w") as snpinfofile:
            json.dump(snpinfo, snpinfofile, indent=2, separators=(',', ': '))
