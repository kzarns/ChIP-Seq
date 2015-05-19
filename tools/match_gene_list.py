#!/usr/bin/python
"""generate fasta for a list of genes"""
from subprocess import call
import sys
import re
#from collections import defaultdict

from optparse import OptionParser
#print "python match_gene_list.py sites_file chrom_dir out range limit"
usage = "USAGE: %prog [options] sites_file chrom_dir out_file"
parser = OptionParser(usage=usage)
parser.add_option("-r", "--range", dest="range", help="The range +-N to retrieve around the tss", default=10, metavar="INT")
parser.add_option("-l", "--limit", dest="limit", help="The limit on the number of lines to read from the input file", default=None, metavar="INT")
parser.add_option("-f", "--filter", dest="filter", help="The name of a 'csv' file to filter by", default=None, metavar="STRING")
comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 't':'a', 'g':'c'}
filter = None

class GeneSite:
    """represent a gene, this may be more overhead than a dict"""
    def __init__(self, line):
        line = line.rstrip()
        #        print line
        words = re.split("[ \t]", line)
#        print words
#        print "line: %s\n" % line
#        print "words: %s\n" % words
#        print "words[1]: %s" % words[1]
        (symbol, sep, tss) = words[1].rpartition(":")
        self.gene = symbol
        self.chromosome = words[2]
        self.start = int(words[3])
        self.end = int(words[4])
        self.strand = words[5]

    def __repr__(self):
        return "<%s %s %s %s %s>" % (self.gene, self.chromosome, self.start, self.end, self.strand)

    def __str__(self):
        return "'%s %s %s %s %s'" % (self.gene, self.chromosome, self.start, self.end, self.strand)


def read_sites(sites_file, chrom_dir, chromes, out, range, limit):
    """read the gene_list file"""

    #skip the header line
    next(sites_file)

    num_sites = 0
    for line in sites_file:

        site = GeneSite(line)
        if filter is not None and site.gene not in filter:
            print "move along: %s" % site.gene
            continue

        num_sites += 1
#        print "limit: %s, site: %s\n" % (limit, num_sites)
        if num_sites % 100 == 0:
            print "limit: %s, site: %s\n" % (limit, num_sites)
        if limit is not None and num_sites > limit:
            break

        if site.chromosome not in chromes:
            try:
                with open("%s/%s.fa" % (chrom_dir, site.chromosome), 'r') as chromosome_file:

                    contents = chromosome_file.read()
                    pos = contents.index('\n')
                    contents = contents[pos+1:]
                    contents = re.sub('[\n\r \t]', '', contents)
                    chromes[site.chromosome] = contents

                print "read chromosome", site.chromosome
            except IOError, ex:
                print "skipping site: %s missing file: %s" % (line, ex)
                continue



        #account for reversed
        if site.start > site.end:
            print "need to swap start/end for:", site
            tmp = site.start
            site.start = site.end
            site.end = tmp

        if site.strand == "plus":
            #print "site before: %s" % site
            site.end = site.start + range
            site.start = site.start - (range + 1)
            #print "site after: %s" % site
        elif site.strand == "minus":
            #print "site before: %s" % site
            site.start = site.end - range
            site.end = site.end + (range + 1)
            #print "site after: %s" % site
        else:
            print "strand error: %s" % site

        process_site(site, out, chromes, site.chromosome)

    print "total sites used: %s" % num_sites

def process_site(site, out, chromes, chromosome):
    """process for one chromosome"""

    snip = chromes[chromosome][site.start:site.end]

    if snip == '':
        print "ERROR: snip was '' for: ", site
    else:
        if site.strand == "minus":
            #print "snip before: %s\n" % snip
            bases = list(snip)
            bases = reversed([comp.get(base, base) for base in bases])
            snip = ''.join(bases)
            #print "snip after: %s\n" % snip

        out.write("> %s\n" % site)
        out.write("%s\n" % snip)
        out.write("\n")

def read_filter(filter_file):
    global filter
    with open(filter_file) as fil_in:
        filter = dict((key.rstrip(), 1) for key in fil_in)
    print filter

def main():
    """keep it all organized"""
    print "python match_gene_list.py sites_file chrom_dir out range limit"
    (options, args) = parser.parse_args()

    print args
    print options
    with open(args[0]) as sites_file:
        
        chrom_dir = args[1]
        chromes = {}

        range = int(options.range)
        limit = None
        if options.limit is not None:
            limit = int(options.limit)

        filter_file = options.filter
        if filter_file:
            filter = read_filter(filter_file)

        with open(args[2], 'w') as out:

            read_sites(sites_file, chrom_dir, chromes, out, range, limit)

main()
