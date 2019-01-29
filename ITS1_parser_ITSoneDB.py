#! /usr/bin/python2
__author__ = 'Bruno Fosso'
__version__ = 2.0

import getopt
import os
import sys
from string import strip
from numpy import mean as calcola_media
from pysam import Samfile


def usage():
    print ('This script performs the Bowtie2 execution:\n'
           'Option:\n'
           '\t-f fasta file containing the ITSoneDB sequences\n'
           '\t-p paired-end sam file\n'
           '\t-s single-end sam file [default = 97]\n'
           '\t-i identity percentage threshold [default = 70]\n'
           '\t-c coverage\n'
           '\t-o output-file path [MANDATORY]\n'
           '\t-h    print this help\n'
           'Usage:\n'
           '\tpython bowtie2-execution_ITSoneDB.py -d bowtie_indexes -v mapping_file\n'
           '\t')


def cigar_parsing(cigar_object):
    cigar_list = list(cigar_object)
    mm = 0
    insertion = 0
    d = 0
    for item in cigar_list:
        if item[0] == 0:
            mm += item[1]
        elif item[0] == 1:
            insertion += item[1]
        elif item[0] == 2:
            d += item[1]
    alen = mm + insertion + d
    qalgn = float(mm + insertion)
    return alen, qalgn


def make_tuple(identifier):
    parts = identifier.split("|")
    return parts[0]


def itsonedb2node(itsonedb_fasta_file):
    diz = {}
    """

    :param itsonedb_fasta_file: a fasta file downladed from ITSoneDB (http://itsonedb.cloud.ba.infn.it)
    :return: a dictionary containing the association between the ITSoneDB accession number and NCBI taxonomy identifier.
    :return: If the input fasta file is None or id doesn't exixt the funciont return None
    """
    import os
    if itsonedb_fasta_file is not None:
        if os.path.exists(itsonedb_fasta_file):
            with open(itsonedb_fasta_file) as b:
                for linea in b:
                    s = map(strip, linea.split("|"))
                    diz[s[0].lstrip(">")] = s[2]
        else:
            diz = None
    else:
        diz = None
    return diz


if __name__ == "__main__":
    fasta_itesondb = None
    paired_sam = None
    single_sam = None
    identity_threshold = 97.0
    coverage = 70.0
    outfile = None
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:p:s:i:c:o:")
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit()
    if len(opts) == 0:
        usage()
        sys.exit()
    for o, a in opts:
        if o == "-h":
            usage()
            sys.exit()
        elif o == "-f":
            fasta_itesondb = a
        elif o == "-p":
            paired_sam = a
        elif o == "-s":
            single_sam = a
        elif o == "-i":
            identity_threshold = float(a)
        elif o == "-c":
            coverage = float(a)
        elif o == "-o":
            outfile = a
        else:
            print "Unhandled option."
            usage()
            sys.exit()

    acc2node = itsonedb2node(fasta_itesondb)

    if outfile is None:
        sys.exit("Output file option is missing")
    match = {}
    # mappiamo prima in modalita' glocal
    if single_sam is not None:
        if os.path.exists(single_sam):
            sam = Samfile(single_sam)
            for align in sam:
                if align.tid != -1:
                    query_name, query_len, ref_name = align.qname, float(align.rlen), sam.getrname(align.tid)
                    if align.cigar is not None:
                        align_len, query_aligned_len = cigar_parsing(align.cigar)
                        nm = -1
                        if (query_aligned_len / query_len)*100 >= coverage:
                            for coppia in align.tags:
                                if coppia[0] == "NM":
                                    nm = float(coppia[1])
                        if align_len != 0 and nm >= 0:
                            paired_perc_id = ((align_len - nm) / align_len) * 100
                            if paired_perc_id >= identity_threshold:
                                match.setdefault(query_name, set())
                                match[query_name].add(ref_name)
            sam.close()
        else:
            print "no mapping data"
            sys.exit()

    if paired_sam is not None:
        if os.path.exists(paired_sam):
            r1_match = {}
            r2_match = {}
            sam = Samfile(paired_sam)
            for align in sam:
                if align.tid != -1:
                    query_name, query_len, ref_name = align.qname, float(align.rlen), sam.getrname(align.tid)
                    if align.cigar is not None:
                        align_len, query_aligned_len = cigar_parsing(align.cigar)
                        nm = -1
                        if (query_aligned_len / query_len)*100 >= coverage:
                            for coppia in align.tags:
                                if coppia[0] == "NM":
                                    nm = float(coppia[1])
                        if align_len != 0 and nm >= 0:
                            paired_perc_id = ((align_len - nm) / align_len) * 100
                            if paired_perc_id >= 90:
                                if align.is_read1:
                                    r1_match.setdefault(query_name, {})
                                    r1_match[query_name].setdefault(ref_name, [])
                                    r1_match[query_name][ref_name].append(paired_perc_id)
                                if align.is_read2:
                                    r2_match.setdefault(query_name, {})
                                    r2_match[query_name].setdefault(ref_name, [])
                                    r2_match[query_name][ref_name].append(paired_perc_id)
            sam.close()
            for query in set(r1_match.keys()).intersection(set(r2_match.keys())):
                for ref in set(r1_match[query].keys()).intersection(r2_match[query].keys()):
                    average_perc_id = calcola_media([max(r1_match[query][ref]), max(r2_match[query][ref])])
                    if average_perc_id >= identity_threshold:
                        match.setdefault(query, set())
                        match[query].add(ref)
        else:
            print "no mapping data"
            sys.exit()

    match_file = open(outfile, "w")
    for acc in match.keys():
        if acc2node is not None:
            match_list = [acc2node[make_tuple(i)] for i in match[acc]]
            match_file.write("%s %s\n" % (acc, " ".join(match_list)))
        else:
            match_list = [make_tuple(i) for i in match[acc]]
            match_file.write("%s %s\n" % (acc, " ".join(match_list)))

    match_file.close()
    print "DONE"
