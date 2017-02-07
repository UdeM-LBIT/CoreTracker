#!/usr/bin/env python

# CoreTracker Copyright (C) 2016  Emmanuel Noutahi
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from coretracker.coreutils import SequenceLoader, CoreFile, makehash, utils
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import generic_protein
import argparse
import json
import os
import sys


class DomainRPR(object):

    def __init__(self, species):
        self.species = species
        self.position = {}
        self.description = {}
        self.evalue = {}
        self.pfam = {}
        self.score = {}

    def add_domain(self, name, pfam, start, end, description="", evalue=None, score=None):
        exist_al = self.pfam.get(name, None)
        if exist_al:
            name += "|%s-%s" % (str(start), str(end))
        self.position[name] = (int(start), int(end))
        self.pfam[name] = pfam
        self.description[name] = description
        self.evalue[name] = float(evalue)
        self.score[name] = float(score)

    def get_domain_evalue_thresh(ethresh):
        for n, ev in evalue.items():
            if ev < ethresh:
                yield n


def parse_json_file(jfile):
    """Load json position file"""
    t_len = 0
    with open(jfile) as INDATA:
        tmpdt = json.load(INDATA)
    data = makehash(2, list)
    for (spec, readict) in tmpdt.items():
        for rea, pos_list in readict.items():
            for info in pos_list:
                gene, pos = info
                t_len += 1
                data[gene][rea][spec].append(pos)
    return data, t_len


def load_alignment(alfile, alpha=generic_protein):
    """Load alignment, expect a fasta file"""
    c_inst = CoreFile(alfile, alpha)
    return c_inst


def core_file_to_dict(corefile):
    """Convert corefile to dict of dict"""
    d = dict((k, SeqIO.to_dict(v)) for (k, v) in corefile.items())
    return d


def cvt_pos_al_to_pos_seq(data, reakey, gene, aldict, gap='-'):
    """Convert position in alignment to position in sequence"""
    gpos = {}
    for species in aldict[gene].keys():
        seq = aldict[gene][species]
        try:
            position = data[gene][reakey][species]
            seq_pos = []
            for pos in position:
                l = len(seq[:pos].seq.ungap(gap))
                seq_pos.append(l - 1)
            gpos[species] = seq_pos
        except Exception, e:
            print e
    return gpos


def get_domain_in_gene(gene, outdir, path_to_db, addargs=""):
    """Get list of domains in a multiple seq"""
    # this function assume that you have downloaded Pfam
    # and run hmmpress on it
    # Also there won't be attempt to check if you have the hmmmer installed
    # work with hmmmer3. You should alias hmmscan and hmmpress to the corresponding
    # binaries
    # this also assume that xport_align was called
    dirname = outdir
    if gene not in dirname:
        dirname = os.path.join(outdir, gene)
    outfile = os.path.join(dirname, "%s-domains.tab" % gene)
    outmp = os.path.join(dirname, "%s-domains.out" % gene)
    seqfile = os.path.join(dirname, "%s.fasta" % gene)
    build_cmd = "hmmscan -o %s --domtblout %s %s %s %s" % (
        outmp, outfile, path_to_db, seqfile, addargs)
    print(build_cmd)
    utils.executeCMD(build_cmd, 'hmmscan')
    cut_cmd = "grep -v '#' %s | sed 's/  */\t/g' |cut -f 1,2,4,7,8,18,19,20,21,23- > %s" % (
        outfile, outfile + ".tmp")
    utils.executeCMD(cut_cmd, 'grep')
    return parse_scan_out(outfile + ".tmp")


def xport_align_to_dir(c_inst, outdir, aligned=False, gap='-'):
    for (gene, align) in c_inst.items():
        dirname = utils.purge_directory(os.path.join(outdir, gene))
        otfile = os.path.join(dirname, gene + ".fasta")
        if aligned:
            AlignIO.write(align, open(otfile, 'w'), format='fasta')
        else:
            for s in align:
                s.seq = s.seq.ungap(gap)
            SeqIO.write(align, open(otfile, 'w'), format='fasta')


def parse_scan_out(outf):
    gene_domains = {}
    with open(outf) as OUT:
        # domain, pfam, seq, eval, score, strt1, end1, strt2, end2 ~description
        for line in OUT:
            line = line.strip()
            if line:
                l = line.split()
                spec = l[2]
                spec_do = gene_domains.get(spec, DomainRPR(spec))
                spec_do.add_domain(l[0], l[1], int(
                    l[7]) - 1, int(l[8]) - 1, " ".join(l[9:]), l[3], l[4])
                gene_domains[spec] = spec_do
    return gene_domains


def in_domain(gdomain, gpos):
    rdata = {}
    for spec in gpos.keys():
        dchecked = []
        havedom = gdomain.get(spec, None)
        if havedom:
            for pos in gpos[spec]:
                doname = None
                start = 0
                end = 0
                for (name, s_e) in gdomain[spec].position.items():
                    start, end = s_e
                    if start <= pos and pos <= end:
                        doname = name
                        break
                if doname:
                    dchecked.append((pos, doname, gdomain[spec].pfam[doname], gdomain[
                                    spec].description[doname], str(start), str(end)))

                else:
                    dchecked.append((pos, "NID", "NID", "NID", '-', '-'))
        else:
            dchecked = [(pos, "N/A", "N/A", "N/A", '-', '-')
                        for pos in gpos[spec]]
        rdata[spec] = dchecked
    return rdata


def domain_checker(posfile, sequencefile, outdir, hmmdb, addargs="", aligned=False):
    corefile = load_alignment(sequencefile)
    aldict = core_file_to_dict(corefile)
    data, t_len = parse_json_file(posfile)
    print("Total of: %d" % t_len)
    xport_align_to_dir(corefile, outdir, aligned)
    new_dt = makehash(2, list)
    for gene in data.keys():
        gdomain = get_domain_in_gene(gene, outdir, hmmdb, addargs)
        rea = data[gene]
        for r in rea.keys():
            gpos = rea[r]
            if not aligned:
                gpos = cvt_pos_al_to_pos_seq(data, r, gene, aldict)
                new_dt[gene][r] = in_domain(gdomain, gpos)
    return new_dt


def export_data_as(data, outfile, format='json'):
    if format == 'csv':
        with open(outfile, 'w') as OUT:
            OUT.write("\t".join(['Species', 'Gene', 'Codon', 'AA', 'Pos',
                                 'Domain', 'Pfam', 'start', 'end', 'description']) + "\n")
            for (gene, rea) in data.items():
                for r, specmap in rea.items():
                    codon, aa = r.split(':')
                    for spec, posmap in specmap.items():
                        for pos in posmap:
                            OUT.write("\t".join([spec, gene, codon, aa, str(pos[0]), pos[
                                      1], pos[2], str(pos[4]), str(pos[5]), pos[3]]) + "\n")
    elif format == 'json':
        with open(outfile, 'w') as OUT:
            json.dump(data, OUT, indent=4)


if __name__ == '__main__':
    curdir = os.getcwd()
    fchoices = ['json', 'csv']
    parser = argparse.ArgumentParser(
        description='Find information about position')
    parser.add_argument('--corefile', '-c', required=True,
                        dest="corefile", help="Protein corefile")
    parser.add_argument('--posfile', '-p', required=True,
                        dest="posfile", help="Position file in json format")
    parser.add_argument('--outdir', dest='outdir',
                        default=curdir, help="Output directory")
    parser.add_argument('--outfile', '--out', '-o',
                        dest='outfile', default='outfile', help="Output file")
    parser.add_argument('--pfamhmm', '--hmm', required=True,
                        dest="hmmdb", help="PFam hmm. ou should run hmmpress before")
    parser.add_argument('--hmmscan_args', '--extra_args', default="", dest="extrarg",
                        help="Position in alignment, separated by comma or a file")
    parser.add_argument('--align', '-a', action='store_true',
                        dest="aligned", help="save data as alignment intead of sequences")
    parser.add_argument('--format', '--fmt', choices=fchoices,
                        dest="format", help="Format to save the data in 9json or csv")

    args = parser.parse_args()
    dt = domain_checker(args.posfile, args.corefile, args.outdir,
                        args.hmmdb, args.extrarg, aligned=args.aligned)
    outfile = args.outfile
    outdirname = os.path.dirname(args.outfile)
    if not outdirname:
        outfile = os.path.join(args.outdir, args.outfile)
    export_data_as(dt, outfile, format=args.format)
