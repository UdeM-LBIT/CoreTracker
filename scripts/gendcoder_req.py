#!/usr/bin/env python
import requests
import sys
import glob
import os
import time
from Bio import SeqIO
from lxml import etree
import argparse

waittime = 5
entropy = {'conserved': "0.0", 'high': "1.0",
           'weak': "2.0", 'variable': "3.0", 'all': "4.3"}
gaps = {'0%': "0.0", '20%': "0.2", '50%': "0.5"}
dataset = {'Metazoa': 'M', 'Nematoda': 'N', 'Platyhelminthes': 'P'}
root = "http://darwin.uvigo.es/cgi-bin/genDecoder.pl"


def parse_gdec_out(cntnt, outfile):
    root = etree.HTML(cntnt.strip())
    a_exp = root.xpath('.//a[text()="Expected"]')[0]
    code = a_exp.get("title", "")
    expected_aa = [x for x in a_exp.tail.split(':')[1].split("\n")[0].strip()]
    a_list = root.findall(".//*b/a[@target='codon']")
    a_parent = a_list[0].getparent()
    parts = [a_parent.text]
    for a in a_list:
        parts += [a.text] + [child.text for child in a.getchildren()] + \
            [a.tail]
    predicted_aa = "".join([x.strip() for x in filter(None, parts)])
    sup_infos = a_parent.tail.strip().split("\n")
    first = sup_infos[0].split(':')[1].strip()
    second = sup_infos[1].split(':')[1].strip()
    third = sup_infos[2].split(':')[1].strip()
    with open(outfile, 'w') as OUT:
        OUT.write("#%s\n" % code)
        OUT.write("#codon\tGenDecoder\tExpected\n")
        print "File ==> ", outfile
        # print "Predicted len ==> ", len(predicted_aa)
        print "Predicted ==>", predicted_aa
        # print(etree.tostring(a_parent, pretty_print=True))
        for i, letter in enumerate(first):
            codon = letter + second[i] + third[i]
            OUT.write("%s\t%s\t%s\n" %
                      (codon, predicted_aa[i], expected_aa[i]))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Predict genetic code from GenDecoder Web server')
    parser.add_argument('--entropy', choices=tuple(entropy.keys()),
                        default="weak", help="Set filtering entropy values")
    parser.add_argument('--gaps', choices=tuple(gaps.keys()),
                        default='20%', help="Gap percent to remove from the alignment")
    parser.add_argument('--dataset', choices=tuple(dataset.keys()),
                        default="Metazoa", help="Dataset to use for pfam")
    parser.add_argument('--accession', nargs="*",
                        help="List of accession number to use")
    parser.add_argument(
        '--gbdir', help="Use a directory that contains genbank files instead")
    parser.add_argument('--outdir', '--wdir', dest="outdir",
                        default="gdecoder", help='Output filename')
    args = parser.parse_args()

    if args.outdir and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    filemode = False
    if args.gbdir:
        naccs = glob.glob(os.path.join(args.gbdir, '*.gbk'))
        filemode = True
    else:
        naccs = [x.strip() for x in args.accession]
    # penality = (len(naccs)//10)
    start = 0
    for f in naccs:
        start += 1
        tryagain = 2
        succeed = False
        nacc = f
        specname = nacc
        if filemode:
            rec = SeqIO.read(f, format="genbank")
            nacc = rec.id
            specname = rec.annotations["organism"]

        while tryagain and not succeed:
            try:
                if filemode:
                    r = requests.post(root, files={'GBFILE': (f, open(f, 'rb'))}, data={'DATASET': dataset[
                                      args.dataset], 'ENTROPY': entropy[args.entropy], 'GAPS': gaps[args.gaps]})
                else:
                    r = requests.post(root, data={'GBACCNUM': nacc, 'ENTROPY': entropy[
                                      args.entropy], 'GAPS': gaps[args.gaps], 'DATASET': dataset[args.dataset]})
                assert r.status_code == requests.codes.ok, 'request to %s failed' % r.url
                succeed = True
                parse_gdec_out(r.text, os.path.join(
                    args.outdir, specname + ".txt"))

            except AssertionError:
                print "*** %s : failed job, retry again" % specname
                tryagain -= 1

            # * (1 if (start%penality) else penality/2.0))
            time.sleep(waittime)
