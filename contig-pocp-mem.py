import sys
import csv
import screed
import argparse
from collections import deque
from sourmash.logging import notify
from sourmash import MinHash, SourmashSignature
from sourmash.index import LinearIndex
from sourmash.distance_utils import containment_to_distance

from typing import NamedTuple


# define NamedTuple for storing fileinfo for each file
class FileInfo(NamedTuple):
    fasta: str
    cidx: LinearIndex
    fidx: LinearIndex
    nskip: int

def read_and_sketch_fromfile(fromfile, moltype='protein', ksize=10, scaled=20, verbose=False):
    """
    Read a CSV file with columns "name", "genome_filename", and "protein_filename".
    Return a dictionary of name: filename, selecting the appropriate filename column
    based on '--moltype' argument passed in via argparse.
    """
    p_skipped = set()
    n=0
    with open(fromfile, 'r') as infile:
        fileinfo = {}
        csvreader = csv.DictReader(infile)
        for n, row in enumerate(csvreader):
            if moltype == 'protein':
                filecol = 'protein_filename'
            else:
                filecol = 'genome_filename'
            name = row["name"]
            fasta = row[filecol]
            fidx,cidx,n_skipped,perc_skipped = read_and_sketch_fasta(fasta, ksize=ksize, moltype=moltype, scaled=scaled, verbose=verbose)
            fileinfo[name] = FileInfo(fasta=fasta,cidx=cidx,fidx=fidx,nskip=n_skipped)
            p_skipped.add(perc_skipped)
             
    avg_p_skipped = sum(p_skipped) / len(p_skipped)
    notify(f"Sketched {n+1} files. On average, skipped {avg_p_skipped:.2f}% of contigs due to scaled value {scaled}.\n" \
           "Decrease scaled value to reduce this (but analysis will take longer).")
    return fileinfo


def read_and_sketch_fasta(fasta_file, ksize, moltype, scaled=1, verbose=False):
    """
    Read and sketch a fasta file.
    Return:
        - contig-signature LinearIndex
        - whole-file signature LinearIndex
        - number of contigs skipped due to 0 hashes
    """
    contigidx = LinearIndex()
    ffidx = LinearIndex()
    n_skipped = 0
    is_protein = (moltype == 'protein')
    mh = MinHash(n=0, ksize=ksize, scaled=scaled, track_abundance=False, is_protein=is_protein)
#    mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled, track_abundance=True, is_protein=(moltype == 'protein'))
    for record in screed.open(fasta_file):
        contig_mh = mh.copy_and_clear()
        if is_protein:
            contig_mh.add_protein(record.sequence)
            mh.add_protein(record.sequence)
        else:
            contig_mh.add_sequence(record.sequence)
            mh.add_sequence(record.sequence)
        if not contig_mh.hashes:
            n_skipped += 1
            continue
        contigidx.insert(SourmashSignature(contig_mh, name=record.name))
    sig = SourmashSignature(mh, name=fasta_file)
    ffidx.insert(sig)
    p_skipped = (n_skipped / len(contigidx)) * 100
    if verbose:
        notify(f"Skipped {p_skipped:.2f}% ({n_skipped} contigs) that were too short due to scaled value")
    return ffidx, contigidx, n_skipped, p_skipped
 
def find_contig_matches(querycontigsdb, searchdb, *, threshold=0.95, ani_threshold=None, verbose=False):
    """Find the contigs in siglist1 that have a containment >= threshold with the db."""
    n_with_match = 0
    percent_with_match = 0
   # match_info = []
    for n, csig in enumerate(querycontigsdb.signatures()):
        # if the query is too short, skip it
        if verbose and n > 0 and n % 100 == 0:
            notify(f"on contig {n} of {len(querycontigsdb)}")
            notify(f"current percent match: {n_with_match / n}")
        if not csig.minhash.hashes:
            raise ValueError(f"query contig '{csig.name}' has no hashes, this should not happen")
        res = searchdb.best_containment(csig, threshold=threshold)
        if not res:
            continue
        # if ani_threshold:
        #     ani = csig.containment_ani(res.signature, containment=res.score).ani
        #     #ani = containment_to_distance(res.containment)
        #     if ani < ani_threshold:
        #         continue
        if res.score < threshold:
            continue
        n_with_match +=1
        if verbose:
            notify(f"found match for '{csig.name}': ({res.score})")
            # optionally, we could store match info. If db is contig-level, we could also store names of matching contigs
            #match_info.append((qcontig.name, res[0].score))
    percent_with_match = (n_with_match / len(querycontigsdb)) * 100
    return percent_with_match, n_with_match


def main(args):
    fileinfo = read_and_sketch_fromfile(args.fromfile, moltype=args.moltype, ksize=args.ksize, scaled=args.scaled, verbose=args.verbose)

    if args.output_csv:
        notify(f"Writing output to '{args.output_csv}'")
        outF = open(args.output_csv, 'w')
        writer = csv.writer(outF)
        writer.writerow(['nameA', 'nameB',  'POCP', 'averagePOCP', 'n_match_fileA', 'n_match_fileB',
                          'ncompared_fileA', 'ncompared_fileB', 'nskipped_fileA', 'nskipped_fileB',
                          'fileA', 'fileB', 'threshold', 'ksize', 'moltype', 'scaled'])
    else:
        notify(f"WARNING: no output file specified. Will not write output.")
        

    names = fileinfo.keys()
    compare_queue = deque(names)
    n_compared = 0
    for name1 in names:
        compare_queue.popleft() # remove name1 from comparison deque
        f1_info = fileinfo[name1]
        # compare all to first file
        for name2 in compare_queue:
            if name1 == name2:
                # this shouldn't happen, but let's notify and fix if it does
                notify('skipping self-comparison')
                continue
            if n_compared % 100 == 0:
                notify(f"starting comparison {n_compared+1} ({name1} x {name2}) ...")
            
            f2_info = fileinfo[name2]
            # find contig->fullfile matches, both directions
            c1_in_f2, c1_matches = find_contig_matches(f1_info.cidx, f2_info.fidx, threshold=args.threshold, ani_threshold=0.95)
            c2_in_f1, c2_matches = find_contig_matches(f2_info.cidx, f1_info.fidx, threshold=args.threshold, ani_threshold=0.95)
            # find contig-contig matches (very slow, better to NOT to do this)
            # c1_in_c2 = find_contig_matches(cidx1, cidx2, args.threshold, verbose=True)
            # c2_in_c1 = find_contig_matches(cidx2, cidx1, args.threshold, verbose=True)
            average_pocp = (c1_in_f2 + c2_in_f1) / 2
            pocp = ((c1_matches + c2_matches) / (len(f1_info.cidx) + len(f2_info.cidx,))*100)
            n_compared += 1

            if args.verbose:
                notify(f"POCP (thresh {args.threshold}): {pocp:.2f} (p1: {c1_in_f2:.2f}, p2: {c2_in_f1:.2f})")
            if args.output_csv:
                writer.writerow([name1, name2, pocp, average_pocp, c1_matches, c2_matches,
                                len(f2_info.cidx,), len(f2_info.cidx,), f1_info.nskip, f2_info.nskip,
                                f1_info.fasta, f2_info.fasta, args.threshold, args.ksize, args.moltype, args.scaled])

    if args.output_csv:
        outF.close()


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser(description='Find the percent of matching contigs between pairs of fasta files.')
    p.add_argument('--fromfile', help='input file with list of fasta files')
    p.add_argument('--ksize', type=int, default=10, help='k-mer size')
    p.add_argument('--scaled', type=int, default=20, help='scaled value')
    p.add_argument('--moltype', default='protein', help='molecule type', choices=['DNA', 'protein'])
    p.add_argument('--threshold', type=float, help='desired threshold', default=0.99)
    p.add_argument('--ani-threshold', '--aai-threshold', type=float, help='desired ANI/AAI threshold', default=0.6)
    p.add_argument('-o', '--output-csv', help='output contig matching details to csv')
    p.add_argument('-v', '--verbose', action='store_true')
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)    
