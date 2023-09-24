from hitparser import parseHit
from functools import partial
from rich.progress import track
from blastparser import HitsDataContainer
import json
from multiprocessing import Pool

#SQLITE_DATABASE_DIR = "../../epitopedia.sqlite3"
#PDB_INPUT="6VXX_A"
#pdb_input_str = "6VXX_A"
#span = 5
#PDB_DATABASE_DIR = "../../mmcif"
#rasafile = "query_pdb_seq.binaryrasa.txt"
#seqnumsfile = "query_pdb_seq.seqnums.txt"
#seqsolvfile = "query_pdb_seq.seqsolv.txt"
#csvfile = "EPI_SEQ_span_filt_acc_hits_6VXX_A.tsv"
#outprefix = "EPI_PDB_fragment_pairs"
import sys
SQLITE_DATABASE_DIR = sys.argv[1]
PDB_INPUT=sys.argv[2]
span = int(sys.argv[3])
PDB_DATABASE_DIR = sys.argv[4]
rasafile = sys.argv[5]
seqnumsfile = sys.argv[6]
seqsolvfile = sys.argv[7]
csvfile = sys.argv[8]
outprefix = sys.argv[9]

#for PDB_INPUT in PDB_INPUTS:
query_pdb_base = PDB_INPUT.split("_")[0].lower()
query_pdb_chain = PDB_INPUT.split("_")[1]

data_m = []

def condInt(x):
    try:
        retval = int(x)
        return retval
    except:
        return x

def condFloat(x):
    try:
        retval = float(x)
        return retval
    except:
        return x


binaryrasa = []
rasa_path = open(rasafile, 'r')
for line in rasa_path:
    line = line.strip()
    binaryrasa.append(condFloat(line))

seqnums = []
seqnums_path = open(seqnumsfile, 'r')
for line in seqnums_path:
    line = line.strip()
    seqnums.append(condInt(line))

seqsolv_path = open(seqsolvfile, 'r')
seqsolv = seqsolv_path.readline().strip()




parseHitP = partial(
                parseHit,
                span=span,
                #pdb_seq=query_pdb_seq,
                seqnums=seqnums,
                seqsolv=seqsolv,
                binaryrasa=binaryrasa,
                query_pdb_base=query_pdb_base,
                query_pdb_chain=query_pdb_chain,
                pdb_input_str=PDB_INPUT,
                SQLITE_DATABASE_DIR=SQLITE_DATABASE_DIR,
                PDB_DATABASE_DIR=PDB_DATABASE_DIR
            )

hits = HitsDataContainer()
hits.fromcsv(f"{csvfile}")

with Pool(12) as p:
                data = list(
                    track(
                        p.imap(parseHitP, hits),
                        total=len(hits),
                        description="Searching and aligning SeqBMM structural representatives",
                    )
                )

data = [datum for datum in data if datum]
data_m.append(data)
data = [data for data in data_m if data]
with open(f"{outprefix}_{PDB_INPUT}.json", "w") as output_handle:
            json.dump({"results": data}, output_handle)
#config.con.close()

