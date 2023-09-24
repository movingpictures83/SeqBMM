# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import sqlite3
from dataclasses import dataclass, field

#from epitopedia.app import config

residue_max_acc = {
    # Miller max acc: Miller et al. 1987 https://doi.org/10.1016/0022-2836(87)90038-6
    # Wilke: Tien et al. 2013 https://doi.org/10.1371/journal.pone.0080635
    # Sander: Sander & Rost 1994 https://doi.org/10.1002/prot.340200303
    "Miller": {
        "A": 113.0,
        "R": 241.0,
        "N": 158.0,
        "D": 151.0,
        "C": 140.0,
        "Q": 189.0,
        "E": 183.0,
        "G": 85.0,
        "H": 194.0,
        "I": 182.0,
        "L": 180.0,
        "K": 211.0,
        "M": 204.0,
        "F": 218.0,
        "P": 143.0,
        "S": 122.0,
        "T": 146.0,
        "W": 259.0,
        "Y": 229.0,
        "V": 160.0,
    },
    "Wilke": {
        "A": 129.0,
        "R": 274.0,
        "N": 195.0,
        "D": 193.0,
        "C": 167.0,
        "Q": 225.0,
        "E": 223.0,
        "G": 104.0,
        "H": 224.0,
        "I": 197.0,
        "L": 201.0,
        "K": 236.0,
        "M": 224.0,
        "F": 240.0,
        "P": 159.0,
        "S": 155.0,
        "T": 172.0,
        "W": 285.0,
        "Y": 263.0,
        "V": 174.0,
    },
    "Sander": {
        "A": 106.0,
        "R": 248.0,
        "N": 157.0,
        "D": 163.0,
        "C": 135.0,
        "Q": 198.0,
        "E": 194.0,
        "G": 84.0,
        "H": 184.0,
        "I": 169.0,
        "L": 164.0,
        "K": 205.0,
        "M": 188.0,
        "F": 197.0,
        "P": 136.0,
        "S": 130.0,
        "T": 142.0,
        "W": 227.0,
        "Y": 222.0,
        "V": 142.0,
    },
}

class AccAgree:
    def __init__(self, query, target, acc_char="A"):
        query_len = len(query)
        target_len = len(target)

        assert query_len == target_len, "Cant compute acc agreement for sequences of different length"

        query_acc_count = 0
        target_acc_count = 0
        agree_count = 0
        for query_char, target_char in zip(query, target):
            if query_char == acc_char:
                query_acc_count += 1
            if target_char == acc_char:
                target_acc_count += 1

            if query_char == target_char:
                agree_count += 1

        self.percentAccQuery = query_acc_count / query_len
        self.percentAccTarget = target_acc_count / target_len
        self.percentAccAgree = agree_count / target_len


class MMCIFSeqs:
    def __init__(self, pdbid, chain, SQLITE_DATABASE_DIR,compute_acc=False, threshold=0.25):
        con = sqlite3.connect(SQLITE_DATABASE_DIR)
        cur = con.cursor()
        cur.execute(f'SELECT seqres, seqsolv, seqnums FROM mmCIF_seqs WHERE pdb_id = "{pdbid}_{chain}"')
        row = cur.fetchone()
        self.seqres = row[0]
        self.seqsolv = row[1]
        self.seqnums = row[2].split(" ")

        if compute_acc:
            cur.execute(f'SELECT acc FROM PDB_DSSP WHERE pdb_id = "{pdbid}_{chain}"')
            row = cur.fetchone()
            if row == None:
                self.acc = None
                self.rasa = None
                self.binaryrasa = None

            else:
                self.acc = row[0].split(" ")

                self.rasa = []
                for val, res in zip(self.acc, self.seqsolv):
                    if res == "?":
                        self.rasa.append("?")
                    elif val == "?":
                        self.rasa.append("?")
                    else:
                        self.rasa.append(float(val) / residue_max_acc["Wilke"][res])

                self.binaryrasa = []
                for val in self.rasa:
                    if val == "?":
                        self.binaryrasa.append("?")
                    elif val >= threshold:
                        self.binaryrasa.append("A")
                    else:
                        self.binaryrasa.append("B")

        else:
            self.acc = None
            self.rasa = None
            self.binaryrasa = None
