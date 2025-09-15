#!/usr/bin/env python3
import sys, re, argparse

def read_ml_sequences(rst):
    nodes = {}
    with open(rst) as f:
        lines = f.readlines()
    i=0
    while i < len(lines):
        m = re.match(r"\s*node\s+#?(\d+)\s*(.*)", lines[i], re.I)
        if m:
            nid = int(m.group(1))
            seq = []
            i+=1
            while i < len(lines) and lines[i].strip():
                chunk = re.sub(r"\s|\d", "", lines[i])
                seq.append(chunk.strip())
                i+=1
            nodes[nid] = ("node_%d"%nid, "".join(seq))
        else:
            i+=1
    return nodes

def write_fasta(d, path, prefix="N"):
    with open(path,"w") as out:
        for nid,(name,seq) in sorted(d.items()):
            out.write(f">{prefix}{nid}_{name}\n")
            for j in range(0,len(seq),70):
                out.write(seq[j:j+70]+"\n")

def make_altanc(rst, p2=0.20):
    with open(rst) as f:
        lines=f.readlines()
    # parse ML
    nodes = read_ml_sequences(rst)
    # parse posterior tables
    alt_nodes = {}
    i=0
    while i < len(lines):
        m = re.match(r"\s*Prob\s+distribution.*node\s+#?(\d+)", lines[i], re.I)
        if m:
            nid = int(m.group(1))
            ml = list(nodes[nid][1])
            pos = 0
            i+=1
            while i < len(lines) and lines[i].strip():
                line = lines[i].strip()
                pairs = re.findall(r"([ACDEFGHIKLMNPQRSTVWY\*])\s+([01]\.\d+|0\.\d+)", line)
                if pairs:
                    probs = [(aa, float(p)) for aa,p in pairs]
                    probs.sort(key=lambda x: x[1], reverse=True)
                    if len(probs)>=2 and probs[1][1] >= p2 and pos < len(ml) and ml[pos] != "-":
                        ml[pos] = probs[1][0]
                    pos += 1
                i+=1
            alt_nodes[nid] = (nodes[nid][0], "".join(ml))
        else:
            i+=1
    return nodes, alt_nodes

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("rst", help="codeml 'rst' file")
    ap.add_argument("--p2", type=float, default=0.20, help="2nd-highest posterior threshold (default 0.20)")
    args=ap.parse_args()
    ml, alt = make_altanc(args.rst, args.p2)
    write_fasta(ml, "anc_ml.fasta")
    if alt:
        write_fasta(alt, f"anc_alt_p{int(args.p2*100)}.fasta")
    print("Written: anc_ml.fasta", "and", f"anc_alt_p{int(args.p2*100)}.fasta")

if __name__ == "__main__":
    main()
