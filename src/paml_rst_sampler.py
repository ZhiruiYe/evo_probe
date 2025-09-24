#!/usr/bin/env python3
import re, random, argparse, sys
from collections import defaultdict

def parse_node_posteriors_from_rst_text(rst_text: str, node_id: int):
    header_pattern = re.compile(rf"Prob distribution at node\s+{node_id}\s*,\s*by site", re.IGNORECASE)
    m = header_pattern.search(rst_text)
    if not m:
        raise ValueError(f"Could not find node {node_id} header in rst text.")
    start_idx = m.end()
    site_line_pattern = re.compile(r"^\s*(\d+)\s+\d+\s+.+?:\s*(.*)$")
    aa_prob_pattern = re.compile(r"\b([A-Z])\(([-+]?\d*\.\d+|\d+)\)")
    lines = rst_text[start_idx:].splitlines()
    posteriors = []
    for line in lines:
        if not line.strip():
            continue
        if header_pattern.search(line):
            break
        msite = site_line_pattern.match(line)
        if msite:
            site_idx = int(msite.group(1))
            tail = msite.group(2)
            pairs = aa_prob_pattern.findall(tail)
            probs = {}
            for aa, val in pairs:
                probs[aa] = float(val)
            # normalize
            total = sum(probs.values())
            if total > 0:
                for k in list(probs.keys()):
                    probs[k] /= total
            posteriors.append((site_idx, probs))
        else:
            # continuation lines with more A(0.xxx) tokens
            pairs = aa_prob_pattern.findall(line)
            if pairs and posteriors:
                site_idx, probs = posteriors[-1]
                for aa, val in pairs:
                    probs[aa] = float(val)
                total = sum(probs.values())
                if total > 0:
                    for k in list(probs.keys()):
                        probs[k] /= total
                posteriors[-1] = (site_idx, probs)
            else:
                # end heuristics
                if line.strip().startswith("Prob distribution at node"):
                    break
                continue
    posteriors.sort(key=lambda x: x[0])
    return [p for (_, p) in posteriors]

def sample_sequences_from_posteriors(posteriors, k=10, seed=None, states=None):
    if seed is not None:
        random.seed(seed)
    if states is None:
        alphabet = sorted(set().union(*[set(d.keys()) for d in posteriors]))
    else:
        alphabet = list(states)
    seqs = []
    for s in range(k):
        chars = []
        for d in posteriors:
            probs = [d.get(a, 0.0) for a in alphabet]
            total = sum(probs)
            if total == 0:
                keys = list(d.keys())
                choice = random.choice(keys) if keys else 'X'
                chars.append(choice)
                continue
            probs = [p/total for p in probs]
            r = random.random()
            cum = 0.0
            choice = alphabet[-1]
            for a, p in zip(alphabet, probs):
                cum += p
                if r <= cum:
                    choice = a
                    break
            chars.append(choice)
        seqs.append("".join(chars))
    return seqs, alphabet

def save_fasta(seqs, name_prefix, out_path, wrap=80):
    with open(out_path, "w") as f:
        for i, s in enumerate(seqs, 1):
            f.write(f">{name_prefix}_sample{i}\n")
            for j in range(0, len(s), wrap):
                f.write(s[j:j+wrap] + "\n")
    return out_path

def main():
    ap = argparse.ArgumentParser(description="Sample ancestor sequences from PAML rst posterior block.")
    ap.add_argument("rst_path", help="Path to PAML rst file")
    ap.add_argument("--node", type=int, required=True, help="Ancestral node id (as shown in rst)")
    ap.add_argument("--k", type=int, default=10, help="Number of sequences to sample")
    ap.add_argument("--seed", type=int, default=None, help="Random seed")
    ap.add_argument("--out", default=None, help="Output FASTA path (default: node<N>_samples.fasta)")
    ap.add_argument("--alphabet", default="AA", choices=["AA","DNA","RNA"], help="Force alphabet order (optional)")
    args = ap.parse_args()

    with open(args.rst_path, "r", errors="ignore") as fh:
        rst_text = fh.read()
    post = parse_node_posteriors_from_rst_text(rst_text, args.node)
    if not post:
        print(f"No posterior entries parsed for node {args.node}.", file=sys.stderr)
        sys.exit(2)
    if args.alphabet == "AA":
        states = list("ACDEFGHIKLMNPQRSTVWY")
    elif args.alphabet == "DNA":
        states = list("ACGT")
    else:
        states = list("ACGU")
    seqs, alpha = sample_sequences_from_posteriors(post, k=args.k, seed=args.seed, states=states)
    out_path = args.out or f"node{args.node}_samples.fasta"
    save_fasta(seqs, f"node{args.node}", out_path)
    print(f"Wrote {len(seqs)} sequences to {out_path} with alphabet {''.join(alpha)}")

if __name__ == "__main__":
    main()

# 为 node 499 抽 100 条序列、固定随机种子，输出到指定 FASTA：
# python paml_rst_sampler.py your_run.rst --node 499 --k 100 --seed 42 --out node499_samples.fasta
# python src/paml_rst_sampler.py data/anc_prob.txt --node 507 --k 1000 --seed 42 --out data/sample/node507_anc4_samples.fasta
