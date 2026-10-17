#!/usr/bin/env python3
"""
build_tf_list.py

Builds a D. melanogaster transcription factor list (FBgn IDs) for Scribe
regulators: genes annotated to GO:0003700 "DNA-binding transcription factor
activity" or any of its is_a descendants (e.g. GO:0000981, GO:0001228,
GO:0001227). Annotations with a NOT qualifier are dropped.

Needs internet (run on the cluster login node) unless --gaf/--obo point to
local copies.

Usage:
    python build_tf_list.py \
        --release FB2025_04 \
        --t2g ../reference/KB_Drosophila_melanogaster_wMel_combined_16S/transcripts_to_genes.txt \
        --out ../cell_culture_system/config/gene_lists/dmel_TF_GO0003700_FBgn.txt
"""
import argparse
import gzip
import os
import urllib.request
from collections import defaultdict

ROOT = "GO:0003700"
OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
GAF_URL = "https://s3ftp.flybase.org/releases/{release}/precomputed_files/go/gene_association.fb.gz"


def fetch(url, dest):
    if not os.path.exists(dest):
        print(f"Downloading {url}")
        urllib.request.urlretrieve(url, dest)
    return dest


def descendants(obo_path, root):
    children = defaultdict(set)
    cur, obsolete = None, False
    with open(obo_path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                cur, obsolete = None, False
            elif line.startswith("id: GO:"):
                cur = line[4:]
            elif line == "is_obsolete: true":
                obsolete = True
            elif line.startswith("is_a: ") and cur and not obsolete:
                children[line[6:16]].add(cur)
    terms, stack = {root}, [root]
    while stack:
        for c in children[stack.pop()]:
            if c not in terms:
                terms.add(c)
                stack.append(c)
    return terms


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--release", default="FB2025_04", help="FlyBase release for the GAF")
    p.add_argument("--gaf", help="local gene_association.fb(.gz); skips download")
    p.add_argument("--obo", help="local go-basic.obo; skips download")
    p.add_argument("--t2g", help="kb transcripts_to_genes.txt, to report coverage in your reference")
    p.add_argument("--exclude-iea", action="store_true", help="drop electronic (IEA) annotations")
    p.add_argument("--out", required=True)
    p.add_argument("--workdir", default=".")
    args = p.parse_args()

    obo = args.obo or fetch(OBO_URL, os.path.join(args.workdir, "go-basic.obo"))
    gaf = args.gaf or fetch(GAF_URL.format(release=args.release),
                            os.path.join(args.workdir, f"gene_association_{args.release}.fb.gz"))

    terms = descendants(obo, ROOT)
    print(f"{ROOT} + descendants: {len(terms)} GO terms")

    genes = {}
    opener = gzip.open if gaf.endswith(".gz") else open
    with opener(gaf, "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or not f[1].startswith("FBgn"):
                continue
            if "NOT" in f[3] or f[4] not in terms:
                continue
            if args.exclude_iea and f[6] == "IEA":
                continue
            genes[f[1]] = f[2]
    print(f"TF genes: {len(genes)}")

    if args.t2g:
        ref = set()
        with open(args.t2g) as fh:
            for line in fh:
                parts = line.split("\t")
                if len(parts) > 1:
                    ref.add(parts[1].strip())
        missing = sorted(set(genes) - ref)
        print(f"In your reference t2g: {len(genes) - len(missing)} / {len(genes)}")
        if missing:
            print(f"  Not in reference (ID changes between releases?): {', '.join(missing[:15])}"
                  + (" ..." if len(missing) > 15 else ""))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(f"# D. melanogaster TFs: {ROOT} + is_a descendants, NOT excluded, "
                 f"IEA {'excluded' if args.exclude_iea else 'included'}; GAF {os.path.basename(gaf)}\n")
        for g in sorted(genes):
            fh.write(g + "\n")
    with open(os.path.splitext(args.out)[0] + "_symbols.tsv", "w") as fh:
        for g in sorted(genes):
            fh.write(f"{g}\t{genes[g]}\n")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
