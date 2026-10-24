#!/usr/bin/env python3
"""Offline MAFFT → IQ-TREE → unrooted taxonomy split assessment."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

from Bio import Phylo, SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

COLORS = {'Primates': '#2563a6', 'Carnivora': '#b44d28', 'Artiodactyla': '#34805c'}


def read_sequences(path):
    records = list(SeqIO.parse(path, 'fasta'))
    ids = [r.id for r in records]
    if not records or len(ids) != len(set(ids)):
        raise ValueError('Empty FASTA or duplicate sequence identifiers')
    return {r.id: str(r.seq) for r in records}


def validate_alignment(raw, aligned):
    if set(raw) != set(aligned):
        raise ValueError('Alignment identifiers differ from inputs')
    if len({len(s) for s in aligned.values()}) != 1:
        raise ValueError('Alignment rows have different lengths')
    if any(aligned[k].replace('-', '') != v for k, v in raw.items()):
        raise ValueError('Alignment changed an ungapped sequence')


def validate_tips(tree, expected):
    tips = [t.name for t in tree.get_terminals()]
    if len(tips) != len(set(tips)) or set(tips) != set(expected):
        raise ValueError('Tree tips differ from inputs or contain duplicates')


def edge_splits(tree):
    """Return nontrivial edge bipartitions without assuming a biological root."""
    all_tips = frozenset(t.name for t in tree.get_terminals())
    splits = {}
    for clade in tree.find_clades():
        subset = frozenset(t.name for t in clade.get_terminals())
        complement = all_tips - subset
        if min(len(subset), len(complement)) < 2:
            continue
        key = frozenset((subset, complement))
        if clade.confidence is not None:
            splits[key] = float(clade.confidence)
        else:
            splits.setdefault(key, None)
    return splits


def assess_orders(tree, metadata):
    splits = edge_splits(tree)
    all_tips = frozenset(t.name for t in tree.get_terminals())
    rows = []
    for order in sorted({m['order'] for m in metadata}):
        tips = frozenset(m['accession'] for m in metadata if m['order'] == order)
        key = frozenset((tips, all_tips - tips))
        rows.append({'order': order, 'taxa': len(tips),
                     'split_present': key in splits, 'support_percent': splits.get(key)})
    return rows


def draw_tree(tree, metadata, output):
    labels = {m['accession']: m for m in metadata}
    fig, ax = plt.subplots(figsize=(13, 8))
    Phylo.draw(tree, axes=ax, do_show=False,
               label_func=lambda c: (f"{labels[c.name]['organism']}  [{c.name}]"
                                     if c.name in labels else None),
               label_colors={f"{v['organism']}  [{k}]": COLORS[v['order']]
                             for k, v in labels.items()},
               show_confidence=True)
    ax.set_title('Mitochondrial cytochrome b: 12 mammals', loc='left', fontsize=17, pad=20)
    ax.set_xlabel('Branch length (expected amino-acid substitutions per site)')
    ax.set_ylabel('')
    ax.set_yticks([])
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.legend(handles=[Line2D([0], [0], color=color, lw=3, label=order)
                       for order, color in COLORS.items()], loc='lower right', frameon=False)
    fig.text(0.04, 0.025, 'IQ-TREE ML · 1,000 ultrafast bootstrap replicates + NNI · seed 22\n'
             'Internal labels: bootstrap support (%). Display orientation is arbitrary; no ancestral root inferred.',
             fontsize=9, color='#444444')
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(output / 'tree.png', dpi=180)
    fig.savefig(output / 'tree.svg')
    plt.close(fig)


def summarize(source, output):
    metadata = json.loads((source / 'metadata.json').read_text())
    raw = read_sequences(source / 'sequences.fasta')
    aligned = read_sequences(output / 'alignment.fasta')
    validate_alignment(raw, aligned)
    if {m['accession'] for m in metadata} != set(raw):
        raise ValueError('Metadata identifiers differ from sequence inputs')
    tree = Phylo.read(output / 'cytb.treefile', 'newick')
    validate_tips(tree, raw)
    bootstrap_count = 0
    for boot in Phylo.parse(output / 'cytb.ufboot', 'newick'):
        validate_tips(boot, raw)
        bootstrap_count += 1
    splits = edge_splits(tree)
    supports = list(splits.values())
    if len(splits) != len(raw) - 3 or any(s is None or not 0 <= s <= 100 for s in supports):
        raise ValueError('Missing or invalid support on unrooted internal edges')
    columns = list(zip(*aligned.values()))
    # Count residues only; a gap is not an amino-acid state.
    informative = sum(sum(n >= 2 for n in Counter(c for c in col if c != '-').values()) >= 2
                      for col in columns)
    variable = sum(len(set(col) - {'-'}) > 1 for col in columns)
    report = (output / 'cytb.iqtree').read_text()
    model = re.search(r'Model of substitution: (.+)', report).group(1).strip()
    orders = assess_orders(tree, metadata)
    lengths = [len(s) for s in raw.values()]
    metrics = {
        'sequence_count': len(raw), 'unique_sequence_count': len(set(raw.values())),
        'accession_count': len({m['accession'] for m in metadata}),
        'taxon_count': len({m['taxon_id'] for m in metadata}),
        'species_count': len({' '.join(m['organism'].split()[:2]) for m in metadata}),
        'order_counts': dict(Counter(m['order'] for m in metadata)),
        'input_length_min': min(lengths), 'input_length_max': max(lengths),
        'alignment_columns': len(columns),
        'gap_cells': sum(s.count('-') for s in aligned.values()),
        'alignment_cells': len(columns) * len(raw),
        'variable_columns': variable, 'parsimony_informative_columns': informative,
        'tree_tip_count': len(tree.get_terminals()),
        'bootstrap_tree_count': bootstrap_count, 'internal_edge_count': len(splits),
        'internal_support_min': min(supports), 'internal_support_max': max(supports),
        'selected_model': model, 'seed': 22, 'threads': 1,
        'orders_with_split': sum(row['split_present'] for row in orders),
        'orders_assessed': len(orders), 'order_assessments': orders,
    }
    (output / 'metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
    with (output / 'order_splits.tsv').open('w') as handle:
        handle.write('order\ttaxa\tsplit_present\tsupport_percent\n')
        for row in orders:
            handle.write('\t'.join(str(row[k]) for k in row) + '\n')
    draw_tree(tree, metadata, output)
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('data/source'))
    parser.add_argument('--output', type=Path, default=Path('results'))
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if not args.summarize_only:
        if (args.output / 'cytb.treefile').exists():
            parser.error('Tree already exists: use a fresh --output directory or --summarize-only')
        commands = [
            ['mafft', '--auto', '--thread', '1', str(args.source / 'sequences.fasta')],
            ['iqtree2', '-s', str(args.output / 'alignment.fasta'), '-st', 'AA',
             '-m', 'MFP', '-msub', 'mitochondrial', '-B', '1000', '-bnni',
             '-wbt', '-seed', '22', '-T', '1', '--prefix', str(args.output / 'cytb')],
        ]
        started = time.perf_counter()
        with (args.output / 'alignment.fasta').open('w') as out, (args.output / 'mafft.log').open('w') as err:
            subprocess.run(commands[0], stdout=out, stderr=err, check=True)
        validate_alignment(read_sequences(args.source / 'sequences.fasta'),
                           read_sequences(args.output / 'alignment.fasta'))
        with (args.output / 'iqtree.stdout.log').open('w') as log:
            subprocess.run(commands[1], stdout=log, stderr=subprocess.STDOUT, check=True)
        versions = {}
        for command in [['mafft', '--version'], ['iqtree2', '--version']]:
            completed = subprocess.run(command, capture_output=True, text=True, check=True)
            versions[command[0]] = (completed.stdout + completed.stderr).strip()
        run = {'commands': commands, 'versions': versions,
               'analysis_wall_seconds': round(time.perf_counter() - started, 3)}
        (args.output / 'run.json').write_text(json.dumps(run, indent=2) + '\n')
    metrics = summarize(args.source, args.output)
    if metrics['bootstrap_tree_count'] != 1000:
        raise ValueError('Expected exactly 1,000 bootstrap trees')
    files = sorted(p for folder in (args.source, args.output) for p in folder.iterdir()
                   if p.is_file() and p.name != 'checksums.sha256')
    (args.output / 'checksums.sha256').write_text(''.join(
        f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.as_posix()}\n' for p in files))
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    main()
