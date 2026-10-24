import copy
import hashlib
import json
from io import StringIO
from pathlib import Path
import subprocess
import sys

from Bio import Phylo
import pytest

from scripts.download_data import AA, TAXA, select_records
from scripts.run_analysis import (
    assess_orders, edge_splits, read_sequences, validate_alignment, validate_tips,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/source'
RESULTS = ROOT / 'results'


@pytest.fixture(scope='module')
def raw_payload():
    return json.loads((SOURCE / 'uniprot_response.json').read_text())


def test_frozen_source_matches_metadata_and_sequences(raw_payload):
    records = select_records(raw_payload)
    metadata = json.loads((SOURCE / 'metadata.json').read_text())
    fasta = read_sequences(SOURCE / 'sequences.fasta')
    assert len(records) == len(metadata) == len(fasta) == 12
    assert len({m['taxon_id'] for m in metadata}) == 12
    assert len({' '.join(m['organism'].split()[:2]) for m in metadata}) == 12
    by_id = {m['accession']: m for m in metadata}
    for entry in records:
        accession = entry['primaryAccession']
        seq = entry['sequence']['value']
        m = by_id[accession]
        assert fasta[accession] == seq
        assert 350 <= len(seq) <= 400 and set(seq) <= AA
        assert m['sequence_sha256'] == hashlib.sha256(seq.encode()).hexdigest()
        assert m['lineage'] == entry['organism']['lineage']
        assert m['order'] in m['lineage']
        assert m['entry_version'] == entry['entryAudit']['entryVersion']
        assert m['sequence_version'] == entry['entryAudit']['sequenceVersion']
    for order in set(TAXA.values()):
        assert sum(m['order'] == order for m in metadata) == 4
    provenance = json.loads((SOURCE / 'provenance.json').read_text())
    assert provenance['source_sha256'] == hashlib.sha256(
        (SOURCE / 'uniprot_response.json').read_bytes()).hexdigest()
    assert provenance['selected_count'] == provenance['returned_count'] == 12
    assert provenance['uniprot_release']


@pytest.mark.parametrize('mutation', ['missing', 'fragment', 'invalid_residue', 'wrong_lineage', 'unreviewed'])
def test_download_rejects_ineligible_source(raw_payload, mutation):
    payload = copy.deepcopy(raw_payload)
    entry = payload['results'][0]
    if mutation == 'missing':
        payload['results'].pop(0)
    elif mutation == 'fragment':
        entry['proteinDescription']['flag'] = 'Fragment'
    elif mutation == 'invalid_residue':
        entry['sequence']['value'] = 'X' + entry['sequence']['value'][1:]
    elif mutation == 'wrong_lineage':
        entry['organism']['lineage'] = []
    else:
        entry['entryType'] = 'UniProtKB unreviewed (TrEMBL)'
    with pytest.raises(ValueError):
        select_records(payload)


@pytest.mark.parametrize('aligned', [
    {'a': 'AC-', 'b': 'ADE'},  # lost a residue
    {'a': 'ACD', 'c': 'ADE'},  # changed identifier
    {'a': 'ACD', 'b': 'A-DE'},  # unequal lengths
])
def test_alignment_detects_corruption(aligned):
    with pytest.raises(ValueError):
        validate_alignment({'a': 'ACD', 'b': 'ADE'}, aligned)


def test_alignment_accepts_gaps_without_changing_residues():
    validate_alignment({'a': 'ACD', 'b': 'ADE'}, {'a': 'A-CD', 'b': 'AD-E'})


def test_duplicate_fasta_rejected(tmp_path):
    path = tmp_path / 'duplicate.fasta'
    path.write_text('>a\nACD\n>a\nACE\n')
    with pytest.raises(ValueError):
        read_sequences(path)


def test_tree_tip_validation_rejects_duplicate():
    tree = Phylo.read(StringIO('(a,a,b);'), 'newick')
    with pytest.raises(ValueError):
        validate_tips(tree, {'a', 'b', 'c'})


def test_order_split_is_independent_of_display_root():
    tree = Phylo.read(StringIO('((a,b)91,c,(d,e)87);'), 'newick')
    metadata = [{'accession': tip, 'order': 'A' if tip in 'ab' else 'B'} for tip in 'abcde']
    before = assess_orders(tree, metadata)
    # Same edge-labeled topology, serialized with a display root inside the target
    # group. Avoid Bio.Phylo.root_with_outgroup: it relocates clades but does not
    # preserve edge confidence labels. Rooted monophyly alone would fail here.
    tree = Phylo.read(StringIO('(a,(b,(c,(d,e)87)91));'), 'newick')
    assert assess_orders(tree, metadata) == before
    assert before[0]['split_present'] and before[0]['support_percent'] == 91
    assert before[1]['split_present']


def test_nonmonophyletic_group_reported_as_absent():
    tree = Phylo.read(StringIO('((a,b)91,c,(d,e)87);'), 'newick')
    metadata = [{'accession': t, 'order': 'A' if t in 'ac' else 'B'} for t in 'abcde']
    assert all(not row['split_present'] and row['support_percent'] is None
               for row in assess_orders(tree, metadata))


def check_outputs(output):
    raw = read_sequences(SOURCE / 'sequences.fasta')
    validate_alignment(raw, read_sequences(output / 'alignment.fasta'))
    tree = Phylo.read(output / 'cytb.treefile', 'newick')
    validate_tips(tree, raw)
    bootstrap = list(Phylo.parse(output / 'cytb.ufboot', 'newick'))
    assert len(bootstrap) == 1000
    for boot in bootstrap:
        validate_tips(boot, raw)
    splits = edge_splits(tree)
    assert len(splits) == 9
    assert all(v is not None and 0 <= v <= 100 for v in splits.values())
    # Independent bootstrap-tree frequency check against IQ-TREE edge labels.
    for split, support in splits.items():
        side, other = tuple(split)
        occurrences = 0
        for boot in bootstrap:
            clade_sets = [{t.name for t in clade.get_terminals()} for clade in boot.find_clades()]
            occurrences += set(side) in clade_sets or set(other) in clade_sets
        assert abs(support - 100 * occurrences / len(bootstrap)) <= 0.51
    metrics = json.loads((output / 'metrics.json').read_text())
    assert metrics['sequence_count'] == metrics['tree_tip_count'] == 12
    assert metrics['accession_count'] == metrics['taxon_count'] == metrics['species_count'] == 12
    assert metrics['bootstrap_tree_count'] == 1000
    assert metrics['orders_assessed'] == len(metrics['order_assessments']) == 3
    assert metrics['internal_edge_count'] == 9
    assert len((output / 'order_splits.tsv').read_text().splitlines()) == 4
    for name in ['tree.png', 'tree.svg']:
        assert (output / name).stat().st_size > 1000
    assert (output / 'tree.png').read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    for line in (output / 'checksums.sha256').read_text().splitlines():
        checksum, name = line.split('  ', 1)
        path = Path(name)
        if not path.is_absolute():
            path = ROOT / path
        assert hashlib.sha256(path.read_bytes()).hexdigest() == checksum


def test_committed_real_outputs():
    check_outputs(RESULTS)


@pytest.mark.integration
def test_live_mafft_iqtree_pipeline(tmp_path):
    output = tmp_path / 'analysis'
    subprocess.run([sys.executable, '-m', 'scripts.run_analysis',
                    '--source', str(SOURCE), '--output', str(output)],
                   cwd=ROOT, check=True, capture_output=True, text=True, timeout=600)
    check_outputs(output)
