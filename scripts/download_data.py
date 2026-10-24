#!/usr/bin/env python3
"""Download and freeze UniProt records; no network access needed for analysis."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

TAXA = {
    9606: 'Primates', 9598: 'Primates', 9595: 'Primates', 9601: 'Primates',
    9685: 'Carnivora', 9615: 'Carnivora', 9644: 'Carnivora', 9720: 'Carnivora',
    9913: 'Artiodactyla', 9940: 'Artiodactyla', 9925: 'Artiodactyla', 9823: 'Artiodactyla',
}
AA = set('ACDEFGHIKLMNPQRSTVWY')


def select_records(payload):
    selected = []
    for taxid, order in TAXA.items():
        candidates = []
        for entry in payload['results']:
            name = entry['proteinDescription'].get('recommendedName', {}).get('fullName', {}).get('value')
            sequence = entry['sequence']['value']
            if (entry['organism']['taxonId'] == taxid and name == 'Cytochrome b'
                    and entry['entryType'] == 'UniProtKB reviewed (Swiss-Prot)'
                    and not entry['proteinDescription'].get('flag')
                    and 350 <= len(sequence) <= 400 and set(sequence) <= AA):
                candidates.append(entry)
        if not candidates:
            raise ValueError(f'No eligible complete CYTB record for taxon {taxid}')
        entry = min(candidates, key=lambda e: e['primaryAccession'])
        if order not in entry['organism']['lineage']:
            raise ValueError(f'Expected order {order} absent for taxon {taxid}')
        selected.append(entry)
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('data/source'))
    args = parser.parse_args()
    query = ('protein_name:"Cytochrome b" AND ('
             + ' OR '.join(f'organism_id:{t}' for t in TAXA)
             + ') AND reviewed:true AND length:[350 TO 400]')
    session = requests.Session()
    session.mount('https://', HTTPAdapter(max_retries=Retry(
        total=4, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])))
    response = session.get('https://rest.uniprot.org/uniprotkb/search',
                           params={'query': query, 'format': 'json', 'size': 500}, timeout=90)
    response.raise_for_status()
    if 'next' in response.links:
        raise ValueError('Unexpected pagination; refusing an incomplete snapshot')
    payload = response.json()
    selected = select_records(payload)  # Validate everything before writing any output.
    metadata = []
    fasta = []
    for entry in selected:
        accession = entry['primaryAccession']
        organism = entry['organism']
        sequence = entry['sequence']['value']
        metadata.append({
            'accession': accession, 'taxon_id': organism['taxonId'],
            'organism': organism['scientificName'], 'order': TAXA[organism['taxonId']],
            'lineage': organism['lineage'], 'length': len(sequence),
            'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest(),
            'entry_version': entry['entryAudit']['entryVersion'],
            'sequence_version': entry['entryAudit']['sequenceVersion'],
            'url': f'https://www.uniprot.org/uniprotkb/{accession}/entry',
        })
        fasta.append(f'>{accession}\n{sequence}\n')
    args.output.mkdir(parents=True, exist_ok=True)
    raw = response.content
    (args.output / 'uniprot_response.json').write_bytes(raw)
    (args.output / 'sequences.fasta').write_text(''.join(fasta))
    (args.output / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    provenance = {
        'retrieved_at_utc': datetime.now(timezone.utc).isoformat(),
        'request_url': response.url, 'query': query,
        'uniprot_release': response.headers.get('X-UniProt-Release'),
        'uniprot_release_date': response.headers.get('X-UniProt-Release-Date'),
        'source_sha256': hashlib.sha256(raw).hexdigest(),
        'selected_count': len(selected), 'returned_count': len(payload['results']),
        'selection': 'Reviewed, exact Cytochrome b name, no fragment flag, exact taxon ID, '
                     '350–400 standard amino acids; lowest accession if multiple candidates.',
    }
    (args.output / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(provenance, indent=2))


if __name__ == '__main__':
    main()
