"""Recompute metadata-only candidate blockers; never open participant files."""
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from research.mortality.freeze30.adapter import strict_json
from research.mortality.freeze30.qualification import audit


def main():
    root = Path(__file__).resolve().parents[1] / 'research/mortality/freeze30'
    candidates = strict_json(root / 'CANDIDATE_MANIFEST.json')
    saved = strict_json(root / 'QUALIFICATION_RESULTS.json')
    actual = [dict(candidate_id=r['candidate_id'], **audit(r['qualification_metadata']))
              for r in candidates['candidates']]
    if (saved['candidates'] != actual or saved['candidate_count'] != len(actual)
            or saved['eligible_external_candidates'] != sum(r['frozen_scoring_allowed'] for r in actual)):
        raise ValueError('Qualification snapshot differs from executable metadata gate')
    print(json.dumps(dict(candidates_checked=len(actual),
        eligible_external_candidates=saved['eligible_external_candidates'],
        new_external_outcome_comparisons=saved['new_external_outcome_comparisons'],
        snapshot_exact=True)))


if __name__ == '__main__':
    main()
