"""Study inclusion decisions, kept separate from measured model validity."""

import csv


REASONS = ('manual_review', 'segmentation', 'registration', 'support_contact',
           'insufficient_coverage', 'acceptable_after_inspection', 'other')
FIELDS = ('model_id', 'site', 'model_path', 'automatic_qc_status', 'automatic_qc_reasons',
          'manual_decision', 'manual_reason', 'manual_note', 'reviewer', 'review_timestamp',
          'final_inclusion', 'solver_eligible', 'qc_policy_version', 'measurement_source')


def identity(row):
    """Return the site and model ID used to match review decisions."""
    return (row['site'], row['model_id'])


def inclusion_row(measurement, decision=None):
    """Combine automatic QC with a manual decision, retaining both for provenance."""
    decision = decision or {}
    manual = decision.get('manual_decision') or 'automatic'
    if manual not in ('automatic', 'include', 'exclude'):
        raise ValueError(f'Invalid manual decision: {manual}')
    row = {key: measurement.get(key, '') for key in FIELDS}
    status = measurement.get('qc_status', 'review')
    row.update(automatic_qc_status=status, automatic_qc_reasons=measurement.get('qc_reasons', ''),
               manual_decision=manual,
               manual_reason=(decision.get('manual_reason') or 'manual_review') if manual != 'automatic' else '',
               manual_note=decision.get('manual_note', '') if manual != 'automatic' else '',
               reviewer=decision.get('reviewer', ''), review_timestamp=decision.get('review_timestamp', ''),
               final_inclusion=manual if manual != 'automatic' else {'pass': 'include', 'fail': 'exclude'}.get(status, 'pending'),
               solver_eligible='False' if measurement.get('generation_status') == 'too_short' else '')
    return row


def load_decisions(path):
    """Read exported decisions, rejecting duplicate or malformed identities."""
    decisions = {}
    with open(path, newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        if not {'model_id', 'site', 'manual_decision'} <= set(reader.fieldnames or []):
            raise ValueError('Review CSV requires model_id, site, and manual_decision')
        for row in reader:
            if not row['model_id'] or row['site'] not in ('hip', 'spine', 'unknown'):
                raise ValueError('Invalid model identity in review CSV')
            key = identity(row)
            if key in decisions:
                raise ValueError(f'Duplicate review identity: {key}')
            inclusion_row(row, row)  # Validate the decision before applying it.
            decisions[key] = row
    return decisions


def write_inclusion(rows, path, decisions=None):
    """Write study decisions separately from the model measurement CSV."""
    decisions = decisions or {}
    with open(path, 'w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(inclusion_row(row, decisions.get(identity(row))) for row in rows)
