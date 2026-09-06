from stiff_medium.dependency_ledger import (
    ClaimTag,
    Risk,
    active_work_queue,
    all_entries,
    entries_by_tag,
    independent_prediction_count,
    tag_summary,
)


def test_ledger_has_required_status_classes():
    tags = tag_summary()
    assert tags[ClaimTag.PRIMITIVE.value] >= 1
    assert tags[ClaimTag.DERIVED.value] >= 1
    assert tags[ClaimTag.INHERITED.value] >= 1
    assert tags[ClaimTag.ANCHORED.value] >= 1
    assert tags[ClaimTag.CALIBRATED.value] >= 1
    assert tags[ClaimTag.FAILED.value] >= 1
    assert tags[ClaimTag.OPEN.value] >= 1


def test_independent_prediction_count_is_conservative():
    derived_entries = entries_by_tag(ClaimTag.DERIVED)
    assert independent_prediction_count() == len(derived_entries)
    assert independent_prediction_count() < len(all_entries())


def test_active_work_queue_prioritizes_blockers():
    queue = active_work_queue()
    assert queue
    first_risks = {entry.risk for entry in queue[:5]}
    assert Risk.BLOCKER in first_risks
    assert any(entry.claim_id == "F002" for entry in queue)
    assert any(entry.claim_id == "O003" for entry in queue)
