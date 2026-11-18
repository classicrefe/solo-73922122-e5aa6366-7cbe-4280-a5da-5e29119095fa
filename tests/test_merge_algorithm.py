"""Unit tests for overlap-finding rules, reverse complement and qualities."""
import random

import pytest

from app import config, merging, processing


def _rec(seq, qual=None):
    if qual is None:
        qual = b"I" * len(seq)
    return processing.Record(b"@r", seq, b"+", qual)


def _rnd(n, seed):
    rng = random.Random(seed)
    return bytes(rng.choice(b"ACGT") for _ in range(n))


def test_reverse_complement_basic_and_n():
    assert merging.reverse_complement(b"ACGTN") == b"NACGT"
    assert merging.reverse_complement(b"AAAA") == b"TTTT"
    # involution
    s = b"ACGTTGCATGCA"
    assert merging.reverse_complement(merging.reverse_complement(s)) == s


def test_mismatch_budget_boundaries():
    # 5% inclusive, floor-based.
    assert merging._max_mismatches(12) == 0
    assert merging._max_mismatches(19) == 0
    assert merging._max_mismatches(20) == 1
    assert merging._max_mismatches(39) == 1
    assert merging._max_mismatches(40) == 2
    assert merging._max_mismatches(100) == 5
    assert merging._max_mismatches(150) == 7


def test_minimum_overlap_is_12():
    assert merging.find_overlap(b"X" * 30 + b"A" * 11,
                                b"A" * 11 + b"Y" * 30) is None
    ok = merging.find_overlap(b"X" * 30 + b"A" * 12,
                              b"A" * 12 + b"Y" * 30)
    # Only one perfect position because the flanking X/Y bases break shorter
    # alignments.
    assert ok == 12


def test_no_overlap_returns_none():
    assert merging.find_overlap(b"A" * 50, b"G" * 50) is None


def test_multiple_qualifying_positions_not_merged():
    # Homopolymer: every length 12..50 is a perfect overlap.
    assert merging.find_overlap(b"A" * 50, b"A" * 50) is None
    result = merging.merge_pair(_rec(b"A" * 50),
                                _rec(merging.reverse_complement(b"A" * 50)))
    assert result.merged is False


def test_unique_overlap_and_merged_sequence_layout():
    r1 = _rnd(50, seed=1)
    r2_rc = r1[-25:] + b"T" * 25          # unique 25 bp overlap, 25 bp tail
    r2 = _rec(merging.reverse_complement(r2_rc))
    result = merging.merge_pair(_rec(r1), r2)
    assert result.merged is True
    assert result.overlap_len == 25
    assert result.seq == r1 + b"T" * 25
    assert len(result.seq) == len(result.qual) == 75


def test_mismatch_inside_tolerance_merges_above_fails():
    r1 = _rnd(40, seed=2)
    # One substitution inside a 40 bp overlap: 1/40 = 2.5% <= 5%.
    bad = b"G" if r1[39] != ord("G") else b"T"
    r2_rc = r1[:39] + bad
    assert merging.find_overlap(r1, r2_rc) == 40
    # Two substitutions: 2/40 = 5% still allowed (inclusive).
    r2_rc2 = bytearray(r1)
    r2_rc2[38] = ord("G") if r1[38] != ord("G") else ord("T")
    r2_rc2[39] = bad[0]
    assert merging.find_overlap(r1, bytes(r2_rc2)) == 40
    # Three substitutions: 7.5% -> rejected everywhere -> no merge.
    r2_rc3 = bytearray(r2_rc2)
    r2_rc3[37] = ord("G") if r1[37] != ord("G") else ord("T")
    assert merging.find_overlap(r1, bytes(r2_rc3)) is None


def test_same_base_qualities_add_and_cap_at_40():
    r1 = _rnd(50, seed=3)
    r2_rc = r1[-25:] + b"T" * 25
    r2_raw = merging.reverse_complement(r2_rc)
    low = chr(33 + 10).encode()  # Phred 10 -> '+'
    result = merging.merge_pair(_rec(r1, low * 50), _rec(r2_raw, low * 50))
    ov = 50 - result.overlap_len
    # 10 + 10 = Phred 20 -> char 53 ('5'), below the cap.
    assert result.qual[ov] == 33 + 20
    high = chr(33 + 30).encode()  # Phred 30; 30+30=60 -> capped 40
    capped = merging.merge_pair(_rec(r1, high * 50), _rec(r2_raw, high * 50))
    assert capped.qual[ov] == 33 + config.MERGE_QUAL_CAP


def test_discordant_base_quality_rules():
    r1 = _rnd(40, seed=4)
    bad = b"G" if r1[39] != ord("G") else b"T"
    r2_rc = r1[:39] + bad
    # R2-prime position 39 (the mismatch) corresponds to raw-R2 index 0.
    # R1 has the higher quality at the mismatch: its base wins, its qual kept.
    q1 = b"I" * 39 + b"K"        # K = Phred 42 on R1 position 39
    q2 = b"C" + b"I" * 39        # C = Phred 34 lands on R2' position 39
    res = merging.merge_pair(
        _rec(r1, q1), _rec(merging.reverse_complement(r2_rc), q2))
    assert res.merged is True
    assert res.seq[-1] == r1[39]
    assert res.qual[-1] == ord("K")
    # R2 higher: reverse-complemented R2 base wins.
    res2 = merging.merge_pair(
        _rec(r1, b"I" * 39 + b"C"),
        _rec(merging.reverse_complement(r2_rc), b"K" + b"I" * 39))
    assert res2.seq[-1] == bad[0]
    assert res2.qual[-1] == ord("K")
    # Equal quality: N, carrying the lower (equal) quality.
    res3 = merging.merge_pair(
        _rec(r1, b"I" * 40),
        _rec(merging.reverse_complement(r2_rc), b"I" * 40))
    assert res3.seq[-1] == ord("N")
    assert res3.qual[-1] == ord("I")


def test_r2_quality_string_is_reversed_during_merge():
    # A marker quality at raw-R2 index 39 must land on R2' position 0.
    r1 = _rnd(40, seed=5)
    r2_rc = r1  # full-length unique overlap
    r2_raw = merging.reverse_complement(r2_rc)
    q2 = b"I" * 39 + b"K"
    res = merging.merge_pair(_rec(r1, b"#" * 40), _rec(r2_raw, q2))
    assert res.merged is True
    # R2' position 0 gets the K; R1 side is '#'=2. Bases agree there, so the
    # merged Phred is min(2 + 42, 40) = 40 -> 'I'.
    assert res.qual[0] == ord("I")


def test_find_overlap_matches_brute_force_fuzz():
    rng = random.Random(99)
    for _ in range(400):
        n = rng.randint(40, 150)
        s1 = bytes(rng.choice(b"ACGT") for _ in range(n))
        L = rng.randint(12, n)
        s2_tail = bytearray(s1[-L:])
        for pos in rng.sample(range(L), (L * 5) // 100):
            s2_tail[pos] = rng.choice(
                [b for b in (ord("A"), ord("C"), ord("G"), ord("T"))
                 if b != s2_tail[pos]])
        s2 = bytes(s2_tail) + bytes(
            rng.choice(b"ACGT") for _ in range(n - L))
        qualifying = []
        for length in range(12, min(len(s1), len(s2)) + 1):
            diffs = sum(1 for a, b in zip(s1[-length:], s2[:length])
                        if a != b)
            if diffs <= (length * 5) // 100:
                qualifying.append(length)
        expected = qualifying[0] if len(qualifying) == 1 else None
        assert merging.find_overlap(s1, s2) == expected
