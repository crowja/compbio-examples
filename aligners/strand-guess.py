#! /usr/bin/env python3



def strand_guess(aligner, q, t, minpct=50):
    """Guess the strand of t (target) that q (query) lies on.

    Given a Bio.Align aligner and two Bio.SeqRecords q and t, guess which strand
    of t that q lies on. The approach is to align both q and q.reverse_complement()
    to t, seeing which scores higher. The score has to be at least minpct of the
    maximum possible score, default 50 percent.

    ARGUMENTS
    aligner:   A Bio.Align aligner.
    q:         Query sequence as a Bio.SeqRecord.
    t:         Target sequence as a Bio.SeqRecord.
    minpct:    Score of best alignment between q and t must be at least minpct
               percent of the maximum possible score.

    RETURNS
    > 0:       q appears to lie on the forward strand of t.
    < 0:       q appears to lie on the reverse strand of t.
    otherwise: unable to determine which strand of t.
    """

    score_max = min(aligner.score(t.seq, t.seq), aligner.score(q.seq, q.seq))
    score_f = aligner.score(q.seq, t.seq)
    score_r = aligner.score(q.reverse_complement().seq, t.seq)

    if score_f > score_r and score_f >= minpct * score_max / 100:
        return 1
    elif score_r > score_f and score_r >= minpct * score_max / 100:
        return -1
    else:
        return 0


if __name__ == "__main__":
    # Do a few checks

    from Bio import Align
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"

    # x aligns with the forward strand of y with a small insertion
    x = SeqRecord(
        Seq("CATCAAGCCCTGTGAGCTGAAAAACTTAGCCGGGGACAATGGAGGTGCTGGGTTTGGTGGATATTCCGACATG")
    )
    y = SeqRecord(Seq("GCTGAACTTAGCCGGGGACAATG"))
    print(f"Test 1: Got {strand_guess(aligner, x, y)}, expecting 1")

    # x aligns with the forward strand of y with a small insertion
    x = SeqRecord(Seq("GCTGAACTTAGCCGGGGACAATG"))
    y = SeqRecord(
        Seq("CATCAAGCCCTGTGAGCTGAAAAACTTAGCCGGGGACAATGGAGGTGCTGGGTTTGGTGGATATTCCGACATG")
    )
    print(f"Test 2: Got {strand_guess(aligner, x, y)}, expecting 1")

    # x aligns with the reverse strand of y with a small insertion
    x = SeqRecord(
        Seq("CATGTCGGAATATCCACCAAACCCAGCACCTCCATTGTCCCCGGCTAAGTTTTTCAGCTCACAGGGCTTGATG")
    )
    y = SeqRecord(Seq("GCTGAACTTAGCCGGGGACAATG"))
    print(f"Test 3: Got {strand_guess(aligner, x, y)}, expecting -1")

    # x aligns with the reverse strand of y with a small insertion
    x = SeqRecord(Seq("GCTGAACTTAGCCGGGGACAATG"))
    y = SeqRecord(
        Seq("CATGTCGGAATATCCACCAAACCCAGCACCTCCATTGTCCCCGGCTAAGTTTTTCAGCTCACAGGGCTTGATG")
    )
    print(f"Test 4: Got {strand_guess(aligner, x, y)}, expecting -1")
