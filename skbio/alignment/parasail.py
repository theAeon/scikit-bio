from skbio.sequence import DNA, RNA, Protein
from skbio.alignment._tabular_msa import TabularMSA

import parasail


class SubstitutionMatrix(object):
    """ Wrapper around a built-in Parasail substitution matrix.
    """

    def __init__(self, parasail_matrix):
        self._matrix = parasail_matrix

    @classmethod
    def from_name(cls, name):
        matrix = getattr(parasail, name)
        return cls(matrix)

    @classmethod
    def from_match_mismatch(cls, match, mismatch, alphabet='ACGTU'):
        matrix = parasail.matrix_create(alphabet, match, mismatch)
        return cls(matrix)

    @classmethod
    def from_dict(cls, d):
        alphabet = str(d.keys())
        matrix = parasail.matrix_create(alphabet, 1, -1)
        for i, x in enumerate(alphabet):
            for j, y in enumerate(alphabet):
                value = d.get(x, {}).get(y)
                if value is not None:
                    matrix.set_value(i, j, value)
        return cls(matrix)


# Local alignment functions

def local_pairwise_align_nucleotide(
        seq1, seq2, gap_open_penalty=5,
        gap_extend_penalty=2,
        match_score=2, mismatch_score=-3,
        substitution_matrix=None):

    # TODO: allow specifying subst. matrix as dict

    _check_seq_types(seq1, seq2)

    if substitution_matrix is None:
        substitution_matrix = SubstitutionMatrix.from_match_mismatch(
            match_score, mismatch_score
        )

    return local_pairwise_align(
        seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix
    )

def local_pairwise_align_protein(seq1, seq2, gap_open_penalty=11,
                                 gap_extend_penalty=1,
                                 substitution_matrix=None):

    _check_seq_types(seq1, seq2, types=(Protein,))
    _check_protein_seq_types(seq1, seq2)

    if substitution_matrix is None:
        substitution_matrix = SubstitutionMatrix.from_name("blosum50")

    return local_pairwise_align(seq1, seq2, gap_open_penalty,
                                gap_extend_penalty, substitution_matrix)


def local_pairwise_align(seq1, seq2, gap_open_penalty,
                         gap_extend_penalty, substitution_matrix):

    if isinstance(substitution_matrix, dict):
        substitution_matrix = SubstitutionMatrix.from_dict(substitution_matrix)
    
    seq1_str = str(seq1)
    seq2_str = str(seq2)

    matrix = substitution_matrix._matrix
    result = parasail.sw_trace(
        seq1_str, seq2_str, gap_open_penalty, gap_extend_penalty, matrix
    )

    cigar = result.cigar
    aligned1, aligned2 = _expand_aligned(cigar, seq1_str, seq2_str)
    msa = TabularMSA([_wrap_aligned(seq1, aligned1),
                      _wrap_aligned(seq2, aligned2)])

    score = result.score
    start_end_positions = [(cigar.beg_query, result.end_query),
                           (cigar.beg_ref, result.end_ref)]

    return msa, score, start_end_positions


# Global alignment functions

def global_pairwise_align_nucleotide(
        seq1, seq2, gap_open_penalty=5,
        gap_extend_penalty=2,
        match_score=2, mismatch_score=-3,
        substitution_matrix=None):

    _check_seq_types(seq1, seq2, types=(DNA, RNA, TabularMSA))
    _check_nucleotide_seq_types(seq1, seq2, types=(DNA, RNA))

    if substitution_matrix is None:
        substitution_matrix = SubstitutionMatrix.from_match_mismatch(
            match_score, mismatch_score
        )

    return global_pairwise_align(
        seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix
    )

def global_pairwise_align_protein(seq1, seq2, gap_open_penalty=11,
                                 gap_extend_penalty=1,
                                 substitution_matrix=None):

    _check_seq_types(seq1, seq2, types=(Protein, TabularMSA))
    _check_protein_seq_types(seq1, seq2)

    if substitution_matrix is None:
        substitution_matrix = SubstitutionMatrix.from_name("blosum50")

    return global_pairwise_align(seq1, seq2, gap_open_penalty,
                                 gap_extend_penalty, substitution_matrix)


def global_pairwise_align(seq1, seq2, gap_open_penalty,
                          gap_extend_penalty, substitution_matrix):

    if isinstance(substitution_matrix, dict):
        substitution_matrix = SubstitutionMatrix.from_dict(substitution_matrix)

    seq1_str = str(seq1)
    seq2_str = str(seq2)

    matrix = substitution_matrix._matrix
    result = parasail.nw_trace(
        seq1_str, seq2_str, gap_open_penalty, gap_extend_penalty, matrix
    )

    cigar = result.cigar
    aligned1, aligned2 = _expand_aligned(cigar, seq1_str, seq2_str)
    msa = TabularMSA([_wrap_aligned(seq1, aligned1),
                      _wrap_aligned(seq2, aligned2)])

    score = result.score
    start_end_positions = [(cigar.beg_query, result.end_query),
                           (cigar.beg_ref, result.end_ref)]

    return msa, score, start_end_positions


# Semiglobal alignment functions

def semiglobal_pairwise_align_nucleotide(
        seq1, seq2, gap_open_penalty=5,
        gap_extend_penalty=2,
        match_score=2, mismatch_score=-3,
        substitution_matrix=None):

    # TODO: allow specifying subst. matrix as dict

    _check_seq_types(seq1, seq2, types=(DNA, RNA, TabularMSA))
    _check_nucleotide_seq_types(seq1, seq2)

    if substitution_matrix is None:
        substitution_matrix = SubstitutionMatrix.from_match_mismatch(
            match_score, mismatch_score
        )

    return semiglobal_pairwise_align(
        seq1, seq2, gap_open_penalty, gap_extend_penalty, substitution_matrix
    )

def semiglobal_pairwise_align_protein(seq1, seq2, gap_open_penalty=11,
                                      gap_extend_penalty=1,
                                      substitution_matrix=None):

    _check_seq_types(seq1, seq2, types=(Protein, TabularMSA))
    _check_protein_seq_types(seq1, seq2)

    if substitution_matrix is None:
        substitution_matrix = SubstitutionMatrix.from_name("blosum50")

    return semiglobal_pairwise_align(seq1, seq2, gap_open_penalty,
                                     gap_extend_penalty, substitution_matrix)


def semiglobal_pairwise_align(seq1, seq2, gap_open_penalty,
                              gap_extend_penalty, substitution_matrix):

    if isinstance(substitution_matrix, dict):
        substitution_matrix = SubstitutionMatrix.from_dict(substitution_matrix)

    seq1_str = str(seq1)
    seq2_str = str(seq2)

    matrix = substitution_matrix._matrix
    result = parasail.sg_trace(
        seq1_str, seq2_str, gap_open_penalty, gap_extend_penalty, matrix
    )

    cigar = result.cigar
    aligned1, aligned2 = _expand_aligned(cigar, seq1_str, seq2_str)
    msa = TabularMSA([_wrap_aligned(seq1, aligned1),
                      _wrap_aligned(seq2, aligned2)])

    score = result.score
    start_end_positions = [(cigar.beg_query, result.end_query),
                           (cigar.beg_ref, result.end_ref)]

    return msa, score, start_end_positions


# Internal helpers

def _expand_aligned(cigar, seq1, seq2):
    """ Expand a parasail cigar sequence into two aligned sequences.
    """
    aligned1 = []
    aligned2 = []
    pos1 = cigar.beg_query
    pos2 = cigar.beg_ref
    for s in cigar.seq:
        op = parasail.Cigar.decode_op(s)
        ln = parasail.Cigar.decode_len(s)
        for j in range(0, ln):
            if op == b'=' or op == b'X':
                c1 = seq1[pos1]
                c2 = seq2[pos2]
                pos1 += 1
                pos2 += 1
            elif op == b'I':
                c1 = seq1[pos1]
                c2 = '-'
                pos1 += 1
            elif op == b'D':
                c1 = '-'
                c2 = seq2[pos2]
                pos2 += 1
            else:
                msg = "Invalid character in cigar string: {!r}".format(op)
                raise ValueError(msg)

            aligned1.append(c1)
            aligned2.append(c2)

    return "".join(aligned1), "".join(aligned2)


def _wrap_aligned(original, aligned):
    """ Wrap aligned string so that it has the same type as the original
    sequence.
    """
    constructor = type(original)
    metadata = None
    if original.has_metadata():
        metadata = original.metadata
    aligned = constructor(aligned, metadata=metadata, validate=False)
    return aligned


def _check_seq_types(*seqs, types=(DNA, RNA)):
    """ Check type of sequences to align.

    Raises
    ------
    TypeError

    """

    if len(seqs) == 0:
        return

    seq_types = set(type(seq) for seq in seqs)
    if len(seq_types) > 1:
        msg = "sequences must be the same type, but got {}"
        raise TypeError(msg.format(
            ", ".join(typ.__name__ for typ in seq_types)
        ))

    seq_type = next(iter(seq_types))
    if not issubclass(seq_type, types):
        msg = "sequences must be one of the following: {}, but got type {!r}"
        raise TypeError(
            msg.format(
                ", ".join(typ.__name__ for typ in types),
                seq_type.__name__
            )
        )


def _check_protein_seq_types(*seqs):
    if len(seqs) == 0:
        return

    for seq in seqs:
        if isinstance(seq, TabularMSA) and not issubclass(seq.dtype, Protein):
            raise TypeError(
                "`seq1` and `seq2` must be TabularMSA with Protein dtype, "
                "not dtype %r" % seq.dtype.__name__
            )        


def _check_nucleotide_seq_types(*seqs, types=(DNA, RNA)):
    if len(seqs) == 0:
        return

    for seq in seqs:
        if isinstance(seq, TabularMSA) and not issubclass(seq.dtype, types):
            raise TypeError(
                "`seq1` and `seq2` must be TabularMSA with DNA or RNA dtype, "
                "not dtype %r" % seq.dtype.__name__
            )        
