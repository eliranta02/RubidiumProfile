from __future__ import division,print_function
import numpy as np
import scipy.sparse as sp




def mat2vec(mat):
    """
    Private function reshaping matrix to vector.
    """
    return mat.T.reshape(np.prod(np.shape(mat)), 1)


def vec2mat(vec):
    """
    Private function reshaping vector to matrix.
    """
    n = int(np.sqrt(len(vec)))
    return vec.reshape((n, n)).T

def vec2mat_index(N, I):
    """
    Convert a vector index to a matrix index pair that is compatible with the
    vector to matrix rearrangement done by the vec2mat function.
    """
    j = int(I / N)
    i = I - N * j
    return i, j


def mat2vec_index(N, i, j):
    """
    Convert a matrix index pair to a vector index that is compatible with the
    matrix to vector rearrangement done by the mat2vec function.
    """
    return i + N * j




