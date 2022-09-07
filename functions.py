import numpy as np
import os
import sys

import random
from datetime import datetime
random.seed(datetime.now())

# Scalar product
# input: - vector 1, vector 2
# output: - scalar product


def scalar(v1, v2):

    scalar_product = np.dot(v1, v2)

    return scalar_product

# Cross product
# input: - vector 1, vector 2
# output: - vector cross product


def cross(v1, v2):

    cross_product = []
    cross_product.append(v1[1] * v2[2] - v1[2] * v2[1])
    cross_product.append(v1[2] * v2[0] - v1[0] * v2[2])
    cross_product.append(v1[0] * v2[1] - v1[1] * v2[0])
    cross_product = np.array(cross_product)

    return cross_product

# Rotation assuming isotropy on a sphere
### input: - vector
# output: - vector rotated


def rotation_iso(vec):

    phi = 2.0 * np.pi * random.random()
    theta = np.arccos(1.0 - 2.0 * random.random())

    result = []
    result.append(vec[0] * np.sin(theta) * np.cos(phi))
    result.append(vec[0] * np.sin(theta) * np.sin(phi))
    result.append(vec[0] * np.cos(theta))
    result = np.array(result)

    return result

# Rotation of a vector
### input: - vector
# output: - vector rotated


def rotation(vec):

    phi = 2.0 * np.pi * random.random()
    theta = np.pi * random.random()

    result = []
    result.append(vec[0] * np.sin(theta) * np.cos(phi))
    result.append(vec[0] * np.sin(theta) * np.sin(phi))
    result.append(vec[0] * np.cos(theta))
    result = np.array(result)

    return result

# Sample assuming power law
# input: - min, max, slope
### output: - sampled value


def sample_powerlaw(mmin, mmax, slope):

    gammaa = 1.0 + slope

    if slope == -1:
        mb = mmin * np.exp(random.random() * np.log(mmax / mmin))
    else:
        mb = (random.random() * (mmax ** gammaa - mmin ** gammaa) + mmin ** gammaa) ** (1.0 / gammaa)

    return mb

def red_formation(z, z_nsc, sigma_nsc):

    return np.exp(-(z - z_nsc) ** 2 / (2.0 * sigma_nsc ** 2))



