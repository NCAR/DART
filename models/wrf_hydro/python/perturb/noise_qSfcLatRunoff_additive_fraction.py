import numpy as np

def noise_model(x, seed: int=None):

    #######################################################
    # Noise model for qSfcLatRunoff
    # 0) Additive noise,
    # 1) Zero-mean,
    # 2) frac: Standard deviation is a fraction of the value,
    # 3) min: Truncated below at min,
    frac = 0.4
    minim = 0.0

    # Do not want variables to be correlated...
    # Adjust the seed per variable using unicode numbers summed for the
    # name of the variable. Make sure it's a 32-bit integer.
    if seed is not None:
        seed = (seed - sum([ord(ii) -96 for ii in list('qsfclatrunoff')])) % (2**32 - 1)
        np.random.seed(seed)
        
    std_dev_vec = np.multiply(frac, x)
    perturb_vec = np.add(x, np.random.normal(0.0, std_dev_vec))
    constrain_perturb_vec = np.fmax(perturb_vec, minim)

    # Dummy for testing that perturbations are zero-mean and uncorrelated.
    # std_dev_vec = np.add(1, np.multiply(0, x))
    # perturb_vec = np.random.normal(0.0, std_dev_vec)
    # constrain_perturb_vec = perturb_vec

    return constrain_perturb_vec
