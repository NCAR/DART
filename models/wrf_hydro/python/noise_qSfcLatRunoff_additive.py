import numpy as np

def noise_model():

    #######################################################
    # Noise model for qSfcLatRunoff
    # 0) Additive noise,
    # 1) Zero-mean,
    # 2) frac: Standard deviation is a fraction of the value,
    # 3) min: Truncated below at min,
    # 4) size: number of samples,
    # 5) Closure takes a single value argument.

    def close_trunc_gauss_sd_pct_value(frac, min):
        def the_closure(x):
            return np.maximum(x+np.random.normal(0.0, frac*x, 1), min)
        return the_closure

    return np.vectorize(close_trunc_gauss_sd_pct_value(.2, 0))
