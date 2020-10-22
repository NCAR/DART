import numpy as np

def noise_model():

    #######################################################
    # Noise model for qBucket
    ## 0) Additive noise,
    ## 1) Zero-mean,
    ## 2) Standard deviation is a fixed at .02
    ## 3) min: Truncated below at min,
    ## 4) size: number of samples,
    ## 5) Closure takes a single value argument.
    def close_trunc_gauss(min: float=0.0):
        def the_closure(x, seed: int=None):
            if seed is not None:
                np.random.seed(seed)
            return np.maximum(x+np.random.normal(0.0, 10.0, 1), min)
        return the_closure

    return np.vectorize(close_trunc_gauss(min=0.0))
