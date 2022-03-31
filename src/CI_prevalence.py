#shebang

from scipy.stats import chi2, beta, randint

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def _beta_confidence(proportion, total_sample, alpha=0.05, integer=True):
    """
    Better confidence intervals for proportions, based on a classic
    contribution by:
    .. [1] Steven A Julious, "Two-sided confidence intervals for the single
    proportion: comparison of seven methods by Robert G. Newcombe, Statistics
    in Medicine 1998; 17:857-872"
    
    Parameters
    ----------
    proportion : float
        A proportion between 0 and 1.
    total_sample : int
        The total sample size the proportion was derived from.
    alpha : float, default 0.05
        A float between 0 and 1, representing the type 1 error rate. Is used
        to define the confidence interval coverage: (1-alpha)*100.
    integer : boolean, default True
        If the function should fail if the proportion times total_sample does
        not result in an integer number of events. Set to `False` to ignore
        the ValueError.
    
    Returns
    -------
    Unpacks the lower and upper bounds
    """
    
    # check input
    # is_type(proportion, float)
    # is_type(total_sample, int)
    # get the number of events
    no_events = proportion * total_sample
    if (not no_events.is_integer()) and (integer==True):
        raise ValueError('proportion * total_sample is not an integer, either \
correct input or set `integer` to `False`')
        
    # lower bound
    lb = 1-beta.ppf(1-alpha/2, total_sample - no_events + 1, no_events)
    # upper bound
    ub = beta.ppf(1-alpha/2, no_events + 1, total_sample - no_events)
    # return
    return lb, ub


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

acm = _beta_confidence(282/200643, 200643, integer = False)
print("ACM prevalence:\t", str(round(200643/282, 0)), " [95% CI", str(round(1/acm[1], 0)), str(round(1/acm[0], 0)) + "]")

dcm = _beta_confidence(694/200643, 200643, integer = False)
print("DCM prevalence:\t", str(round(200643/694, 0)), " [95% CI", str(round(1/dcm[1], 0)), str(round(1/dcm[0], 0)) + "]")

dhm = _beta_confidence(566/200643, 200643, integer = False)
print("DCM prevalence:\t", str(round(200643/566, 0)), " [95% CI", str(round(1/dhm[1], 0)), str(round(1/dhm[0], 0)) + "]")

hcm = _beta_confidence(772/200643, 200643, integer = False)
print("HCM prevalence:\t", str(round(200643/772, 0)), " [95% CI", str(round(1/hcm[1], 0)), str(round(1/hcm[0], 0)) + "]")
