"""
Specific warning categories to more easily test whether specific warnings
have been emitted.
"""


class NonAutosomalWarning(UserWarning):
    pass


class QCMissingWarning(UserWarning):
    pass


class SLiMScalingFactorWarning(UserWarning):
    pass


class SLiMOddSampleWarning(UserWarning):
    pass


class UnspecifiedSLiMWarning(UserWarning):
    pass
