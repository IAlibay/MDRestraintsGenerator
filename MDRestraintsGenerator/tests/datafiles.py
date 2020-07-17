"""
Location of datafiles for tests
"""

__all__ = [
        "T4_TPR", "T4_XTC", "T4_NC",
        "T4_OGRO", "T4_OTOP",
]

from pkg_resources import resource_filename

T4_TPR = resource_filename(__name__, '../data/4RBN/npt_prod1.tpr')
T4_XTC = resource_filename(__name__,
                           '../data/4RBN/npt_prod1_center.skip200.xtc')
T4_NC = resource_filename(__name__,
                          '../data/4RN/npt_prod1_center.skip200.nc')
T4_OGRO = resource_filename(__name__,
                            '../data/4RBN/ClosestRestraintFrame.gro')
T4_OTOP = resource_filename(__name__,
                            '../data/4RBN/BoreschRestraint.top')
