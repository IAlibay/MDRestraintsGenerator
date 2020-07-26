"""
Location of datafiles for tests
"""

__all__ = [
        "T4_TPR", "T4_XTC",
        "T4_OGRO", "T4_OTOP",
        "CB8_TPR", "CB8_XTC",
        "CB8_OGRO", "CB8_OTOP",
]

from pkg_resources import resource_filename

T4_TPR = resource_filename(__name__, '../data/4RBN/npt_prod1.tpr')
T4_XTC = resource_filename(__name__,
                           '../data/4RBN/npt_prod1_center.skip200.xtc')
T4_OGRO = resource_filename(__name__,
                            '../data/4RBN/ClosestRestraintFrame.gro')
T4_OTOP = resource_filename(__name__,
                            '../data/4RBN/BoreschRestraint.top')
CB8_TPR = resource_filename(__name__,
                            '../data/CB8-G3/npt_prod1.tpr')
CB8_XTC = resource_filename(__name__,
                            '../data/CB8-G3/npt_prod_center.skip200.xtc')
CB8_OGRO = resource_filename(__name__,
                             '../data/CB8-G3/ClosestRestraintFrame.gro')
CB8_OTOP = resource_filename(__name__,
                             '../data/CB8-G3/BoreschRestraint.top')
