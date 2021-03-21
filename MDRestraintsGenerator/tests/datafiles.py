"""
Location of datafiles for tests
"""

__all__ = [
        "T4_TPR", "T4_XTC", "T4_OGRO", "T4_OTOP",
        "T4_FB_NDX", "T4_FB_MDP", "T4_FB_OGRO", "T4_FB_MDP_DEBUG",
        "T4_H_MDP", "T4_H_MDP_DEBUG",
        "CB8_TPR", "CB8_XTC", "CB8_OGRO", "CB8_OTOP",
]

from pkg_resources import resource_filename

T4_TPR = resource_filename(__name__, '../data/4RBN/npt_prod1.tpr')
T4_XTC = resource_filename(__name__,
                           '../data/4RBN/npt_prod1_center.skip200.xtc')
T4_OGRO = resource_filename(__name__, '../data/4RBN/ClosestRestraintFrame.gro')
T4_OTOP = resource_filename(__name__, '../data/4RBN/BoreschRestraint.top')
T4_FB_NDX = resource_filename(__name__, '../data/4RBN/flatbottom_index.ndx')
T4_FB_MDP = resource_filename(__name__, '../data/4RBN/flatbottom.mdp')
T4_FB_OGRO = resource_filename(
        __name__, '../data/4RBN/flatbottom_ClosestRestraintFrame.gro')
T4_FB_MDP_DEBUG = resource_filename(__name__,
                                    '../data/4RBN/flatbottom_debug.mdp')
T4_H_MDP = resource_filename(__name__, '../data/4RBN/harmonic.mdp')
T4_H_MDP_DEBUG = resource_filename(__name__, '../data/4RBN/harmonic_debug.mdp')
CB8_TPR = resource_filename(__name__, '../data/CB8-G3/npt_prod1.tpr')
CB8_XTC = resource_filename(__name__,
                            '../data/CB8-G3/npt_prod_center.skip200.xtc')
CB8_OGRO = resource_filename(__name__,
                             '../data/CB8-G3/ClosestRestraintFrame.gro')
CB8_OTOP = resource_filename(__name__,
                             '../data/CB8-G3/BoreschRestraint.top')
