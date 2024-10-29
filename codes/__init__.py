from .gs2 import gs2
from .cgyro import cgyro
from .tglf import tglf



codes = {'GS2': gs2(), 'CGYRO': cgyro(), 'TGLF': tglf()}

__all__ = ['codes']
