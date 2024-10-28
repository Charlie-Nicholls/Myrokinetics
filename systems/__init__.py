from .archer2 import archer2
from .viking import viking
from .ypi_server import ypi_server

systems = {'viking': viking, 'archer2': archer2, 'ypi_server': ypi_server}

__all__ = ['systems']
