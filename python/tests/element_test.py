import os, sys

#tracy_dir = os.getenv('TRACY_LIB')
#sys.path.append(tracy_dir+'/tracy/lib')

import unittest
import thor_scsi.lib as tslib
from thor_scsi import pyflame
import pytest
import abc


class _ElementTestBasis(unittest.TestCase, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def setUp(self):
        self.elem = None
        raise ValueError('Use derived class instead')

    def test00_name(self):
        '''test if element name can be read
        '''
        elem = self.elem
        name = elem.name

    @pytest.mark.skip
    def test01_reverse(self):
        '''Can element be reverted
        '''
        self.assertEqual(self.elem.reverse, False)
        self.elem.Reverse = True
        self.assertEqual(self.elem.Reverse, True)

    def test02_repr(self):
        repr(self.elem)

    def test10_print(self):
        '''test if element can be printed
        '''
        print(self.elem)
        print(repr(self.elem))


class DriftTest(_ElementTestBasis):
    def setUp(self):
        config = pyflame.Config()
        config.setAny("name", "a_drift")
        config.setAny("length", 0.5)
        self.elem = tslib.Drift(config)


del _ElementTestBasis


if __name__ == '__main__':
    unittest.main()
