import os
import unittest
import uuid

import medaka.common

class TestMkdir(unittest.TestCase):

    def test_000_directory_created(self):
        dirname = str(uuid.uuid4())
        medaka.common.mkdir_p(dirname)
        self.assertTrue(os.path.exists(dirname))
        os.rmdir(dirname)

    def test_001_directory_not_overwritten(self):
        dirname = str(uuid.uuid4())
        medaka.common.mkdir_p(dirname)
        t0 = os.path.getmtime(dirname)
        medaka.common.mkdir_p(dirname)
        t1 = os.path.getmtime(dirname)
        # time stamp of directory unchanged after calling mkdir_p twice
        self.assertEqual(t0, t1)
        os.rmdir(dirname)

class TestBase2IndexDict(unittest.TestCase):
    def testb2iu(self):
        b2i_dict = medaka.common.base2index
        for idx,character in enumerate('acgtACGTdD'):
            assert b2i_dict[character] == idx
