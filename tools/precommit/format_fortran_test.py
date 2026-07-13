#!/usr/bin/env python3

import unittest
import os
import tempfile
import shutil
import sys

from format_fortran import main
from prettify_cp2k import normalizeFortranFile
from prettify_cp2k import selftest


class TestSingleFileFolder(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.mkdtemp()
        self.fname = os.path.join(self.tempdir, "prettify_selftest.F")

        # create temporary file with example code
        with open(self.fname, "w", encoding="utf8") as fhandle:
            fhandle.write(selftest.content)

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_prettify(self):
        # call prettify, the return value should be 0 (OK)
        self.assertEqual(main([self.fname]), 0)

        # check if file was altered (it shouldn't)
        with open(self.fname, encoding="utf8") as fhandle:
            result = fhandle.read()

        self.assertEqual(result.splitlines(), selftest.content.splitlines())

    def test_large_declaration_dependency_pass(self):
        declarations = []
        for idx in range(500):
            declarations.append(
                {
                    "attributes": [],
                    "parameters": None,
                    "vars": [f"var_{idx:03d}"],
                    "normalizedType": f"TYPE_{idx:03d}",
                }
            )

        normalizeFortranFile.enforceDeclDependecies(declarations)
        self.assertEqual(sum(len(d["vars"]) for d in declarations), 500)


if __name__ == "__main__":
    unittest.main()
