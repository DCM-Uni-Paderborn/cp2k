#!/usr/bin/env python3
"""Test Wannier90 overlap and projection-export diagnostics.

These tests cover file parsing, overlap reciprocity and projection subspaces.
The Fortran wannier90_projections_unittest instead exercises AO trial selection,
rank checks, and projection covariance.
"""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import numpy as np

from compare_wannier90_mmn import check_reciprocity, compare_mmn, read_mmn
from validate_wannier90_export import compare_amn, main, read_amn, read_eig


class ReciprocityTests(unittest.TestCase):
    def setUp(self):
        self.forward = (1, 2, 1, 0, -1)
        self.reverse = (2, 1, -1, 0, 1)
        self.matrix = np.array([[0.9, 0.2j], [0.1 - 0.2j, 0.8j]])

    def data(self, reverse=None):
        if reverse is None:
            reverse = self.matrix.conj().T
        return (2, 2, 1, [(self.forward, self.matrix), (self.reverse, reverse)])

    def test_nonhermitian_block_and_reverse(self):
        self.assertEqual(check_reciprocity(self.data()), 0.0)

    def test_complex_unitary_gauge(self):
        left = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)
        right = np.diag([1j, -1])
        forward = left.conj().T @ self.matrix @ right
        reverse = right.conj().T @ self.matrix.conj().T @ left
        data = (2, 2, 1, [(self.forward, forward), (self.reverse, reverse)])
        self.assertLess(check_reciprocity(data), 1e-14)
        raw, singular, _, _ = compare_mmn(self.data(), data)
        self.assertGreater(raw, 0.1)
        self.assertLess(singular, 1e-14)

    def test_wrong_sign(self):
        self.assertGreater(check_reciprocity(self.data(-self.matrix.conj().T)), 1.0)

    def test_missing_reverse(self):
        data = (2, 2, 1, [(self.forward, self.matrix), ((2, 1, 0, 0, 0), self.matrix)])
        with self.assertRaisesRegex(ValueError, "Missing reverse"):
            check_reciprocity(data)

    def test_duplicate_header(self):
        data = (2, 2, 1, [(self.forward, self.matrix), (self.forward, self.matrix)])
        with self.assertRaisesRegex(ValueError, "Repeated or missing"):
            check_reciprocity(data)

    def test_nonfinite_overlap(self):
        reverse = self.matrix.conj().T.copy()
        reverse[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "Invalid overlap"):
            check_reciprocity(self.data(reverse))
        with self.assertRaisesRegex(ValueError, "Invalid overlap"):
            compare_mmn(self.data(), self.data(reverse))


class ExportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.matrix = np.array([[0.9, 0.2j], [0.1 - 0.2j, 0.8j]])
        self.mmn = "test export\n2 1 1\n1 1 1 0 0\n" + "".join(
            f"{z.real:.16e} {z.imag:.16e}\n" for z in self.matrix.ravel(order="F")
        )
        self.amn = "test export\n2 1 1\n1 1 1 1.0 0.0\n2 1 1 0.0 0.0\n"
        self.eig = "2 1 2.0\n1 1 1.0\n"
        for seed in ("left", "right"):
            for suffix, text in (
                ("mmn", self.mmn),
                ("amn", self.amn),
                ("eig", self.eig),
            ):
                (self.root / f"{seed}.{suffix}").write_text(text, encoding="utf-8")

    def test_mmn_file_validation(self):
        path = self.root / "left.mmn"
        np.testing.assert_array_equal(read_mmn(str(path))[3][0][1], self.matrix)
        for text in (
            self.mmn.replace("2 1 1", "-2 1 1", 1),
            self.mmn.rsplit("\n", 2)[0],
            self.mmn + "1 1 1 0 0\n",
            self.mmn.replace("9.0000000000000002e-01", "nan"),
            self.mmn.replace("1 1 1 0 0", "2 1 1 0 0"),
        ):
            with self.subTest(text=text):
                path.write_text(text, encoding="utf-8")
                with self.assertRaises(ValueError):
                    read_mmn(str(path))

    def test_indexed_eig_and_amn_validation(self):
        self.assertEqual(read_eig(self.root / "left.eig", 2, 1), [1.0, 2.0])
        np.testing.assert_array_equal(
            read_amn(self.root / "left.amn", 2, 1)[0], [[1.0], [0.0]]
        )
        for suffix, malformed, reader in (
            ("eig", "1 1 1\n1 1 2\n", read_eig),
            ("eig", "1 1 nan\n2 1 2\n", read_eig),
            ("eig", "1 1 1\n", read_eig),
            ("amn", self.amn.replace("2 1 1 0.0", "1 1 1 0.0"), read_amn),
            ("amn", self.amn.replace("1.0", "nan"), read_amn),
            ("amn", self.amn.rsplit("\n", 2)[0], read_amn),
        ):
            with self.subTest(suffix=suffix, malformed=malformed):
                path = self.root / f"bad.{suffix}"
                path.write_text(malformed, encoding="utf-8")
                with self.assertRaises(ValueError):
                    reader(path, 2, 1)

    def test_rectangular_projection_subspaces(self):
        a = np.eye(3, 2, dtype=complex)
        gauge = np.array([[1, 1j], [1j, 1]]) / np.sqrt(2)
        metrics = compare_amn([a], [a @ gauge])
        self.assertGreater(metrics["amn_max_abs_delta"], 0.1)
        self.assertLess(metrics["amn_max_singular_delta"], 1e-14)
        self.assertLess(metrics["amn_subspace_deviation"], 1e-14)
        changed = a.copy()
        changed[:, 0] = [1 / np.sqrt(2), 0, 1 / np.sqrt(2)]
        metrics = compare_amn([a], [changed])
        self.assertLess(metrics["amn_max_singular_delta"], 1e-14)
        self.assertGreater(metrics["amn_subspace_deviation"], 0.7)
        with self.assertRaisesRegex(ValueError, "rank-deficient"):
            compare_amn([a], [np.zeros_like(a)])

    def test_export_cli_thresholds(self):
        arguments = [
            "compare-export",
            "left",
            "right",
            "--left-workdir",
            str(self.root),
            "--right-workdir",
            str(self.root),
            "--compare-amn",
            "--eig-tol",
            "1e-12",
            "--mmn-tol",
            "1e-12",
            "--mmn-singular-tol",
            "1e-12",
            "--amn-tol",
            "1e-12",
            "--amn-singular-tol",
            "1e-12",
            "--amn-subspace-tol",
            "1e-12",
        ]
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(arguments), 0)
            (self.root / "right.eig").write_text("1 1 1.1\n2 1 2.0\n", encoding="utf-8")
            self.assertEqual(main(arguments), 1)


if __name__ == "__main__":
    unittest.main()
