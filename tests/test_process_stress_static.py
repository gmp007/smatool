import ast
import pathlib
import unittest


PROCESS_STRESS = pathlib.Path(__file__).resolve().parents[1] / "src" / "process_stress.py"


class ProcessStressStaticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = PROCESS_STRESS.read_text(encoding="utf-8")
        cls.tree = ast.parse(cls.source)

    def test_uses_scipy_trapezoid_not_removed_trapz_import(self):
        scipy_integrate_imports = [
            alias.name
            for node in ast.walk(self.tree)
            if isinstance(node, ast.ImportFrom) and node.module == "scipy.integrate"
            for alias in node.names
        ]

        self.assertIn("trapezoid", scipy_integrate_imports)
        self.assertNotIn("trapz", scipy_integrate_imports)

    def test_bulk_modulus_fits_validate_data_points(self):
        self.assertIn("def _validate_curve_fit_data", self.source)
        self.assertIn("Birch-Murnaghan energy-volume fit", self.source)
        self.assertIn("Birch-Murnaghan pressure-volume fit", self.source)


if __name__ == "__main__":
    unittest.main()
