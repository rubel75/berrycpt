import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class TestReconstruction(unittest.TestCase):

    @staticmethod
    def normalize_log(text):
        lines = text.splitlines()

        for index, line in enumerate(lines):
            if line.startswith("Wrote reconstructed full-BZ data to "):
                lines[index] = (
                    "Wrote reconstructed full-BZ data to <OUTPUT_FILE>"
                )

        return "\n".join(lines).strip()

    def run_reference_test(
        self,
        input_name,
        expected_output_name,
        expected_log_name,
    ):
        utility_dir = Path(__file__).resolve().parent.parent
        input_dir = utility_dir / "tests" / "input"
        expected_dir = utility_dir / "tests" / "expected_output"

        expected_output = expected_dir / expected_output_name
        expected_log = expected_dir / expected_log_name

        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "reconstructed_bcurv_fbz.dat"

            result = subprocess.run(
                [
                    sys.executable,
                    str(utility_dir / "berrycpt_reconstruct.py"),
                    str(input_dir / "OUTCAR"),
                    str(input_dir / "IBZKPT"),
                    str(input_dir / input_name),
                    "--output",
                    str(output_file),
                ],
                check=True,
                cwd=utility_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
            )

            self.assertEqual(
                output_file.read_text(),
                expected_output.read_text(),
            )

            self.assertEqual(
                self.normalize_log(result.stdout),
                self.normalize_log(expected_log.read_text()),
            )

    def test_current_format(self):
        self.run_reference_test(
            "bcurv_ij.dat",
            "reconstructed_bcurv_fbz.dat",
            "reconstruction.log",
        )

    def test_legacy_format(self):
        self.run_reference_test(
            "bcurv_ij_legacy.dat",
            "reconstructed_bcurv_fbz_legacy.dat",
            "reconstruction_legacy.log",
        )


if __name__ == "__main__":
    unittest.main()