import pathlib
import tempfile
import unittest
from unittest import mock

from typer.testing import CliRunner

from primalscheme3.cli import app

MSA = "tests/core/test_mismatch.fasta"

DOWNSAMPLE_ARGS = {
    "--downsample-target": 0.95,
    "--downsample-max-dct": -1.5,
    "--downsample-always-add-prop": 0.3,
    "--primer-scorer-model": "q5_dct_v1",
}


class TestDownsampleOptions(unittest.TestCase):
    """scheme-create and panel-create expose the downsampling settings."""

    def test_shown_in_help(self):
        for command in ("scheme-create", "panel-create"):
            result = CliRunner().invoke(
                app, [command, "--help"], env={"COLUMNS": "200"}
            )
            self.assertEqual(result.exit_code, 0, result.output)
            for option in ("--downsample", *DOWNSAMPLE_ARGS):
                with self.subTest(command=command, option=option):
                    self.assertIn(option, result.output)

    def test_passed_to_config(self):
        args = ["--downsample"]
        for option, value in DOWNSAMPLE_ARGS.items():
            args += [option, str(value)]

        for command, target in (
            ("scheme-create", "schemecreate"),
            ("panel-create", "panelcreate"),
        ):
            with (
                self.subTest(command=command),
                tempfile.TemporaryDirectory() as d,
                mock.patch(f"primalscheme3.cli.{target}") as run,
            ):
                out = str(pathlib.Path(d) / "out")
                result = CliRunner().invoke(
                    app, [command, "--msa", MSA, "--output", out, *args]
                )
                self.assertEqual(result.exit_code, 0, result.output)
                config = run.call_args.kwargs["config"]
                self.assertTrue(config.downsample)
                self.assertEqual(config.downsample_target, 0.95)
                self.assertEqual(config.downsample_max_dct, -1.5)
                self.assertEqual(config.downsample_always_add_prop, 0.3)
                self.assertEqual(config.primer_scorer_model, "q5_dct_v1")

    def test_out_of_range_rejected(self):
        for option, value in (
            ("--downsample-max-dct", "1"),
            ("--downsample-target", "1.5"),
            ("--downsample-always-add-prop", "-0.1"),
        ):
            with self.subTest(option=option):
                result = CliRunner().invoke(
                    app,
                    ["scheme-create", "--msa", MSA, "--output", "x", option, value],
                )
                self.assertNotEqual(result.exit_code, 0)


if __name__ == "__main__":
    unittest.main()
