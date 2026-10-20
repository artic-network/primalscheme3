import pathlib
import tempfile
import unittest

from primalscheme3.core.logger import setup_rich_logger


def flush(logger):
    for handler in logger.handlers:
        handler.flush()


class TestSetupRichLogger(unittest.TestCase):
    def test_rerun_with_same_logfile_replaces_handlers(self):
        with tempfile.TemporaryDirectory() as d:
            logfile = pathlib.Path(d) / "file.log"
            first = setup_rich_logger(str(logfile))
            first.debug("first run")
            flush(first)

            second = setup_rich_logger(str(logfile))
            self.assertIs(second, first)
            self.assertEqual(len(second.handlers), 2)
            second.debug("second run")
            flush(second)
            text = logfile.read_text()
            for handler in second.handlers:
                handler.close()
        self.assertEqual(text.count("second run"), 1)
        self.assertNotIn("first run", text)

    def test_console_only_without_logfile(self):
        logger = setup_rich_logger(None)
        self.assertEqual(len(logger.handlers), 1)


if __name__ == "__main__":
    unittest.main()
