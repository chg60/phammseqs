"""Test the code used for sub-processing."""

import unittest

from phamerate.subprocess import run


class MyTestCase(unittest.TestCase):
    def test_valid_command_ok(self):
        """Test a command that should work on any macOS or Linux system."""
        command = f"uname -a"
        stdout, stderr = run(command)
        with self.subTest():
            self.assertEqual(stderr, "")
        with self.subTest():
            self.assertGreater(stdout, "")

    def test_invalid_command_bad(self):
        """Test a command that should not work on any macOS or Linux system."""
        command = f"fake_command --fake-argument"
        with self.assertRaises(FileNotFoundError):
            run(command)


if __name__ == '__main__':
    unittest.main()
