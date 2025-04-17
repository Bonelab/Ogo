import unittest
from unittest import mock
from ogo.cli import GenerateFEM


class TestGenerateFEM(unittest.TestCase):
    def setUp(self) -> None:
        """ Set up before each test. """
        pass

    def tearDown(self) -> None:
        """ Clean up after each test. """
        pass

    @mock.patch("ogo.cli.GenerateFEM.spine_compression_main")
    def test_dispatch_to_vertebra(self, mock_spine):
        """ Test that 'vertebra' dispatches to spine_compression_main """
        test_args = ["program", "--model_type", "vertebra", "--arg1", "value"]
        with mock.patch("sys.argv", test_args):
            GenerateFEM.main()
            mock_spine.assert_called_once()

    @mock.patch("ogo.cli.GenerateFEM.sideways_fall_main")
    def test_dispatch_to_femur(self, mock_femur):
        """ Test that 'femur' dispatches to sideways_fall_main """
        test_args = ["program", "--model_type", "femur", "--arg1", "value"]
        with mock.patch("sys.argv", test_args):
            GenerateFEM.main()
            mock_femur.assert_called_once()


if __name__ == "__main__":
    unittest.main()
