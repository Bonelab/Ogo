import unittest
from unittest import mock
from ogo.cli import GenerateFEM


class TestGenerateFEM(unittest.TestCase):
    def setUp(self) -> None:
        pass

    def tearDown(self) -> None:
        pass

    @mock.patch("ogo.cli.GenerateFEM.spine_compression_main")
    def test_dispatch_to_vertebra(self, mock_spine):
        """Test that 'vertebra' dispatches to spine_compression_main"""
        test_args = [
            "program",
            "--model_type", "vertebra",
            "--calibrated_image", "dummy_image.nii.gz",
            "--bone_mask", "dummy_mask.nii.gz",
            "--vertebral_body_label", "1",
            "--vertebral_process_label", "2"
        ]
        with mock.patch("sys.argv", test_args):
            GenerateFEM.main()
            mock_spine.assert_called_once()

    @mock.patch("ogo.cli.GenerateFEM.sideways_fall_main")
    def test_dispatch_to_femur(self, mock_femur):
        """Test that 'femur' dispatches to sideways_fall_main"""
        test_args = [
            "program",
            "--model_type", "femur",
            "--calibrated_image", "dummy_image.nii.gz",
            "--bone_mask", "dummy_mask.nii.gz",
            "--femur_side", "1"
        ]
        with mock.patch("sys.argv", test_args):
            GenerateFEM.main()
            mock_femur.assert_called_once()


if __name__ == "__main__":
    unittest.main()
