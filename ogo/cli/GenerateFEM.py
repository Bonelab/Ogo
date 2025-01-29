import argparse
import sys

from ogo.cli.ref.SidewaysFallFe import main as sideways_fall_main # Adjust import based on file structure
from ogo.cli.ref.SpineCompressionFe import main as spine_compression_main # Adjust import based on file structure

def main():
    """
    Wrapper function to run the vertebra or femur FE model setup.
    'One day we will have a better way to do this.' - Matthias Walle
    
    - vertebra (OgoSpineCompressionFe):
        Required Arguments:
            - calibrated_image: Path to the *_K2HPO4.nii image file.
            - bone_mask: Path to the *_MASK.nii bone mask image.
            - --vertebral_body_label: Label for the vertebral body in the mask.
            - --vertebral_process_label: Label for the vertebral process in the mask.
        Optional Arguments:
            - Refer to OgoSpineCompressionFe documentation for additional parameters.

    - femur (ogoSidewaysFallFe):
        Required Arguments:
            - calibrated_image: Path to the *_K2HPO4.nii image file.
            - bone_mask: Path to the *_MASK.nii bone mask image.
            - --femur_side: 1 for left femur, 2 for right femur.
        Optional Arguments:
            - Refer to ogoSidewaysFallFe documentation for additional parameters.
    """
    # Wrapper argument parsing
    parser = argparse.ArgumentParser(
        description=(
            "Wrapper to select vertebra or femur model setup.\n\n"
            "Model Types:\n"
            "  vertebra: Run the OgoSpineCompressionFe script to set up an L4 vertebral compression FE model.\n"
            "  femur: Run the ogoSidewaysFallFe script to set up a sideways fall FE model for the femur."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--model_type",
        choices=["vertebra", "femur"],
        required=True,
        help="Specify the model type to run: 'vertebra' or 'femur'.",
    )
    parser.add_argument(
        "remaining_args",
        nargs=argparse.REMAINDER,
        help=(
            "Additional arguments for the selected model.\n"
            "For vertebra, refer to OgoSpineCompressionFe.\n"
            "For femur, refer to ogoSidewaysFallFe."
        ),
    )
    args = parser.parse_args()

    # Forward to the selected script
    if args.model_type == "vertebra":
        # Simulate calling OgoSpineCompressionFe directly
        sys.argv = ["OgoSpineCompressionFe"] + args.remaining_args
        spine_compression_main()
    elif args.model_type == "femur":
        # Simulate calling ogoSidewaysFallFe directly
        sys.argv = ["OgoSidewaysFallFe"] + args.remaining_args
        sideways_fall_main()


if __name__ == "__main__":
    GenerateFE()