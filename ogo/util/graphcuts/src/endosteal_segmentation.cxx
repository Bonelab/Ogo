#include <iostream>

#include "itkEndostealSegmentationImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"

constexpr unsigned int ImageDimension = 3;
using InputPixelType 	= float;
using MaskPixelType 	= unsigned char;

using InputImageType	= itk::Image< InputPixelType, ImageDimension >;
using MaskImageType		= itk::Image< MaskPixelType, ImageDimension >;

using InputReaderType 	= itk::ImageFileReader< InputImageType >;
using MaskWReaderType		= itk::ImageFileReader< MaskImageType >;
using OutputWriterType	= itk::ImageFileWriter< MaskImageType >;

using EndostealSegmentationFilterType = itk::EndostealSegmentationImageFilter< InputImageType, MaskImageType, MaskImageType >;
using BinaryThresholdFilterType = itk::BinaryThresholdImageFilter< MaskImageType, MaskImageType >;

int main(int argc, char** argv) {
  if( argc != 13 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <MaskFileName> <OutputSegmentation> ";
		std::cerr << " <Lambda> <Sigma>";
		std::cerr << " <LowerThresh> <UpperThresh>";
    std::cerr << " <CortcialLabel> <CancellousLabel> <BackgroundLabel>";
		std::cerr << " <MinDistance> <MaxDistance>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

	/* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string maskFileName = argv[2];
  std::string outputFileName = argv[3];

  double lambda = atof(argv[4]);
	double sigma = atof(argv[5]);

	double lowerThresh = atof(argv[6]);
	double upperThresh = atof(argv[7]);

  int corticalLabel = atoi(argv[8]);
  int cancellousLabel = atoi(argv[9]);
	int backgroundLabel = atoi(argv[10]);

	double minDistance = atof(argv[11]);
	double maxDistance = atof(argv[12]);

	std::cout << "Parameters:" << std::endl;
  std::cout << "  InputFilePath:    " << inputFileName << std::endl;
  std::cout << "  MaskFilePath:     " << maskFileName << std::endl;
  std::cout << "  OutputFilePath:   " << outputFileName << std::endl;
  std::cout << "  Lambda:           " << lambda << std::endl;
  std::cout << "  Sigma:            " << sigma << std::endl;
	std::cout << "  Lower Thresh:     " << lowerThresh << std::endl;
	std::cout << "  Upper Thresh:     " << upperThresh << std::endl;
  std::cout << "  CorticalLabel:    " << corticalLabel << std::endl;
  std::cout << "  CancellousLabel:  " << cancellousLabel << std::endl;
	std::cout << "  BackgroundLabel:  " << backgroundLabel << std::endl;
	std::cout << "  Min Distance:     " << minDistance << std::endl;
	std::cout << "  Max Distance:     " << maxDistance << std::endl;
  std::cout << std::endl;

	std::cout << "Reading input " << inputFileName << std::endl;
	InputReaderType::Pointer input_reader = InputReaderType::New();
	input_reader->SetFileName(inputFileName);
	input_reader->Update();

	std::cout << "Reading mask " << maskFileName << std::endl;
	MaskWReaderType::Pointer mask_reader = MaskWReaderType::New();
	mask_reader->SetFileName(maskFileName);
	mask_reader->Update();

	std::cout << "Thresholding..." << std::endl;
	BinaryThresholdFilterType::Pointer thresh = BinaryThresholdFilterType::New();
	thresh->SetLowerThreshold( lowerThresh );
  thresh->SetUpperThreshold( upperThresh );
  thresh->SetOutsideValue( 0 );
  thresh->SetInsideValue( 1 );
	thresh->SetInput( mask_reader->GetOutput());
	thresh->Update();

	std::cout << "Running graph cut filter" << std::endl;
	EndostealSegmentationFilterType::Pointer filter = EndostealSegmentationFilterType::New();
	filter->SetLambda(lambda);
	filter->SetSigma(sigma);
	filter->SetCorticalLabel(corticalLabel);
  filter->SetCancellousLabel(cancellousLabel);
	filter->SetBackgroundLabel(backgroundLabel);
	filter->SetMinDistance(minDistance);
	filter->SetMaxDistance(maxDistance);
	filter->SetInput(input_reader->GetOutput());
	filter->SetMask(thresh->GetOutput());
	filter->Update();

  std::cout << "  Max Flow: " << filter->GetMaxFlow() << std::endl;

	std::cout << "Writing result to " << outputFileName << std::endl;
	OutputWriterType::Pointer writer = OutputWriterType::New();
	writer->SetFileName(outputFileName);
	writer->SetInput(filter->GetOutput());
	writer->Update();

	std::cout << "Finished!" << std::endl;

	return EXIT_SUCCESS;
}
