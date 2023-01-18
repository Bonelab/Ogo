#include <iostream>

#include "itkPeriostealSegmentationImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

/* Type definitions */
constexpr unsigned int ImageDimension = 3;
using InputPixelType 	= float;
using MaskPixelType 	= unsigned long;
using OutputPixelType = unsigned char;

using InputImageType	= itk::Image< InputPixelType, ImageDimension >;
using MaskImageType		= itk::Image< MaskPixelType, ImageDimension >;
using OutputImageType	= itk::Image< OutputPixelType, ImageDimension >;

using InputReaderType 	= itk::ImageFileReader< InputImageType >;
using MaskWReaderType		= itk::ImageFileReader< MaskImageType >;
using OutputWriterType	= itk::ImageFileWriter< OutputImageType >;

using PeriostealSegmentationFilterType = itk::PeriostealSegmentationImageFilter< InputImageType, MaskImageType, OutputImageType >;

using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter< OutputImageType, MaskImageType >;
using LabelShapeKeepNObjectsImageFilterType = itk::LabelShapeKeepNObjectsImageFilter< MaskImageType >;
using ThresholdFilterType = itk::BinaryThresholdImageFilter< MaskImageType, OutputImageType >;

int main(int argc, char** argv) {
  if( argc != 8 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <MaskFileName> <OutputSegmentation> ";
		std::cerr << " <Lambda> <Sigma> <Label> <ConnFilter>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

	/* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string maskFileName = argv[2];
  std::string outputFileName = argv[3];

  double lambda = atof(argv[4]);
	double sigma = atof(argv[5]);
	int label = atoi(argv[6]);
	int connFilter = atoi(argv[7]);

	std::cout << "Parameters:" << std::endl;
  std::cout << "  InputFilePath:    " << inputFileName << std::endl;
  std::cout << "  MaskFilePath:     " << maskFileName << std::endl;
  std::cout << "  OutputFilePath:   " << outputFileName << std::endl;
  std::cout << "  Lambda:           " << lambda << std::endl;
  std::cout << "  Sigma:            " << sigma << std::endl;
	std::cout << "  Label:            " << label << std::endl;
	std::cout << "  ConnFilter:       " << connFilter << std::endl;
  std::cout << std::endl;

	std::cout << "Reading input " << inputFileName << std::endl;
	InputReaderType::Pointer input_reader = InputReaderType::New();
	input_reader->SetFileName(inputFileName);
	input_reader->Update();

	std::cout << "Reading mask " << maskFileName << std::endl;
	MaskWReaderType::Pointer mask_reader = MaskWReaderType::New();
	mask_reader->SetFileName(maskFileName);
	mask_reader->Update();

	std::cout << "Running graph cut filter" << std::endl;
	PeriostealSegmentationFilterType::Pointer filter = PeriostealSegmentationFilterType::New();
	filter->SetLambda(lambda);
	filter->SetSigma(sigma);
	filter->SetForegroundLabel(label);
	filter->SetBackgroundLabel(0);
	filter->SetInput(input_reader->GetOutput());
	filter->SetMask(mask_reader->GetOutput());
	filter->Update();

	std::cout << "  Max Flow: " << filter->GetMaxFlow() << std::endl;

	std::cout << "Running connectivity filter" << std::endl;
  ConnectedComponentImageFilterType::Pointer fgConnected = ConnectedComponentImageFilterType::New ();
  fgConnected->SetInput(filter->GetOutput());

  LabelShapeKeepNObjectsImageFilterType::Pointer fgKeeper = LabelShapeKeepNObjectsImageFilterType::New();
  fgKeeper->SetInput( fgConnected->GetOutput() );
  fgKeeper->SetBackgroundValue( 0 );
  fgKeeper->SetNumberOfObjects( connFilter );
  fgKeeper->SetAttribute( LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
	fgKeeper->Update();

	std::cout << "  Found " << fgConnected->GetObjectCount() << " foreground objects" << std::endl;

	ThresholdFilterType::Pointer thresh = ThresholdFilterType::New();
	thresh->SetInput(fgKeeper->GetOutput());
	thresh->SetLowerThreshold(1);
	thresh->SetInsideValue (0);
	thresh->SetOutsideValue(1);
	thresh->Update();

  ConnectedComponentImageFilterType::Pointer bkgConnected = ConnectedComponentImageFilterType::New ();
  bkgConnected->SetInput(thresh->GetOutput());

  LabelShapeKeepNObjectsImageFilterType::Pointer bkgKeeper = LabelShapeKeepNObjectsImageFilterType::New();
  bkgKeeper->SetInput( bkgConnected->GetOutput() );
  bkgKeeper->SetBackgroundValue( 0 );
  bkgKeeper->SetNumberOfObjects( 1 );
  bkgKeeper->SetAttribute( LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
	bkgKeeper->Update();

	std::cout << "  Found " << bkgConnected->GetObjectCount() << " background objects" << std::endl;

	ThresholdFilterType::Pointer thresh2 = ThresholdFilterType::New();
	thresh2->SetInput(bkgKeeper->GetOutput());
	thresh2->SetLowerThreshold(1);
	thresh2->SetInsideValue (0);
	thresh2->SetOutsideValue(1);
	thresh2->Update();

	std::cout << "Writing result to " << outputFileName << std::endl;
	OutputWriterType::Pointer writer = OutputWriterType::New();
	writer->SetFileName(outputFileName);
	writer->SetInput(thresh2->GetOutput());
	writer->Update();

	std::cout << "Finished!" << std::endl;

	return EXIT_SUCCESS;
}
