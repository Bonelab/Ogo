#include <iostream>

#include "itkArray.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMultiScaleHessianEnhancementImageFilter.h"
#include "itkCalgaryEigenToMeasureImageFilter.h"
#include "itkCalgaryEigenToMeasureParameterEstimationFilter.h"
#include "itkCommand.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkFlatStructuringElement.h"

/* Setup Types */
constexpr unsigned int ImageDimension = 3;
using InputPixelType = short;
using InputImageType = itk::Image<InputPixelType, ImageDimension>;
using MaskPixelType = unsigned char;
using MaskImageType = itk::Image<MaskPixelType, ImageDimension>;
using CCPixelType = unsigned long;
using CCImageType = itk::Image< CCPixelType, ImageDimension >;
using OutputPixelType = float;
using OutputImageType = itk::Image<OutputPixelType, ImageDimension>;

using ReaderType = itk::ImageFileReader< InputImageType >;
using MaskWriterType = itk::ImageFileWriter< MaskImageType >;
using MeasureWriterType = itk::ImageFileWriter< OutputImageType >;
using MaskSpatialObject = itk::ImageMaskSpatialObject< ImageDimension >;

using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter< MaskImageType, CCImageType >;
using LabelShapeKeepNObjectsImageFilterType = itk::LabelShapeKeepNObjectsImageFilter< CCImageType >;
using MaskThesholdFilter = itk::BinaryThresholdImageFilter< CCImageType, MaskImageType >;
using MaskThesholdFilter2 = itk::BinaryThresholdImageFilter< MaskImageType, MaskImageType >;

using StructuringElementType = itk::FlatStructuringElement< ImageDimension  >;
using ErodeFilterType = itk::BinaryErodeImageFilter< MaskImageType, MaskImageType, StructuringElementType >;

using BinaryThesholdFilter = itk::BinaryThresholdImageFilter< InputImageType, MaskImageType >;
using MultiScaleHessianFilterType = itk::MultiScaleHessianEnhancementImageFilter< InputImageType, OutputImageType >;
using CalgaryEigenToMeasureImageFilterType = itk::CalgaryEigenToMeasureImageFilter< MultiScaleHessianFilterType::EigenValueImageType, OutputImageType >;
using CalgaryEigenToMeasureParameterEstimationFilterType = itk::CalgaryEigenToMeasureParameterEstimationFilter< MultiScaleHessianFilterType::EigenValueImageType >;

int main(int argc, char * argv[])
{
  if( argc != 11 )
  {
    std::cerr << "Usage: "<< std::endl;
    std::cerr << argv[0];
    std::cerr << " <InputFileName> <OutputMask> <OutputMeasure> ";
    std::cerr << " <SetEnhanceBrightObjects[0,1]> ";
    std::cerr << " <NumberOfSigma> <MinSigma> <MaxSigma> ";
    std::cerr << " <LowThreshold> <HighThreshold> <Weight>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  /* Read input Parameters */
  std::string inputFileName = argv[1];
  std::string maskFileName = argv[2];
  std::string outputMeasureFileName = argv[3];

  int enhanceBrightObjects = atoi(argv[4]);
  int numberOfSigma = atoi(argv[5]);
  double minSigma = atof(argv[6]);
  double maxSigma = atof(argv[7]);
  double lowThreshold = atof(argv[8]);
  double highThreshold = atof(argv[9]);
  double weight = atof(argv[10]);

  std::cout << "Read in the following parameters:" << std::endl;
  std::cout << "  InputFilePath:               " << inputFileName << std::endl;
  std::cout << "  MaskFilePath:                " << maskFileName << std::endl;
  std::cout << "  OutputMeasure:               " << outputMeasureFileName << std::endl;
  if (enhanceBrightObjects == 1) {
    std::cout << "  SetEnhanceBrightObjects:     " << "Enhancing bright objects" << std::endl;
  } else {
    std::cout << "  SetEnhanceBrightObjects:     " << "Enhancing dark objects" << std::endl;
  }
  std::cout << "  NumberOfSigma:               " << numberOfSigma << std::endl;
  std::cout << "  Minimum Sigma:               " << minSigma << std::endl;
  std::cout << "  Maximum Sigma:               " << maxSigma << std::endl;
  std::cout << "  Low Threshold:               " << lowThreshold << std::endl;
  std::cout << "  High Threshold:              " << highThreshold << std::endl;
  std::cout << "  Weight:                      " << weight << std::endl;
  std::cout << std::endl;

  /* Do preprocessing */
  std::cout << "Reading in " << inputFileName << std::endl;
  ReaderType::Pointer  reader = ReaderType::New();
  reader->SetFileName(inputFileName);
  reader->Update();

  std::cout << "Creating mask by thresholding outside [" << lowThreshold << "-" << highThreshold << "]" << std::endl;
  BinaryThesholdFilter::Pointer thresholder = BinaryThesholdFilter::New();
  thresholder->SetInput(reader->GetOutput());
  thresholder->SetLowerThreshold(lowThreshold);
  thresholder->SetUpperThreshold(highThreshold);
	thresholder->SetInsideValue(0);
	thresholder->SetOutsideValue(1);
  thresholder->Update();

	std::cout << "Connected components filter on background" << std::endl;
  ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New ();
  connected->SetInput(thresholder->GetOutput());

  LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter = LabelShapeKeepNObjectsImageFilterType::New();
  labelShapeKeepNObjectsImageFilter->SetInput( connected->GetOutput() );
  labelShapeKeepNObjectsImageFilter->SetBackgroundValue( 0 );
  labelShapeKeepNObjectsImageFilter->SetNumberOfObjects( 1 );
  labelShapeKeepNObjectsImageFilter->SetAttribute( LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);

  MaskThesholdFilter::Pointer inverter = MaskThesholdFilter::New();
  inverter->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());
	inverter->SetUpperThreshold(0);
	inverter->Update();

  std::cout << "Writing mask to " << maskFileName << std::endl;
  MaskWriterType::Pointer maskWriter = MaskWriterType::New();
  maskWriter->SetInput(inverter->GetOutput());
  maskWriter->SetFileName(maskFileName);
  maskWriter->Write();

  std::cout << "Eroding skin image" << std::endl;
  std::cout << "  Radius: ";
  StructuringElementType::RadiusType radius;
  radius.Fill( 1 );
  for (unsigned int i = 0; i < ImageDimension; ++i)
  {
    radius[i] = static_cast< unsigned int >(maxSigma / inverter->GetOutput()->GetSpacing()[i]);
    std::cout << radius[i] << " ";
  }
  std::cout << std::endl;

  MaskThesholdFilter2::Pointer inverter2 = MaskThesholdFilter2::New();
  inverter2->SetInput(thresholder->GetOutput());
	inverter2->SetUpperThreshold(0);
	inverter2->Update();

  StructuringElementType structuringElement = StructuringElementType::Ball( radius );
  ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
  erodeFilter->SetInput( inverter2->GetOutput() );
  erodeFilter->SetKernel( structuringElement );
  erodeFilter->Update();

	MaskSpatialObject::Pointer skinMaskSpatialObject = MaskSpatialObject::New();
	skinMaskSpatialObject->SetImage(inverter->GetOutput());

	MaskSpatialObject::Pointer erodedMaskSpatialObject = MaskSpatialObject::New();
	erodedMaskSpatialObject->SetImage(erodeFilter->GetOutput());

  /* Multiscale measure */
  MultiScaleHessianFilterType::Pointer multiScaleFilter = MultiScaleHessianFilterType::New();

  MultiScaleHessianFilterType::SigmaArrayType sigmaArray = multiScaleFilter->GenerateLogarithmicSigmaArray(minSigma, maxSigma, numberOfSigma);
  std::cout << "Sigma Array: " << sigmaArray << std::endl;

  CalgaryEigenToMeasureParameterEstimationFilterType::Pointer estimationFilter = CalgaryEigenToMeasureParameterEstimationFilterType::New();
  CalgaryEigenToMeasureImageFilterType::Pointer calgaryFilter = CalgaryEigenToMeasureImageFilterType::New();

  estimationFilter->SetMask(erodedMaskSpatialObject);
	calgaryFilter->SetMask(skinMaskSpatialObject);

  estimationFilter->SetFrobeniusNormWeight(weight);

  std::cout << "Running multiScaleFilter..." << std::endl;
  multiScaleFilter->SetInput(reader->GetOutput());
  multiScaleFilter->SetEigenToMeasureImageFilter(calgaryFilter);
  multiScaleFilter->SetEigenToMeasureParameterEstimationFilter(estimationFilter);
  multiScaleFilter->SetSigmaArray(sigmaArray);
  multiScaleFilter->Update();

  std::cout << "Writing results to " << outputMeasureFileName << std::endl;
  MeasureWriterType::Pointer measureWriter = MeasureWriterType::New();
  measureWriter->SetInput(multiScaleFilter->GetOutput());
  measureWriter->SetFileName(outputMeasureFileName);
  measureWriter->Write();

  return EXIT_SUCCESS;
}
