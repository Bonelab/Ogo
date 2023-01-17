 /*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef itkCalgaryEigenToMeasureParameterEstimationFilter_hxx
#define itkCalgaryEigenToMeasureParameterEstimationFilter_hxx

#include "itkCalgaryEigenToMeasureParameterEstimationFilter.h"

namespace itk {

template< typename TInputImage, typename TOutputImage >
CalgaryEigenToMeasureParameterEstimationFilter< TInputImage, TOutputImage >
::CalgaryEigenToMeasureParameterEstimationFilter() :
  Superclass(),
  m_FrobeniusNormWeight(0.5)
{
  /* Set parameter size to 3 */
  ParameterArrayType parameters = this->GetParametersOutput()->Get();
  parameters.SetSize(2);
  parameters[0] = 0.5;
  parameters[1] = 1.0;
  this->GetParametersOutput()->Set(parameters);
}

template< typename TInputImage, typename TOutputImage >
void
CalgaryEigenToMeasureParameterEstimationFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  m_MaxFrobeniusNorm = NumericTraits< RealType >::NonpositiveMin();
}

template< typename TInputImage, typename TOutputImage >
void
CalgaryEigenToMeasureParameterEstimationFilter< TInputImage, TOutputImage >
::AfterThreadedGenerateData()
{
  /* Determine default parameters */
  RealType alpha, c;
  alpha = 0.5f;

  /* Scale c */
  c = m_FrobeniusNormWeight * m_MaxFrobeniusNorm;

  /* Assign outputs parameters */
  ParameterArrayType parameters;
  parameters.SetSize(2);
  parameters[0] = alpha;
  parameters[1] = c;
  this->GetParametersOutput()->Set(parameters);
}

template< typename TInputImage, typename TOutputImage >
void
CalgaryEigenToMeasureParameterEstimationFilter< TInputImage, TOutputImage >
::DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread)
{
  /* If size is zero, return */
  const SizeValueType size0 = outputRegionForThread.GetSize(0);
  if (size0 == 0)
  {
    return;
  }

  /* Keep track of the current max */
  RealType max = NumericTraits< RealType >::NonpositiveMin();

  /* Get input and mask pointer */
  InputImageConstPointer inputPointer = this->GetInput();
  MaskSpatialObjectTypeConstPointer maskPointer = this->GetMask();
  typename InputImageType::PointType point;

  OutputImageType      *outputPtr = this->GetOutput(0);
  
  // Define the portion of the input to walk for this thread, using
  // the CallCopyOutputRegionToInputRegion method allows for the input
  // and output images to be different dimensions
  InputImageRegionType inputRegionForThread;

  this->CallCopyOutputRegionToInputRegion(inputRegionForThread, outputRegionForThread);

  /* Setup iterator */
  ImageRegionConstIteratorWithIndex< TInputImage > inputIt(inputPointer, inputRegionForThread);
  ImageRegionIterator< OutputImageType > outputIt(outputPtr, outputRegionForThread);

  /* Iterate and count */
  inputIt.GoToBegin();
  outputIt.GoToBegin();
  while ( !inputIt.IsAtEnd() )
  {
    // Process point
    inputPointer->TransformIndexToPhysicalPoint(inputIt.GetIndex(), point);
    if ( (!maskPointer) ||  (maskPointer->IsInsideInWorldSpace(point)) )
    {
      /* Compute max norm */
      max = std::max( max, this->CalculateFrobeniusNorm(inputIt.Get()) );
    }

    // Set 
    outputIt.Set( static_cast< OutputImagePixelType >( inputIt.Get() ) );

    // Increment
    ++inputIt;
    ++outputIt;
  }

  /* Block and store */
  std::lock_guard<std::mutex> guard(m_Mutex);
  m_MaxFrobeniusNorm = std::max( m_MaxFrobeniusNorm, max );
}

template< typename TInputImage, typename TOutputImage >
typename CalgaryEigenToMeasureParameterEstimationFilter< TInputImage, TOutputImage >::RealType
CalgaryEigenToMeasureParameterEstimationFilter< TInputImage, TOutputImage >
::CalculateFrobeniusNorm(const InputImagePixelType& pixel) const
{
  RealType norm = 0;
  for( unsigned int i = 0; i < pixel.Length; ++i) {
    norm += pixel[i]*pixel[i];
  }
  return sqrt(norm);
}

template< typename TInputImage, typename TOutputImage >
void
CalgaryEigenToMeasureParameterEstimationFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "FrobeniusNormWeight: " << GetFrobeniusNormWeight() << std::endl;
}

} /* end namespace */

#endif /* itkCalgaryEigenToMeasureParameterEstimationFilter_hxx */
