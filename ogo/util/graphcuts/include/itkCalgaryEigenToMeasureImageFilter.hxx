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

#ifndef itkCalgaryEigenToMeasureImageFilter_hxx
#define itkCalgaryEigenToMeasureImageFilter_hxx

#include "itkCalgaryEigenToMeasureImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"

namespace itk {
template< typename TInputImage, typename TOutputImage >
CalgaryEigenToMeasureImageFilter< TInputImage, TOutputImage >
::CalgaryEigenToMeasureImageFilter() :
  m_EnhanceType(-1.0)
{}

template< typename TInputImage, typename TOutputImage >
void
CalgaryEigenToMeasureImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  ParameterArrayType parameters = this->GetParametersInput()->Get();
  if (parameters.GetSize() != 2)
  {
    itkExceptionMacro(<< "Parameters must have size 2. Given array of size " << parameters.GetSize());
  }
}

template< typename TInputImage, typename TOutputImage >
typename CalgaryEigenToMeasureImageFilter< TInputImage, TOutputImage >::OutputImagePixelType
CalgaryEigenToMeasureImageFilter< TInputImage, TOutputImage >
::ProcessPixel(const InputImagePixelType& pixel)
{
  /* Grab parameters */
  ParameterArrayType parameters = this->GetParametersInput()->Get();
  RealType alpha = parameters[0];
  RealType c = parameters[1];

  /* Grab pixel values */
  double sheetness = 0.0;
  double a1 = static_cast<double>( pixel[0] );
  double a2 = static_cast<double>( pixel[1] );
  double a3 = static_cast<double>( pixel[2] );
  double l1 = Math::abs(a1);
  double l2 = Math::abs(a2);
  double l3 = Math::abs(a3);

  /* Avoid divisions by zero (or close to zero) */
  if ( static_cast<double>( l3 ) < Math::eps ) {
      return static_cast<OutputImagePixelType>( sheetness );
  }

  /**
   * Compute sheet, noise, and tube like measures.
   */
  const double Rbone = (l1*l2)/(l3*l3);
  const double Rnoise = sqrt(l1*l1 + l2*l2 + l3*l3);

  /* Multiply together to get sheetness */
  sheetness = (m_EnhanceType*a3/l3);
  sheetness *= std::exp(-(Rbone * Rbone) / (2 * alpha * alpha));
  sheetness *= (1.0 - std::exp(-(Rnoise * Rnoise) / (2 * c * c)));

  return static_cast<OutputImagePixelType>( sheetness );
}

template< typename TInputImage, typename TOutputImage >
void
CalgaryEigenToMeasureImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Direction: " << GetEnhanceType() << std::endl;
}

} /* end namespace */

#endif /* itkCalgaryEigenToMeasureImageFilter_hxx */
