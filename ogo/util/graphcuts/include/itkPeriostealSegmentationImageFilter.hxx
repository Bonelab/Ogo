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

#ifndef itkPeriostealSegmentationImageFilter_hxx
#define itkPeriostealSegmentationImageFilter_hxx

#include "itkPeriostealSegmentationImageFilter.h"

namespace itk {
template< typename TInputImage, typename TMaskImage, typename TOutputImage >
PeriostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::PeriostealSegmentationImageFilter() :
  m_BackgroundLabel(0.0),
  m_ForegroundLabel(1.0),
  m_Lambda(5.0),
  m_Sigma(0.2)
{
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
typename PeriostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >::CostType
PeriostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::ComputeDataTerm(const InputPixelType p, const LabelType l, const MaskPixelType m) 
{
  DistanceType weight = 0;

  /*
  std::cout << "----------------:   ComputeDataTerm" << std::endl;
  std::cout << "  InputPixelType:   " << p << std::endl;
  std::cout << "       LabelType:   " << l << std::endl;
  std::cout << "   MaskPixelType:   " << m << std::endl;
  std::cout << std::endl;
  */
  
  switch (l) {
    case 0:
      // {p,S}
      if (m == this->m_ForegroundLabel) {
        weight = this->m_Lambda * this->GetnNeighbours() + 1;
      } else if (m != this->m_BackgroundLabel && m != this->m_ForegroundLabel ) {
        weight = 0;
      } else if (p > 0) {
        weight = 1;
      } else {
        weight = 0;
      }
      break;
    case 1:
      // {p,T}
      if (m == this->m_ForegroundLabel) {
        weight = 0;
      } else if (m != this->m_BackgroundLabel && m != this->m_ForegroundLabel ) {
        weight = this->m_Lambda * this->GetnNeighbours() + 1;
      } else {
        weight = 1;
      }
      break;
    default:
      weight = 0;
  }

  assert(weight >= 0);
  return static_cast< CostType > (this->GetWeightScale() * weight);
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
typename PeriostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >::CostType
PeriostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::ComputeSmoothnessTerm(const InputPixelType p, const InputPixelType q, const DistanceType d, const MaskPixelType m_p, const MaskPixelType m_q)
{
  DistanceType weight = 0;

  /*
  std::cout << "----------------:   ComputeSmoothnessTerm" << std::endl;
  std::cout << "  InputPixelType:   " << p << std::endl;
  std::cout << "  InputPixelType:   " << q << std::endl;
  std::cout << "    DistanceType:   " << d << std::endl;
  std::cout << "   MaskPixelType:   " << m_p << std::endl;
  std::cout << "   MaskPixelType:   " << m_q << std::endl;
  std::cout << std::endl;
  */
  
  if ( (p > q) )
  {
    weight = std::exp(-1.0 * (std::pow(p - q, 2)) / (2.0* std::pow(this->m_Sigma, 2) ) );
  }
  else
  {
    weight = 1;
  }
  weight *= this->m_Lambda;

  assert(weight >= 0);
  return static_cast< CostType > (this->GetWeightScale() * weight);
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
PeriostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Lambda: " << this->m_Lambda << std::endl;
  os << indent << "Sigma: " << this->m_Sigma << std::endl;
  os << indent << "Background label: " << this->m_BackgroundLabel << std::endl;
  os << indent << "Foreground label: " << this->m_ForegroundLabel << std::endl;
}

} /* end namespace */

#endif /* itkPeriostealSegmentationImageFilter_hxx */
