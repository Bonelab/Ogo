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

#ifndef itkEndostealSegmentationImageFilter_hxx
#define itkEndostealSegmentationImageFilter_hxx

#include "itkEndostealSegmentationImageFilter.h"

namespace itk {
template< typename TInputImage, typename TMaskImage, typename TOutputImage >
EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::EndostealSegmentationImageFilter() :
  m_BackgroundLabel(0.0),
  m_CorticalLabel(1.0),
  m_CancellousLabel(2.0),
  m_Lambda(5.0),
  m_Sigma(0.2),
  m_MaxDistance(2.0),
	m_MinDistance(1.0)
{
  m_DistanceFilter = nullptr;
  m_MaskingFilter = nullptr;
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  /* Set everything up */
  Superclass::BeforeThreadedGenerateData();

  /* Create the distance image */
  m_DistanceFilter = DistanceFilterType::New();
  m_DistanceFilter->SetInput(this->GetMask());
  m_DistanceFilter->SquaredDistanceOff();
  m_DistanceFilter->UseImageSpacingOn();
  m_DistanceFilter->InsideIsPositiveOn();
  m_DistanceFilter->Update();
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::AfterThreadedGenerateData()
{
  /* Set everything up */
  Superclass::AfterThreadedGenerateData();

  /* Mask by input */
  m_MaskingFilter = MaskFilterType::New();
  m_MaskingFilter->SetInput(this->GetOutput(0));
  m_MaskingFilter->SetMaskImage(this->GetMask());
  m_MaskingFilter->SetMaskingValue(0);
  m_MaskingFilter->SetOutsideValue(0);
  m_MaskingFilter->Update();
  this->GraftOutput(m_MaskingFilter->GetOutput());

  /* Clean up filters */
  m_DistanceFilter = nullptr;
  m_MaskingFilter = nullptr;
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread)
{
  /* Get Inputs */
  InputImageConstPointer input = this->GetInput(0);
  MaskImageConstPointer mask = this->GetMask();
  InputImageRegionType image_region = input->GetLargestPossibleRegion();
  DistanceConstImageType dist = m_DistanceFilter->GetOutput();

  /* Setup regions */
  DistanceImageRegionType distanceRegionForThread;
  InputImageRegionType inputRegionForThread;
  MaskImageRegionType maskRegionForThread;
  this->CallCopyOutputRegionToInputRegion(inputRegionForThread, outputRegionForThread);
  this->CallCopyOutputRegionToInputRegion(maskRegionForThread, outputRegionForThread);
  this->CallCopyOutputRegionToInputRegion(distanceRegionForThread, outputRegionForThread);

  /* Setup Iterator */
  typename ShapedIteratorType::RadiusType radius;
  radius.Fill(1);
  ShapedIteratorType it(radius, input, inputRegionForThread);

  typename ShapedIteratorType::OffsetType center = {{0, 0, 0}};

  MaskIterator mi(radius, mask, maskRegionForThread);
  DistanceIteratorType di(radius, dist, distanceRegionForThread);

  /* Iterate and perform logic */
  bool pixelIsValid;
  for( it.GoToBegin(), mi.GoToBegin(), di.GoToBegin(); !it.IsAtEnd(); ++it, ++mi, ++di )
  {
    /* Get index and iterator */
    auto p_value = it.GetPixel(center);
    auto p = it.GetIndex(center);
    auto d_value = di.GetPixel(center);

    auto m_p_value = mi.GetPixel(center);

    /* Process data term */
    for (LabelType l = 0; l < this->GetnLabels(); ++l)
    {
      this->SetDataTerm(p, l, this->ComputeDataTerm(p_value, d_value, l, m_p_value));
    }

    /* Process smooth term */
    for (unsigned int i = 0; i < this->GetnNeighbours(); ++i)
    {
      /* Get q index, see if inside */
      auto q = it.GetIndex( this->GetNeighbors()[i] );
      auto q_value = it.GetPixel( this->GetNeighbors()[i], pixelIsValid );
      auto m_q_value = mi.GetPixel( this->GetNeighbors()[i] );
      if (!pixelIsValid)
      {
        continue;
      }

      /* Compute distance */
      typename InputImageType::PointType p_point, q_point;
      input->TransformIndexToPhysicalPoint(p, p_point);
      input->TransformIndexToPhysicalPoint(q, q_point);
      DistanceType d = p_point.EuclideanDistanceTo(q_point);

      /* Compute smooth term */
      this->SetSmoothTerm(p, i, this->ComputeSmoothnessTerm(p_value, q_value, d, m_p_value, m_q_value));
    }
  }
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
typename EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >::CostType
EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::ComputeDataTerm(const InputPixelType p, const RealType d, const LabelType l, const MaskPixelType m)
{
  DistanceType weight = 0;

  switch (l) {
    case 0:
      // {p,S}
      if (m == this->m_BackgroundLabel) {
        weight = 0;
      } else if (d > m_MaxDistance) {
        weight = 0;
			} else if (d < m_MinDistance) {
				weight = this->m_Lambda * this->GetnNeighbours() + 1;
      } else {
        weight = 1;
      }
      break;
    case 1:
      // {p,T}
      if (m == this->m_BackgroundLabel) {
        weight = this->m_Lambda * this->GetnNeighbours() + 1;
      } else if (d > m_MaxDistance) {
        weight = 1;
			} else if (d < m_MinDistance) {
				weight = 0;
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
typename EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >::CostType
EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::ComputeSmoothnessTerm(const InputPixelType p, const InputPixelType q, const DistanceType d, const MaskPixelType m_p, const MaskPixelType m_q)
{
  DistanceType weight = 0;

  if ( (m_p == this->m_BackgroundLabel) || (m_q == this->m_BackgroundLabel) ) {
    weight = 0;
  } else if ( (p > q) )
  {
    weight = std::exp(-1.0 * (std::pow(p - q, 2)) / (2.0 * std::pow(this->m_Sigma, 2) ) );
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
typename EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >::OutputImagePixelType
EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::GetLabel(const LabelType l) const
{
  switch(l) {
    case 0:
      return this->m_CorticalLabel;
    case 1:
      return this->m_CancellousLabel;
    default:
      return this->m_BackgroundLabel;
  }
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
EndostealSegmentationImageFilter< TInputImage, TMaskImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Lambda: " << this->m_Lambda << std::endl;
  os << indent << "Sigma: " << this->m_Sigma << std::endl;
  os << indent << "Background label: " << this->m_BackgroundLabel << std::endl;
  os << indent << "Cortical label: " << this->m_CorticalLabel << std::endl;
  os << indent << "Cancellous label: " << this->m_CancellousLabel << std::endl;
}

} /* end namespace */

#endif /* itkEndostealSegmentationImageFilter_hxx */
