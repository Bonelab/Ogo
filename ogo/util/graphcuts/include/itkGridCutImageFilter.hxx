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

#ifndef itkGridCutImageFilter_hxx
#define itkGridCutImageFilter_hxx

#include "itkGridCutImageFilter.h"


namespace itk {
template< typename TInputImage, typename TMaskImage, typename TOutputImage >
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::GridCutImageFilter() :
  m_nLabels(2),
  m_nNeighbours(0),
  m_nVoxels(0),
  m_MaxFlow(0.0),
  m_BlockSize(100),
  m_Grid(nullptr),
  m_WeightScale(1000.0)
{
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  /* May have specific inputs here */
  this->SetupNeighbourhood();

  /* Setup Dimensions */
  this->m_Dimensions = this->GetInput(0)->GetLargestPossibleRegion().GetSize();
  this->m_nVoxels = this->m_Dimensions[0]*this->m_Dimensions[1]*this->m_Dimensions[2];

  /* Create graph */
  this->m_Grid = std::unique_ptr<Grid>(new Grid(
    this->m_Dimensions[0], this->m_Dimensions[1], this->m_Dimensions[2],
    this->GetMultiThreader()->GetMaximumNumberOfThreads(), this->m_BlockSize
  ));

  /* Create arrays */
  this->m_tLinks.resize(this->m_nLabels);
  for (unsigned int i = 0; i < this->m_nLabels; ++i)
  {
    this->m_tLinks[i].resize(this->m_nVoxels, 0);
  }

  this->m_nLinks.resize(this->m_nNeighbours);
  for (unsigned int i = 0; i < this->m_nNeighbours; ++i)
  {
    this->m_nLinks[i].resize(this->m_nVoxels, 0);
  }

  /* Reset max flow */
  this->m_MaxFlow = 0.0;
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::AfterThreadedGenerateData()
{
  /* Solve the graph */
  m_Grid->set_caps(
    this->m_tLinks[0].data(),  // cap_source
    this->m_tLinks[1].data(),  // cap_sink

    this->m_nLinks[0].data(),  // [-1, 0, 0]
    this->m_nLinks[1].data(),  // [+1, 0, 0]
    this->m_nLinks[2].data(),  // [ 0,-1, 0]
    this->m_nLinks[3].data(),  // [ 0,+1, 0]
    this->m_nLinks[4].data(),  // [ 0, 0,-1]
    this->m_nLinks[5].data()   // [ 0, 0,+1]
  );
  m_Grid->compute_maxflow();

  /* Read out */
  this->m_MaxFlow = m_Grid->get_flow() / this->m_WeightScale;

  OutputImagePointer          output = this->GetOutput(0);
  OutputIteratorType          ot(output, output->GetLargestPossibleRegion());

  ot.GoToBegin();
  while ( !ot.IsAtEnd() )
  {
    auto p = ot.GetIndex();
    auto id = this->m_Grid->node_id(p[0], p[1], p[2]);
    ot.Set( this->GetLabel( this->m_Grid->get_segment(id) ) );

    ++ot;
  }

  /* Free memory */
  for (unsigned int i = 0; i < this->m_nLabels; ++i)
  {
    this->m_tLinks[i].resize(0);
  }
  this->m_tLinks.resize(0);

  for (unsigned int i = 0; i < this->m_nNeighbours; ++i)
  {
    this->m_nLinks[i].resize(0);
  }
  this->m_nLinks.resize(0);

  this->m_Grid.release();
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread)
{
  /* Get Inputs */
  InputImageConstPointer input = this->GetInput(0);
  MaskImageConstPointer mask = this->GetMask();
  InputImageRegionType image_region = input->GetLargestPossibleRegion();

  /* Setup regions */
  InputImageRegionType inputRegionForThread;
  MaskImageRegionType maskRegionForThread;
  this->CallCopyOutputRegionToInputRegion(inputRegionForThread, outputRegionForThread);
  this->CallCopyOutputRegionToInputRegion(maskRegionForThread, outputRegionForThread);

  /* Setup Iterator */
  typename ShapedIteratorType::RadiusType radius;
  radius.Fill(1);
  ShapedIteratorType it(radius, input, inputRegionForThread);

  typename ShapedIteratorType::OffsetType center = {{0, 0, 0}};

  MaskIterator mi(radius, mask, maskRegionForThread);

  /* Iterate and perform logic */
  bool pixelIsValid;
  for( it.GoToBegin(), mi.GoToBegin(); !it.IsAtEnd(); ++it, ++mi )
  {
    /* Get index and iterator */
    auto p_value = it.GetPixel(center);
    auto p = it.GetIndex(center);

    auto m_p_value = mi.GetPixel(center);

    /* Process data term */
    for (LabelType l = 0; l < this->m_nLabels; ++l)
    {
      this->SetDataTerm(p, l, this->ComputeDataTerm(p_value, l, m_p_value));
    }

    /* Process smooth term */
    for (unsigned int i = 0; i < this->m_nNeighbours; ++i)
    {
      /* Get q index, see if inside */
      auto q = it.GetIndex( this->m_Neighbors[i] );
      auto q_value = it.GetPixel( this->m_Neighbors[i], pixelIsValid );
      auto m_q_value = mi.GetPixel( this->m_Neighbors[i] );
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
void
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::SetupNeighbourhood()
{
  this->m_nNeighbours = 6;
  this->m_Neighbors.resize(this->m_nNeighbours);

  this->m_Neighbors[0] = {{-1,0,0}};
  this->m_Neighbors[1] = {{+1,0,0}};
  this->m_Neighbors[2] = {{0,-1,0}};
  this->m_Neighbors[3] = {{0,+1,0}};
  this->m_Neighbors[4] = {{0,0,-1}};
  this->m_Neighbors[5] = {{0,0,+1}};
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
typename GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >::OutputImagePixelType
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::GetLabel(const LabelType l) const
{
  return static_cast< OutputImagePixelType > ( l == 0 );
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::SetDataTerm(const IndexType p, const LabelType l, const CostType cost)
{
  auto id = this->GetIndex(p);
  this->m_tLinks[l][id] = cost;
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::SetSmoothTerm(const IndexType p, const LabelType n_i, const CostType cost)
{
  auto id = this->GetIndex(p);
  this->m_nLinks[n_i][id] = cost;
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
unsigned long
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::GetIndex(const IndexType p)
{
  unsigned int id = p[0] + 
    p[1]*this->m_Dimensions[0] +
    p[2]*this->m_Dimensions[1] * this->m_Dimensions[0];

  return id;
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
typename GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >::CostType
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::ComputeDataTerm(const InputPixelType p, const LabelType l, const MaskPixelType m)
{
  itkExceptionMacro(<< "ComputeDataTerm not overwritten");
}


template< typename TInputImage, typename TMaskImage, typename TOutputImage >
typename GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >::CostType
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::ComputeSmoothnessTerm(const InputPixelType p, const InputPixelType q, const DistanceType d, const MaskPixelType m_p, const MaskPixelType m_q)
{
  itkExceptionMacro(<< "ComputeSmoothnessTerm not overwritten");
}

template< typename TInputImage, typename TMaskImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Number of labels: " << this->m_nLabels << std::endl;
  os << indent << "Neighbourhood size: " << this->m_nNeighbours << std::endl;
  os << indent << "Block size: " << this->m_BlockSize << std::endl;
  os << indent << "Max flow: " << this->m_MaxFlow << std::endl;
  os << indent << "Weight scale: " << this->m_WeightScale << std::endl;
}

} /* end namespace */

#endif /* itkGridCutImageFilter_hxx */
