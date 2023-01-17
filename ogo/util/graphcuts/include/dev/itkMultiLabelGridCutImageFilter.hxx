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
#include "itkImageRegionIteratorWithIndex.h"

#include <iostream>

namespace itk {
template< typename TInputImage, typename TOutputImage >
GridCutImageFilter< TInputImage, TOutputImage >
::GridCutImageFilter() :
  m_nLabels(0),
  m_nNeighbours(0),
  m_nVoxels(0),
  m_BackgroundLabel(0)
  m_BlockSize(100)
{
}

template< typename TInputImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  /* May have specific inputs here */
  this->SetupNeighbourhood();
  this->SetupLabels();

  /* Setup Dimensions */
  SizeType dimensions = this->GetInput(0)->GetLargestPossibleRegion().GetSize();
  this->m_nVoxels = dimensions[0]*dimensions[1]*dimensions[2];

  /* Setup arrays */
  this->m_Data.resize(this->m_nVoxels * this->m_nLabels, 0);
  this->m_Smooth.resize(this->m_nVoxels * this->m_nNeighbours / 2);
  for (auto&& table :  this->m_Smooth){
   	table.resize(this->m_nLabels * this->m_nLabels, 0);
  }
}

template< typename TInputImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TOutputImage >
::AfterThreadedGenerateData()
{
  /* Solve the graph */
  auto grid = std::make_unique<Grid>(width, height, depth, nLabels, this->m_Data, this->m_Smooth, this->GetMultiThreader()->GetMaximumNumberOfThreads(), block_size);
  grid->set_labels(this->m_BackgroundLabel);
  grid->perform();

  /* Read out */


  /* Free memory */
  for (auto&& table :  this->m_Smooth){
   	table.resize(0);
  }
  this->m_Smooth.resize(0);
  this->m_Data.resize(0);
}

template< typename TInputImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TOutputImage >
::DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread)
{
  /* Get Inputs */
  InputImageConstPointer input = this->GetInput(0);
  InputImageRegionType image_region = input->GetLargestPossibleRegion();

  /* Setup Iterator */
  InputImageRegionType inputRegionForThread;
  this->CallCopyOutputRegionToInputRegion(inputRegionForThread, outputRegionForThread);
  ImageRegionIteratorWithIndex< TInputImage > it(input, inputRegionForThread);

  /* Iterate and perform logic */
  it.GoToBegin();
  while ( !it.IsAtEnd() )
  {
    auto p = it.GetIndex();
    auto p_i = input.ComputeOffset(p);

    /* Process data term */
    for (LabelType l = 0; l < this->m_nLabels; ++l)
    {
      this->SetDataTerm(p_i, l,
        this->ComputeDataTerm(p, l);
      );
    }

    /* Process smooth term */
    for (auto&& offset: this->m_Neighbors)
    {
      /* Get q index, see if inside */
      auto q = p + offset;
      if (!image_region.IsInside(q))
      {
        continue;
      }
      auto q_i = input.ComputeOffset(q);

      for (LabelType l_p = 0; l_p < this->m_nLabels; ++l_p)
      {
        for (LabelType l_q = 0; l_q < this->m_nLabels; ++l_q)
        {
          /* Skip if labels are the same */
          if (l_p == l_q)
          {
            continue;
          }

          /* Compute smooth term */
          this->SetDataTerm(p_i, q_i, l_p, l_q,
            this->ComputeSmoothnessTerm(p, q, l_p, l_q);
          );
        }
      }
    }
  }
}

template< typename TInputImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TOutputImage >
::SetupNeighbourhood()
{
  this->m_nNeighbours = 6;
  this->m_Neighbors.resize(this->m_nNeighbours);
  this->m_Neighbors[0] = {{+1,0,0}};
  this->m_Neighbors[1] = {{0,+1,0}};
  this->m_Neighbors[2] = {{0,0,+1}};
  this->m_Neighbors[3] = {{-1,0,0}};
  this->m_Neighbors[4] = {{0,-1,0}};
  this->m_Neighbors[5] = {{0,0,-1}};

  /* TODO: Setup which graph function we use here*/
}

template< typename TInputImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TOutputImage >
::SetupLabels()
{
  /* There is where you would loop over an input image to determine the number of inputs */
  this->m_nLabels = 2;
}

template< typename TInputImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TOutputImage >
::SetDataTerm(const LabelType p_i, const LabelType l, const EnergyType cost)
{
  this->m_Data[p_i * this->m_nLabels + l] = cost;
}

template< typename TInputImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TOutputImage >
::SetSmoothTerm(const LabelType p_i, const LabelType q_i, const LabelType l_p, const LabelType l_q, const EnergyType cost)
{

  this->m_Smooth[p_i * this->m_nNeighbours / 2 + ][] = cost;
}

template< typename TInputImage, typename TOutputImage >
void
GridCutImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent)
{
  Superclass::PrintSelf(os, indent);
  // os << indent << "Direction: " << GetEnhanceType() << std::endl;
}

} /* end namespace */

#endif /* itkGridCutImageFilter_hxx */
