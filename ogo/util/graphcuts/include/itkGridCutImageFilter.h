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

#ifndef itkGridCutImageFilter_h
#define itkGridCutImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "GridGraph_3D_6C_MT.h"
#include <vector>

namespace itk {
/** \class GridCutImageFilter
 * \brief Abstract class for performing multi label graph cut using grid cut
 * 
 * \author: Bryce Besler
 * \ingroup BoneEnhancement
 */
template< typename TInputImage, typename TMaskImage, typename TOutputImage >
class ITK_TEMPLATE_EXPORT GridCutImageFilter
  : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(GridCutImageFilter);

  /** Standard Self typedef */
  using Self          = GridCutImageFilter;
  using Superclass    = ImageToImageFilter< TInputImage, TOutputImage >;
  using Pointer       = SmartPointer< Self >;
  using ConstPointer  = SmartPointer< const Self >;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GridCutImageFilter, ImageToImageFilter);

  /** Input Image typedefs. */
  using InputImageType          = TInputImage;
  using InputImagePointer       = typename InputImageType::Pointer;
  using InputImageConstPointer  = typename InputImageType::ConstPointer;
  using InputImageRegionType    = typename InputImageType::RegionType;
  using InputPixelType          = typename InputImageType::PixelType;
  using SizeType                = typename InputImageType::SizeType;
  itkStaticConstMacro(ImageDimension, unsigned int,  TInputImage::ImageDimension);

  /** Mask image typedefs. */
  using MaskImageType         = TMaskImage;
  using MaskImageConstPointer = typename MaskImageType::ConstPointer;
  using MaskPixelType         = typename MaskImageType::PixelType;
  using MaskImageRegionType   = typename MaskImageType::RegionType;

  /** Output image typedefs. */
  using OutputImageType       = TOutputImage;
  using OutputImagePointer    = typename OutputImageType::Pointer;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using OutputImagePixelType  = typename OutputImageType::PixelType;

  /** Grid cut definitions */
  using LabelType     = int;
  using CostType      = int;
  using EnergyType    = typename NumericTraits< InputPixelType >::RealType;
  using Grid          = GridGraph_3D_6C_MT< CostType, CostType, EnergyType >;
  using VectorType    = std::vector< CostType >;
  using NLinkType     = std::vector< VectorType >;
  using TLinkType     = std::vector< VectorType >;
  using RealType      = typename NumericTraits< InputPixelType >::RealType;
  using DistanceType  = typename NumericTraits< InputPixelType >::RealType;

  /** Iterator types */
  using ShapedIteratorType          = ConstShapedNeighborhoodIterator< InputImageType >;
  using MaskIterator                = ConstShapedNeighborhoodIterator< MaskImageType >;
  using OutputIteratorType          = ImageRegionIteratorWithIndex< OutputImageType >;
  using IndexType                   = typename OutputIteratorType::IndexType;
  using OffsetType                  = typename OutputIteratorType::OffsetType;
  using NeighboursType              = std::vector< OffsetType >; 

  /** Methods to set/get the mask image */
  itkSetInputMacro(Mask, MaskImageType);
  itkGetInputMacro(Mask, MaskImageType);

  /** Set/Get macros for BlockSize */
  itkSetMacro(BlockSize, LabelType);
  itkGetConstMacro(BlockSize, LabelType);

  /** Get max flow */
  itkGetConstMacro(MaxFlow, EnergyType);

  /** Set/Get macros for WeightScale */
  itkSetMacro(WeightScale, DistanceType);
  itkGetConstMacro(WeightScale, DistanceType);

  /** Get number of nNeighbours */
  itkGetConstMacro(nNeighbours, LabelType);

  /** Get number of Neighbors */
  itkGetConstMacro(Neighbors, NeighboursType);

  /** Get number of nLabels */
  itkGetConstMacro(nLabels, LabelType);

protected:
  GridCutImageFilter();
  virtual ~GridCutImageFilter() {}

  void PrintSelf(std::ostream & os, Indent indent) const override;

  /** Functions to be overwritten by inheritance */
  virtual CostType ComputeDataTerm(const InputPixelType p, const LabelType l, const MaskPixelType m);
  virtual CostType ComputeSmoothnessTerm(const InputPixelType p, const InputPixelType q, const DistanceType d, const MaskPixelType m_p, const MaskPixelType m_q);

  /* Accessing grid cut data */
  virtual void SetupNeighbourhood();
  unsigned long GetIndex(const IndexType p);
  virtual OutputImagePixelType GetLabel(const LabelType l) const;

  /* Helper functions */
  void SetDataTerm(const IndexType p, const LabelType l, const CostType cost);
  void SetSmoothTerm(const IndexType p, const LabelType n_i, const CostType cost);

  /** Multi-threading. */
  void BeforeThreadedGenerateData() override;
  void AfterThreadedGenerateData() override;
  void DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread) override;

private:
  /** Grid cut terms */
  NLinkType             m_nLinks;
  TLinkType             m_tLinks;
  LabelType             m_nLabels;
  LabelType             m_nNeighbours;
  LabelType             m_nVoxels;
  NeighboursType        m_Neighbors;
  LabelType             m_BlockSize;
  EnergyType            m_MaxFlow;
  std::unique_ptr<Grid> m_Grid;
  DistanceType          m_WeightScale;
  SizeType              m_Dimensions;
}; // end class
} /* end namespace */

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGridCutImageFilter.hxx"
#endif

#endif /* itkGridCutImageFilter_h */
