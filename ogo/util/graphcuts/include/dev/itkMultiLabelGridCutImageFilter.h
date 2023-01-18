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

#include "AlphaExpansion_3D_6C_MT.h"
#include <vector>

namespace itk {
/** \class GridCutImageFilter
 * \brief Abstract class for performing multi label graph cut using grid cut
 * 
 * \author: Bryce Besler
 * \ingroup BoneEnhancement
 */
template< typename TInputImage, typename TOutputImage >
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
  using InputImagePixelType     = typename InputImageType::PixelType;
  using PixelValueType          = typename InputImagePixelType::ValueType;
  using SizeType                = typename InputImageType::SizeType;
  itkStaticConstMacro(ImageDimension, unsigned int,  TInputImage::ImageDimension);

  /** Output image typedefs. */
  using OutputImageType       = TOutputImage;
  using OutputImagePointer    = typename OutputImageType::Pointer;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using OutputImagePixelType  = typename OutputImageType::PixelType;

  /** Grid cut definitions */
  using LabelType       = int;
  using CostType        = int;
  using EnergyType      = typename NumericTraits< PixelValueType >::RealType;
  using Grid            = AlphaExpansion_3D_6C_MT< LabelType, CostType, EnergyType >;
  using DataCostType    = std::vector< CostType >;
  using SmoothCostType  = std::vector< DataCostType >;

  /** Iterator types */
  using InputIteratorWithIndexType  = typename ImageRegionIteratorWithIndex< InputImageType >;
  using InputIteratorType           = typename ImageRegionIterator< InputImageType >;
  using OutputIteratorType          = typename ImageRegionIterator< OutputImageType >;
  using IndexType                   = typename InputIteratorWithIndexType::IndexType;
  using OffsetType                  = typename InputIteratorWithIndexType::OffsetType;
  using NeighboursType              = std::vector< OffsetType >; 

  /** Set/Get macros for BlockSize */
  itkSetMacro(BlockSize, LabelType);
  itkGetConstMacro(BlockSize, LabelType);

protected:
  GridCutImageFilter() {};
  virtual ~GridCutImageFilter() {}

  void PrintSelf(std::ostream & os, Indent indent) const override;

  /** Functions to be overwritten by inheritance */
  virtual OutputImagePixelType ComputeDataTerm(const IndexType p, const LabelType l) = 0;
  virtual OutputImagePixelType ComputeSmoothnessTerm(const IndexType p, const IndexType q, const LabelType l_p, const LabelType l_q) = 0;

  /* Accessing grid cut data */
  virtual void SetupNeighbourhood();
  virtual void SetupLabels();

  /* Helper functions */
  void SetDataTerm(const LabelType p_i, const LabelType l, const EnergyType cost);
  void SetSmoothTerm(const LabelType p_i, const LabelType q_i, const LabelType l_p, const LabelType l_q, const EnergyType cost);

  /** Multi-threading. */
  void BeforeThreadedGenerateData() override;
  void AfterThreadedGenerateData() override;
  void DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread) override;

private:
  /** Grid cut terms */
  DataCostType          m_Data;
  SmoothCostType        m_Smooth;
  LabelType             m_nLabels;
  LabelType             m_nNeighbours;
  LabelType             m_nVoxels;
  OutputImagePixelType  m_BackgroundLabel;
  NeighboursType        m_Neighbors;
  LabelType             m_BlockSize;
}; // end class
} /* end namespace */

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGridCutImageFilter.hxx"
#endif

#endif /* itkGridCutImageFilter_h */
