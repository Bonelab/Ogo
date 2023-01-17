
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

#ifndef itkEndostealSegmentationImageFilter_h
#define itkEndostealSegmentationImageFilter_h

#include "itkGridCutImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkMaskImageFilter.h"

namespace itk {
/** \class EndostealSegmentationImageFilter
 * \brief Perform femur segmentation
 * 
 * \author: Bryce Besler
 * \ingroup BoneEnhancement
 */
template< typename TInputImage, typename TMaskImage, typename TOutputImage >
class ITK_TEMPLATE_EXPORT EndostealSegmentationImageFilter
  : public GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(EndostealSegmentationImageFilter);

  /** Standard Self typedef */
  using Self          = EndostealSegmentationImageFilter;
  using Superclass    = GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >;
  using Pointer       = SmartPointer< Self >;
  using ConstPointer  = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(EndostealSegmentationImageFilter, GridCutImageFilter);

  /** Grid cut definitions */
  using InputPixelType        = typename Superclass::InputPixelType;
  using MaskPixelType         = typename Superclass::MaskPixelType;
  using CostType              = typename Superclass::CostType;
  using LabelType             = typename Superclass::LabelType;
  using DistanceType          = typename Superclass::DistanceType;
  using RealType              = typename Superclass::RealType;

  /** Iterators */
  using OutputImagePixelType  = typename Superclass::OutputImagePixelType;
  using InputImageConstPointer = typename Superclass::InputImageConstPointer;
  using MaskImageConstPointer = typename Superclass::MaskImageConstPointer;
  using InputImageRegionType = typename Superclass::InputImageRegionType;
  using MaskImageRegionType = typename Superclass::MaskImageRegionType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;
  using ShapedIteratorType = typename Superclass::ShapedIteratorType;
  using OffsetType = typename Superclass::OffsetType;
  using MaskIterator = typename Superclass::MaskIterator;
  using InputImageType = typename Superclass::InputImageType;

  /** Signed distance filter */
  using MaskImageType           = typename Superclass::MaskImageType;
  using DistanceImageType       = Image< RealType, MaskImageType::ImageDimension >;
  using DistanceConstImageType  = typename DistanceImageType::ConstPointer; 
  using DistanceFilterType      = SignedMaurerDistanceMapImageFilter< MaskImageType, DistanceImageType >;
  using DistanceIteratorType    = ConstShapedNeighborhoodIterator< DistanceImageType >;
  using DistanceImageRegionType = typename DistanceImageType::RegionType;

  /** Masking filter */
  using MaskFilterType = MaskImageFilter< MaskImageType, MaskImageType, MaskImageType >;

  /** Set/Get macros for BackgroundLabel */
  itkSetMacro(BackgroundLabel, MaskPixelType);
  itkGetConstMacro(BackgroundLabel, MaskPixelType);

  /** Set/Get macros for CorticalLabel */
  itkSetMacro(CorticalLabel, MaskPixelType);
  itkGetConstMacro(CorticalLabel, MaskPixelType);

  /** Set/Get macros for CancellousLabel */
  itkSetMacro(CancellousLabel, MaskPixelType);
  itkGetConstMacro(CancellousLabel, MaskPixelType);

  /** Set/Get macros for Lambda */
  itkSetMacro(Lambda, RealType);
  itkGetConstMacro(Lambda, RealType);

  /** Set/Get macros for Sigma */
  itkSetMacro(Sigma, RealType);
  itkGetConstMacro(Sigma, RealType);

  /** Set/Get macros for MaxDistance */
  itkSetMacro(MaxDistance, DistanceType);
  itkGetConstMacro(MaxDistance, DistanceType);

  /** Set/Get macros for MaxDistance */
  itkSetMacro(MinDistance, DistanceType);
  itkGetConstMacro(MinDistance, DistanceType);

protected:
  EndostealSegmentationImageFilter();
  virtual ~EndostealSegmentationImageFilter() {}

  void PrintSelf(std::ostream & os, Indent indent) const override;

  /** Functions to be overwritten by inheritance */
  CostType ComputeDataTerm(const InputPixelType p, const RealType d, const LabelType l, const MaskPixelType m);
  CostType ComputeSmoothnessTerm(const InputPixelType p, const InputPixelType q, const DistanceType d, const MaskPixelType m_p, const MaskPixelType m_q) override;

  OutputImagePixelType GetLabel(const LabelType l) const override;

  /** Multi-threading. */
  void BeforeThreadedGenerateData() override;
  void AfterThreadedGenerateData() override;
  void DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread) override;

private:
  /** Grid cut terms */
  MaskPixelType               m_BackgroundLabel;
  MaskPixelType               m_CorticalLabel;
  MaskPixelType               m_CancellousLabel;
  RealType                    m_Lambda;
  RealType                    m_Sigma;
  RealType                    m_MaxDistance;
	RealType										m_MinDistance;

  typename DistanceFilterType::Pointer m_DistanceFilter;
  typename MaskFilterType::Pointer     m_MaskingFilter;
}; // end class
} /* end namespace */

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEndostealSegmentationImageFilter.hxx"
#endif

#endif /* itkEndostealSegmentationImageFilter_h */
