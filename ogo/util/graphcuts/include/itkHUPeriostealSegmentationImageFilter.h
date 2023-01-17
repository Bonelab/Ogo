
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

#ifndef itkHUPeriostealSegmentationImageFilter_h
#define itkHUPeriostealSegmentationImageFilter_h

#include "itkGridCutImageFilter.h"

namespace itk {
/** \class HUPeriostealSegmentationImageFilter
 * \brief Perform femur segmentation
 * 
 * \author: Bryce Besler
 * \ingroup BoneEnhancement
 */
template< typename TInputImage, typename TMaskImage, typename TOutputImage >
class ITK_TEMPLATE_EXPORT HUPeriostealSegmentationImageFilter
  : public GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(HUPeriostealSegmentationImageFilter);

  /** Standard Self typedef */
  using Self          = HUPeriostealSegmentationImageFilter;
  using Superclass    = GridCutImageFilter< TInputImage, TMaskImage, TOutputImage >;
  using Pointer       = SmartPointer< Self >;
  using ConstPointer  = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(HUPeriostealSegmentationImageFilter, GridCutImageFilter);

  /** Grid cut definitions */
  using InputPixelType  = typename Superclass::InputPixelType;
  using MaskPixelType   = typename Superclass::MaskPixelType;
  using CostType        = typename Superclass::CostType;
  using LabelType       = typename Superclass::LabelType;
  using DistanceType    = typename Superclass::DistanceType;
  using RealType        = typename Superclass::RealType;


  /** Set/Get macros for BackgroundLabel */
  itkSetMacro(BackgroundLabel, MaskPixelType);
  itkGetConstMacro(BackgroundLabel, MaskPixelType);

  /** Set/Get macros for ForegroundLabel */
  itkSetMacro(ForegroundLabel, MaskPixelType);
  itkGetConstMacro(ForegroundLabel, MaskPixelType);

  /** Set/Get macros for Lambda */
  itkSetMacro(Lambda, DistanceType);
  itkGetConstMacro(Lambda, DistanceType);

  /** Set/Get macros for Sigma */
  itkSetMacro(Sigma, DistanceType);
  itkGetConstMacro(Sigma, DistanceType);

protected:
  HUPeriostealSegmentationImageFilter();
  virtual ~HUPeriostealSegmentationImageFilter() {}

  void PrintSelf(std::ostream & os, Indent indent) const override;

  /** Functions to be overwritten by inheritance */
  CostType ComputeDataTerm(const InputPixelType p, const LabelType l, const MaskPixelType m) override;
  CostType ComputeSmoothnessTerm(const InputPixelType p, const InputPixelType q, const DistanceType d, const MaskPixelType m_p, const MaskPixelType m_q) override;

private:
  /** Grid cut terms */
  MaskPixelType   m_BackgroundLabel;
  MaskPixelType   m_ForegroundLabel;
  RealType        m_Lambda;
  RealType        m_Sigma;
}; // end class
} /* end namespace */

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHUPeriostealSegmentationImageFilter.hxx"
#endif

#endif /* itkHUPeriostealSegmentationImageFilter_h */
