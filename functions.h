#pragma once

#include "xpm/pore_network_model.hpp"

#include <vtkActor.h>
#include <vtkLookupTable.h>

namespace xpm
{
  using v3i = dpl::vector3i;
  using v3d = dpl::vector3d;
  
  vtkSmartPointer<vtkActor> CreateNodeActor(const pore_network_model& pnm, vtkLookupTable* lut);
  vtkSmartPointer<vtkActor> CreateThroatActor(const pore_network_model& pnm, vtkLookupTable* lut);
}
