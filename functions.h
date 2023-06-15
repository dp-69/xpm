#pragma once

#include "xpm/pore_network_model.hpp"

#include <vtkActor.h>
#include <vtkLookupTable.h>
#include <vtkCylinderSource.h>
#include <vtkFloatArray.h>
#include <vtkGlyph3DMapper.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>

namespace xpm
{
  using v3i = dpl::vector3i;
  using v3d = dpl::vector3d;

  

  




  // vtkSmartPointer<vtkActor> CreateNodeActor(const pore_network_model& pnm, vtkLookupTable* lut);
  // vtkSmartPointer<vtkActor> CreateThroatActor(const pore_network_model& pnm, vtkLookupTable* lut);

  inline dpl::vector3d angles_for_j_norm(const dpl::vector3d& dir) {
    auto norm = dir.normalise();
    
    return {
      asin(norm.z())/3.141592654*180.0, 
      0.0f,                                                  
      (atan2(-norm.x(), norm.y())/3.141592654*180.0)      
    };
  }
  
  vtkSmartPointer<vtkActor> CreateNodeActor(const pore_network_model& pnm, vtkLookupTable* lut, const auto& color_map) {
    using namespace attribs;

    
    vtkNew<vtkPolyData> polydata;
      
    vtkNew<vtkSphereSource> cylinder;
    cylinder->SetPhiResolution(10);
    cylinder->SetThetaResolution(10);
    cylinder->SetCenter(0, 0, 0);
    cylinder->SetRadius(1);
      
    vtkNew<vtkGlyph3DMapper> node_glyphs;
    node_glyphs->SetSourceConnection(cylinder->GetOutputPort());
    node_glyphs->SetInputData(polydata);
    node_glyphs->OrientOff();

    vtkNew<vtkFloatArray> scale_array;
    scale_array->SetName("scale");
    scale_array->SetNumberOfComponents(1);

    polydata->GetPointData()->AddArray(scale_array);
    node_glyphs->SetScaleArray(scale_array->GetName());
    node_glyphs->SetScaleModeToScaleByMagnitude();



    vtkNew<vtkFloatArray> color_array;
    color_array->SetName("color");
    color_array->SetNumberOfComponents(1);
    polydata->GetPointData()->SetScalars(color_array);

    node_glyphs->SetLookupTable(lut);
    node_glyphs->SetColorModeToMapScalars();
    node_glyphs->UseLookupTableScalarRangeOn();
    node_glyphs->SetScalarModeToUsePointData();
      
    vtkNew<vtkPoints> points;
      
    for (pnm_idx i = 0, count = pnm.node_count_; i < count; ++i) {
      points->InsertNextPoint(pnm.node_[pos][i]);
      scale_array->InsertNextTuple1(pnm.node_[r_ins][i]);
      color_array->InsertNextTuple1(color_map(i));
    }
      
    polydata->SetPoints(points);
      
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(node_glyphs);

    actor->GetProperty()->SetSpecularColor(1, 1, 1);
    actor->GetProperty()->SetSpecular(0);
    actor->GetProperty()->SetSpecularPower(75);
    actor->GetProperty()->SetAmbient(0.15);
    actor->GetProperty()->SetDiffuse(0.9);
      
    vtkNew<vtkNamedColors> colors;
    actor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());
    return actor;
  }

  vtkSmartPointer<vtkActor> CreateThroatActor(const pore_network_model& pnm, vtkLookupTable* lut, const auto& color_map) {
    using namespace attribs;

    vtkNew<vtkPolyData> polydata;
      
    vtkNew<vtkCylinderSource> cylinder;
    cylinder->SetResolution(20);
    cylinder->SetCenter(0.0, 0.5, 0.0);
    cylinder->SetRadius(1);
    cylinder->SetHeight(1);
      
    vtkNew<vtkGlyph3DMapper> throat_glyphs;
    throat_glyphs->SetSourceConnection(cylinder->GetOutputPort());
    throat_glyphs->SetInputData(polydata);
    
                     

    vtkNew<vtkFloatArray> orient_array;
    orient_array->SetName("orient");
    orient_array->SetNumberOfComponents(3);

    polydata->GetPointData()->AddArray(orient_array);
    throat_glyphs->SetOrientationArray(orient_array->GetName());
    throat_glyphs->SetOrientationModeToRotation();

    vtkNew<vtkFloatArray> scale_array;
    scale_array->SetName("scale");
    scale_array->SetNumberOfComponents(3);

    polydata->GetPointData()->AddArray(scale_array);
    throat_glyphs->SetScaleArray(scale_array->GetName());
    throat_glyphs->SetScaleModeToScaleByVectorComponents();



      
      

    vtkNew<vtkFloatArray> color_array;
    color_array->SetName("color");
    color_array->SetNumberOfComponents(1);
    polydata->GetPointData()->SetScalars(color_array);

    throat_glyphs->SetLookupTable(lut);
    throat_glyphs->SetColorModeToMapScalars();
    throat_glyphs->UseLookupTableScalarRangeOn();
    throat_glyphs->SetScalarModeToUsePointData();

    vtkNew<vtkPoints> points;

    for (pnm_idx i = 0, count = pnm.throat_count_; i < count; ++i)
      if (auto [n0, n1] = pnm.throat_[adj][i];
        pnm.inner_node(n0) && pnm.inner_node(n1)) {
        auto& n0_pos = pnm.node_[pos][n0];

        points->InsertNextPoint(n0_pos);
        orient_array->InsertNextTuple(angles_for_j_norm(pnm.node_[pos][n1] - n0_pos));

        scale_array->InsertNextTuple(v3d{
          pnm.throat_[r_ins][i],
          (pnm.node_[pos][n1] - n0_pos).length(),
          pnm.throat_[r_ins][i]
        });

        color_array->InsertNextTuple1(color_map(i));
      }

    polydata->SetPoints(points);
      
    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(throat_glyphs);
      

    actor->GetProperty()->SetSpecularColor(1, 1, 1);
    actor->GetProperty()->SetSpecular(0);
    actor->GetProperty()->SetSpecularPower(75);
    actor->GetProperty()->SetAmbient(0.15);
    actor->GetProperty()->SetDiffuse(0.9);

      
      
    vtkNew<vtkNamedColors> colors;
    actor->GetProperty()->SetColor(colors->GetColor3d("Salmon").GetData());
    return actor;        
  }
}
