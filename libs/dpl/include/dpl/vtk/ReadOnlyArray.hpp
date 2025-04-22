/*
 * This file is part of Dmytro Petrovskyy Library (dpl).
 *
 * Copyright (c) 2023
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * dpl is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * dpl is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with dpl. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <vtkAOSDataArrayTemplate.h>
#include <vtkGenericDataArray.h>


namespace dpl::vtk
{
  template <typename Functor>
  class ReadOnlyArray : public vtkGenericDataArray<ReadOnlyArray<Functor>, std::invoke_result_t<Functor, vtkIdType>>
  {
    using ValueType = std::invoke_result_t<Functor, vtkIdType>;
    using BaseType = vtkGenericDataArray<ReadOnlyArray, ValueType>;
    friend class vtkGenericDataArray<ReadOnlyArray, ValueType>;
    
    std::unique_ptr<Functor> func_;
    
  protected:
    bool AllocateTuples(vtkIdType) { return true; }
    bool ReallocateTuples(vtkIdType) { return true; }
    
    vtkObjectBase* NewInstanceInternal() const override { return vtkAOSDataArrayTemplate<ValueType>::New(); }

  public:
    template <typename... Args>
    void ConstructFunctor(Args&&... args) {
      func_ = std::make_unique<Functor>(std::forward<Args>(args)...);
    }
    
    auto& GetFunctor() { return *func_; }
    
    static auto* New() { VTK_STANDARD_NEW_BODY(ReadOnlyArray) }
    vtkAbstractTemplateTypeMacro(ReadOnlyArray, BaseType)
    
    ValueType GetValue(vtkIdType) const { throw std::logic_error("not implemented"); }
    void SetValue(vtkIdType, ValueType) { throw std::logic_error("not implemented"); }
    void GetTypedTuple(vtkIdType, ValueType*) const { throw std::logic_error("not implemented"); }
    void SetTypedTuple(vtkIdType, const ValueType*) { throw std::logic_error("not implemented"); }
    void SetTypedComponent(vtkIdType, int, ValueType) { throw std::logic_error("not implemented"); }
    
    ValueType GetTypedComponent(vtkIdType tupleIdx, int) const {
      return (*func_)(tupleIdx);
    }
  };

  template <typename Func>
  auto CreateReadOnlyArray(Func f) {
    auto result = vtkSmartPointer<ReadOnlyArray<Func>>::New();
    result->ConstructFunctor(f);
    return result;
  }
}
