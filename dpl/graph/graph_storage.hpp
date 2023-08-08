#pragma once
#include "HW/dynamic_connectivity/dynamic_connectivity_graph.hpp"

namespace HW { namespace dynamic_connectivity
{
  class graph_storage : non_copyable_movable
  {
    size_t _vertexCount;
    size_t _edgeCount;

    unique_ptr<vertex[]> _vertexBuffer;
    unique_ptr<directed_edge[]> _directedEdgeBuffer;
    
    vertex_ptr _v0Ptr;
    directed_edge_ptr _de0Ptr;

  public:    
        
    vertex_ptr inletVertex;
    vertex_ptr outletVertex;    

    size_t total_idx(const vertex_ptr v) const {
      return _edgeCount + (v - _v0Ptr);
    }

    size_t vertex_idx(const vertex_ptr v) const {
      return v - _v0Ptr;      
    }

    size_t edge_idx(const directed_edge_ptr de) const {
      return (de - _de0Ptr)/2;      
    }

    size_t total_idx(const directed_edge_ptr de) const {
      return edge_idx(de);      
    }
        
    vertex_ptr vertex_by_total_idx(size_t totalIdx) const {
      return _v0Ptr + (totalIdx - _edgeCount);
    }

    vertex_ptr vertex_by_vertex_idx(size_t vertexIdx) const {
      return _v0Ptr + vertexIdx;
    }

    directed_edge_ptr directed_edge_by_edge_idx(size_t edgeIdx) const {
      return _de0Ptr + 2*edgeIdx;
    }

    directed_edge_ptr directed_edge_by_total_idx(size_t totalIdx) const {
      return directed_edge_by_edge_idx(totalIdx);
    }
    
    

    dynamic_connectivity_graph g;
    

//    graph_storage(size_t vertexCount, size_t edgeCount, const pair<size_t, size_t>* edges)
//      : _vertexCount(vertexCount),
//        _edgeCount(edgeCount),
//        _vertexBuffer(make_unique<vertex[]>(vertexCount)),
//        _directedEdgeBuffer(make_unique<directed_edge[]>(2*edgeCount)),
//        _v0Ptr(_vertexBuffer.get()),
//        _de0Ptr(_directedEdgeBuffer.get()),
//        inletVertex(_v0Ptr + vertexCount - 2),
//        outletVertex(inletVertex + 1) 
//    {            
//
//      for (auto i = 0; i < vertexCount; i++) {
//        auto v = &_vertexBuffer[i];
//        add_vertex(v, g);        
//      }
//
//      for (auto i = 0; i < edgeCount; i++) {
//        auto ab = _directedEdgeBuffer.get() + 2 * i;
//        add_edge(_vertexBuffer[edges[i].first], _vertexBuffer[edges[i].second], *ab, *(ab + 1), g);
//      }
//    }

    void init_topology(size_t vertexCount, size_t edgeCount, const pair<size_t, size_t>* edges) {
      _vertexCount = vertexCount;
      _edgeCount = edgeCount;
      _vertexBuffer = make_unique<vertex[]>(vertexCount);
      _directedEdgeBuffer = make_unique<directed_edge[]>(2*edgeCount),
      _v0Ptr = _vertexBuffer.get(),
      _de0Ptr =_directedEdgeBuffer.get();
      inletVertex = _v0Ptr + vertexCount - 2;
      outletVertex = inletVertex + 1;             

      for (size_t i = 0; i < vertexCount; i++) {
        auto v = &_vertexBuffer[i];
        add_vertex(v, g);    
      }

      for (size_t i = 0; i < edgeCount; i++) {
        auto ab = _directedEdgeBuffer.get() + 2 * i;

        add_edge(_vertexBuffer[edges[i].first], _vertexBuffer[edges[i].second], *ab, *(ab + 1), g);
      }
    }

//    void clear_et_entries() const {
//      for (size_t i = 0; i < _vertexCount; ++i)
//        _v0Ptr[i].et_entry_ = nullptr;
//
//      for (size_t i = 0; i < 2*_edgeCount; ++i)
//        directed_edge::set_null_et_entry(_de0Ptr + i);        
//    }
  };  
}}
