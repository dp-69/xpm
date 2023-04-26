#pragma once


#include "dpl/static_vector.hpp"

#include <boost/iostreams/device/mapped_file.hpp>

#include <concepts>
#include <filesystem>
// #include <fstream>
#include <vector>
#include <cstring>

#include <stdint.h>



namespace xpm
{



  
  using pnm_idx = int64_t;
  
  class pore_network_model
  {
    static void skip_until(char* &ptr, char val) {   
      while (*ptr != '\0' && *ptr++ != val) {}    
    }

    static void skip_line(char* &ptr) {   
      skip_until(ptr, '\n');
    }
    
    static void skip_word(char* &ptr) {       
      while (isspace(*ptr))
        ++ptr;
      while (*ptr != '\0' && !isspace(*ptr))
        ++ptr;
    }
    
    static void parse_text(char* &ptr, int& val) {
      val = strtol(ptr, &ptr, 10);
    }
    
    static void parse_text(char* &ptr, long long& val) {
      val = strtoll(ptr, &ptr, 10);
    }

    static void parse_text(char* &ptr, double& val) {
      val = strtod(ptr, &ptr);
    }

    static void parse_text(char* &ptr, dpl::vector3d& val) {
      parse_text(ptr, val.x());
      parse_text(ptr, val.y());
      parse_text(ptr, val.z());
    }

    static void parse_bin(char* &ptr, auto& val) {
      std::memcpy(&val, ptr, sizeof val);
      ptr += sizeof val;
    }

    // template <typename T>
    // static void parse_bin(char* &ptr, std::vector<T>& vec, size_t count) {
    //   vec.resize(count);
    //   auto total_size = sizeof(T)*count;
    //   std::memcpy(vec.data(), ptr, total_size);
    //   ptr += total_size;
    // }

    template <typename T>
    static void parse_bin(char* &ptr, T* ptr_out, size_t count) {
      auto total_size = sizeof(T)*count;
      std::memcpy(ptr_out, ptr, total_size);
      ptr += total_size;
    }
    

    void reset() {      // TODO: Not good function
      throats.resize(throat_count_);
      lengthThroat.resize(throat_count_);
      _length0.resize(throat_count_, 0);
      _length1.resize(throat_count_, 0);


      auto inlet_total_idx = throat_count_ + node_count_ - 2;
      auto outlet_total_idx = inlet_total_idx + 1;
      r_ins_.resize(throat_count_ + node_count_);            
      r_ins_[inlet_total_idx] = 1e30;
      r_ins_[outlet_total_idx] = 1e-30;

      volume.resize(throat_count_ + node_count_);
      volume[inlet_total_idx] = 0.0;
      volume[outlet_total_idx] = 0.0;

      shapeFactor.resize(throat_count_ + node_count_);
      shapeFactor[inlet_total_idx] = -9999.0;
      shapeFactor[outlet_total_idx] = -9999.0;
    }
    
    
    pnm_idx parse_statoil_text_idx(std::integral auto idx) const {
      if (idx == -1)
        return node_count_ - 2;
  
      if (idx == 0)
        return node_count_ - 1;
    
      return idx - 1;
    }

    pnm_idx node_count_ = 2; // includes inlet and outlet nodes
    pnm_idx throat_count_ = 0;
    
  public:
    dpl::vector3d physical_size{1};
    
    std::vector<std::pair<pnm_idx, pnm_idx>> throats;        

    std::vector<dpl::vector3d> node_pos;

    
    // Throat values go first, then node values.      
    std::vector<double> r_ins_; 
    std::vector<double> volume;
    std::vector<double> shapeFactor;
    std::vector<double> lengthThroat;
    std::vector<double> _length0;
    std::vector<double> _length1;

    enum class file_format
    {
      statoil,
      binary_dp69
    };

    pore_network_model(const pore_network_model& other) = default;
    pore_network_model(pore_network_model&& other) noexcept = default;
    pore_network_model& operator=(const pore_network_model& other) = default;
    pore_network_model& operator=(pore_network_model&& other) noexcept = default;

    pore_network_model(const std::filesystem::path& p, file_format ff) {
      if (ff == file_format::statoil)
        read_from_text_file(p);
      else if (ff == file_format::binary_dp69)
        read_from_binary_file(p);
    }


    
    
    auto inner_node_count() const {
      return node_count_ - 2;
    }

    auto throat_count() const {
      return throat_count_;
    }

    auto inner_count() const {
      return inner_node_count() + throat_count();
    }

    void scale(double x) {
      for (auto& p : node_pos)
        p = p*x;

      for (pnm_idx i = 0, count = inner_count(); i < count; ++i) {
        r_ins_[i] *= x;
      }
    }
    
    /**
     * \brief inlet and outlet are outer nodes, other nodes are inner.
     */
    bool inner_node(pnm_idx i) const {
      return i < inner_node_count();
    }

    auto node_r_ins(pnm_idx i) const {
      return r_ins_[throat_count_ + i];
    }

    void read_from_binary_file(const std::filesystem::path& network_path) {
      boost::iostreams::mapped_file_source file(network_path.string());
      auto* ptr = const_cast<char*>(file.data());


      
      // parse(ptr, k);
      
      parse_bin(ptr, node_count_);
      parse_bin(ptr, throat_count_);
      
      node_pos.resize(node_count_);
      r_ins_.resize(node_count_ + throat_count_);
      parse_bin(ptr, node_pos.data(), node_count_);
      parse_bin(ptr, r_ins_.data() + throat_count_, node_count_);
      
      for (pnm_idx i = 0; i < throat_count_; ++i) {
        dpl::vector2i pair;
        parse_bin(ptr, pair);
        throats.emplace_back(pair.x(), pair.y());
      }

      parse_bin(ptr, r_ins_.data(), throat_count_);

      for (auto& r : r_ins_)
        r = r/2;

      // for (pnm_idx i = 0; i < throat_count_ + node_count_; ++i)
      //   r_ins_[i] /= 2;

      node_count_ += 2;


      reset();
    }
    
    void read_from_text_file(const std::filesystem::path& network_path) {
      auto network_path_str = network_path.string();
      
      {
        boost::iostreams::mapped_file_source node1_file(network_path_str + "_node1.dat");
        auto* node1_ptr = const_cast<char*>(node1_file.data());
        parse_text(node1_ptr, node_count_);  
        parse_text(node1_ptr, physical_size);

        for (pnm_idx i = 0; i < node_count_; ++i) {
          dpl::vector3d p;
          skip_word(node1_ptr);
          parse_text(node1_ptr, p);
          skip_line(node1_ptr);
          node_pos.push_back(p);
        }


        node_count_ += 2;
      }

      {
        boost::iostreams::mapped_file_source throat1_file(network_path_str + "_link1.dat");
        auto* throat1_ptr = const_cast<char*>(throat1_file.data());
        parse_text(throat1_ptr, throat_count_);   

        reset();

        int param_int;
        
        for (pnm_idx i = 0; i < throat_count_; ++i) {       
          skip_word(throat1_ptr);

          parse_text(throat1_ptr, param_int);
          throats[i].first = parse_statoil_text_idx(param_int); 

          parse_text(throat1_ptr, param_int);
          throats[i].second = parse_statoil_text_idx(param_int);        

          parse_text(throat1_ptr, r_ins_[i]);
          parse_text(throat1_ptr, shapeFactor[i]);

          skip_line(throat1_ptr);
        }        
      }      

      {
        boost::iostreams::mapped_file_source throat2_file(network_path_str + "_link2.dat");
        auto* throat2_ptr = const_cast<char*>(throat2_file.data());

        for (pnm_idx i = 0; i < throat_count_; ++i) {
          dpl::sfor<3>([&throat2_ptr] {
            skip_word(throat2_ptr);  
          });
                    
          parse_text(throat2_ptr, _length0[i]);
          parse_text(throat2_ptr, _length1[i]);
          parse_text(throat2_ptr, lengthThroat[i]);
          parse_text(throat2_ptr, volume[i]);
          skip_line(throat2_ptr);
        }  
      }

      {
        boost::iostreams::mapped_file_source node2_file(network_path_str + "_node2.dat");
        auto* node2_ptr = const_cast<char*>(node2_file.data());        

        auto r_ins_node = r_ins_.begin() + throat_count_;
        auto volume_node = volume.begin() + throat_count_;
        auto shape_factor_node = shapeFactor.begin() + throat_count_;

        for (pnm_idx i = 0; i < node_count_ - 2; i++) {          
          skip_word(node2_ptr);
          parse_text(node2_ptr, volume_node[i]);
          parse_text(node2_ptr, r_ins_node[i]);
          parse_text(node2_ptr, shape_factor_node[i]);
          skip_line(node2_ptr);
        }        
      }        
    }
  };
}
