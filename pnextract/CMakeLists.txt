cmake_minimum_required (VERSION 3.3.0)
project (pnextract_project)

add_executable(pnextract # WIN32
  "src/pnm/pnextract/nextract.cpp"  
  "src/pnm/pnextract/blockNet.cpp"  
  "src/pnm/pnextract/medialSurf.cpp"  
  "src/pnm/pnextract/writers_vxl.cpp"  
  "src/pnm/pnextract/writers_vtk.cpp"  
  "src/libvoxel/voxelImage.cpp"    
)

target_compile_features(pnextract PRIVATE cxx_std_17)

target_include_directories(pnextract PRIVATE "src/include")
target_include_directories(pnextract PRIVATE "src/libvoxel")