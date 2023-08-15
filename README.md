Extensive Pore Modelling (xpm)
===

## Downloads

The latest release can be found at

| Platform                | Files          |
|-------------------------|----------------|
| Windows ZIP (portable)  | *coming soon*  |




To run xpm, unpack a zip archive and run the
executable `xpm.exe`.


## Build instructions





#### Windows

All commands are issued from a **x64 Native Tools Command Prompt for VS 2022** command prompt (which is a part of
Microsoft's Visual Studio toolset). It is assumed that Step 1 is conducted in `C:\` drive. Note that Step 1 may take at least an hour to execute and occupy around 100 GB of free storage.

- Installing MPI

Download installer from https://www.microsoft.com/en-us/download/details.aspx?id=57467 and run it.

- Step 1 - Installing xpm's dependencies in vcpkg
```cmd
> git clone https://github.com/microsoft/vcpkg
> cd vcpkg
> bootstrap-vcpkg.bat
> vcpkg.exe install vtk[qt]:x64-windows hypre:x64-windows boost-interprocess:x64-windows boost-iostreams:x64-windows boost-graph:x64-windows fmt:x64-windows
```

- Step 2 - xpm compilation
```cmd
> git clone https://github.com/dp-69/xpm.git
> cd xpm
> cmake --preset=win-rel
> cmake --build --preset=win-rel
```

- Step 3 - Executing xpm
```cmd
> cd build/Release/bin
> xpm.exe
```

#### Ubuntu

The following instructions assume that vcpkg's repository will be cloned into the user's home folder (`cd ~`), and that all commands are issued from a terminal. Note that Step 1 may take at least an hour to execute and occupy around 100 GB of free storage.

- Required Ubuntu packages (install with `sudo apt install`)

```cmd
build-essential
git
tar curl zip unzip                                                                                                
pkg-config
meson
libxi-dev libgl1-mesa-dev libglu1-mesa-dev mesa-common-dev libxrandr-dev libxxf86vm-dev
gfortran
autoconf autoconf-archive
libtool
'^libxcb.*-dev' libx11-xcb-dev libglu1-mesa-dev libxrender-dev libxi-dev libxkbcommon-dev libxkbcommon-x11-dev
cmake 
openmpi-bin
```

- Step 1 - Installing xpm's dependencies in vcpkg
```cmd
> git clone https://github.com/microsoft/vcpkg
> cd vcpkg
> ./bootstrap-vcpkg.sh
> ./vcpkg install vtk[qt] hypre boost-interprocess boost-iostreams boost-graph fmt
```

- Step 2 - xpm compilation
```cmd
> git clone https://github.com/dp-69/xpm.git
> cd xpm
> git checkout ubuntu-compilation
> cmake --preset=lin-rel-cfg
> cmake --build --preset=lin-rel-build
```

- Step 3 - Executing xpm
```cmd
> cd build/Release/bin
> ./xpm
```

## Copyright

The xpm project's research results are a copyright of the respective members and associates.

### Heriot-Watt University, UK
[Dmytro Petrovskyy](https://www.linkedin.com/in/dmytro-petrovskyy/)<br/>
Julien Maes<br/>
Hannah Menke<br/>
Kamaljit Singh<br/>

### Ghent University, Belgium
Tom Bultreys<br/>




## License

xpm is Free Software released under the 
[GPLv3](https://www.gnu.org/licenses/gpl.html) license.  It includes 
[third party libraries](https://bitbucket.org/rapidreservoirmodelling/rrm2/src/main/rrm_3rd_party_libraries.md)
with other licenses, which we understand to be compatible. An overview of these libraries including
brief descriptions, the last known URL the code was made available at, the date
it was retrieved, the license(s) in which the code was made available by their
authors at the retrieving date, and URL links to publicly available copies of
such licenses is available
[here](https://bitbucket.org/rapidreservoirmodelling/rrm2/src/main/rrm_3rd_party_libraries.md).
Please, check the source code for details on the licenses that apply to each
piece of software.


