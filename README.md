# Extensive Pore Modelling (xpm)

xpm is a software for predicting flow properties of the multi-scale pore space.

xpm uses [pnextract](https://github.com/ImperialCollegeLondon/pnextract) to acquire a network model from an image.

<!-- 
## Downloads

The latest release can be found at

| Platform                | Files          |
|-------------------------|----------------|
| Windows ZIP (portable)  | *coming soon*  |

To run xpm, unpack a zip archive and run the executable `xpm.exe`. -->

## Publications

- [Poster](https://doi.org/10.6084/m9.figshare.25902862.v1) presented at Interpore 2024, Qingdao, China

## Build instructions

### Windows

Commands are issued from a **x64 Native Tools Command Prompt for VS 2022** command prompt (which is a part of Microsoft's Visual Studio toolset). It is assumed that Step 1 is conducted in `C:\` drive. Note that Step 1 may take at least an hour to execute.

- Installing MPI

Download installer from https://www.microsoft.com/en-us/download/details.aspx?id=57467 and run it.

- Step 1 - Installing xpm's dependencies in vcpkg
```cmd
> git clone https://github.com/microsoft/vcpkg
> cd vcpkg
> git clone -b xpm https://github.com/dp-69/vcpkg-ports ports-xpm
> bootstrap-vcpkg.bat
> vcpkg.exe install vtk[qt] qtcharts hypre boost-intrusive boost-iostreams boost-graph fmt argh --overlay-ports=ports-xpm --clean-after-build
```

- Step 2 - xpm compilation
```cmd
> git clone https://github.com/dp-69/xpm
> cd xpm
> cmake --preset=win-rel
> cmake --build --preset=win-rel
```

- Step 3 - Executing xpm
```cmd
> cd build/Release/bin
> xpm.exe
```

### Ubuntu

The following instructions assume that vcpkg's repository will be cloned into the user's home folder (`cd ~`), and that all commands are issued from a terminal. Note that Step 1 may take at least an hour to execute.

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
linux-libc-dev
python3-jinja2
openmpi-bin
```

- Step 1 - Installing xpm's dependencies in vcpkg
```cmd
> git clone https://github.com/microsoft/vcpkg
> cd vcpkg
> git clone -b xpm https://github.com/dp-69/vcpkg-ports ports-xpm
> ./bootstrap-vcpkg.sh
> ./vcpkg install vtk[qt] qtcharts hypre boost-intrusive boost-iostreams boost-graph fmt argh --overlay-ports=ports-xpm --clean-after-build
```

- Step 2 - xpm compilation
```cmd
> git clone https://github.com/dp-69/xpm
> cd xpm
> cmake --preset=lin-rel
> cmake --build --preset=lin-rel
```

- Step 3 - Executing xpm
```cmd
> cd build/Release/bin
> ./xpm
```

## Copyright

The xpm project's research results are a copyright of the respective members and associates.

[Dmytro Petrovskyy](https://www.linkedin.com/in/dmytro-petrovskyy/) <sup>1</sup><br/>
Julien Maes <sup>2</sup><br/>
Hannah Menke <sup>2</sup><br/>
Kamaljit Singh <sup>2</sup><br/>

<sup>1</sup> Independant Consultant, Ivano-Frankivsk, Ukraine<br/>
<sup>2</sup> Heriot-Watt University, Edinburgh, UK

*Tom Bultreys (Ghent University, Belgium) is thanked for their advisory in multi-scale pore-network modelling.*

## License

xpm is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

xpm is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with xpm. If not, see <http://www.gnu.org/licenses/>.