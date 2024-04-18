/*
 * This file is part of Dmytro Petrovskyy Library (DPL).
 *
 * Copyright (c) 2024
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * DPL is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * DPL is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with RRM. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include "units.hpp"

#ifdef _WIN32

#include <windows.h>
#undef max
#undef min

namespace dpl
{
  template <typename T = units::byte>
  T get_memory_consumption() {
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);

    GlobalMemoryStatusEx(&status);

    return units::byte{status.ullTotalPhys - status.ullAvailPhys};

    // _tprintf (TEXT("There is  %*ld percent of memory in use.\n"),WIDTH, mem_status.dwMemoryLoad);
    // _tprintf (TEXT("There are %*I64d total Gbytes of physical memory.\n"),WIDTH,mem_status.ullTotalPhys/DIV);
    // _tprintf (TEXT("There are %*I64d free GBytes of physical memory.\n"),WIDTH, mem_status.ullAvailPhys/DIV);
    // _tprintf (TEXT("There are %*I64d total GBytes of paging file.\n"),WIDTH, statex.ullTotalPageFile/DIV);
    // _tprintf (TEXT("There are %*I64d free GBytes of paging file.\n"),WIDTH, statex.ullAvailPageFile/DIV);
    // _tprintf (TEXT("There are %*I64d total GBytes of virtual memory.\n"),WIDTH, statex.ullTotalVirtual/DIV);
    // _tprintf (TEXT("There are %*I64d free GBytes of virtual memory.\n"),WIDTH, statex.ullAvailVirtual/DIV);
    // _tprintf (TEXT("There are %*I64d free GBytes of extended memory.\n"),WIDTH, statex.ullAvailExtendedVirtual/DIV);
  }
}

#endif