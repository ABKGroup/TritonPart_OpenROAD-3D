#############################################################################
##
## Copyright (c) 2019, OpenROAD
## All rights reserved.
##
## BSD 3-Clause License
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
## * Redistributions of source code must retain the above copyright notice, this
##   list of conditions and the following disclaimer.
##
## * Redistributions in binary form must reproduce the above copyright notice,
##   this list of conditions and the following disclaimer in the documentation
##   and/or other materials provided with the distribution.
##
## * Neither the name of the copyright holder nor the names of its
##   contributors may be used to endorse or promote products derived from
##   this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
## ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
## LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
## INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.

include("openroad")

# The location used by etc/DependencyInstaller.sh

add_library(mlsa_lib
  src/rtl_bp.cpp
  src/hier_mlsa.cpp
  src/object.cpp
  src/MultilevelSACore.cpp
  src/MultilevelSACore_New.cpp
)

target_include_directories(mlsa_lib
  PUBLIC
    include
)

target_link_libraries(mlsa_lib
  PRIVATE
    utl_lib
    dbSta_lib
    OpenSTA
    odb
    dl
    par_lib
)

swig_lib(NAME      mlsa
         NAMESPACE mlsa
         I_FILE    src/mlsa.i
         SCRIPTS   src/mlsa.tcl
)

target_sources(mlsa
  PRIVATE
    src/MakeBlobPlacer.cpp
)

messages(
  TARGET mlsa
)

target_include_directories(mlsa
  PUBLIC
    include
    # src
)

target_link_libraries(mlsa
  PRIVATE
    mlsa_lib
    odb
    OpenSTA
    dbSta_lib
    par_lib
)
