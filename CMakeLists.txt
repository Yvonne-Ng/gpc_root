
atlas_subdir( GPC )

## As a Fortran compiler may not have been set up yet when this file is included.
#if( ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" ) OR
#    ( "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" ) )
#   _add_flag( CMAKE_Fortran_FLAGS "-ffixed-line-length-132" )
#   _add_flag( CMAKE_Fortran_FLAGS "-DFVOIDP=INTEGER*${CMAKE_SIZEOF_VOID_P}" )
#   _add_flag( CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops" )
#   _add_flag( CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-funroll-all-loops" )
#endif()
#
# no include directory
#include_directories(${CMAKE_SOURCE_DIR}/stuff/include) # -I flags for compiler
find_package( ROOT COMPONENTS Tree Hist Gpad )

#===gp====
# adding gp.cpp as an executable
#atlas_add_executable( gp gp.cpp gp.h
#INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS}
#LINK_LIBRARIES  ${ROOT_LIBRARIES} ${BOOST_LIBRARIES} CClctrl CGp CMatrix ndlfortran CNoise ndlutil ndlstrutil CTransform COptimisable CKern CDist ndlassert )
##target_link_library( gp CClctrl CGp CMatrix ndlfortran CNoise ndlutil ndlstrutil CTransform COptimisable CKern CDist ndlassert )
#gp.o: gp.cpp gp.h ndlexceptions.h ndlstrutil.h CMatrix.h ndlassert.h \
#  CNdlInterfaces.h ndlutil.h ndlfortran.h lapack.h CKern.h CTransform.h \
#  CDataModel.h CDist.h CGp.h CMltools.h COptimisable.h CNoise.h CClctrl.h

# linking gp with its dependence
#======


#===ndlfortran===
atlas_add_library( ndlfortranlib ndlfortran.f )
#SET_TARGET_PROPERTIES(ndlfortran PROPERTIES LINKER_LANGUAGE FORTRAN)

##===ndlassert==
## adding a library ndlassert
#atlas_add_library(ndlassert ndlassert.cpp ndlassert.h ndlexceptions.h
#    PUBLIC_HEADER GPc)
## no linking of libs required
##===============
##=================

##===ndlstrutil
#atlas_add_library(ndlstrutil ndlstrutil.cpp ndlstrutil.h ndlexceptions.h )
##=================
#
##===ndlutil===
#atlas_add_library(ndlutil ndlutil.cpp ndlutil.h ndlassert.h ndlexceptions.h ndlfortran.h ndlfortranlib )
##=================
#
##===CIvm===
#atlas_add_library(CIvm CIvm.cpp CIvm.h CMltools.h ndlassert.h ndlexceptions.h ndlstrutil.h COptimisable.h CMatrix.h CNdlInterfaces.h ndlutil.h ndlfortran.h lapack.h CKern.h CTransform.h CDataModel.h CDist.h CNoise.h ndlfortranlib )
##=================
#
##=====CGP======
#atlas_add_library(CGp CGp.cpp CGp.h CMltools.h ndlassert.h ndlexceptions.h ndlstrutil.h COptimisable.h CMatrix.h CNdlInterfaces.h ndlutil.h ndlfortran.h lapack.h CKern.h CTransform.h CDataModel.h CDist.h CNoise.h ndlfortranlib )
##=================
#
##=====CMatrix======
#atlas_add_library(CMatrix CMatrix.cpp CMatrix.h ndlassert.h ndlexceptions.h CNdlInterfaces.h ndlstrutil.h ndlutil.h ndlfortran.h lapack.h ndlfortranlib)
##=================
#
#
##=====CDist======
#atlas_add_library(CDist CDist.cpp CDist.h CMatrix.h ndlassert.h ndlexceptions.h CNdlInterfaces.h ndlstrutil.h ndlutil.h ndlfortran.h lapack.h CTransform.h ndlfortranlib )
##=================
#
#
##=====COptimisable===
#atlas_add_library(COptimisable COptimisable.cpp COptimisable.h CMatrix.h ndlassert.h ndlexceptions.h CNdlInterfaces.h ndlstrutil.h ndlutil.h ndlfortran.h lapack.h ndlfortranlib )
##=================
#
##=====CTransform-===
#atlas_add_library(CTransform CTransform.cpp CTransform.h CMatrix.h ndlassert.h ndlexceptions.h CNdlInterfaces.h ndlstrutil.h ndlutil.h ndlfortran.h lapack.h)
##=================
#
#
##=====CKern-===
#atlas_add_library(CKern CKern.cpp CKern.h ndlassert.h ndlexceptions.h CTransform.h CMatrix.h CNdlInterfaces.h ndlstrutil.h ndlutil.h ndlfortran.h lapack.h CDataModel.h CDist.h)
##=================
#
#
##=====CNoise-===
#atlas_add_library(CNoise CNoise.cpp CNoise.h ndlexceptions.h ndlutil.h ndlassert.h ndlfortran.h ndlstrutil.h CMatrix.h CNdlInterfaces.h lapack.h
#  CTransform.h COptimisable.h CDist.h CKern.h CDataModel.h)
##=================
#
#
#atlas_add_library(CClctrl CClctrl.cpp CClctrl.h ndlstrutil.h ndlexceptions.h ndlutil.h ndlassert.h ndlfortran.h CMatrix.h CNdlInterfaces.h lapack.h)
#
#atlas_add_library(CGplvm CGplvm.cpp CGplvm.h CMltools.h ndlassert.h ndlexceptions.h ndlstrutil.h COptimisable.h CMatrix.h CNdlInterfaces.h ndlutil.h ndlfortran.h lapack.h CKern.h CTransform.h CDataModel.h CDist.h CNoise.h)
#
##====
##suspect CClgplvm is not being used, skipping
## skipping test that uses matlab
