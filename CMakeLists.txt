atlas_subdir( GPC )

find_package( ROOT COMPONENTS Tree Hist Gpad TreePlayer RooStats RooFit RooFitCore TFile )

#===ndlfortran===
atlas_add_library( ndlfortranlib ndlfortran.f )


##===ndlassert==
## adding a library ndlassert
atlas_add_library( ndlassert ndlassert.cpp ndlassert.h ndlexceptions.h )
## no linking of libs required
##===============
##=================

##===ndlstrutil
atlas_add_library( ndlstrutil ndlstrutil.cpp ndlstrutil.h ndlexceptions.h )
##=================
#
##===ndlutil===
atlas_add_library( ndlutil ndlutil.cpp ndlutil.h ndlassert.h ndlexceptions.h ndlfortran.h LINK_LIBRARIES ndlfortranlib ndlstrutil ndlassert )
##=================
#

##=====CMatrix======
atlas_add_library(CMatrix CMatrix.cpp CMatrix.h ndlassert.h ndlexceptions.h CNdlInterfaces.h ndlstrutil.h ndlutil.h ndlfortran.h lapack.h LINK_LIBRARIES ndlfortranlib ndlstrutil ndlassert ndlutil )

#=====CTransform-===
atlas_add_library(CTransform CTransform.cpp CTransform.h CMatrix.h ndlassert.h ndlexceptions.h CNdlInterfaces.h ndlstrutil.h ndlutil.h ndlfortran.h lapack.h
LINK_LIBRARIES ndlfortranlib ndlstrutil ndlassert ndlutil CMatrix )
###=================

###=====CDist======
atlas_add_library(CDist CDist.cpp CDist.h CMatrix.h ndlassert.h ndlexceptions.h CNdlInterfaces.h ndlstrutil.h ndlutil.h ndlfortran.h lapack.h CTransform.h
LINK_LIBRARIES ndlfortranlib ndlstrutil ndlassert ndlutil CTransform CMatrix )

###=====COptimisable===
atlas_add_library(COptimisable COptimisable.cpp COptimisable.h CMatrix.h ndlassert.h ndlexceptions.h CNdlInterfaces.h ndlstrutil.h ndlutil.h ndlfortran.h lapack.h
LINK_LIBRARIES ndlfortranlib ndlstrutil ndlassert ndlutil CMatrix )
##=================
###=================

#=====CKern-===
atlas_add_library(CKern CKern.cpp CKern.h ndlassert.h ndlexceptions.h CTransform.h CMatrix.h CNdlInterfaces.h ndlstrutil.h ndlutil.h ndlfortran.h lapack.h CDataModel.h CDist.h LINK_LIBRARIES ndlutil ndlassert ndlfortranlib CTransform CMatrix ndlstrutil ndlutil CDist )
#=================


#=====CNoise-===
atlas_add_library(CNoise CNoise.cpp CNoise.h ndlexceptions.h ndlutil.h ndlassert.h ndlfortran.h ndlstrutil.h CMatrix.h CNdlInterfaces.h lapack.h
    CTransform.h COptimisable.h CDist.h CKern.h CDataModel.h LINK_LIBRARIES ndlutil ndlassert ndlfortranlib ndlstrutil CMatrix CTransform COptimisable CKern CDist )
#=================
#
#
atlas_add_library(CClctrl CClctrl.cpp CClctrl.h ndlstrutil.h ndlexceptions.h ndlutil.h ndlassert.h ndlfortran.h CMatrix.h CNdlInterfaces.h lapack.h
LINK_LIBRARIES ndlstrutil ndlutil ndlassert ndlfortranlib CMatrix )
#
atlas_add_library(CGplvm CGplvm.cpp CGplvm.h CMltools.h
ndlassert.h ndlexceptions.h ndlstrutil.h COptimisable.h CMatrix.h CNdlInterfaces.h ndlutil.h ndlfortran.h lapack.h CKern.h CTransform.h CDataModel.h CDist.h CNoise.h
LINK_LIBRARIES ndlassert ndlstrutil COptimisable CMatrix  ndlutil ndlfortranlib CKern CTransform CDist CNoise )


#
##=====CGP======
atlas_add_library(CGp CGp.cpp CGp.h CMltools.h ndlassert.h ndlexceptions.h ndlstrutil.h COptimisable.h CMatrix.h CNdlInterfaces.h ndlutil.h ndlfortran.h lapack.h CKern.h CTransform.h CDataModel.h CDist.h CNoise.h
LINK_LIBRARIES ndlassert ndlstrutil COptimisable CMatrix ndlutil ndlfortranlib CKern CTransform CDist CNoise)
###=================
##
##
##
###=================
##
##
##
#

###===CIvm===
atlas_add_library(CIvm CIvm.cpp CIvm.h CMltools.h ndlassert.h ndlexceptions.h ndlstrutil.h COptimisable.h CMatrix.h CNdlInterfaces.h ndlutil.h ndlfortran.h lapack.h CKern.h CTransform.h CDataModel.h CDist.h CNoise.h
LINK_LIBRARIES ndlassert ndlstrutil COptimisable CMatrix ndlutil ndlfortranlib CKern CTransform CDist CNoise )
##=================

#===gp====
# adding gp.cpp as an executable
atlas_add_executable( gp gp.cpp gp.h
    #INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS}
    #LINK_LIBRARIES  ${ROOT_LIBRARIES} ${BOOST_LIBRARIES} CClctrl CGp CMatrix ndlfortran CNoise ndlutil ndlstrutil CTransform COptimisable CKern CDist ndlassert )
LINK_LIBRARIES  CClctrl CGp CMatrix ndlfortranlib CNoise ndlutil ndlstrutil CTransform COptimisable CKern CDist ndlassert )
##target_link_library( gp CClctrl CGp CMatrix ndlfortran CNoise ndlutil ndlstrutil CTransform COptimisable CKern CDist ndlassert )

atlas_add_library( gp_TH1 gp_TH1.cpp gp_TH1.h gp_config.h ndlexceptions.h CMatrix.h CKern.h ndlexceptions.h ndlstrutil.h CGp.h ndlutil.h
        INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
        LINK_LIBRARIES  ${ROOT_LIBRARIES} CMatrix CKern ndlstrutil CGp ndlutil)

atlas_add_executable( dummy_main dummy_main.cpp dummy_main.h gp_TH1.h gp_config.h CMatrix.h
INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
LINK_LIBRARIES  ${ROOT_LIBRARIES} gp_TH1 CMatrix )

#gp.o: gp.cpp gp.h ndlexceptions.h ndlstrutil.h CMatrix.h ndlassert.h \
#  CNdlInterfaces.h ndlutil.h ndlfortran.h lapack.h CKern.h CTransform.h \
#  CDataModel.h CDist.h CGp.h CMltools.h COptimisable.h CNoise.h CClctrl.h
#
# linking gp with its dependence
#======
#
#
##====
##suspect CClgplvm is not being used, skipping
## skipping test that uses matlab
