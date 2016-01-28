# basic source code CMakeLists.txt

cet_find_library( BOOST_SYSTEM NAMES boost_system PATHS ENV BOOST_LIB NO_DEFAULT_PATH )

art_make_library( LIBRARY_NAME  ub_TrkEffTagAlgs
	SOURCE TagCreatorAlg.cxx TrkEffLinearityCheck.cxx
	LIBRARIES
		Geometry
		RecoBase
		RecoAlg
		#Utilities
		DetectorProperties_service
		${FHICLCPP}
                ${ROOT_BASIC_LIB_LIST}
                ${CETLIB}
		${BOOST_SYSTEM} #why the hell do I need this?
)

install_headers()
install_source()
