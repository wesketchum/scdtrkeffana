# basic source code CMakeLists.txt

art_make_library( LIBRARY_NAME  ub_TrkEffTagAlgs
	SOURCE TagCreatorAlg.cxx
	LIBRARIES
		Geometry
		RecoBase
		${FHICLCPP}
                ${ROOT_BASIC_LIB_LIST}
                ${CETLIB}

)

install_headers()
install_source()
