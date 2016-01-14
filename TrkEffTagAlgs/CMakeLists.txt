# basic source code CMakeLists.txt

art_make_library( LIBRARY_NAME  ub_TrkEffTagAlgs
	SOURCE TagCreatorAlg.cxx
	LIBRARIES
		 RecoBase
		 ${FHICLCPP}
                 ${ROOT_BASIC_LIB_LIST}

)

install_headers()
install_source()
