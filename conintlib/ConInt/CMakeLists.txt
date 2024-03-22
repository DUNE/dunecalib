
FILE (GLOB src_files *.cxx)

art_make_library(  LIBRARY_NAME dunecalib::condb
		        SOURCE ${src_files}
			LIBRARIES
	                wda::wda
	                nuevdb::IFDatabase
                        cetlib::cetlib 
			cetlib_except::cetlib_except
                        messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
                        ROOT_BASIC_LIB_LIST
			
        )

install_headers()
install_fhicl()
install_source()
