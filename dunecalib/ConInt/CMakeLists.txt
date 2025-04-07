
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
                        ROOT::Tree
			
        )

art_make_exec( NAME getRunConditionsPDUNE
               SOURCE getRunConditionsPDUNE.cc
               LIBRARIES
                        dunecalib::condb
               )

install_headers()
install_fhicl()
install_source()
