
FILE (GLOB src_files *.cxx)

art_make_library(  LIBRARY_NAME  dunecalib::Calib
		        SOURCE ${src_files}
			LIBRARIES
	                wda::wda
	                nuevdb::IFDatabase
                        cetlib::cetlib 
			cetlib_except::cetlib_except
                        messagefacility::MF_MessageLogger
			fhiclcpp::fhiclcpp
			ROOT::Core ROOT::Hist
			
        )
 
cet_make_exec( NAME getXYZCalibProtoDUNE
               SOURCE getXYZCalibProtoDUNE.cc                 
               LIBRARIES
			dunecalib::Calib
               ) ### MIGRATE-CHANGED (migrate-3.22.02) - art_make_exec() -> cet_make_exec()

cet_make_exec( NAME getLinConstsProtoDUNE
               SOURCE getLinConstsProtoDUNE.cc                 
               LIBRARIES
			dunecalib::Calib
               ) ### MIGRATE-CHANGED (migrate-3.22.02) - art_make_exec() -> cet_make_exec()

install_headers()
install_fhicl()
install_source()
