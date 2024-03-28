cet_build_plugin(XYZCalibServiceProtoDUNE art::service
              dunecalib::Calib
              art::Framework_Core
              art::Framework_Principal
              art::Framework_Services_Registry
	      art::Persistency_Common
	      art::Persistency_Provenance
	      canvas::canvas
	      art_root_io::tfile_support
	      ROOT::Core
	      art_root_io::TFileService_service
	      art::Framework_Services_System_FileCatalogMetadata_service
	      art::Utilities
	      messagefacility::MF_MessageLogger
  	      cetlib::cetlib 
	      cetlib_except::cetlib_except
   	      ROOT_BASIC_LIB_LIST
              
)

cet_build_plugin(LifetimeCalibServiceProtoDUNE art::service
      dunecalib::Calib
      art::Framework_Core
      art::Framework_Principal
      art::Framework_Services_Registry
      art::Persistency_Common
      art::Persistency_Provenance
      art_root_io::tfile_support
      ROOT::Core
      art_root_io::TFileService_service
      art::Framework_Services_System_FileCatalogMetadata_service
      art::Utilities
      messagefacility::MF_MessageLogger
      canvas::canvas
      cetlib::cetlib 
      cetlib_except::cetlib_except
      ROOT_BASIC_LIB_LIST

) 


install_headers()
install_fhicl()
install_source()
