cet_build_plugin(RunConditionsServicePDUNE art::service
              dunecalib::Calib
              dunecalib::condb
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
              ROOT::Tree
              larcore::ServiceUtil
	    )



install_headers()
install_fhicl()
install_source()
