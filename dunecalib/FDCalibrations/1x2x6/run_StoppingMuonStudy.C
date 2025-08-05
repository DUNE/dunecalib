{
  gROOT->ProcessLine(".L $CALIBDIR/srcs/Setup.cpp+");
  gROOT->ProcessLine(".L $CALIBDIR/srcs/Plane.cpp+");
  gROOT->ProcessLine(".L $CALIBDIR/srcs/Geometry.cpp+");
  gROOT->ProcessLine(".L $CALIBDIR/srcs/Utilities.cpp+");
  gROOT->ProcessLine(".L $CALIBDIR/srcs/ReadFiles.cpp+");
  gROOT->ProcessLine(".L $CALIBDIR/srcs/ConfigReader.cpp+");
  gROOT->ProcessLine(".L $CALIBDIR/srcs/EventProcessor.cpp+");
  gROOT->ProcessLine(".L $CALIBDIR/srcs/StoppingMuonStudy.cpp+");
}
