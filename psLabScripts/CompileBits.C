{    
    //Using the gROOT pointer one has access to basically every object created in a ROOT based program
//     gROOT->Macro("${LAB_MAIN_DIR}/llhTimeDep/loadlibs.C");
    gROOT->Macro("${LAB_MAIN_DIR}/psLab/llhTimeDep/loadlibs.C");
    
    gSystem->SetBuildDir("$LAB_MAIN_DIR/psLabScripts/",true);
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/psLabScripts/ArkTime.C+");
    gROOT->ProcessLine(".L $LAB_MAIN_DIR/psLabScripts/RescaledSigma_IC86_II_III_IV_SplineMPE.C+");
}