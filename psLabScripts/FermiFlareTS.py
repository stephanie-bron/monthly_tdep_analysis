#!/usr/bin/env python
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
gROOT       =ROOT.gROOT
gSystem     =ROOT.gSystem
TCanvas     =ROOT.TCanvas
TGaxis      =ROOT.TGaxis
TH1D        =ROOT.TH1D
gPad        =ROOT.gPad
TF1         =ROOT.TF1

import argparse
from math import log10,exp
import numpy as np
from scipy.special import gamma

def chi2(x,p):
    return p[0] * pow(x[0], p[1]/2 - 1)*exp(-x[0]*0.5)/(pow(2,p[1]/2)*gamma(p[1]/2))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run psLab for given time period and source')
    parser.add_argument('-b', metavar='path to lc_source_ra_dec_1day.rawBlock.dat',  required=True,  type=str,  default="", help='the blocks file')
    parser.add_argument('-ra', required=True,  type=float, help='RA of the source in degrees')
    parser.add_argument('-dec',required=True,  type=float, help='DEC of the source in degrees')
    parser.add_argument('-lc', metavar='out_file_base_name',  required=True,  type=str,  default="", help='out file base name')
    parser.add_argument('-s',required=True,  type=float, help='analysis period start')
    parser.add_argument('-e',required=True,  type=float, help='analysis period end')
    parser.add_argument('-nts',required=True,type=int, help='TS0 trials')

    args = parser.parse_args()
    
    blocksFile=args.__dict__["b"]
    source_ra=args.__dict__["ra"]
    source_dec=args.__dict__["dec"]
    outfile=args.__dict__["lc"]
    Ttmin=args.__dict__["s"]
    Ttmax=args.__dict__["e"]
    ntrials=args.__dict__["nts"]

    print "Opening blockfile",blocksFile
    print "With coordinates",source_ra,source_dec
    outfiledat=outfile+"_TS0fit.dat"
    outfilepng=outfile+"_TS0.png"
    print "outfile",outfilepng
    print "Ttmin",Ttmin,"Ttmax",Ttmax
    
    ROOT.OPT_USEREALDATA=False
    
    gROOT.SetBatch()
    gROOT.Macro("$LAB_MAIN_DIR/psLab/llhTimeDep/loadlibs.C")
    ROOT.initialize_ran1(-156)
    gSystem.Load("$LAB_MAIN_DIR/psLabScripts/ArkTime_C.so")
    gSystem.Load("$LAB_MAIN_DIR/psLabScripts/RescaledSigma_IC86_II_III_IV_SplineMPE_C.so");
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/psLabScripts/TreeLoader_IC86_II_III_IV_ic.C");
    gROOT.ProcessLine("I3Ark ark_month;")
    gROOT.ProcessLine(".x $LAB_MAIN_DIR/psLabScripts/load_ark_ic86_II_III_IV_tdep_ic.C(ark_month, OPT_USEREALDATA,\"SplineMPE\")")

    nSrcEvents = 0

    newllh_month=ROOT.NewLlhBlockTime()
    newllh_month.SetUseEnergy(True)
    newllh_month.SetOptimizeTolerance(0.01)
    newllh_month.SetMonitorLevel(0)
    newllh_month.SetEMaxRatioWarnOnlyOnce(1)
    newllh_month.JimsTerm_ = True
    newllh_month.SpectralPenalty_ = True
    newllh_month.ndof_ = 2.
    newllh_month.laglimit_ = 0.5
    newllh_month.SetLivetime(ROOT.ark_month.livetime/86400.)
    newllh_month.SetBlocks(blocksFile,0.)
    newllh_month.SetAnalysisSet(ROOT.ark_month.psData)

    newllh_month.SetTimeBounds(Ttmin, Ttmax)
    print "IC_month: arkTime : " ,ROOT.ark_month.tmin , ROOT.ark_month.tmax

    testSearch=ROOT.EquatorialDeg(source_ra, source_dec)
    newllh_month.SetSearchCoord( testSearch )
    tPdf =ROOT.BlockTimePdf1()
    tPdf.SetBlockLevels(blocksFile,0.)
    tPdf.fillHisto(1000)

    ROOT.ark_month.SetPointSource( testSearch, ROOT.PowerLawFlux(1.,-2.), tPdf)

    TS0=np.zeros(ntrials)
    for trial in range(ntrials):
        ROOT.ark_month.psData.GenerateDataSet_with_nSrcEvents(0);
        newllh_month.MaximizeLlh()
        TS0[trial]=newllh_month.GetTestStatistic()

    TS_histo=TH1D("TS0","the Null test statistic", 100,0,max(TS0)*1.01)
    TS_histo.FillN(ntrials,TS0,np.ones(ntrials))
    
    cTS=TCanvas()
    TS_histo.Draw()
    gPad.SetLogy()
    fchi2 = TF1("fchi2",chi2,TS_histo.GetBinLowEdge(2),max(TS0)*1.01,2)
    fchi2.SetParameter(0,20)
    fchi2.SetParLimits(0,1E-2,ntrials)
    fchi2.SetParameter(1,2)
    fchi2.SetParLimits(1,1,5)
    fchi2.SetLineColor(2)
    ROOT.gStyle.SetOptFit()
    TS_histo.Fit(fchi2,"R")
    # in /public_html/monthly_tdep_web/source/ana_results/month_57158.0_57218.0/selectedSources
    cTS.SaveAs(outfilepng)
    
    
    rezfile=open(outfiledat,"w")
    rezfile.write("ndof "+str(fchi2.GetParameter(1))+" norm "+str(fchi2.GetParameter(0)) +"\n")
    rezfile.close()

    
