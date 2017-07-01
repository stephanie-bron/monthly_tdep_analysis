#!/usr/bin/env python
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
gROOT       =ROOT.gROOT
gSystem     =ROOT.gSystem
TCanvas     =ROOT.TCanvas
TLine       =ROOT.TLine
TGaxis      =ROOT.TGaxis
TLatex      =ROOT.TLatex
TArrow      =ROOT.TArrow
import argparse
from math import log10
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run psLab for given time period and source')
    parser.add_argument('-b', metavar='path to lc_source_ra_dec_1day.rawBlock.dat',  required=True,  type=str,  default="", help='the blocks file')
    parser.add_argument('-ra', required=True,  type=float, help='RA of the source in degrees')
    parser.add_argument('-dec',required=True,  type=float, help='DEC of the source in degrees')
    parser.add_argument('-lc', metavar='out_file_base_name',  required=True,  type=str,  default="", help='out file base name')
    parser.add_argument('-s',required=True,  type=float, help='analysis period start')
    parser.add_argument('-e',required=True,  type=float, help='analysis period end')
    
    args = parser.parse_args()
    
    blocksFile=args.__dict__["b"]
    source_ra=args.__dict__["ra"]
    source_dec=args.__dict__["dec"]
    outfile=args.__dict__["lc"]
    Ttmin=args.__dict__["s"]
    Ttmax=args.__dict__["e"]

    print "Opening blockfile",blocksFile
    print "With coordinates",source_ra,source_dec
    outfiledat=outfile+"_unblRez.dat"
    print "outfile",outfiledat
    print "Ttmin",Ttmin,"Ttmax",Ttmax
    
    ROOT.OPT_USEREALDATA=False
    
    gROOT.SetBatch()
    gROOT.Macro("$LAB_MAIN_DIR/psLab/llhTimeDep/loadlibs.C")
    ROOT.initialize_ran1(-156) # in psLab/trunk/rootExt/public/randomfunctions.C
    gSystem.Load("$LAB_MAIN_DIR/psLabScripts/ArkTime_C.so")
    gSystem.Load("$LAB_MAIN_DIR/psLabScripts/RescaledSigma_IC86_II_III_IV_SplineMPE_C.so"); #correction for the angular uncertainty bias: difference between angular resolution and paraboloid sigma (signal PDF: ), energy dependant correction based on a polynomial fit
    gROOT.ProcessLine(".L $LAB_MAIN_DIR/psLabScripts/TreeLoader_IC86_II_III_IV_ic.C");  # loading the data
    gROOT.ProcessLine("I3Ark ark_month;")
    gROOT.ProcessLine(".x $LAB_MAIN_DIR/psLabScripts/load_ark_ic86_II_III_IV_tdep_ic.C(ark_month, OPT_USEREALDATA,\"SplineMPE\")") #configuring the event loader with background/real data sample

    nSrcEvents = 0

    #loads blocks from Fermi denoised lightcurves files
    newllh_month=ROOT.NewLlhBlockTime()
    newllh_month.SetUseEnergy(True)
    newllh_month.SetOptimizeTolerance(0.01)
    newllh_month.SetMonitorLevel(0)
    newllh_month.SetEMaxRatioWarnOnlyOnce(1)
    newllh_month.JimsTerm_ = True
    newllh_month.SpectralPenalty_ = True
    newllh_month.ndof_ = 2.
    newllh_month.laglimit_ = 0.5
    newllh_month.SetLivetime(ROOT.ark_month.livetime/86400.) # 86400sec = 2 months
    newllh_month.SetBlocks(blocksFile,0.)
    newllh_month.SetAnalysisSet(ROOT.ark_month.psData)

    newllh_month.SetTimeBounds(Ttmin, Ttmax) #Ttmin, Ttmax given in the parameters of FermiFlareAna.py

    #ark_month = instance of I3Ark class
    print "IC_month: arkTime : " ,ROOT.ark_month.tmin , ROOT.ark_month.tmax # comes from load_ark_ic86... : ark.tmin = GetValueFromTree(ark.evLoader.GetBkgTree(),"tmin");
                                                                            # ark.tmax = GetValueFromTree(ark.evLoader.GetBkgTree(),"tmax");

    testSearch=ROOT.EquatorialDeg(source_ra, source_dec)
    
    newllh_month.SetSearchCoord( testSearch ) # in psLab/llhTimeDep/public/MultiGaussAnalysisFn.h:  srcCoord_ = testSearch

    #Get denoised lightcurves files in order to set the timePDF
    tPdf =ROOT.BlockTimePdf1() # in psLab/llhTimeDep/public/BlockLevel.h

    # blocksFile are the denoised lightcurves files that are stored in FERMI/lightcurves/<date>/*1day.rawBlock.dat
    # structure of those files is BeginTime(MJD) FluxLevel Duration(day)
    tPdf.SetBlockLevels(blocksFile,0.) # in llhTimeDep/public/BlockLevel.h, will create the tPDF, 0 when block is under threashold

    tPdf.fillHisto(1000)
                
    ROOT.ark_month.SetPointSource( testSearch, ROOT.PowerLawFlux(1.,-2.), tPdf) # in ArkTime.C, sets the source at the given location
    
    newllh_month.MaximizeLlh()
    print "prob , ts , p0 , p1 , p2 , p3"
    print newllh_month.GetEstProb(),newllh_month.GetTestStatistic(),newllh_month.GetPar(0),newllh_month.GetPar(1),newllh_month.GetPar(2),newllh_month.GetPar(3)
    
    c1 = TCanvas("c1","c1",1000,500)
    c1.SetGridx()
    tPdf_histo=tPdf.GetHisto()

    tPdf_histo.Scale(tPdf.GetNorm())
    tPdf_histo.SetLineWidth(2)
    tPdf_histo.SetLineColor(2)
    tPdf_histo.GetXaxis().SetTitle("MJD")
    tPdf_histo.GetYaxis().SetTitle("#Phi [ph s^{-1} cm^{-2}]")
    tPdf_histo.SetTitle("")
    tPdf_histo.GetYaxis().SetLabelColor(2)
    tPdf_histo.GetYaxis().SetTitleColor(2)
    tPdf_histo.GetYaxis().SetAxisColor(2)
    tPdf_histo.SetMinimum(0.)
    tPdf_histo.Draw()
        
    gamma = newllh_month.GetPar(1) # GetPar(1) = gammaBest_ in llhTimeDep/public/NewllhBlockTime.C
    eventVector = newllh_month.GetEventVector()


    max_ev_weight=1E-10
    evw=[]
    for j in range(eventVector.size()):

        if eventVector[j].GetMJD()<Ttmin or eventVector[j].GetMJD()>Ttmax:
            continue
        #event weights, blue lines in plot
        sRatio = eventVector[j].ProbFrom(testSearch) / eventVector[j].GetBkgSpaceProbFn().GetBkgProbDensity(eventVector[j])
        eweight = eventVector[j].GetEnergyProbFn().GetEnergyProbGamma(eventVector[j],gamma)/eventVector[j].GetEnergyProbFn().GetEnergyProbBkg(eventVector[j])
        logw=log10(sRatio*eweight)
        evw=evw+[[eventVector[j].GetMJD(),logw]]
        if sRatio*eweight > max_ev_weight: 
            max_ev_weight=sRatio*eweight
            
    min_draw_log_weight=-3.
    print max_ev_weight
    scale=tPdf_histo.GetMaximum()/(log10(max_ev_weight)-min_draw_log_weight)
    l=[]
    
    for ev in evw:
        if min_draw_log_weight<ev[1]:
            l.append(TLine(ev[0], 0.,ev[0] , (ev[1]-min_draw_log_weight)*scale))
            l[-1].SetLineWidth(1)
            l[-1].SetLineColor(4)
            l[-1].Draw("same")

    newllh_month.GetPar(3)
    l.append(TLine(tPdf_histo.GetXaxis().GetXmin(), newllh_month.GetPar(3), tPdf_histo.GetXaxis().GetXmax(), newllh_month.GetPar(3)))
    l[-1].SetLineWidth(2)
    l[-1].SetLineColor(2)
    l[-1].Draw("same")
    
    ax=TGaxis(tPdf_histo.GetXaxis().GetXmax(),0,tPdf_histo.GetXaxis().GetXmax(),tPdf_histo.GetMaximum(),min_draw_log_weight,log10(max_ev_weight),510,"+L")
    ax.SetLabelColor(4)
    ax.SetTitleColor(4)
    ax.SetLineColor(4)
    ax.SetTitle("log_{10}(time independent event weight)")
    ax.Draw("same")
    lat=TLatex()
    lat.SetTextColor(2)
    
    txtX=tPdf_histo.GetXaxis().GetXmin()+(tPdf_histo.GetXaxis().GetXmax()-tPdf_histo.GetXaxis().GetXmin())*3./4.
    txtY=tPdf_histo.GetMaximum()*1.07
    lat.DrawLatex(txtX,txtY,"Fitted Threshold")

    ar=TArrow(txtX-0.5,txtY,txtX-3,newllh_month.GetPar(3),0.02,">")
    ar.SetLineColor(2)
    ar.SetFillColor(2)
    ar.SetLineWidth(2)
    ar.Draw();
    c1.Update()
    c1.SaveAs(outfile+"_tPDF.root")
    c1.SaveAs(outfile+"_tPDF.png")
    
    rezfile=open(outfiledat,"w")
    rezfile.write("ndof 2: prob , ts , p0 , p1 , p2 , p3\n")
    rezfile.write(str(newllh_month.GetEstProb())+" "+str(newllh_month.GetTestStatistic())+" "+str(newllh_month.GetPar(0))+" "+str(newllh_month.GetPar(1))+" "+str(newllh_month.GetPar(2))+" "+str(newllh_month.GetPar(3))+"\n")
    rezfile.write("Ttmin "+str(Ttmin)+"\n")
    rezfile.write("Ttmax "+str(Ttmax)+"\n")
    rezfile.write("Time pdf norm "+str(tPdf.GetNorm())+"\n")
    rezfile.write("Time pdf bins (low_edge bincont.)\n")
    for i in range(tPdf_histo.GetNbinsX()):
        rezfile.write(str(tPdf_histo.GetBinLowEdge(i))+" "+str(tPdf_histo.GetBinContent(i))+"\n")
    rezfile.write("mjd log10(ew)\n")
    for ev in evw:
        rezfile.write(str(ev[0])+" "+str(ev[1])+"\n")
    rezfile.close()
 