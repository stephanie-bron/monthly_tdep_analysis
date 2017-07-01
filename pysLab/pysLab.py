#!/usr/bin/env python
from subprocess import Popen,PIPE,STDOUT
import uuid,shutil,os
import pickle

emin = 100.0
mjd_start = 57158
mjd_end = 57218

def RunAnaSingleSource(s_name, coord, lcFinal_dir):
        
        print "performing the analysis for source", s_name
        
        s_dir = "/data/user/sbron/monthly_tdep_analysis_test/AnaSources/" + s_name
        if not os.path.isdir(s_dir):
                os.makedirs(s_dir)

        block_file_path = lcFinal_dir + "/lc_" + s_name + "_" + str(coord[0]) + "_" + str(
                coord[1]) + "_1day.rawBlock.dat"

        ana_py_pars = " -b " + block_file_path + " -ra " + str(coord[0]) + " -dec " + str(
                coord[1]) + " -lc lc_" + s_name + "_" + str(coord[0]) + "_" + str(coord[1]) + "_1day" + " -s " + \
                      str(mjd_start) + " -e " + str(mjd_end)

        log_file_path = s_dir + "/" + s_name + "_Ana.log"
        logfile = open(log_file_path, "w")
       
        print "calling", "./start.sh python " + "/data/user/sbron/monthly_tdep_analysis_test/" + "/psLabScripts/FermiFlareAna.py" + ana_py_pars
      
        pro = Popen(". /data/user/sbron/monthly_tdep_analysis_test/pysLab/start.sh python " + "/data/user/sbron/monthly_tdep_analysis_test/" + \
                    "/psLabScripts/FermiFlareAna.py " + ana_py_pars, shell=True, stdout=logfile, stderr=logfile, close_fds=True, \
                    cwd="/data/user/sbron/monthly_tdep_analysis_test/")

        logfile.close()
        pro.wait()

        #For now commenting out the final part of the analysis which calls 'FermiFlareTS.py'

        # log_file_path=s_dir+"/"+s_name+"_TS0.log"
        # logfile=open(log_file_path, "w")
        # #self.par.lab['svn_revision']
        # #ana_py_pars=ana_py_pars+" -nts "+ self.par.lab['TS_trials']
        # ana_py_pars = ana_py_pars + " -nts " + "1000"
        # #print "calling","./start.sh python "+self.co_dir+"/psLabScripts/FermiFlareTS.py"+ana_py_pars
        # print "calling", "./start.sh python " + "/data/user/sbron/monthly_tdep_analysis_test/" + "/psLabScripts/FermiFlareTS.py" + ana_py_pars
        # #pro = Popen("./start.sh python "+self.co_dir+"/psLabScripts/FermiFlareTS.py "+ana_py_pars, shell=True,stdout=logfile, stderr=logfile, close_fds=True,cwd=self.co_dir)
        # pro = Popen(". /data/user/sbron/monthly_tdep_analysis_test/pysLab/start.sh python " + "/data/user/sbron/monthly_tdep_analysis_test/" + "/psLabScripts/FermiFlareTS.py " + ana_py_pars, shell=True,
        #             stdout=logfile, stderr=logfile, close_fds=True, cwd="/data/user/sbron/monthly_tdep_analysis_test/")
        #
        # logfile.close()
        # pro.wait()


lcFinal_dir="/data/user/sbron/FERMI/lightcurves/lcFinalFMonitored_from_"+str(mjd_start)+"_to_"+str(mjd_end)+"_EMin_"+str(emin)+"MeV"

#file containing informations about selected sources
filePath = '/data/user/sbron/monthly_tdep_analysis_test/fileSelectedSources.pickle'
pickleFile = pickle.load(open(filePath))

#for each source, start the analysis
for s_name, values in pickleFile.items():
        print s_name, values[:2]
        if s_name == '3C_279':
                runAna = RunAnaSingleSource(s_name,values[:2],lcFinal_dir)
        