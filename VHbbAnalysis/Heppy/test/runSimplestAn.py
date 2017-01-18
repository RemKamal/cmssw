#!/usr/bin/env python

import optparse
import os,sys
import json
import pickle
import ROOT
from subprocess import Popen, PIPE

import sys
#sys.path.insert(0, r'/afs/cern.ch/user/r/rkamalie/workspace/private/tdr-style/tdrstyles')
sys.path.append('/afs/cern.ch/user/r/rkamalie/workspace/private/tdr-style')
import tdrstyle

import copy

#from PIL import Image                                                    
#import matplotlib.pyplot as plt 
#import matplotlib.image as mpimg                     

inGlobalFile="/afs/cern.ch/user/r/rkamalie/workspace/public/HHtree_from_Loop_4copy.root"

#several debug functions                                                                                                                       
debugMode = False
def bazinga (mes):
    if debugMode:
        print mes


def whoami():
    return inspect.stack()[1][3]

def whoisdaddy():
    return inspect.stack()[2][3]

def printme():
    print("Calling within the fcn '{me}', caller fcn is '{dad}'".format( me=whoami(), dad=whoisdaddy()) )



"""
Perform the analysis on a single file
"""
def runSimplestAn(inFileURL, outFileURL, xsec=None):

    print '...analysing %s' % inFileURL

    binLabels = ['no cut', '2jets>20', '2b-jets', '90<HCSV<150', 'Zto ee/mm', '60<Z<120', 'met>30']

    #book some histograms
    histos={
        
        
        # Book analysis independent histograms
        #*************************************
        'cutFlowLogScale'       :ROOT.TH1F('cutFlowLogScale', ';; Events', 7, 0, 7),
        'cutFlow'       :ROOT.TH1F('cutFlow', ';; Events', 7, 0, 7),
        
        
        'metpt_raw'     :ROOT.TH1F('metpt_raw'      ,';p_{T}(MET) [GeV]; Events'                , 35,0,350),
        
        'njets_raw'	:ROOT.TH1F('njets_raw'	,';Jet multiplicity; Events'				, 9,0,9),
        # 'nbtags_raw':ROOT.TH1F('nbtags_raw'	,';b-tag multiplicity; Events'				, 7,0,7),
        # 'nleps_raw'		:ROOT.TH1F('nleps_raw'		,';Lepton multiplicity; Events'				, 4,0,4),
        
        # #'jeten_raw'	:ROOT.TH1F('jeten_raw'		,';Energy [GeV]; Jets'						, 35,0,350),
        # 'jetpt_raw'     :ROOT.TH1F('jetpt_raw'      ,';p_{T}(jet) [GeV]; Events'                , 35,0,350),
        # 'jeteta_raw'    :ROOT.TH1F('jeteta_raw'     ,';#eta(jet); Events'                       , 24,-2.4,2.4),
        
        # #'bjetenls'	:ROOT.TH1F('bjetenls'	,';log(E);	1/E dN_{b jets}/dlog(E)'		, 20,3.,7.,),
        # #'bmjeteta'  :ROOT.TH1F('bmjeteta'   ,';#eta(b matched jet); Events'             , 24,-2.4,2.4),
        # #'bjeten_raw'	:ROOT.TH1F('bjeten_raw'		,';Energy [GeV]; bJets'						, 35,0,350),
        # 'bjetpt_raw'    :ROOT.TH1F('bjetpt_raw'     ,';p_{T}(b jet) [GeV]; Events'              , 35,0,350),
        # 'bjeteta_raw'   :ROOT.TH1F('bjeteta_raw'    ,';#eta(b jet); Events'                     , 24,-2.4,2.4),


        # #'lepen_raw'	:ROOT.TH1F('lepen_raw'		,';Energy(lepton) [GeV]; Events'						, 35,0,350),        
        # 'leppt_raw'     :ROOT.TH1F('leppt_raw'      ,';p_{T}(lepton) [GeV]; Events'             , 25,0,250),
        # 'lepeta_raw'    :ROOT.TH1F('lepeta_raw'     ,';#eta(lepton); Events'                    , 24,-2.4,2.4),

        
        # 'jet0pt_raw'	:ROOT.TH1F('jet0pt_raw'	        ,';Leading Jet p_{T} [GeV]; Events'	    , 35, 0., 350.),
        # 'jet1pt_raw'	:ROOT.TH1F('jet1pt_raw'	        ,';Subleading Jet p_{T} [GeV]; Events'	, 25, 0., 250.),
        # 'jet0eta_raw'	:ROOT.TH1F('jet0eta_raw'	,';Leading Jet #eta; Events'				, 24, -2.4, 2.4),
        # 'jet1eta_raw'	:ROOT.TH1F('jet1eta_raw'	,';Subleading Jet #eta; Events'			, 24, -2.4, 2.4),
        
        # 'bjet0pt_raw'	:ROOT.TH1F('bjet0pt_raw'	,';Leading b-Jet p_{T} [GeV]; Events'	    , 35, 0., 350.),
        # 'bjet1pt_raw'	:ROOT.TH1F('bjet1pt_raw'	,';Subleading b-Jet p_{T} [GeV]; Events'	, 25, 0., 250.),
        # 'bjet0eta_raw'	:ROOT.TH1F('bjet0eta_raw'	,';Leading b-Jet #eta; Events'				, 24, -2.4, 2.4),
        # 'bjet1eta_raw'	:ROOT.TH1F('bjet1eta_raw'	,';Subleading b-Jet #eta; Events'			, 24, -2.4, 2.4),
        
        # 'lep0pt_raw'	:ROOT.TH1F('lep0pt_raw'		,';Leading Lepton p_{T} [GeV]; Events'		, 25, 0., 250.),
        # 'lep1pt_raw'	:ROOT.TH1F('lep1pt_raw'		,';Subleading Lepton p_{T} [GeV]; Events'	, 20, 0., 200.),
        # 'lep0eta_raw'	:ROOT.TH1F('lep0eta_raw'	,';Leading Lepton #eta ; Events'		    , 24, -2.4, 2.4),
        # 'lep1eta_raw'	:ROOT.TH1F('lep1eta_raw'	,';Subleading Lepton #eta ; Events'	        , 24, -2.4, 2.4),

        
        # #'elpt'      :ROOT.TH1F('elpt'       ,';p_{T}(e) [GeV]; Events'                  , 25,0,250),
        # #'eleta'     :ROOT.TH1F('eleta'      ,';#eta(e) ; Events'                        , 24,-2.4,2.4),
        # #'mupt'      :ROOT.TH1F('mupt'       ,';p_{T}(#mu) [GeV]; Events'                , 25,0,250),
        # #'mueta'     :ROOT.TH1F('mueta'      ,';#eta(#mu) ; Events'                      , 24,-2.4,2.4),
        

        
        # #   Book histograms for HH specific signature of 2b, 2l, 2nu(met)
        # #****************************************************************
        
        # 'metpt'     :ROOT.TH1F('metpt'      ,';p_{T}(MET) [GeV]; Events'                , 35,0,350),
        
        # 'njets'	:ROOT.TH1F('njets'	,';Jet multiplicity; Events'				, 9,0,9),
        # 'nbtags':ROOT.TH1F('nbtags'	,';b-tag multiplicity; Events'				, 7,0,7),
        # 'nleps'		:ROOT.TH1F('nleps'		,';Lepton multiplicity; Events'				, 4,0,4),
        
        # #'jeten'	:ROOT.TH1F('jeten'		,';Energy [GeV]; Jets'						, 35,0,350),
        # 'jetpt'     :ROOT.TH1F('jetpt'      ,';p_{T}(jet) [GeV]; Events'                , 35,0,350),
        # 'jeteta'    :ROOT.TH1F('jeteta'     ,';#eta(jet); Events'                       , 24,-2.4,2.4),
        
        # #'bjetenls'	:ROOT.TH1F('bjetenls'	,';log(E);	1/E dN_{b jets}/dlog(E)'		, 20,3.,7.,),
        # #'bmjeteta'  :ROOT.TH1F('bmjeteta'   ,';#eta(b matched jet); Events'             , 24,-2.4,2.4),
        # #'bjeten'	:ROOT.TH1F('bjeten'		,';Energy [GeV]; bJets'						, 35,0,350),
        # 'bjetpt'    :ROOT.TH1F('bjetpt'     ,';p_{T}(b jet) [GeV]; Events'              , 35,0,350),
        # 'bjeteta'   :ROOT.TH1F('bjeteta'    ,';#eta(b jet); Events'                     , 24,-2.4,2.4),

        
        # #'lepen'	:ROOT.TH1F('lepen'		,';Energy(lepton) [GeV]; Events'						, 35,0,350),  
        # 'leppt'     :ROOT.TH1F('leppt'      ,';p_{T}(lepton) [GeV]; Events'             , 25,0,250),
        # 'lepeta'    :ROOT.TH1F('lepeta'     ,';#eta(lepton); Events'                    , 24,-2.4,2.4),
        
        
        # 'jet0pt'	:ROOT.TH1F('jet0pt'	        ,';Leading Jet p_{T} [GeV]; Events'	    , 35, 0., 350.),
        # 'jet1pt'	:ROOT.TH1F('jet1pt'	        ,';Subleading Jet p_{T} [GeV]; Events'	, 25, 0., 250.),
        # 'jet0eta'	:ROOT.TH1F('jet0eta'	,';Leading Jet #eta; Events'				, 24, -2.4, 2.4),
        # 'jet1eta'	:ROOT.TH1F('jet1eta'	,';Subleading Jet #eta; Events'			, 24, -2.4, 2.4),
        
        # 'bjet0pt'	:ROOT.TH1F('bjet0pt'	,';Leading b-Jet p_{T} [GeV]; Events'	    , 35, 0., 350.),
        # 'bjet1pt'	:ROOT.TH1F('bjet1pt'	,';Subleading b-Jet p_{T} [GeV]; Events'	, 25, 0., 250.),
        # 'bjet0eta'	:ROOT.TH1F('bjet0eta'	,';Leading b-Jet #eta; Events'				, 24, -2.4, 2.4),
        # 'bjet1eta'	:ROOT.TH1F('bjet1eta'	,';Subleading b-Jet #eta; Events'			, 24, -2.4, 2.4),
        
        # 'lep0pt'	:ROOT.TH1F('lep0pt'		,';Leading Lepton p_{T} [GeV]; Events'		, 25, 0., 250.),
        # 'lep1pt'	:ROOT.TH1F('lep1pt'		,';Subleading Lepton p_{T} [GeV]; Events'	, 20, 0., 200.),
        # 'lep0eta'	:ROOT.TH1F('lep0eta'	,';Leading Lepton #eta ; Events'		    , 24, -2.4, 2.4),
        # 'lep1eta'	:ROOT.TH1F('lep1eta'	,';Subleading Lepton #eta ; Events'	        , 24, -2.4, 2.4),

        



        }
    for key in histos:
        histos[key].Sumw2()
        histos[key].SetDirectory(0)

    #open file and loop over events tree
    fIn=ROOT.TFile.Open(inFileURL)
    tree=fIn.Get('tree')
    totalEntries= 100 if debugMode else tree.GetEntriesFast()
    for i in xrange(0,totalEntries):

        tree.GetEntry(i)
        if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )

        #generator level weight only for MC
        evWgt=1.0
        if xsec				 : evWgt  = xsec*tree.genWeight#[0]*tree.PUWeights[0]
        if tree.genWeight>0 : evWgt *= tree.genWeight

        histos['cutFlow'].Fill(0,evWgt)    

        nJets = nBJets = 0
        jetsP4 = []
        bJetsP4 = []
        # independent loop over jets in the event
        for ij in xrange(0,tree.nJet):

            #get the kinematics and select the jet
            jp4=ROOT.TLorentzVector()
            jp4.SetPtEtaPhiM(tree.Jet_pt[ij],tree.Jet_eta[ij],tree.Jet_phi[ij],tree.Jet_mass[ij])
            if jp4.Pt()<20 or ROOT.TMath.Abs(jp4.Eta())>2.4 : continue
            jetsP4.append(jp4)

            #count selected jet
            nJets +=1

            #save P4 for b-tagged jet
            if tree.Jet_btagCSV[ij]>0.4: # very loose CSV
                nBJets += 1
                bJetsP4.append(jp4)
                #if abs(tree.Jet_id[ij]) == 5: # is it id() fcn the right one???
                #    matchedBJetsP4.append(jp4)
        
        
        nLeps = 0
        lepsP4 = []
        
        # independent loop over leptons in the event
        for emu in xrange(0, tree.nvLeptons):
            if tree.genHiggsDecayMode%10000 == 23:  # work only with HtoZZ decay
                if tree.Vtype ==0 or tree.Vtype==1: # limit to Ztoee or Ztomm
                    #save P4 for leptons
                    lp4 = ROOT.TLorentzVector()
                    lp4.SetPtEtaPhiM(tree.vLeptons_pt[emu],tree.vLeptons_eta[emu],tree.vLeptons_phi[emu],tree.vLeptons_mass[emu])
                    if lp4.Pt()<20 or ROOT.TMath.Abs(lp4.Eta())>2.4 : continue #attention, Eta is limited to 2.4
                    lepsP4.append(lp4)
                    nLeps += 1
            #apply WaP from EG POG
            #nGoodLep +=1
                    
        if nJets != 2: continue
        histos['cutFlow'].Fill(1,evWgt)    

        if nBJets != 2: continue
        histos['cutFlow'].Fill(2,evWgt)    
        
        if tree.HCSV_mass < 90 or tree.HCSV_mass > 150: continue 
        #nGoodHiggs += 1
        histos['cutFlow'].Fill(3,evWgt)    

        if nLeps != 2: continue
        histos['cutFlow'].Fill(4,evWgt)    

        if tree.genHiggsDecayMode%10000 != 23: continue
        if not (tree.Vtype ==0 or tree.Vtype==1): continue # limit to Ztoee or Ztomm
        if tree.V_mass <60 or tree.V_mass > 120 : continue
        histos['cutFlow'].Fill(5,evWgt)    

        if tree.met_pt < 30: continue
        histos['cutFlow'].Fill(6,evWgt)    

        

        # histos['metpt_raw'].Fill(tree.met_pt,evWgt)    
        # histos['njets_raw'].Fill(nJet,evWgt)
        # histos['nbtags_raw'].Fill(nJet,evWgt) # attention I do not know which one bJet is right to use
        # histos['nleps_raw'].Fill(tree.nGenLep,evWgt)
        
        # #histos['jeten_raw'] .Fill(Jet_taggedJetsP4[0].Pt() ,evWgt)
        # histos['jetpt_raw'] .Fill(tree.Jet_pt ,evWgt)
        # histos['jeteta_raw'].Fill(tree.Jet_eta,evWgt)
        
        # histos['bjetpt_raw'] .Fill(tree.Jet_pt ,evWgt)# attention I do not know which one bJet is right to use
        # histos['bjeteta_raw'].Fill(tree.Jet_eta,evWgt)# attention I do not know which one bJet is right to use

        # histos['leppt_raw']    .Fill(tree.GenLep_pt ,evWgt)
        # histos['lepeta_raw']   .Fill(tree.GenLep_eta,evWgt)

        # # twelve lines below - check how to access leading Pt objects and what to use for btags
        # histos['jet0pt_raw'] .Fill(tree.Jet_pt[0] ,evWgt)
        # histos['jet0eta_raw'].Fill(tree.Jet_eta[0],evWgt)
        # histos['jet1pt_raw'] .Fill(tree.Jet_pt[1] ,evWgt)
        # histos['jet1eta_raw'].Fill(tree.Jet_eta[1],evWgt)
        # histos['bjet0pt_raw'] .Fill(tree.Jet_pt[0] ,evWgt)
        # histos['bjet0eta_raw'].Fill(tree.Jet_eta[0],evWgt)
        # histos['bjet1pt_raw'] .Fill(tree.Jet_pt[1] ,evWgt)
        # histos['bjet1eta_raw'].Fill(tree.Jet_eta[1],evWgt)
        # histos['lep0pt_raw']    .Fill(tree.GenLep_pt[0] ,evWgt)
        # histos['lep1pt_raw']    .Fill(tree.GenLep_pt[1] ,evWgt)
        # histos['lep0eta_raw']   .Fill(tree.GenLep_eta[0],evWgt)
        # histos['lep1eta_raw']   .Fill(tree.GenLep_eta[1],evWgt)






















        # # now start applying cuts
        # if nJets<2 : continue
        # histos['nbtags_raw'].Fill(nBtags,evWgt)
        
        # if nBtags !=2 : continue
        

        # #ready to fill the histograms
        # histos['nbtags'].Fill(nBtags,evWgt)
        # histos['njets'].Fill(nJets,evWgt)
        

        

        # # lepton pt and eta; also checks pt-ordering
        # if tree.Lepton_pt[0] > tree.Lepton_pt[1]:
        #     histos['lep0pt']	.Fill(tree.Lepton_pt[0]	,evWgt)
        #     histos['lep1pt']	.Fill(tree.Lepton_pt[1]	,evWgt)
        #     histos['lep0eta']	.Fill(tree.Lepton_eta[0],evWgt)
        #     histos['lep1eta']	.Fill(tree.Lepton_eta[1],evWgt)
        # else:
        #     histos['lep0pt']	.Fill(tree.Lepton_pt[1]	,evWgt)
        #     histos['lep1pt']	.Fill(tree.Lepton_pt[0]	,evWgt)
        #     histos['lep0eta']	.Fill(tree.Lepton_eta[1],evWgt)
        #     histos['lep1eta']	.Fill(tree.Lepton_eta[0],evWgt)

        # # tagged jets pt and eta
        # if len(taggedJetsP4) > 0:
        #     histos['bjet0pt'] .Fill(taggedJetsP4[0].Pt() ,evWgt)
        #     histos['bjet0eta'].Fill(taggedJetsP4[0].Eta(),evWgt)

        # if len(taggedJetsP4) > 1:
        #     histos['bjet1pt'] .Fill(taggedJetsP4[1].Pt() ,evWgt)
        #     histos['bjet1eta'].Fill(taggedJetsP4[1].Eta(),evWgt)

        # #use up to two leading b-tagged jets
        # for ij in xrange(0,len(taggedJetsP4)):
        #     histos['bjeten'].Fill(taggedJetsP4[ij].E(),evWgt)
        #     histos['bjetenls'].Fill(ROOT.TMath.Log(taggedJetsP4[ij].E()),evWgt/taggedJetsP4[ij].E())
        #     histos['bjetpt'].Fill(taggedJetsP4[ij].Pt(),evWgt)
        #     histos['bjeteta'].Fill(taggedJetsP4[ij].Eta(),evWgt)
        # for ij in xrange(0,len(matchedJetsP4)):
        #     histos['bmjeteta'].Fill(matchedJetsP4[ij].Eta(),evWgt)
        # #use up to two leading b-tagged jets
        # for ij in xrange(0,len(jetsP4)):
        #     histos['jetpt'].Fill(jetsP4[ij].Pt(),evWgt)
        #     histos['jeteta'].Fill(jetsP4[ij].Eta(),evWgt)
        # for ij in xrange(0,tree.nLepton):
        #     histos['leppt'].Fill(tree.Lepton_pt[ij],evWgt)
        #     histos['lepeta'].Fill(tree.Lepton_eta[ij],evWgt)
        #     if abs(tree.Lepton_id[ij]) == 11:
        #         histos['elpt'].Fill(tree.Lepton_pt[ij],evWgt)
        #         histos['eleta'].Fill(tree.Lepton_eta[ij],evWgt)
        #     if abs(tree.Lepton_id[ij]) == 13:
        #         histos['mupt'].Fill(tree.Lepton_pt[ij],evWgt)
        #         histos['mueta'].Fill(tree.Lepton_eta[ij],evWgt)

        
    
    #histos['cutFlowLogScale'] = copy.deepcopy(histos['cutFlow']) 
    histos['cutFlowLogScale'] = histos['cutFlow'].Clone() 

    for bin, label in enumerate(binLabels):
        histos['cutFlow'].GetXaxis().SetBinLabel(bin+1,label)
        histos['cutFlowLogScale'].GetXaxis().SetBinLabel(bin+1,label)

    
        

        
    #all done with this file
    fIn.Close()
    
    tdrstyle.setTDRStyle() 
    c1 = ROOT.TCanvas("c1", "c1")
    c1.cd()
    #ROOT.gPad.SetLogy()
    histos['cutFlowLogScale'].Draw()
    #histos['cutFlowLogScale'].Draw()
    c1.SaveAs("cutFlow.pdf") 
    c1.SaveAs("cutFlow.png")

    #c1.SaveAs("cutFlowLogScale.pdf") 
    #c1.SaveAs("cutFlowLogScale.png")
    #raw_input("press some key")

    #save histograms to file
    fOut=ROOT.TFile.Open(outFileURL,'RECREATE')
    for key in histos: histos[key].Write()
    fOut.Close()


"""
Wrapper to be used when run in parallel
"""
def runSimplestAnPacked(args):

    try:
        return runSimplestAn(inFileURL=args[0],
                             outFileURL=args[1],
                             xsec=args[2])
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],args[0])
        print 50*'<'
        return False


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-j', '--json',		 dest='json'  ,		 help='json with list of files',	  default=None,		   type='string')
    parser.add_option('-i', '--inDir',		 dest='inDir',		 help='input directory with files',   default=None,		   type='string')
    parser.add_option('-o', '--outDir',		 dest='outDir',		 help='output directory',			  default='analysis',  type='string')
    parser.add_option('-n', '--njobs',		 dest='njobs',		 help='# jobs to run in parallel',	  default=0,		   type='int')
    (opt, args) = parser.parse_args()

    #read list of samples
    if len(args) > 1:
        jsonFile = open(opt.json,'r')
        samplesList=json.load(jsonFile,encoding='utf-8').items()
        jsonFile.close()

    #prepare output
    if len(opt.outDir)==0	 : opt.outDir='./'
    os.system('mkdir -p %s' % opt.outDir) #overwrite if exists

    #create the analysis jobs
    taskList = []
    if len(args) > 1 and sampleList:
        for sample, sampleInfo in samplesList: 
            inFileURL  = '/%s/%s.root' % (opt.inDir,sample)
        #if not os.path.isfile(inFileURL): continue
            xsec=sampleInfo[0] if sampleInfo[1]==0 else None		
            outFileURL = '%s/%s.root' % (opt.outDir,sample)
            taskList.append( (inFileURL,outFileURL,xsec) )
    else:
        inFileURL  = inGlobalFile
        xsec = None
        outFileURL = '%s/%s.root' % (opt.outDir,"plots")
        taskList.append( (inFileURL,outFileURL,xsec) )
  
    #run the analysis jobs
    if opt.njobs == 0:
        for inFileURL, outFileURL, xsec in taskList:
            runSimplestAn(inFileURL=inFileURL, outFileURL=outFileURL, xsec=xsec)
    else:
        from multiprocessing import Pool
        pool = Pool(opt.njobs)
        pool.map(runSimplestAnPacked,taskList)

    #all done here
    print '\nAnalysis results are available in %s' % opt.outDir
    #exit(0)
    


"""
for execution from another script
"""
if __name__ == "__main__":
    #sys.exit(main())
    main()
    #img = mpimg.imread('cutFlow.png')
    #plt.imshow(img)
    #img = Image.open('cutFlow.png')
    #img.show() 


