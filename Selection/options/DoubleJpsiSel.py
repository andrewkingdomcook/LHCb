
from Gaudi.Configuration             import *
from Configurables                   import DaVinci
from DoubleJpsiSel.DoubleJpsiSelConf import DoubleJpsiSel
from Configurables                   import GaudiSequencer
from PhysConf.Filters                import LoKi_Filters
from GaudiKernel.SystemOfUnits       import MeV,GeV

JpsiSelectionSeq = GaudiSequencer("JpsiSelectionSeq")
MyJpsi2MuMu  = DoubleJpsiSel("Jpsi2MuMuTruth");

#set options
dataType      = "2011"
isSimulation = True # Data when False
is2012 = False
is2011sim08 = True


inputType     = "MDST"
stripping20r1 = False # set True for reStripping for MCtruth. Do not want on when using 20r1 latest tuple
if isSimulation:
    inputType = "DST"
    MyJpsi2MuMu.MCTruth = False #False for MC/ True for truth
    MyJpsi2MuMu.fillMCTruthTuple = False # controls if fills test tuple in MCTrue, want false for big run
    if MyJpsi2MuMu.MCTruth == True:
         MyJpsi2MuMu.NoAccCutMC = False # set to true or false to check effect of gen level acc in MC sample
    if dataType == "2011":
        stripping13b2011 = False #Set True want 13b starts from muons to avoid 3GeV cut on 2011 DimuonMDST

DaVinci().Simulation    = isSimulation           

inputFromHere = False
if inputFromHere:
    lfns = ["/lhcb/MC/2011/ALLSTREAMS.DST/00028319/0000/00028319_00000001_1.allstreams.dst",
            "/lhcb/MC/2011/ALLSTREAMS.DST/00028319/0000/00028319_00000002_1.allstreams.dst",
            "/lhcb/MC/2011/ALLSTREAMS.DST/00028319/0000/00028319_00000003_1.allstreams.dst",
            ]
    pre = "DATAFILE='PFN:root://eoslhcb.cern.ch//eos/lhcb/grid/prod" 
    post = "?svcClass=lhcbdisk' TYP='ROOT' OPT='READ'"
    complete = []
    for lfn in lfns:
        complete.append(pre + lfn + post)
    print complete
    EventSelector().Input = complete

inputFromHereDoubleJpsiUp = False
if inputFromHereDoubleJpsiUp:
    lfns = [
            "/lhcb/MC/2011/ALLSTREAMS.DST/00040646/0000/00040646_00000001_2.AllStreams.dst",
            "/lhcb/MC/2011/ALLSTREAMS.DST/00040646/0000/00040646_00000002_2.AllStreams.dst",
            "/lhcb/MC/2011/ALLSTREAMS.DST/00040646/0000/00040646_00000003_2.AllStreams.dst",
            "/lhcb/MC/2011/ALLSTREAMS.DST/00040646/0000/00040646_00000004_2.AllStreams.dst",
            ]
    pre = "DATAFILE='PFN:root://eoslhcb.cern.ch//eos/lhcb/grid/prod" 
    post = "?svcClass=lhcbdisk' TYP='ROOT' OPT='READ'"
    complete = []
    for lfn in lfns:
        complete.append(pre + lfn + post)
    print complete
    EventSelector().Input = complete

inputFromHereBuToJpsi = False
if inputFromHereBuToJpsi:
    lfns = [
        "/lhcb/MC/2011/ALLSTREAMS.DST/00031954/0000/00031954_00000001_1.allstreams.dst",
        "lhcb/MC/2011/ALLSTREAMS.DST/00031954/0000/00031954_00000002_1.allstreams.dst",
        "lhcb/MC/2011/ALLSTREAMS.DST/00031954/0000/00031954_00000004_1.allstreams.dst",
            ]
    pre = "DATAFILE='PFN:root://eoslhcb.cern.ch//eos/lhcb/grid/prod" 
    post = "?svcClass=lhcbdisk' TYP='ROOT' OPT='READ'"
    complete = []
    for lfn in lfns:
        complete.append(pre + lfn + post)
    print complete
    EventSelector().Input = complete

MyJpsi2MuMu.MCTen = (isSimulation and (dataType == "2010")) or (isSimulation and stripping13b2011) #restripping
MyJpsi2MuMu.MaxJpsiChi2 = 20.0
MyJpsi2MuMu.MinJpsiMass = 3.0 * GeV
MyJpsi2MuMu.MaxJpsiMass = 3.2 * GeV
MyJpsi2MuMu.MaxMuChi2 = 20.0
MyJpsi2MuMu.MinMuDll  = 0.0
MyJpsi2MuMu.MinMuPt   = 650 * MeV
MyJpsi2MuMu.FillVerbose            = False
MyJpsi2MuMu.FillAllJpsi            = True
DaVinci().DataType      = dataType
DaVinci().TupleFile     = "Tuple.root"
DaVinci().EvtMax        = -1 # 10000 # Number of events
DaVinci().PrintFreq     = 10000               # Reduces how often print success reading file message
print DaVinci().TupleFile


#DB Tags and input locations
if isSimulation == False:
    if dataType == "2011":
        ddDB     = "head-20110914" 
        condDB   = "head-20110914"
        lineName = 'MicroDSTDiMuonDiMuonIncLine'
        stream        = "Leptonic"
        rootInTES     = '/Event/' + stream + "/"
        inputLocation = ['Phys/' + lineName + '/Particles']       
else:
    import LoKiPhysMC.Track2MC_Configuration    
    if dataType == "2011":
        ddDB =  "MC11-20111102"   #2011 MC      
        condDB =  "sim-20111111-vc-md100" #2011 MC (vc-mu100 for magup)
        stream = "AllStreams"
        rootInTES = '/Event/'
        if not MyJpsi2MuMu.MCTruth and not stripping13b2011:
            lineName = 'AllStreams/MicroDSTDiMuonDimuonIncLine'
            inputLocation = ['AllStreams/Phys/MicroDSTDiMuonDiMuonIncLine/Particles']
        else:
            lineName = 'AllStreams/StdAllNoPIDsMuons'
            inputLocation = ['Phys/StdAllNoPIDsMuons/Particles']   
    elif dataType == "2010":
        ddDB =  "head-20101206"   #2010 MC      
        condDB =  "sim-20101210-vc-md100" #2010 MC
        stream = "AllStreams"
        rootInTES = '/Event/'
        if not MyJpsi2MuMu.MCTruth:
            #inputLocation = ['Phys/StdAllLooseMuons/Particles'] #changed to below for test
            lineName = 'AllStreams/StdAllNoPIDsMuons'
            inputLocation = ['Phys/StdAllNoPIDsMuons/Particles']
        else:
            lineName = 'AllStreams/StdAllNoPIDsMuons'
            inputLocation = ['Phys/StdAllNoPIDsMuons/Particles']
            
if is2012 == True:
    print "****** 2012 ******"
    inputType     = "DST"
    dataType = "2012"
    ddDB =  "dddb-20120831"   #2012 MC      
    condDB =  "sim-20121025-vc-mu100" #2012 MC (vc-mu100 for magup)
    #stream = "AllStreams"
    #rootInTES = '/Event/'
    lineName = 'AllStreams/StdLooseDiMuon'
    inputLocation = ['Phys/StdLooseDiMuon/Particles']
if is2011sim08 == True:
    ddDB = "Sim08-20130503-1"
    condDB = "Sim08-20130503-1-vc-mu100"
if not stripping20r1:
    MyJpsi2MuMu.RootInTES = rootInTES ##?????
    MyJpsi2MuMu.Inputs    = inputLocation #??????
    DaVinci().InputType     = inputType #???????


# carries out re-stripping - (a test that is generaly not needed)
if stripping20r1:
    import StrippingArchive.Stripping20.StrippingDiMuonNew            as DiMuon
    import StrippingSettings.Stripping20.LineConfigDictionaries_BandQ as LineSettings
    config = LineSettings.MicroDSTDiMuon['CONFIG']
    name   = 'MicroDST'
    builder = DiMuon.DiMuonConf( name , config )
    the_line  = builder.DiMuonLine
    selection = the_line.selection()
    from  PhysSelPython.Wrappers import Selection, SelectionSequence
    mySel = Selection('JPSI', Algorithm = MyJpsi2MuMu,RequiredSelections = [selection])
    mySeq = SelectionSequence('JPSI2MUMU', TopSelection = mySel)
    JpsiSelectionSeq.Members += [mySeq.sequence()]
    

if not isSimulation:
    from Configurables import L0TriggerTisTos
    from Configurables import TriggerTisTos
    MyJpsi2MuMu.addTool(TriggerTisTos())
    MyJpsi2MuMu.addTool(L0TriggerTisTos())

if not MyJpsi2MuMu.MCTruth:
    if dataType == "2011":
        fltrs = LoKi_Filters (
            STRIP_Code = """
            HLT_PASS ('StrippingMicroDSTDiMuonDiMuonIncLineDecision')  
            """
            )
    elif dataType == "2010" and isSimulation == False:
        fltrs = LoKi_Filters (
            STRIP_Code = """
            HLT_PASS ('StrippingMicroDSTDiMuonIncLineDecision')  
            """
            )


if not stripping20r1:
    JpsiSelectionSeq.Members += [MyJpsi2MuMu] 

DaVinci().DDDBtag       = ddDB
DaVinci().CondDBtag     = condDB

print DaVinci().TupleFile
DaVinci().UserAlgorithms = [ JpsiSelectionSeq ]


#NOTE ADDITION BELOW OF STRIPPING!!!!!!!!!
if dataType == "2011" and not stripping20r1:
    rawEventLoc = rootInTES + "DAQ/RawEvent"
    from Configurables import RawEventSelectiveCopy
    rawCopy = RawEventSelectiveCopy("CopyRawEvent")
    rawCopy.InputRawEventLocation = "/Event/Trigger/RawEvent"
    rawCopy.OutputRawEventLocation = rawEventLoc
    rawCopy.RawBanksToCopy = ["ODIN",
                              "HltSelReports" ,
                              "HltDecReports",
                              "L0Calo",
                              "L0CaloFull",
                              "L0DU",
                              "L0Muon",
                              "L0MuonProcCand",
                              "L0PU",
                              "HltRoutingBits"]
    from Configurables import CopyODIN
    CopyODIN().OutputPrefix = stream 
    from Configurables import (OdinTimeDecoder,ODINDecodeTool)
    OdinTimeDecoder().addTool( ODINDecodeTool, 'ODINDecodeTool' )
    #removed 2014.11.03 as does not have that property
    #OdinTimeDecoder().ODINDecodeTool.RawEventLocation = "/Event/Trigger/RawEvent"
    from Configurables import EventNodeKiller
    eventNodeKiller = EventNodeKiller('DAQkiller')
    eventNodeKiller.Nodes = ['DAQ']
    JpsiSelectionSeq.Members += [eventNodeKiller, CopyODIN(), rawCopy] 

#from Configurables import StoreExplorerAlg
#DaVinci().appendToMainSequence( [ StoreExplorerAlg(PrintEvt=1, ExploreRelations=
#True, Load=True) ] )

