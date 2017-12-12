import ROOT
from DataFormats.FWLite import Events, Handle
import sys
from array import array

doByParticle=True

events = Events(sys.argv[1])

outf = ROOT.TFile(sys.argv[2],'recreate')
tree = ROOT.TTree('tree','tree')

class AllBr:
    def __init__(self):
        pass
br  = AllBr()

def addCollection(name,quantities,maxN):
    setattr(br,'max_n%s'%name,maxN)
    setattr(br,'n%s'%name,array('i',[0]))
    tree.Branch('n%s'%name, getattr(br,'n%s'%name), 'n%s/I'%name)
    for quant in quantities:
        setattr(br,'%s_%s'%(name,quant),array('i' if 'id_' in quant else 'f',getattr(br,'max_n%s'%name)*[0]))
        tree.Branch('%s_%s'%(name,quant),getattr(br,'%s_%s'%(name,quant)),'%s_%s[n%s]/%s'%(name,quant,name,'I' if 'id_' in quant else 'F'))

def zside(id_):
    return 1 if (id_&0x10000) else -1
def ietaAbs(id_):
    return (id_>>9)&0x7F
def ieta(id_):
    return zside(id_)*ietaAbs(id_)
def iphi(id_):
    return id_&0x1FF
import math
def dphi(phi1,phi2):
    x = phi1 - phi2
    while x > math.pi:
        x -= 2*math.pi
    while x < -math.pi:
        x += 2*math.pi
    return x
def dieta(eta1,eta2):
    return eta1-eta2 if (eta1*eta2>0) else eta1-eta2-1
def diphi(phi1,phi2):
    x = phi1 - phi2
    while x > 180:
        x -= 360
    while x < -180:
        x += 360
    return x
def zsideHCAL(id_):
    return 1 if (id_&0x2000) else -1
def ietaAbsHCAL(id_):
    return (id_>>7)&0x3F 
def ietaHCAL(id_):
    return zsideHCAL(id_)*ietaAbsHCAL(id_)     
def iphiHCAL(id_):
    return id_&0x7F       

h_ECALxtals = Handle('vector<l1tpf::Particle>')
h_ECALclusters = Handle('vector<l1tpf::Particle>')
h_HCALtowers = Handle('vector<l1tpf::Particle>')
h_genparticles = Handle('vector<reco::GenParticle>')

addCollection('ECALxtal',['pt','eta','ieta','iphi'],10000)
addCollection('ECALclus',['pt','eta','phi'],10000)
addCollection('HCALtowers',['pt','eta','phi','dr','ieta','iphi'],10000)
addCollection('HCALclus',['pt', 'eta', 'phi'],10000)
addCollection('GenPart',['pt','eta','phi','pdgId','status'],100)
if doByParticle:
    addCollection('ECALxtalInClus',['pt','eta','dieta','diphi','ieta','iphi'],1000)
    addCollection('HCALtowInClus',['pt','dieta','diphi','ieta','iphi'],1000)

def Reset():
    br.nECALxtal[0] = 0
    br.nECALclus[0] = 0
    br.nGenPart[0] = 0
    br.nHCALtowers[0] = 0
    br.nHCALclus[0] = 0
    if hasattr(br,'nECALxtalInClus'): br.nECALxtalInClus[0] = 0
    if hasattr(br,'nHCALtowInClus'): br.nHCALtowInClus[0] = 0

def FillGenPart(part,i):
    br.GenPart_pt[i] = part.pt()
    br.GenPart_eta[i] = part.eta()
    br.GenPart_phi[i] = part.phi()
    br.GenPart_pdgId[i] = part.pdgId()
    br.GenPart_status[i] = part.status()
    br.nGenPart[0] += 1
    #print "FOUND gen particle id=", part.pdgId(), " status=", part.status(), " pt=", part.pt(), " eta=",part.eta(), " phi=", part.phi()

'''
################### MARCO ECAL stuff ######################## 
def FillECALxtal(xtal,i,relativeto=None):
    br.ECALxtal_pt[i] = xtal.pt()
    id_ = (xtal.getSeed())()
    br.ECALxtal_ieta[i] = ieta(id_) #if relativeto==None else dieta(ieta(id_),ieta(relativeto))
    br.ECALxtal_iphi[i] = iphi(id_) #if relativeto==None else diphi(iphi(id_),iphi(relativeto))
    br.nECALxtal[0] += 1
    #print i,xtal.iEta(),ieta(id_),br.ECALxtal_ieta[i], ieta(id_) if relativeto==None else dieta(ieta(id_),ieta(relativeto)),xtal.pt()
def FillECALclus(clus,i):
    br.ECALclus_pt[i] = clus.pt()
    br.ECALclus_eta[i] = clus.eta()
    br.ECALclus_phi[i] = clus.phi()
    br.nECALclus[0] += 1
def FillECALxtalInClus(xtal,i,seedid):
    br.ECALxtalInClus_pt[i] = xtal.pt()
    id_ = (xtal.getSeed())()
    br.ECALxtalInClus_dieta[i] = dieta(ieta(id_),ieta(seedid))
    br.ECALxtalInClus_diphi[i] = diphi(iphi(id_),iphi(seedid))
    br.ECALxtalInClus_ieta[i] = ieta(id_)
    br.ECALxtalInClus_iphi[i] = iphi(id_)
    br.nECALxtalInClus[0] += 1
'''

################### MY ECAL stuff ######################## 
def FillECALclus(clusters,genp,xtals):

    best_clus = -1
    minpt = 0.5
    #print "nclust = ", len(clusters)
    for i in xrange(len(clusters)):
     clus = clusters[i]
     deta = part.eta()-clus.eta()
     dphi = part.phi()-clus.phi()
     dr = math.sqrt(deta*deta+dphi*dphi)
     if clus.pt() > minpt and dr < 0.5:
      best_clus = i
      minpt = clus.pt()
      #print "* found cluster ", i, " pt ", clus.pt(), " eta ", clus.eta(), " phi ", clus.phi(), " gen pt ", part.pt(), " gen eta ", part.eta(), " gen phi ", part.phi(), " dr ", dr
     
    if best_clus == -1: return
    #print "-----------> FOUND BEST CLUS=", best_clus  
    
    clus = clusters[best_clus]
    br.ECALclus_pt[0] = clus.pt()
    br.ECALclus_eta[0] = clus.eta()
    br.ECALclus_phi[0] = clus.phi()
    br.nECALclus[0] += 1
    
    seedid = (clus.getSeed())()
    for i in xrange(clus.getDetIdCollection().size()):
     FillECALxtalInClus( filter(lambda xtal: (xtal.getSeed())()==clus.getDetId(i)(), xtals)[0],i,seedid)
     #j = 0
     #for i in xrange(len(xtals)):
      #id_ = (xtals[i].getSeed())()
      #if abs(dieta(ieta(id_),ieta(seedid)))<10 and abs(diphi(iphi(id_),iphi(seedid)))<10:
      #FillECALxtal(xtals[i],j,seedid)
      #j+=1
	     
def FillECALxtalInClus(xtal,i,seedid):
    br.ECALxtalInClus_pt[i] = xtal.pt()
    br.ECALxtalInClus_eta[i] = xtal.eta()
    id_ = (xtal.getSeed())()
    br.ECALxtalInClus_dieta[i] = dieta(ieta(id_),ieta(seedid))
    br.ECALxtalInClus_diphi[i] = diphi(iphi(id_),iphi(seedid))
    br.ECALxtalInClus_ieta[i] = ieta(id_)
    br.ECALxtalInClus_iphi[i] = iphi(id_)
    br.nECALxtalInClus[0] += 1     		      

def FillECALxtal(xtal,i,relativeto=None):
    br.ECALxtal_pt[i] = xtal.pt()
    br.ECALxtal_eta[i] = xtal.eta()
    id_ = (xtal.getSeed())()
    br.ECALxtal_ieta[i] = ieta(id_) #if relativeto==None else dieta(ieta(id_),ieta(relativeto))
    br.ECALxtal_iphi[i] = iphi(id_) #if relativeto==None else diphi(iphi(id_),iphi(relativeto))
    br.nECALxtal[0] += 1
    #print i,xtal.iEta(),ieta(id_),br.ECALxtal_ieta[i], ieta(id_) if relativeto==None else dieta(ieta(id_),ieta(relativeto)),xtal.pt()
    
################### HCAL stuff ########################    
def FillHCALtower(tow,i,dr):
    id_ = (tow.getDetId(0))()
    #print i,tow.pt(),tow.eta(),tow.phi()," ieta tow=", tow.iEta()," iphi tow=", tow.iPhi()
    br.HCALtowers_pt[i] = tow.pt() 
    br.HCALtowers_eta[i] = tow.eta()
    br.HCALtowers_phi[i] = tow.phi()
    br.HCALtowers_ieta[i] = ietaHCAL(id_)
    br.HCALtowers_iphi[i] = iphiHCAL(id_)
    br.HCALtowers_dr[i] = dr   
    br.nHCALtowers[0] += 1
def FillHCALclus(towers,genp):
    minpt = 0.5
    seed = -1
    for i in xrange(len(towers)):
     tow = towers[i]
     deta = part.eta()-tow.eta()
     dphi = part.phi()-tow.phi()
     dr = math.sqrt(deta*deta+dphi*dphi)
     if tow.pt() > minpt and dr < 0.5:
      seed = i
      minpt = tow.pt()    
    #print "-----------> FOUND SEED=", seed,
    
    if seed == -1: return
    
    FillHCALtowInClus(towers[seed],0,(towers[seed].getDetId(0))())
    
    ptclus = towers[seed].pt()
    etaclus = towers[seed].eta()
    phiclus = towers[seed].phi()
    ntows = 1
    for i in xrange(len(HCALtowers)):
     towSeed = towers[seed]
     if i == seed: continue
     tow = towers[i]
     if tow.pt() <= 0.5: continue
     for ie in xrange(-1,1):
      for ip in xrange(-1,1):
       if tow.iEta() == towSeed.iEta()+ie and tow.iPhi() == towSeed.iPhi()+ip:
        #print " : ie=",ie," ip=",ip," i=",i," ieta=",tow.iEta()," iphi=",tow.iPhi()," pt=", tow.pt(), 
	ptclus+=tow.pt()
	FillHCALtowInClus(tow,ntows,(towSeed.getDetId(0))())
	ntows += 1
    #print " ===== ptclus = ", ptclus
    br.HCALclus_pt[0] = ptclus
    br.HCALclus_eta[0] = etaclus
    br.HCALclus_phi[0] = phiclus
    br.nHCALclus[0] += 1 	 
def FillHCALtowInClus(tow,i,seedid):
    br.HCALtowInClus_pt[i] = tow.pt()
    id_ = (tow.getDetId(0))()
    br.HCALtowInClus_ieta[i] = ietaHCAL(id_)
    br.HCALtowInClus_iphi[i] = iphiHCAL(id_)
    br.HCALtowInClus_dieta[i] = dieta(ietaHCAL(id_),ietaHCAL(seedid))
    br.HCALtowInClus_diphi[i] = diphi(iphiHCAL(id_),iphiHCAL(seedid))
    br.nHCALtowInClus[0] += 1
    
################### main ######################## 
for iev,ev in enumerate(events):
    if iev%1000==0: print 'Processing event %d'%(iev+1,)
    if iev > 2000: break
    #print "*************** ", iev
    ev.getByLabel(("l1tPFEcalProducerFromTPDigis","crystals"),h_ECALxtals)
    ev.getByLabel(("l1tPFEcalProducerFromL1EGCrystalClusters",""),h_ECALclusters)
    ev.getByLabel(("l1tPFHcalProducerFromTPDigis",""),h_HCALtowers)
    ev.getByLabel(("genParticles",""),h_genparticles)
    ECALxtals = h_ECALxtals.product()
    ECALclusters = h_ECALclusters.product()
    HCALtowers = h_HCALtowers.product()
    genparticles = h_genparticles.product()

    if not doByParticle:
        Reset()
        j = 0
	print len(HCALtowers)
        for i in xrange(min(br.max_nGenPart,len(genparticles))):
            part = genparticles[i]
            if len(genparticles)>2 and part.status()!=23: continue
            FillGenPart(part,j)
            j += 1
        for i in xrange(min(br.max_nECALclus,len(ECALclusters))):
            FillECALclus(ECALclusters[i],i)
        for i in xrange(min(br.max_nECALxtal,len(ECALxtals))):
            FillECALxtal(ECALxtals[i],i)
        for i in xrange(len(HCALtowers)):
	    FillHCALtower(HCALtowers[i],i,0) 
	tree.Fill()
    else:
        for k in xrange(min(br.max_nGenPart,len(genparticles))):
            Reset()
            part = genparticles[k]
            if len(genparticles)>2 and part.status()!=23: continue
            FillGenPart(part,0)

            '''
            ### MARCO ECAL TREE ####
            best = 0            
	    for i in xrange(len(ECALclusters)):
                clus = ECALclusters[-1-i] 
                if abs(clus.pt()/part.pt()-1)>0.5: continue
                if abs(clus.eta()-part.eta())>0.2: continue
                if abs(dphi(clus.phi(),part.phi()))>0.2: continue
                best = -1-i	        
	    if best!=0:
                clus = ECALclusters[best]
                seedid = (clus.getSeed())()
                FillECALclus(clus,0)
                for i in xrange(clus.getDetIdCollection().size()):
                    FillECALxtalInClus(filter(lambda xtal: (xtal.getSeed())()==clus.getDetId(i)(), ECALxtals)[0],i,seedid)
                j = 0
                for i in xrange(min(br.max_nECALxtal,len(ECALxtals))):
                    id_ = (ECALxtals[i].getSeed())()
                    if abs(dieta(ieta(id_),ieta(seedid)))<10 and abs(diphi(iphi(id_),iphi(seedid)))<10:
                        FillECALxtal(ECALxtals[i],j,seedid)
                        j+=1
            '''
	    ### MY ECAL TREE ####
	    FillECALclus(ECALclusters,part,ECALxtals)
	    
	    
	    ### HCAL TREE ####
	    FillHCALclus(HCALtowers,part)
	    for i in xrange(len(HCALtowers)):
	        tow = HCALtowers[i]
		detaN = tow.eta()-part.eta()
		dphiN = tow.phi()-part.phi()
		dr = math.sqrt(detaN*detaN+dphiN*dphiN)
		FillHCALtower(tow,i,dr)
						
	    tree.Fill()
            
outf.Write()
outf.Close()
