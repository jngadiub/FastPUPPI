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


h_xtals = Handle('vector<pair<DetId,float> >')
h_clusters = Handle('vector<l1slhc::L1EGCrystalCluster>')
h_genparticles = Handle('vector<reco::GenParticle>')

addCollection('Xtal',['pt','ieta','iphi'],10000)
addCollection('Clus',['pt','eta','phi'],10000)
addCollection('GenPart',['pt','eta','phi','pdgId','status'],100)
if doByParticle:
    addCollection('XtalInClus',['pt','dieta','diphi'],1000)

def Reset():
    br.nXtal[0] = 0
    br.nClus[0] = 0
    br.nGenPart[0] = 0
    if hasattr(br,'nXtalInClus'): br.nXtalInClus[0] = 0

def FillGenPart(part,i):
    br.GenPart_pt[i] = part.pt()
    br.GenPart_eta[i] = part.eta()
    br.GenPart_phi[i] = part.phi()
    br.GenPart_pdgId[i] = part.pdgId()
    br.GenPart_status[i] = part.status()
    br.nGenPart[0] += 1
def FillXtal(xtal,i,relativeto=None):
    br.Xtal_pt[i] = xtal.second
    id_ = (xtal.first)()
    br.Xtal_ieta[i] = ieta(id_) if relativeto==None else dieta(ieta(id_),ieta(relativeto))
    br.Xtal_iphi[i] = iphi(id_) if relativeto==None else diphi(iphi(id_),iphi(relativeto))
    br.nXtal[0] += 1
def FillClus(clus,i):
    br.Clus_pt[i] = clus.pt()
    br.Clus_eta[i] = clus.eta()
    br.Clus_phi[i] = clus.phi()
    br.nClus[0] += 1
def FillXtalInClus(xtal,i,seedid):
    br.XtalInClus_pt[i] = xtal.second
    id_ = (xtal.first)()
    br.XtalInClus_dieta[i] = dieta(ieta(id_),ieta(seedid))
    br.XtalInClus_diphi[i] = diphi(iphi(id_),iphi(seedid))
    br.nXtalInClus[0] += 1
    

for iev,ev in enumerate(events):
    if iev%100==0: print 'Processing event %d'%(iev+1,)
    ev.getByLabel(("L1EGammaCrystalsProducer","L1EGXtals"),h_xtals)
    ev.getByLabel(("l1tPFEcalProducerFromL1EGCrystalClusters","FilteredL1EGXtalClusterNoCuts"),h_clusters)
    ev.getByLabel(("genParticles",""),h_genparticles)
    xtals = h_xtals.product()
    clusters = h_clusters.product()
    genparticles = h_genparticles.product()

    if not doByParticle:
        Reset()
        j = 0
        for i in xrange(min(br.max_nGenPart,len(genparticles))):
            part = genparticles[i]
            if len(genparticles)>2 and part.status()!=23: continue
            FillGenPart(part,j)
            j += 1
        for i in xrange(min(br.max_nClus,len(clusters))):
            FillClus(clusters[i],i)
        for i in xrange(min(br.max_nXtal,len(xtals))):
            FillXtal(xtals[i],i)
        tree.Fill()
    else:
        for k in xrange(min(br.max_nGenPart,len(genparticles))):
            Reset()
            part = genparticles[k]
            if len(genparticles)>2 and part.status()!=23: continue
            FillGenPart(part,0)
            best = 0
            for i in xrange(len(clusters)):
                clus = clusters[-1-i] 
                if abs(clus.pt()/part.pt()-1)>0.5: continue
                if abs(clus.eta()-part.eta())>0.2: continue
                if abs(dphi(clus.phi(),part.phi()))>0.2: continue
                best = -1-i
            if best!=0:
                clus = clusters[best]
                seedid = (clus.seedCrystal())()
                FillClus(clus,0)
                for i in xrange(clus.GetNXtals()):
                    FillXtalInClus(filter(lambda xtal: (xtal.first)()==clus.GetDetId(i)(), xtals)[0],i,seedid)
                j = 0
                for i in xrange(min(br.max_nXtal,len(xtals))):
                    id_ = (xtals[i].first)()
                    if abs(dieta(ieta(id_),ieta(seedid)))<10 and abs(diphi(iphi(id_),iphi(seedid)))<10:
                        FillXtal(xtals[i],j,seedid)
                        j+=1
            tree.Fill()
            
outf.Write()
outf.Close()
