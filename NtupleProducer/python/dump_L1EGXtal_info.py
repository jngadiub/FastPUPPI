import ROOT
from DataFormats.FWLite import Events, Handle

events = Events('inputs_17D.root')

h_xtals = Handle('vector<pair<DetId,float> >')
h_clusters = Handle('vector<l1slhc::L1EGCrystalCluster>')
h_particles = Handle('vector<l1tpf::Particle>')
for ev in events:
    ev.getByLabel(("L1EGammaCrystalsProducer","L1EGXtals"),h_xtals)
    ev.getByLabel(("l1tPFEcalProducerFromL1EGCrystalClusters","FilteredL1EGXtalClusterNoCuts"),h_clusters)
    ev.getByLabel(("l1tPFEcalProducerFromL1EGCrystalClusters",""),h_particles)
    xtals = h_xtals.product()
    clusters = h_clusters.product()
    particles = h_particles.product()
    print len(xtals),'crystals,',len(clusters),'clusters,',len(particles),'particles'
    print 'crystals:',[((xtal.first)(),xtal.second) for xtal in xtals]
    print 'clusters:',[(clus.GetNXtals(),clus.pt(),clus.eta(),clus.phi(),tuple([clus.GetDetId(i)() for i in xrange(clus.GetNXtals())])) for clus in clusters]
    print 'particles:',[(part.pt(),part.eta(),part.phi()) for part in particles]
    print '---'
