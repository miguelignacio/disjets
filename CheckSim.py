import numpy as np
bins = np.logspace(np.log10(10),np.log10(100),7)

#minimal comparison
from dataloader import get_Dataframe, applyCut, applyCutsJets
mc_nominal = get_Dataframe('/clusterfs/ml4hep/bpnachman/H1/hera/', name='Rapgap', tag='nominal')
print("HERE")
mc_nominal   = applyCutsJets(mc_nominal, isMC=True,verbose=False)
mc_sys7 = get_Dataframe('/clusterfs/ml4hep/bpnachman/H1/hera/', name='Rapgap', tag='sys_7')
mc_sys7   = applyCutsJets(mc_sys7, isMC=True,verbose=False)
binspt = np.logspace(np.log10(10),np.log10(100),7)
mc_nominal_truth = np.array(mc_nominal['pass_truth'])
mc_nominal_pT = mc_nominal[['genjet_pt']].to_numpy()
mc_nominal_wgt = mc_nominal['wgt'].to_numpy()
nom,_= np.histogram(mc_nominal_pT[mc_nominal_truth==1][:,0],bins=bins,weights=mc_nominal_wgt[mc_nominal_truth==1],density=True)
mc_sys7_truth = np.array(mc_sys7['pass_truth'])
mc_sys7_pT = mc_sys7[['genjet_pt']].to_numpy()
mc_sys7_wgt = mc_sys7['wgt'].to_numpy()
sysval7,_= np.histogram(mc_sys7_pT[mc_sys7_truth==1][:,0],bins=bins,weights=mc_sys7_wgt[mc_sys7_truth==1],density=True)
print(sysval7/nom)
