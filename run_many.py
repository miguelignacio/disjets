import os
import numpy as np

mycounter = 0
mycounter_submit = 0
for i in range(1,100):
    if ("stat"+str(i)+"_losshistory_Rapgap_nominal_b4_i4_s1.json") in os.listdir("/clusterfs/ml4hep/bpnachman/H1/"):
        mycounter += 1
        pass
    elif (mycounter_submit < 10):
        myrand = np.random.randint(0,10000)
        print("python strap.py "+str(i)+" "+str(myrand)+" &> mylog"+str(i)+".txt &")
        os.system("python strap.py "+str(i)+" "+str(myrand)+" &> mylog"+str(i)+".txt &")
        mycounter_submit += 1
    else:
        continue

print(mycounter,mycounter_submit)
