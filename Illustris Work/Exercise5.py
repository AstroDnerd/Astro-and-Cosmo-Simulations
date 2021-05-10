import illustris_python as il
import matplotlib.pyplot as plt
import h5py
import numpy as np
import matplotlib
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)
h = 0.704
def exe3_FoF_Subfind_group():
    #Opening halocatalog
    basePath = './Illustris-3 L75n455FP/postprocessing/halocatalogs'
    HCpath = basePath+'/groups_135.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Fields
        #Group center
        field_choice = 'GroupPos'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum/h
        Gposx = xnum[0]
        Gposy = xnum[1]
        Gposz = xnum[2]
        #Subhalo Center
        field_choice = 'SubhaloPos'
        shape = f['Subhalo'][field_choice].shape[0]
        result['field'] = f['Subhalo'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum/h
        Sposx = xnum[0]
        Sposy = xnum[1]
        Sposz = xnum[2]
        #FOF Parent
        field_choice = 'SubhaloGrNr'
        shape = f['Subhalo'][field_choice].shape[0]
        result['field'] = f['Subhalo'][field_choice][0:shape]
        GrNr = result['field']
        #Subfind Parent
        field_choice = 'SubhaloParent'
        shape = f['Subhalo'][field_choice].shape[0]
        result['field'] = f['Subhalo'][field_choice][0:shape]
        SubParent = result['field']

        #Choosing the relevant entities
        Gposx = Gposx[GrNr.astype(int)]     #Choses indices in Gposx dictated by GrNr
        Gposy = Gposy[GrNr.astype(int)]
        Gposz = Gposz[GrNr.astype(int)]
        Sposx = Sposx[SubParent.astype(int)]
        Sposy = Sposy[SubParent.astype(int)]
        Sposz = Sposz[SubParent.astype(int)]
        #left-shifting by 10,000ckpc
        Gposx = np.mod(Gposx+10000,np.max(Gposx))
        Sposx = np.mod(Sposx+10000,np.max(Sposx))
        plt.scatter(Gposx,Gposy,s=0.1,label='FoF Parent Group')
        plt.scatter(Sposx,Sposy,s=0.1,c='red',label='Subfind Parent Group')
        plt.xlabel('x-coordinate (ckpc)')
        plt.ylabel('y-coordinate (ckpc)')
        plt.legend(loc='upper left')
        plt.show()


exe3_FoF_Subfind_group()  
