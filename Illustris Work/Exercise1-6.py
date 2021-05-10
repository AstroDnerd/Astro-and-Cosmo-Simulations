import illustris_python as il
import matplotlib.pyplot as plt
import h5py
import numpy as np
h = 0.704
#This File Contains Exercise 1-7 divided into different fucntions.
#Please Call the respective function at the of the file to see the output

#Functions

def exe1():
    #Finding CM and radius of the Subhalos
    basePath = './Illustris-3 L75n455FP/Output/groups_135'
    CM = []     #Stores all the centre of masses of subhalos
    radius = [] #Stores all the radii (R_Crit200)
    for i in range(3):
        HCpath = basePath+'/fof_subhalo_tab_135.'+str(i)+'.hdf5'
        result = {}
        with h5py.File(HCpath,'r') as f:
            #Storing Fields
            fields = list(f['Group'].keys())
            #Storing data specific to a Field
            field_choice = 'GroupCM'
            shape = f['Group'][field_choice].shape[0]
            result['field'] = f['Group'][field_choice][0:shape]
            xnum = result['field']
            SubhaloCM = (xnum)/h
            SubhaloCM = SubhaloCM.tolist()
            field_choice = 'Group_R_Crit200'
            shape = f['Group'][field_choice].shape[0]
            result['field'] = f['Group'][field_choice][0:shape]
            xnum = result['field']
            R_Crit200 = (xnum)/h
            R_Crit200 = R_Crit200.tolist()
            CM+=SubhaloCM
            radius+=R_Crit200
    #Finding the respective particles
    basePath = './Illustris-3 L75n455FP/Output/Snapdir_135'
    Stellar_mass=[]
    for k in range(140):
        print("Loading Group"+str(k))
        mass_snap = []
        for i in range(15):
            HCpath = basePath+'/snap_135.'+str(i)+'.hdf5'
            result = {}
            with h5py.File(HCpath,'r') as f:
                PT = ['PartType0','PartType1','PartType4','PartType5']
                PTType = ['Gas','Dark Matter','Stars/Wind','Black Holes']
                arr = np.array(f[PT[2]]['Coordinates'])
                mass = np.array(f[PT[2]]['Masses'])
                arr = arr/h
                arr1 = arr-CM[k]
                arrf=[]
                for j in range(len(arr1)):
                    rad = np.linalg.norm(arr1[j])
                    if rad>radius[k]:
                        arrf.append(j)
                mass = np.delete(mass,arrf)
                mass_snap.append(np.sum(mass))
            if np.sum(mass_snap)==np.sum(mass_snap[0:-1]) and np.sum(mass_snap)!=0.0:
                break
        Stellar_mass.append(np.sum(mass_snap))
        print("Mass is: "+str(np.sum(mass_snap))+" 1e10*Msun/h")
    rad = radius[0:len(Stellar_mass)]
    Stellar_mass = np.array(Stellar_mass)*1e10/h
    plt.scatter(rad,np.log10(Stellar_mass),s=4)
    plt.xlabel("Radius of the Group [ckpc]")
    plt.ylabel("Stellar Mass contained in the group [log($M_\odot$)]")
    plt.show()

def exe2():
    #Finding CM and radius of the Subhalos
    basePath = './Illustris-3 L75n455FP/Output/groups_135'
    CM = []     #Stores all the centre of masses of subhalos
    radius = [] #Stores all the radii (R_Crit200)
    for i in range(3):
        HCpath = basePath+'/fof_subhalo_tab_135.'+str(i)+'.hdf5'
        result = {}
        with h5py.File(HCpath,'r') as f:
            #Storing Fields
            fields = list(f['Group'].keys())
            #Storing data specific to a Field
            field_choice = 'GroupCM'
            shape = f['Group'][field_choice].shape[0]
            result['field'] = f['Group'][field_choice][0:shape]
            xnum = result['field']
            SubhaloCM = (xnum)/h
            SubhaloCM = SubhaloCM.tolist()
            field_choice = 'Group_R_Crit200'
            shape = f['Group'][field_choice].shape[0]
            result['field'] = f['Group'][field_choice][0:shape]
            xnum = result['field']
            R_Crit200 = (xnum)/h
            R_Crit200 = R_Crit200.tolist()
            CM+=SubhaloCM
            radius+=R_Crit200
    #Finding the respective particles
    basePath = './Illustris-3 L75n455FP/Output/Snapdir_135'
    Gas_mass=[]
    Total_mass = []
    for k in range(50):
        print("Loading Group"+str(k))
        mass_snap = []
        mass_snap_gas = []
        for i in range(15):
            HCpath = basePath+'/snap_135.'+str(i)+'.hdf5'
            result = {}
            with h5py.File(HCpath,'r') as f:
                header = dict(f['Header'].attrs.items())
                result['count'] = header['NumPart_ThisFile']
                PT = ['PartType0','PartType1','PartType4','PartType5']
                PTType = ['Gas','Dark Matter','Stars/Wind','Black Holes']
                arr_gas = np.array(f[PT[0]]['Coordinates'])/h
                mass_gas = np.array(f[PT[0]]['Masses'])
                arr_dm = np.array(f[PT[1]]['Coordinates'])/h
                mass_dm = (header['MassTable'][1])*np.ones((len(f[PT[1]]['ParticleIDs'])))
                arr_stellar = np.array(f[PT[2]]['Coordinates'])/h
                mass_stellar = np.array(f[PT[2]]['Masses'])
                arr_bh = np.array(f[PT[3]]['Coordinates'])/h
                mass_bh = np.array(f[PT[3]]['Masses'])
                arr1_gas = arr_gas-CM[k]
                arr1_dm = arr_dm-CM[k]
                arr1_stellar = arr_stellar-CM[k]
                arr1_bh = arr_bh-CM[k]
                arrf_gas=[]
                arrf_dm=[]
                arrf_stellar=[]
                arrf_bh=[]
                for j in range(len(arr_gas)):
                    rad = np.linalg.norm(arr1_gas[j])
                    if rad>radius[k]:
                        arrf_gas.append(j)
                for j in range(len(arr_dm)):
                    rad = np.linalg.norm(arr1_dm[j])
                    if rad>radius[k]:
                        arrf_dm.append(j)
                for j in range(len(arr_stellar)):
                    rad = np.linalg.norm(arr1_stellar[j])
                    if rad>radius[k]:
                        arrf_stellar.append(j)
                for j in range(len(arr_bh)):
                    rad = np.linalg.norm(arr1_bh[j])
                    if rad>radius[k]:
                        arrf_bh.append(j)
                mass_gas = np.delete(mass_gas,arrf_gas)
                mass_dm = np.delete(mass_dm,arrf_dm)
                mass_stellar = np.delete(mass_stellar,arrf_stellar)
                mass_bh = np.delete(mass_bh,arrf_bh)
                mass_snap_gas.append(np.sum(mass_gas))
                mass_snap.append(np.sum(mass_gas)+np.sum(mass_dm)+np.sum(mass_stellar)
                                 +np.sum(mass_dm))
            if np.sum(mass_snap)==np.sum(mass_snap[0:-1]) and np.sum(mass_snap)!=0.0:
                break
        gm = np.sum(mass_snap_gas)
        tm = np.sum(mass_snap)
        Gas_mass.append(gm)
        Total_mass.append(tm)
        print("Gas Mass: "+str(gm))
        print("Total Mass: "+str(tm))
    Gas_mass = np.array(Gas_mass)*1e10/h
    Total_mass = np.array(Total_mass)*1e10/h
    plt.scatter(np.log10(Total_mass),np.log10(Gas_mass),s=5)
    plt.xlabel("Total Mass contained in the group [log($M_\odot$)]")
    plt.ylabel("Gas Mass contained in the group [log($M_\odot$)]")
    plt.show()


def exeA1():
    #Finding CM and radius of the Central SUBFIND Subhalo
    basePath = './Illustris-3 L75n455FP/Output/groups_135'
    #define validHalos which contains the fof_subhalo number as well as the index of the
    #Group which belongs to the criteria defined in the exercise
    validHalos = []         
    for j in range(32):
        HCpath = basePath+'/fof_subhalo_tab_135.'+str(j)+'.hdf5'
        result = {}
        with h5py.File(HCpath,'r') as f:
            #Storing Fields
            fields = list(f['Group'].keys())
            #Storing data specific to a Field
            field_choice = 'Group_M_Crit200'
            shape = f['Group'][field_choice].shape[0]
            result['field'] = f['Group'][field_choice][0:shape]
            xnum = result['field']
            Mtot = 1e10*(xnum)/h
            for i in range(len(Mtot)):
                if Mtot[i]>10**13 and Mtot[i]<5*10**13:
                    validHalos.append((j,i))
                    print((j,i))
    #This analysis shows that only all such Subhalo belong to
    #fof_subhalo_tab_135_2.hdf5. We can now load that file separately.
    basePath = './Illustris-3 L75n455FP/Output/groups_135'
    HCpath = basePath+'/fof_subhalo_tab_135.2.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Fields
        fields = list(f['Group'].keys())
        #Storing data specific to a Field
        field_choice = 'GroupCM'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        SubhaloCM = (xnum)/h
        field_choice = 'Group_R_Crit200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        R_Crit200 = (xnum)/h
        arr = np.transpose(SubhaloCM)
        x = arr[0][:]
        y = arr[1][:]
        xnum = []
        ynum = []
        radius = []
        fig = plt.gcf()
        fig.set_size_inches(11,11)
        ax = fig.gca()
        for i in range(len(validHalos)):
            xnum.append(x[i])
            ynum.append(y[i])
            radius.append(R_Crit200[i])
            circle = plt.Circle((xnum[i],ynum[i]),radius[i],fill = False,alpha=1)
            ax.add_artist(circle)
        plt.scatter(xnum,ynum,s=1)
        plt.xlabel('x-coordinate (ckpc)')
        plt.ylabel('y-coordinate (ckpc)')
        plt.title('2-D plot of All selected halos and their radii')
        plt.show()
          
def exe3():
    #Finding CM and radius of the Central SUBFIND Subhalo
    basePath = './Illustris-3 L75n455FP/Output/groups_135'
    HCpath = basePath+'/fof_subhalo_tab_135.0.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Fields
        fields = list(f['Group'].keys())
        #Storing data specific to a Field
        field_choice = 'GroupCM'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        SubhaloCM = (xnum[0][:])/h
        field_choice = 'Group_R_Crit200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        R_Crit200 = (xnum[0])/h
    #Finding the respective particles
    basePath = './Illustris-3 L75n455FP/Output/Snapdir_135'
    HCpath = basePath+'/snap_135.0.hdf5'
    result = {}
    print('Loading . . . ')
    with h5py.File(HCpath,'r') as f:
        #Storing Header
        header = dict(f['Header'].attrs.items())
        result['count'] = header['NumPart_ThisFile']
        PT = ['PartType0','PartType1','PartType4','PartType5']
        color = ['yellow','blue','red','black']
        PTType = ['Gas','Dark Matter','Stars/Wind','Black Holes']
        for i in range(4):
            arr = np.array(f[PT[i]]['Coordinates'])
            arr = arr/h
            arr1 = arr-SubhaloCM
            arrf=[]
            for j in range(len(arr1)):
                rad = np.linalg.norm(arr1[j])
                if rad>R_Crit200:
                    arrf.append(j)
            arr = np.delete(arr,arrf,axis=0)
            arr = np.transpose(arr)
            xnum = arr[0][:]
            ynum = arr[1][:]
            plt.scatter(xnum,ynum,s=1,c=color[i],label = PTType[i], alpha = 0.5)
        plt.xlabel('x-coordinate (ckpc)')
        plt.ylabel('y-coordinate (ckpc)')
        plt.title('2-D plot of All Particles of Central SUBFIND halo')
        plt.legend(loc='upper right')
        plt.show()
        
def exe4():
    #Opening halocatalog
    basePath = './Illustris-3 L75n455FP/postprocessing/halocatalogs'
    HCpath = basePath+'/groups_135.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Header
        header = dict(f['Header'].attrs.items())
        result['count'] = f['Header']['Ngroups_Total'][0]
        #Plotting the Histogram
        arr = np.array(f['Subhalo']['SubhaloCM'])
        arr = arr/h
        xnum = arr[0][:]
        ynum = arr[1][:]
        plt.scatter(xnum,ynum,s=1,c = 'green',label = 'Halo')
    #Opening Snap
    basePath = './Illustris-3 L75n455FP/Output/Snapdir_135'
    HCpath = basePath+'/snap_135.0.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Header
        header = dict(f['Header'].attrs.items())
        result['count'] = header['NumPart_ThisFile']
        PT = ['PartType0','PartType1','PartType4','PartType5']
        color = ['yellow','blue','red','black']
        PTType = ['Gas','Dark Matter','Stars/Wind','Black Holes']
        for i in range(4):
            arr = np.array(f[PT[i]]['Coordinates'])
            arr = arr/h
            arr = np.transpose(arr)
            xnum = arr[0][:]
            ynum = arr[1][:]
            plt.scatter(xnum,ynum,s=1,c=color[i],label = PTType[i])
        plt.xlabel('x-coordinate (ckpc)')
        plt.ylabel('y-coordinate (ckpc)')
        plt.title('2-D plot of All Particles + Halo')
        plt.legend(loc='upper right')
        plt.show()
    #It would be fun to overplot the centre of mass and the correspoding R_Crit200

def exe5():
    #Opening Snap
    basePath = './Illustris-3 L75n455FP/Output/Snapdir_135'
    HCpath = basePath+'/snap_135.0.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Header
        header = dict(f['Header'].attrs.items())
        result['count'] = header['NumPart_ThisFile']
        PT = ['PartType0','PartType1','PartType4','PartType5']
        color = ['yellow','blue','red','black']
        PTType = ['Gas','Dark Matter','Stars/Wind','Black Holes']
        for i in range(4):
            arr = np.array(f[PT[i]]['Coordinates'])
            arr = arr/h
            arr = np.transpose(arr)
            xnum = arr[0][:]
            ynum = arr[1][:]
            plt.scatter(xnum,ynum,s=1,c=color[i],label = PTType[i])
        fig = plt.gcf()
        fig.set_size_inches(11,11)
        ax = fig.gca()
    #Finding CM and radius of the SUBFIND Subhalo
    basePath = './Illustris-3 L75n455FP/Output/groups_135'
    HCpath = basePath+'/fof_subhalo_tab_135.0.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Fields
        fields = list(f['Group'].keys())
        #Storing data specific to a Field
        field_choice = 'GroupCM'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        cnt = len(xnum)
        xnum = np.transpose(xnum)
        SubhaloCM = (xnum)/h
        field_choice = 'Group_R_Crit200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        R_Crit200 = (xnum)/h
        rcolor = ['orange','green']
        plt.scatter(SubhaloCM[0][:],SubhaloCM[1][:],s=2,color='green',label = 'Subhalo CM',
                    alpha = 1)
        for j in range(cnt):
            circle = plt.Circle((SubhaloCM[0][j],SubhaloCM[1][j]),R_Crit200[j],
                                color=rcolor[j],fill = False,alpha=1)
            ax.add_artist(circle)
        plt.xlabel('x-coordinate (ckpc)')
        plt.ylabel('y-coordinate (ckpc)')
        plt.title('2-D plot of All Particles + SUBFIND HALO CENTRE')
        plt.legend(loc='upper right')

    plt.show()

def exe6():
    #Opening Snap
    basePath = './Illustris-3 L75n455FP/Output/Snapdir_135'
    mass_array=[]
    for i in range(15):
        HCpath = basePath+'/snap_135.'+str(i)+'.hdf5'
        result = {}
        with h5py.File(HCpath,'r') as f:
            #Storing Header
            header = dict(f['Header'].attrs.items())
            result['count'] = header['NumPart_ThisFile']
            PT = ['PartType0','PartType1','PartType4','PartType5']
            color = ['yellow','blue','red','black']
            PTType = ['Gas','Dark Matter','Stars/Wind','Black Holes']
            mass_gas = 1e10*np.sum(np.array(f[PT[0]]['Masses']))/h
            mass_dm = 1e10*(header['MassTable'][1])*(len(f[PT[1]]['ParticleIDs']))/h
            mass_stellar = 1e10*np.sum(np.array(f[PT[2]]['Masses']))/h
            mass_bh = 1e10*np.sum(np.array(f[PT[3]]['Masses']))/h
            total_mass = mass_gas+mass_dm+mass_stellar+mass_bh
            print('Total mass contained in the snapshot '+str(i)+ ' is: '+
                  "{:e}".format(total_mass)+ ' * Msun')
            mass_array.append(total_mass)
    print('\n')
    print("Total mass from all Snapshots combined is (in multiples of Mass_sun): " +
          "{:e}".format(np.sum(mass_array)))
    #Opening halocatalog
    basePath = './Illustris-3 L75n455FP/postprocessing/halocatalogs'
    HCpath = basePath+'/groups_135.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Header
        header = dict(f['Header'].attrs.items())
        result['count'] = f['Header']['Ngroups_Total'][0]
        grp_mass = 1e10*np.sum(np.array(f['Group']['GroupMass']))/h
        print('Total Group mass in Halocatalog is: ' + str(grp_mass)+ ' * Msun')

            
#Call Functions here
