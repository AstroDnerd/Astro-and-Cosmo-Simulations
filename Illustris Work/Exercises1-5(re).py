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

#Functions

def exe1_r_compare():
    basePath = './Illustris-3 L75n455FP/postprocessing/halocatalogs'
    HCpath = basePath+'/groups_135.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Fields
        fields = list(f['Group'].keys())
        #Storing data specific to a Field
        field_choice = 'GroupPos'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        GroupPos = (xnum)/h
        xloc = GroupPos[0]
        yloc = GroupPos[1]
        field_choice = 'Group_R_Crit200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        R_Crit200 = (xnum)/h
        field_choice = 'Group_R_Crit500'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        R_Crit500 = (xnum)/h
        field_choice = 'Group_R_Mean200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        R_Mean200 = (xnum)/h
        field_choice = 'Group_R_TopHat200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        R_TH200 = (xnum)/h
        lim = 5
        xloc = xloc[:lim]
        yloc = yloc[:lim]
        plt.scatter(xloc,yloc,s=2,c='blue')
        plt.xlabel('x-coordinate (ckpc)')
        plt.ylabel('y-coordinate (ckpc)')
        plt.title('2-D plot of Group Positions')
        fig = plt.gcf()
        ax = fig.gca()
        RType = ['R_Crit200','R_Crit500','R_Mean200','R_TopHat200']
        for i in range(lim):
            circle1 = plt.Circle((xloc[i],yloc[i]),R_Crit200[i],fill = False,alpha=1,
                                 color='green')
            ax.add_artist(circle1)
            circle2 = plt.Circle((xloc[i],yloc[i]),R_Crit500[i],fill = False,alpha=1,
                                 color='black')
            ax.add_artist(circle2)
            circle3 = plt.Circle((xloc[i],yloc[i]),R_Mean200[i],fill = False,alpha=1,
                                 color='red')
            ax.add_artist(circle3)
            circle4 = plt.Circle((xloc[i],yloc[i]),R_TH200[i],fill = False,alpha=1,
                                 color='yellow')
            ax.add_artist(circle4)
        plt.legend((circle1,circle2,circle3,circle4),
                   ('R_Crit200','R_Crit500','R_Mean200','R_TopHat200'),loc='upper left')
        plt.show()

def exe1_mass():
    basePath = './Illustris-3 L75n455FP/postprocessing/halocatalogs'
    HCpath = basePath+'/groups_135.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Fields
        fields = list(f['Group'].keys())
        #Radii
        field_choice = 'Group_R_Crit200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum[xnum!= 0.0]         #removing zeroes
        R_Crit200 = np.log10(1e-16+ (xnum)/h)
        field_choice = 'Group_R_Crit500'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum[xnum!= 0.0]
        R_Crit500 = np.log10(1e-16+(xnum)/h)
        field_choice = 'Group_R_Mean200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum[xnum!= 0.0]
        R_Mean200 = np.log10(1e-16+(xnum)/h)
        field_choice = 'Group_R_TopHat200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum[xnum!= 0.0]
        R_TH200 = np.log10(1e-16+(xnum)/h)
        #Masses
        field_choice = 'Group_M_Crit200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum[xnum!= 0.0]
        M_Crit200 = np.log10(1e-16+ 1e10*(xnum)/h)
        field_choice = 'Group_M_Crit500'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum[xnum!= 0.0]
        M_Crit500 = np.log10(1e-16+ 1e10*(xnum)/h)
        field_choice = 'Group_M_Mean200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum[xnum!= 0.0]
        M_Mean200 = np.log10(1e-16+ 1e10*(xnum)/h)
        field_choice = 'Group_M_TopHat200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum[xnum!= 0.0]
        M_TH200 = np.log10(1e-16+ 1e10*(xnum)/h)
        plt.scatter(R_Crit200,M_Crit200,s=1,c='green',label='Crit_200')
        plt.scatter(R_Crit500,M_Crit500,s=1,c='black',label='Crit_500')
        plt.scatter(R_Mean200,M_Mean200,s=1,c='red',label='Mean_200')
        plt.scatter(R_TH200,M_TH200,s=1,c='yellow',label='Crit_$\Delta_c$')
        plt.xlabel("Radii of the halo [log10(ckpc)]")
        plt.ylabel("Mass of the halo [log10($M_\odot$)]")
        plt.legend(loc='upper left')
        plt.show()

def exe1_DM_vs_stellar():
    basePath = './Illustris-3 L75n455FP/postprocessing/halocatalogs'
    HCpath = basePath+'/groups_135.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Fields
        fields = list(f['Group'].keys())
        #Position
        field_choice = 'SubhaloPos'
        shape = f['Subhalo'][field_choice].shape[0]
        result['field'] = f['Subhalo'][field_choice][0:shape]
        xnum = result['field']
        Subpos = (xnum)/h
        xloc = Subpos[0]
        yloc = Subpos[1]
        print(len(xloc))
        print(len(yloc))
        dist = (xloc**2+yloc**2)**0.5
        #Mass
        field_choice = 'SubhaloMassType'
        shape = f['Subhalo'][field_choice].shape[0]
        result['field'] = f['Subhalo'][field_choice][0:shape]
        xnum = result['field']
        Mass = 1e10*(xnum)/h
        Stellar = Mass[4]
        #Stellar = Stellar[Stellar!= 0.0]
        Stellar = np.log10(1e10*(Stellar)/h)
        DM = Mass[1]
        #DM = DM[DM!= 0.0]
        DM = np.log10(1e10*(DM)/h)
        #subplot for Dark matter and stellar mass separately
        plt.subplot(2,1,1)
        plt.scatter(dist,Stellar,s=1,c='red',label='Stellar Mass')
        plt.scatter(dist,DM,s=1,c='blue',label='Dark Matter Mass')
        plt.xlabel("Distance from origin of periodic box (ckpc)")
        plt.ylabel("Mass [log10($M_\odot$)]")
        plt.legend(loc='upper left')
        #subplot for ratio between Dark Matter and Stellar mass
        plt.subplot(2,1,2)
        plt.scatter(dist,DM-Stellar,s=1,c='cyan',label='DM mass/Stellar Mass')
        plt.xlabel("Distance from origin of periodic box (ckpc)")
        plt.ylabel("Mass [log10($M_\odot$)]")
        plt.legend(loc='upper left')
        plt.show()

def exe2_re():
    basePath = './Illustris-3 L75n455FP/postprocessing/halocatalogs'
    HCpath = basePath+'/groups_135.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Fields
        #Masshalo
        field_choice = 'GroupMass'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        Mhalo = np.log10(1e10*(xnum)/h)
        #Mgashalo
        field_choice = 'GroupMassType'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        Mgashalo = np.log10(1e10*(xnum[0])/h)
        #MassSubgas
        field_choice = 'SubhaloMassType'
        shape = f['Subhalo'][field_choice].shape[0]
        result['field'] = f['Subhalo'][field_choice][0:shape]
        xnum = result['field']
        Mgassub = np.log10(1e10*(xnum[0])/h)
        #MassSub
        field_choice = 'SubhaloMass'
        shape = f['Subhalo'][field_choice].shape[0]
        result['field'] = f['Subhalo'][field_choice][0:shape]
        xnum = result['field']
        Msub = np.log10(1e10*(xnum)/h)
        #Subplot for Halo/Group
        plt.subplot(2,1,1)
        plt.scatter(Mhalo,Mgashalo,s=1)
        plt.xlabel("Mass of halo [log10($M_\odot$)]")
        plt.ylabel("Mass of gas in halo [log10($M_\odot$)]")
        #Subplot for Subhalo
        plt.subplot(2,1,2)
        plt.scatter(Msub,Mgassub,s=1)
        plt.xlabel("Mass of Subhalo [log10($M_\odot$)]")
        plt.ylabel("Mass of gas in Subhalo [log10($M_\odot$)]")
        plt.show()

def exeA1_re():
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
        Gposx = xnum[0][0]
        Gposy = xnum[1][0]
        Gposz = xnum[2][0]
        field_choice = 'Group_R_Crit200'
        shape = f['Group'][field_choice].shape[0]
        result['field'] = f['Group'][field_choice][0:shape]
        xnum = result['field']
        Radius = xnum[0]/h
        fig = plt.gcf()
        ax = fig.gca()
        circle1 = plt.Circle((Gposx,Gposy),Radius,fill = False,alpha=1,color='red')
        ax.add_artist(circle1)
        plt.xlabel('x-coordinate (ckpc)')
        plt.ylabel('y-coordinate (ckpc)')
        plt.xlim(-1500+Gposx,1500+Gposx)
        plt.ylim(-1500+Gposy,1500+Gposy)
        #Subhalos
        field_choice = 'SubhaloPos'
        shape = f['Subhalo'][field_choice].shape[0]
        result['field'] = f['Subhalo'][field_choice][0:shape]
        xnum = result['field']
        xnum = xnum/h
        Sx = xnum[0]
        Sy = xnum[1]
        Sz = xnum[2]
        Sd = ((Sx-Gposx)**2+(Sy-Gposy)**2+(Sz-Gposz)**2)**0.5
        Sx = Sx[Sd<=Radius]     #Selects those x coordinates of SubhaloPos which lie inside
        Sy = Sy[Sd<=Radius]
        Sz = Sz[Sd<=Radius]
        plt.scatter(Sx,Sy,s=4,label='Subhalos')
        ax.scatter(Gposx,Gposy,c='black',marker="+",label='Halo Center')
        plt.legend(loc='upper right')
        plt.show()
        plt.subplot(2,1,1)
        plt.hist(Sd,bins=20)
        plt.xlabel('Distance from the center of the Group (ckpc)')
        plt.ylabel('Number of Subhalos')
        plt.subplot(2,1,2)
        plt.hist(Sd,bins=20,cumulative='True')
        plt.xlabel('Distance from the center of the Group (ckpc)')
        plt.ylabel('Total Number of Subhalos')
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
            plt.subplot(2,2,i+1)
            s = 0.1
            if i==3:
                s=4
            plt.scatter(xnum,ynum,s=s,c=color[i],label = PTType[i], alpha = 0.5)
            plt.xlabel('x-coordinate (ckpc)')
            plt.ylabel('y-coordinate (ckpc)')
            plt.legend(loc='upper right')
        plt.show()

def exe4_re():
    #Opening halocatalog
    #plt.style.use('dark_background')
    basePath = './Illustris-3 L75n455FP/postprocessing/halocatalogs'
    HCpath = basePath+'/groups_135.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Header
        header = dict(f['Header'].attrs.items())
        result['count'] = f['Header']['Ngroups_Total'][0]
        #Plotting the Histogram
        arr = np.array(f['Subhalo']['SubhaloPos'])
        arr = arr/h
        xnum = arr[0][:]
        ynum = arr[1][:]
        xmin = np.min(xnum)
        xmax = np.max(xnum)
        ymin = np.min(ynum)
        ymax = np.max(ynum)
        #Displacement from original by xmax-80000 ckpc in x and +19000ckpc in y
        halox = np.mod(xnum+(xmax-80000),xmax)
        haloy = np.mod(ynum+19000,ymax)
    #Opening Snap
    basePath = './Illustris-3 L75n455FP/Output/Snapdir_135'
    HCpath = basePath+'/snap_135.0.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Header
        header = dict(f['Header'].attrs.items())
        result['count'] = header['NumPart_ThisFile']
        PT = ['PartType0','PartType1','PartType4','PartType5']
        color = ['orange','blue','red','black']
        PTType = ['Gas','Dark Matter','Stars/Wind','Black Holes']
        for i in range(4):
            arr = np.array(f[PT[i]]['Coordinates'])
            arr = arr/h
            arr = np.transpose(arr)
            xnum = arr[0][:]
            ynum = arr[1][:]
            xmin = np.min(xnum)
            xmax = np.max(xnum)
            ymin = np.min(ynum)
            ymax = np.max(ynum)
            xnum = np.mod(xnum+(xmax-80000),xmax)
            ynum = np.mod(ynum+19000,ymax)
            plt.subplot(2,2,i+1)
            plt.scatter(halox,haloy,s=0.1,label = 'Halo')
            plt.scatter(xnum,ynum,s=0.1,c=color[i],label = PTType[i])
            plt.xlabel('x-coordinate (ckpc)')
            plt.ylabel('y-coordinate (ckpc)')
            plt.legend(loc='upper right')
        plt.show()

    

    
exe4_re()
