import illustris_python as il
import matplotlib.pyplot as plt
import h5py
import numpy as np
h = 0.704
#This File Contains Exercise 0 divided into its three subparts
#as different fucntions. Please Call the respective function at the end
#of the file to see the output

#Functions
def exe0a():
    #Opening halocatalog
    basePath = './Illustris-3 L75n455FP/postprocessing/halocatalogs'
    HCpath = basePath+'/groups_135.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Header
        header = dict(f['Header'].attrs.items())
        result['count'] = f['Header']['Ngroups_Total'][0]
        #Storing Fields
        fields = list(f['Subhalo'].keys())
        print("Fields in Subhalo are:")
        print(fields)
        #Storing data specific to a Field
        field_choice = input("Enter the Field name (for Histogram):")
        shape = f['Subhalo'][field_choice].shape[0]
        result['field'] = f['Subhalo'][field_choice][0:shape]
        #Plotting the Histogram
        xnum = result['field']+1e-16
        unit = ''
        if field_choice == 'SubhaloBHMdot':
            xnum = xnum *(1e10/h)*(0.978/h)
            unit = '[$M_\odot Gyr$]'
        elif field_choice in ['SubhaloBHMass','SubhaloMass','SubhaloMassInHalfRad',
                              'SubhaloMassInHalfRadType','SubhaloMassInMaxRad',
                              'SubhaloMassInMaxRadType','SubhaloMassInRad',
                              'SubhaloMassInRadType','SubhaloMassType',
                              'SubhaloStellarPhotometricsMassInRad','SubhaloWindMass']:
            xnum = xnum*(1e10/h)
            unit = '[$M_\odot$]'
        elif field_choice in ['SubhaloCM','SubhaloHalfmassRad','SubhaloHalfmassRadType',
                              'SubhaloPos','SubhaloStellarPhotometricsRad','SubhaloVmaxRad']:
            xnum = xnum/h
            unit = '[ckpc]'
        xnum = np.log10(xnum)
        print("The number of values in the field are: "+str(result['count']))
        print("The max and min values in the field are: "+str(max(xnum))+" and "
              + str(min(xnum)) + "respectively.")
        num_bins = int(input("The number of bins in the histogram: "))
        n,bins,patches = plt.hist(xnum,num_bins,facecolor='blue',alpha=0.5)
        for i in range(num_bins):
            plt.text(bins[i],n[i],str(n[i]))
        plt.xlabel(field_choice+unit+'$(log_{10})$')
        plt.ylabel('Frequency')
        plt.title('Histogram plot of '+ field_choice)
        plt.show()
        
def exe0b():
    #Opening Snap
    basePath = './Illustris-3 L75n455FP/Output/Snapdir_135'
    HCpath = basePath+'/snap_135.0.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Header
        header = dict(f['Header'].attrs.items())
        result['count'] = header['NumPart_ThisFile']
        Ptype= '''Particle Types in Snap are:
PartType0: Gas
PartType1: Dark Matter
PartType3: Tracer Particles
PartType4: Stars/Wind Particles
PartType5: Black Holes'''
        print(Ptype)
        PT = int(input("Enter Particle Type (integer):"))
        gpname = 'PartType'+str(PT)
        #Storing Fields
        fields = list(f[gpname].keys())
        print("Fields in "+ gpname+ " are:")
        print(fields)
        #Storing data specific to a Field
        field_choice = input("Enter the Field name (for Histogram):")
        shape = f[gpname][field_choice].shape[0]
        result['field'] = f[gpname][field_choice][0:shape]
        #Plotting the Histogram 
        xnum = result['field']+1e-16
        unit = ''
        if field_choice in ['Density','SubfindDensity','BH_Density']:
            xnum = xnum *(1e10/h)*(1/(h**3))
            unit = '[$M_\odot Gyr^3$]'
        elif field_choice in ['Masses','GFM_InitialMass','BH_CumMassGrowth_QM',
                              'BH_Mass','BH_Mass_bubbles','BH_Mass_ini','HostHaloMass']:
            xnum = xnum*(1e10/h)
            unit = '[$M_\odot$]'
        elif field_choice in ['Coordinates','SmoothingLength','SubfindHsml','BH_Hsml']:
            xnum = xnum/h
            unit = '[ckpc]'
        elif field_choice in ['Volume']:
            xnum = xnum/(h**3)
            unit = '[$ckpc^3$]'
        elif field_choice == 'BH_CumEgyInjection_QM':
            xnum = xnum*(1e10.h)*(1/(h**2))*((0.978/h)**2)
            unit = '[$M_\odot ckpc^2 Gyr^2$]'
        elif field_choice == 'BH_Pressure':
            xnum = xnum*(1e10.h)*(1/h)*((0.978/h)**2)
            unit = '[$M_\odot ckpc Gyr^2$]'
        elif field_choice == 'BH_Mdot':
            xnum = xnum*(1e10.h)*(0.978/h)
            unit = '[$M_\odot Gyr$]'
        xnum = np.log10(xnum)
        print("The number of values in the field are: "+str(result['count']))
        print("The max and min values in the field are: "+str(max(xnum))+" and "
              + str(min(xnum)) + "respectively.")
        num_bins = int(input("The number of bins in the histogram: "))
        n,bins,patches = plt.hist(xnum,num_bins,facecolor='blue',alpha=0.5)
        for i in range(num_bins):
            plt.text(bins[i],n[i],str(n[i]))
        plt.xlabel(field_choice + "("+gpname+")" +unit+'$(log_{10})$' )
        plt.ylabel('Frequency')
        plt.title('Histogram plot of '+ field_choice + "("+gpname+")")
        plt.show()

def exe0c():
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
        plt.scatter(xnum,ynum,s=1)
        plt.xlabel('x-coordinate (ckpc)')
        plt.ylabel('y-coordinate (ckpc)')
        plt.title('2-D plot of Halos')
        plt.show()
    #Opening Snap
    basePath = './Illustris-3 L75n455FP/Output/Snapdir_135'
    HCpath = basePath+'/snap_135.0.hdf5'
    result = {}
    with h5py.File(HCpath,'r') as f:
        #Storing Header
        header = dict(f['Header'].attrs.items())
        result['count'] = header['NumPart_ThisFile']
        Ptype= '''Particle Types in Snap are:
PartType0: Gas
PartType1: Dark Matter
PartType3: Tracer Particles
PartType4: Stars/Wind Particles
PartType5: Black Holes'''
        print(Ptype)
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
        plt.title('2-D plot of All Particles')
        plt.legend(loc='upper right')
        plt.show()
            
#Function Call
exe0a()
exe0b()
exe0c()

