import netCDF4 as nc
import numpy as np
from os import listdir as l
import matplotlib.pyplot as plt
from tqdm import tqdm
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs 
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import pandas as pd
from tabulate import tabulate

fig = plt.figure(figsize=(30,20))
#ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([60.0, 100.0, 0.0, 40.0], crs=ccrs.PlateCarree())
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0, color='gray', alpha=0.5,linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE,linewidth=0.9)
ax.add_feature(cfeature.BORDERS)#, linestyle=':')
#ax.add_feature(cfeature.LAKES, alpha=0.5)
#ax.add_feature(cfeature.RIVERS)
#l=['1','2','3','4','5','6','7','8','9','10']
lats=[26.9124 ,28.7041,   29.3919, 26.4499,25.3176, 25.7631,  25.5941, 22.8046, 22.5726, 23.8103]
lons=[75.7873, 77.1025,   79.4542, 80.3319, 82.9739, 84.1496,  85.1376, 86.2029,88.3639, 90.4125]
d=[48.6130,74.4700,46.0700,65.9200, 65.2000,71.0200,72.4400,59.3600,59.5900,55.3000]
cm = plt.cm.get_cmap('rainbow')
xy = range(10)
z = xy
sc = plt.scatter(lons, lats, c=d,s=500, vmin=40, vmax=75, cmap=cm)

for i in range(0,10):
    plt.text(lons[i]+0.1, lats[i],str(i+1),fontsize=20,horizontalalignment='right',transform=ccrs.Geodetic())

plt.colorbar(sc,location = 'right',extend='both',aspect=35,pad=0.01)#,shrink=0.80)
#plt.show()
NAINITAL=[]
KANPUR=[]
DELHI=[]
PATNA=[]
KOLKATA=[]
VARANASI=[]
BALLIA=[]
DHAKA=[]
JAMSEDHPUR=[]
JAIPUR=[]


print('processing please wait')
site=['NAINITAL','KANPUR','DELHI','PATNA','KOLKATA','VARANASI','BALLIA','DHAKA','JAMSEDHPUR','JAIPUR']
for u in range(0,10):
    if u==0:
        #NAINITAL
        latbounds = [ 29 , 30 ]
        lonbounds = [ 79 , 80 ]
    elif u==1:
        #KANPUR
        latbounds = [ 26 , 27 ]
        lonbounds = [ 80 , 81 ]
    elif u==2:
        #DELHI
        latbounds = [ 28 , 29 ]
        lonbounds = [ 77 , 78 ]
    elif u==3:
        #PATNA
        latbounds = [ 25 , 26 ]
        lonbounds = [ 85 , 86 ]
    elif u==4:
        #KOLKATA
        latbounds = [ 22 , 23 ]
        lonbounds = [ 88 , 89 ]
    elif u==5:
        #VARANASI
        latbounds = [ 25 , 26 ]
        lonbounds = [ 82 , 83 ]
    elif u==6:
        #BALLIA
        latbounds = [ 25 , 27 ]
        lonbounds = [ 84 , 86 ]
    elif u==7:
        #DHAKA
        latbounds = [ 23 , 24 ]
        lonbounds = [ 89 , 90 ]
    elif u==8:
        #JAMSEDHPUR
        latbounds = [ 22 , 23 ]
        lonbounds = [ 86 , 87 ]
    elif u==9:
        #JAIPUR
        latbounds = [ 26 , 27 ]
        lonbounds = [ 75 , 76 ]
    else:
        pass
        
    for i in tqdm(range(1,13)):
          year = "2021"
          if i<10 :
              month = "0"+str(i)
          else:
                month = str(i)
          color = ['#8B8378', '#00FFFF', '#458B74', '#E3CF57', '#00008B',  '#000000', '#8A2BE2', '#9C661F', '#FFD39B', '#FF6103', '#7FFF00', '#FF1493', '#E066FF', '#00FA9A', '#8B8B00','#4B0082']
          MAH = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG','SEP','OCT','NOV','DEC']
    
          content = ['time', 'lon', 'lat', 'lev', 'AIRDENS',  'BCPHILIC', 'BCPHOBIC', 'DU001', 'DU002', 'OCPHILIC', 'OCPHOBIC', 'SO4', 'SS001', 'SS002', 'SS003']   # ds.variables.keys()
          # fn = f"verticalprofile/MERRA2_400.inst3_3d_aer_Nv.{year}{month}{date}.SUB.nc"
         #THESE ARE BOUNDED LAT AND LONS OF PATICULAR LOCATIO
          
          
    
    
    #opening file
          
          if i in (6,7,8,9):
                files = ["vp_2021/"+x for x in l("vp_2021/") if f"MERRA2_401.inst3_3d_aer_Nv.{year}{month}" in x]
          else:
                files = ["vp_2021/"+x for x in l("vp_2021/") if f"MERRA2_400.inst3_3d_aer_Nv.{year}{month}" in x]
                
          dayAvgs = []
          for file in files:
              ds = nc.Dataset(file)
              #extracting paticular lat lon from india map
              lons = np.array(ds['lon'][:])
              lats = np.array(ds['lat'][:])
              # latitude lower and upper index
              latli = np.argmin( np.abs( lats - latbounds[0] ) )
              latui = np.argmin( np.abs( lats - latbounds[1] ) )
              # longitude lower and upper index
              lonli = np.argmin( np.abs( lons - lonbounds[0] ) )
              lonui = np.argmin( np.abs( lons - lonbounds[1] ) )
    
    
              #extracting data with paticular  lat lon from india map
              SO4 = np.array(ds["SO4"][:, : , latli:latui , lonli:lonui ])
              AIRDENS = np.array(ds['AIRDENS'][:, : , latli:latui , lonli:lonui ])
              BCPHILIC = np.array(ds['BCPHILIC'][:, : , latli:latui , lonli:lonui ])
              BCPHOBIC = np.array(ds['BCPHOBIC'][:, : , latli:latui , lonli:lonui ])
              OCPHILIC = np.array(ds['OCPHILIC'][:, : , latli:latui , lonli:lonui ])
              OCPHOBIC = np.array(ds['OCPHOBIC'][:, : , latli:latui , lonli:lonui ])
              DU001 = np.array(ds['DU001'][:, : , latli:latui , lonli:lonui ])
              DU002 = np.array(ds['DU002'][:, : , latli:latui , lonli:lonui ])
              SS001 = np.array(ds['SS001'][:, : , latli:latui , lonli:lonui ])
              SS002 = np.array(ds['SS002'][:, : , latli:latui , lonli:lonui ])
              SS003 = np.array(ds['SS003'][:, : , latli:latui , lonli:lonui ])
              levs = np.array(ds['lev'][:])
    
              PM25 = (1.375*(SO4) + 1.6*( OCPHILIC + OCPHOBIC) +  ( BCPHILIC + BCPHOBIC ) + ( DU001 + (DU002)*0.38) + (SS001) + (SS002) + ((SS003)*0.83))*(AIRDENS)
              dayAvgs.append(np.mean(PM25, axis=0))#DAYAVG
    
          dayAvgs = np.array(dayAvgs)
          moAvgs = np.mean(dayAvgs, axis = 0)
          moAvgs = moAvgs*(10**9)#MONTHLY AVG
          #print(np.max(moAvgs),np.min(moAvgs))
          
          levs=[0.0100,0.0200,0.0327,0.0476,0.0660,0.0893,0.1197,0.1595,0.2113,0.2785,0.3650,0.4758,0.6168,0.7951,1.0194,1.3005,1.6508,2.0850,2.6202,3.2764,4.0766,5.0468,6.2168,7.6198,9.2929,11.2769,13.6434,16.4571,19.7916,23.7304,28.3678,33.8100,40.1754,47.6439,56.3879,66.6034,78.5123,92.3657,108.6630,127.8370,150.3930,176.9300,208.1520,244.8750,288.0830,337.5000,375.0000,412.5000,450.0000,487.5000,525.0000,562.5000,600.0000,637.5000,675.0000,700.0000,725.0000,750.0000,775.0000,800.0000,820.0000,835.0000,850.0000,865.0000,880.0000,895.0000,910.0000,925.0000,940.0000,955.0000,970.0000,985.0000]
          p=moAvgs[:,0,0]
          if u==0:
              NAINITAL.append(p)
          elif u==1:
              KANPUR.append(p)
          elif u==2:
              DELHI.append(p)
          elif u==3:
              PATNA.append(p)
          elif u==4:
              KOLKATA.append(p)
          elif u==5:
              VARANASI.append(p)
          elif u==6:
              BALLIA.append(p)
          elif u==7:
              DHAKA.append(p)
          elif u==8:
              JAMSEDHPUR.append(p)
          elif u==9:
              JAIPUR.append(p)
          else:
              pass
          

NAINITAL=np.array(NAINITAL)
KANPUR=np.array(KANPUR)
DELHI=np.array(DELHI)
PATNA=np.array(PATNA)
KOLKATA=np.array(KOLKATA)
VARANASI=np.array(VARANASI)
BALLIA=np.array(BALLIA)
DHAKA=np.array(DHAKA)
JAMSEDHPUR=np.array(JAMSEDHPUR)
JAIPUR=np.array(JAIPUR)
head=[0.0100,0.0200,0.0327,0.0476,0.0660,0.0893,0.1197,0.1595,0.2113,0.2785,0.3650,0.4758,0.6168,0.7951,1.0194,1.3005,1.6508,2.0850,2.6202,3.2764,4.0766,5.0468,6.2168,7.6198,9.2929,11.2769,13.6434,16.4571,19.7916,23.7304,28.3678,33.8100,40.1754,47.6439,56.3879,66.6034,78.5123,92.3657,108.6630,127.8370,150.3930,176.9300,208.1520,244.8750,288.0830,337.5000,375.0000,412.5000,450.0000,487.5000,525.0000,562.5000,600.0000,637.5000,675.0000,700.0000,725.0000,750.0000,775.0000,800.0000,820.0000,835.0000,850.0000,865.0000,880.0000,895.0000,910.0000,925.0000,940.0000,955.0000,970.0000,985.0000]
n=tabulate(NAINITAL, headers=head)
k=tabulate(KANPUR, headers=head)
d=tabulate(DELHI, headers=head)
p=tabulate(PATNA, headers=head)
v=tabulate(VARANASI, headers=head)
b=tabulate(BALLIA, headers=head)
dh=tabulate(DHAKA, headers=head)
jh=tabulate(JAMSEDHPUR, headers=head)
j=tabulate(JAIPUR, headers=head)

with open("NAINITAL.txt", "a") as o:
    o.write(n)
    
with open("KANPUR.txt", "a") as o:
    o.write(k)
    
with open("DELHI.txt", "a") as o:
    o.write(d)

with open("PATNA.txt", "a") as o:
    o.write(p)

with open("VARANASI.txt", "a") as o:
    o.write(v)
    
with open("BALLIA.txt", "a") as o:
    o.write(b)
    
with open("DHAKA.txt", "a") as o:
    o.write(dh)
    
with open("JAMSEDHPUR.txt", "a") as o:
    o.write(jh)

with open("JAIPUR.txt", "a") as o:
    o.write(j)

PL=tabulate(KOLKATA, headers=head)
with open("KOLKATA.txt", "a") as o:
    o.write(PL)
 
        


          
