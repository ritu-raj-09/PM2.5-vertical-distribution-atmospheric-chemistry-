import netCDF4 as nc
import numpy as np
from os import listdir as l
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

print('processing please wait')


for i in range(1,13):
      year = "2021"
      if i<10 :
          month = "0"+str(i)
      else:
            month = str(i)
      color = ['#8B8378', '#00FFFF', '#458B74', '#E3CF57', '#00008B',  '#000000', '#8A2BE2', '#9C661F', '#FFD39B', '#FF6103', '#7FFF00', '#FF1493', '#E066FF', '#00FA9A', '#8B8B00','#4B0082']
      MAH = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG','SEP','OCT','NOV','DEC']

      content = ['time', 'lon', 'lat', 'lev', 'AIRDENS',  'BCPHILIC', 'BCPHOBIC', 'DU001', 'DU002', 'OCPHILIC', 'OCPHOBIC', 'SO4', 'SS001', 'SS002', 'SS003']   # ds.variables.keys()
      # fn = f"verticalprofile/MERRA2_400.inst3_3d_aer_Nv.{year}{month}{date}.SUB.nc"

      #NAINITAL
      #latbounds = [ 29 , 31 ]
      #lonbounds = [ 79 , 81 ]
      
      #KANPUR
      latbounds = [ 26 , 27 ]
      lonbounds = [ 80 , 81 ]
      
      '''#DELHI
      latbounds = [ 28 , 29 ]
      lonbounds = [ 77 , 78 ]

      
      #PATNA
      latbounds = [ 25 , 26 ]
      lonbounds = [ 85 , 86 ]

      
      #KOLKATA
      latbounds = [ 22 , 23 ]
      lonbounds = [ 88 , 89 ]

      

      #VARANASI
      latbounds = [ 25 , 26 ]
      lonbounds = [ 82 , 83 ]


      
      #BALLIA
      latbounds = [ 25 , 27 ]
      lonbounds = [ 84 , 86 ]

      
      #DHAKA
      latbounds = [ 23 , 24 ]
      lonbounds = [ 90 , 91 ]

      
      #KARACHI
      latbounds = [ 24 , 26 ]
      lonbounds = [ 67 , 69 ]

      #JAIPUR
      latbounds = [ 26 , 27 ]
      lonbounds = [ 75 , 76 ]'''
      


#opening file
      
      if i in (6,7,8,9):
            files = ["verticalprofile/"+x for x in l("verticalprofile/") if f"MERRA2_401.inst3_3d_aer_Nv.{year}{month}" in x]
            k=0
      else:
            files = ["verticalprofile/"+x for x in l("verticalprofile/") if f"MERRA2_400.inst3_3d_aer_Nv.{year}{month}" in x]
            k=0
            
      dayAvgs = []
      for file in files:
          ds = nc.Dataset(file)
          k=k+1
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
          PM25 = np.mean(PM25, axis=0)#DAY AVGs
          PM25 = PM25*(10**9)
          levs=[108.663,208.152,337.5,450,562.5,637.5,725,820,910,985]
          p=PM25[:,0,0]
          ax = plt.plot(p,levs)
          plt.xlim(0, 100)
          if i==1:
                plt.gca().invert_yaxis()
          else:
                pass
          plt.ylim(1000, 100)
          plt.xlabel('CONC.PM2.5(\u03BCg per meter cube)')
          plt.ylabel('Pressure hPa')
          plt.title('kanpur VERTICAL PROFILE (AVERAGE FOR MONTH) '+str(k)+'-'+ MAH[i-1]+ '-2021' )##note: please change the name according to city
          plt.show()
          #plt.savefig(r'C:\Users\tryri\Desktop\cook'+'\\'+'vp'+str(k)+'-'+ MAH[i-1]+'.png')
