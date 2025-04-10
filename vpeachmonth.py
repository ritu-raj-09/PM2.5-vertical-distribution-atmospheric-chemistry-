#each month diffrent plot 


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
      #color = ['#8B8378', '#00FFFF', '#458B74', '#E3CF57', '#00008B',  '#000000', '#8A2BE2', '#9C661F', '#FFD39B', '#FF6103', '#7FFF00', '#FF1493', '#E066FF', '#00FA9A', '#8B8B00','#4B0082']
      MAH = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG','SEP','OCT','NOV','DEC']

      content = ['time', 'lon', 'lat', 'lev', 'AIRDENS',  'BCPHILIC', 'BCPHOBIC', 'DU001', 'DU002', 'OCPHILIC', 'OCPHOBIC', 'SO4', 'SS001', 'SS002', 'SS003']   # ds.variables.keys()
      # fn = f"verticalprofile/MERRA2_400.inst3_3d_aer_Nv.{year}{month}{date}.SUB.nc"

      #NAINITAL
      latbounds = [ 20 , 21 ]
      lonbounds = [ 85 , 86 ]
      
      '''#KANPUR
      latbounds = [ 26 , 27 ]
      lonbounds = [ 80 , 81 ]
      
      #DELHI
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
          dayAvgs.append(np.mean(PM25, axis=0))

      dayAvgs = np.array(dayAvgs)
      moAvgs = np.mean(dayAvgs, axis = 0)
      moAvgs = moAvgs*(10**9)#print(moAvgs)
      levs=[0.0100,0.0200,0.0327,0.0476,0.0660,0.0893,0.1197,0.1595,0.2113,0.2785,0.3650,0.4758,0.6168,0.7951,1.0194,1.3005,1.6508,2.0850,2.6202,3.2764,4.0766,5.0468,6.2168,7.6198,9.2929,11.2769,13.6434,16.4571,19.7916,23.7304,28.3678,33.8100,40.1754,47.6439,56.3879,66.6034,78.5123,92.3657,108.6630,127.8370,150.3930,176.9300,208.1520,244.8750,288.0830,337.5000,375.0000,412.5000,450.0000,487.5000,525.0000,562.5000,600.0000,637.5000,675.0000,700.0000,725.0000,750.0000,775.0000,800.0000,820.0000,835.0000,850.0000,865.0000,880.0000,895.0000,910.0000,925.0000,940.0000,955.0000,970.0000,985.0000]
      p=moAvgs[:,0,0]
      ax = plt.plot(p,levs)
      plt.xlim(0, 100)
      if i==1:
            plt.gca().invert_yaxis()
      else:
            pass
      plt.ylim(1000, 100)
      plt.xlabel('CONC.PM2.5(\u03BC g per meter cube)')
      plt.ylabel('Pressure hPa')
      plt.title('Bhubaneswar Vertical Profile'+MAH[i-1]+ '2021' )#note: please change the name according to city
      plt.show()
