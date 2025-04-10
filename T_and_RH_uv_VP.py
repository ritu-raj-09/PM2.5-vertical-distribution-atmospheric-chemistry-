# estimation of PM2.5 mass concentration using the re-analysis data over the Indian subcontinent.

import netCDF4 as nc
import numpy as np
from os import listdir as l
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import warnings
warnings.filterwarnings('ignore')

latbounds = [ 29 , 30 ]
lonbounds = [ 79 , 80 ]

MAH = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG','SEP','OCT','NOV','DEC']  # month list
content = ['time', 'lon', 'lat', 'BCSMASS', 'DUSMASS25', 'OCSMASS', 'SO4SMASS', 'SSSMASS25']   # ds.variables.keys()
#MAKING LIST FOR EACH MONTH TO STORING DATA
jan = []
wjan = []

feb = []
wfeb = []

mar = []
wmar = []

apr = []
wapr = []

may = []
wmay = []

jun = []
wjun = []

jul = []
wjul = []

aug = []
waug = []

sep = []
wsep = []

octo = []
wocto = []

nov = []
wnov = []

dec = []
wdec = []


for i in range(1,13):
    if i<10 :
          month = "0"+str(i)
    else:
          month = str(i)
    

    # fn = f"data/MERRA2_400.tavg1_2d_aer_Nx.{year}{month}{date}.SUB.nc" formate of file
    
    for year in range(2021,2022):
        k=f"MERRA2_400.inst3_3d_aer_Nv.{year}{month}"
        y=f"MERRA2_401.inst3_3d_aer_Nv.{year}{month}"
        se=f"MERRA2_400.tavg3_3d_asm_Nv.{year}{month}"
        sel=f"MERRA2_401.tavg3_3d_asm_Nv.{year}{month}"
        
        if i == 1:
              f1 = ["vp_2021/"+x for x in l("vp_2021/") if k in x] # PM2.5 DATA FILES
              w1 = ["wb_2021/"+x for x in l("wb_2021/") if se in x] #WIND DATA FILES
              jan=jan+f1#APPENDING ALL FILES OF JAN 
              wjan=wjan+w1#APPENDING ALL FILES OF WIND JAN
        elif i == 2:
              f2 = ["vp_2021/"+x for x in l("vp_2021/") if k in x]
              w2 = ["wb_2021/"+x for x in l("wb_2021/") if se in x]
              feb=feb+f2
              wfeb=wfeb+w2
        elif i == 3:
              f3 = ["vp_2021/"+x for x in l("vp_2021/") if k in x]
              w3 = ["wb_2021/"+x for x in l("wb_2021/") if se in x]
              mar=mar+f3
              wmar=wmar+w3
        elif i == 4:
              f4 = ["vp_2021/"+x for x in l("vp_2021/") if k in x]
              w4 = ["wb_2021/"+x for x in l("wb_2021/") if se in x]
              apr=apr+f4
              wapr=wapr+w4
        elif i == 5:
              f5 = ["vp_2021/"+x for x in l("vp_2021/") if k in x]
              w5 = ["wb_2021/"+x for x in l("wb_2021/") if se in x]
              may=may+f5
              wmay=wmay+w5
        elif i == 6:
              f6 = ["vp_2021/"+x for x in l("vp_2021/") if y in x]
              w6 = ["wb_2021/"+x for x in l("wb_2021/") if sel in x]
              jun=jun+f6
              wjun=wjun+w6

        elif i == 7:
              f7 = ["vp_2021/"+x for x in l("vp_2021/") if y in x]
              w7 = ["wb_2021/"+x for x in l("wb_2021/") if sel in x]
              jul=jul+f7
              wjul=wjul+w7
        elif i == 8:
              f8 = ["vp_2021/"+x for x in l("vp_2021/") if y in x]
              w8 = ["wb_2021/"+x for x in l("wb_2021/") if sel in x]
              aug=aug+f8
              waug=waug+w8

        elif i == 9:
              f9 = ["vp_2021/"+x for x in l("vp_2021/") if y in x]
              w9 = ["wb_2021/"+x for x in l("wb_2021/") if sel in x]
              sep=sep+f9
              wsep=wsep+w9

        elif i == 10:
              f10 = ["vp_2021/"+x for x in l("vp_2021/") if k in x]
              w10 = ["wb_2021/"+x for x in l("wb_2021/") if se in x]
              octo=octo+f10
              wocto=wocto+w10

        elif i == 11:
              f11 = ["vp_2021/"+x for x in l("vp_2021/") if k in x]
              w11 = ["wb_2021/"+x for x in l("wb_2021/") if se in x]
              nov=nov+f11
              wnov=wnov+w11
        else:
              f12 = ["vp_2021/"+x for x in l("vp_2021/") if k in x]
              w12 = ["wb_2021/"+x for x in l("wb_2021/") if se in x]
              dec=dec+f12
              wdec=wdec+w12



#MAKING HAVING ALL YAER DATA BY MONTH AND STROING IN LIST M
m = [jan, feb, mar, apr, may, jun, jul, aug, sep, octo, nov, dec]
a = [wjan, wfeb, wmar, wapr, wmay, wjun, wjul, waug, wsep, wocto, wnov, wdec]

for i in range(0,12):
      dayAvgs = []
      for file in m[i]:
          #EXTRACTING PM2.5 DATA
          ds = nc.Dataset(file)
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
      print(np.max(moAvgs),np.min(moAvgs))


    #EXTRACTING WIND  DATA
      wU500 = []
      wV500 = []
      HJ =[]
      for file in a[i]:
            dwws = nc.Dataset(file)
            U500 = np.array(dwws["U"][:, : , latli:latui , lonli:lonui])
            T = np.array(dwws["T"][:, : , latli:latui , lonli:lonui])
            RH=np.array(dwws["RH"][:, : , latli:latui , lonli:lonui])
            V500 = np.array(dwws["V"][:, : , latli:latui , lonli:lonui])
            wU500.append(np.mean(U500, axis=0))
            wV500.append(np.mean(V500, axis=0))
            HJ.append(np.mean(T, axis=0))#IF YOU WANT TO OTHER PARAMETER CHANGE HERE 'T' WITH SOMETHING ELSE RH
            
      HJ = np.array(HJ)
      KO = np.mean(HJ, axis = 0)

      wU500 = np.array(wU500)
      wV500 = np.array(wV500)

      mwU500 = np.mean(wU500, axis=0)
      mwV500 = np.mean(wV500, axis=0)
      fig = plt.figure()
      ax = fig.add_subplot(111)
      p=moAvgs[:,0,0]
      p=np.array(p)
      p = np.round_(p, decimals = 2)
      levs=[0.0100,0.0200,0.0327,0.0476,0.0660,0.0893,0.1197,0.1595,0.2113,0.2785,0.3650,0.4758,0.6168,0.7951,1.0194,1.3005,1.6508,2.0850,2.6202,3.2764,4.0766,5.0468,6.2168,7.6198,9.2929,11.2769,13.6434,16.4571,19.7916,23.7304,28.3678,33.8100,40.1754,47.6439,56.3879,66.6034,78.5123,92.3657,108.6630,127.8370,150.3930,176.9300,208.1520,244.8750,288.0830,337.5000,375.0000,412.5000,450.0000,487.5000,525.0000,562.5000,600.0000,637.5000,675.0000,700.0000,725.0000,750.0000,775.0000,800.0000,820.0000,835.0000,850.0000,865.0000,880.0000,895.0000,910.0000,925.0000,940.0000,955.0000,970.0000,985.0000]
      levs=np.array(levs)
      levs = np.round_(levs, decimals = 2)
      plt.plot(p,levs,label = MAH[i],marker='o')
      plt.gca().invert_yaxis()
      #plt.xlim(0, 200)
      for i,j in zip(p,levs):
          ax.annotate('%s)' %j, xy=(i,j), xytext=(30,0), textcoords='offset points')
          ax.annotate('(%s,' %i, xy=(i,j))




     # plt.ylim(1000, 100)
      plt.xlabel('air temp. (k)')
      plt.ylabel('PM2.5 Concentration (\u03BCg/m\u00b3)')
      plt.title('nanital VERTICAL PROFILE (AVERAGE) 2021' )##note: please change the name according to city
      plt.legend()
      plt.grid()
      plt.show()


            

        
