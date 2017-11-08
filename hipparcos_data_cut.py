"""
Geometric and magnitude cuts to obtain population of tracer A- and F-stars from the new reduction of the Hipparcos catalog. For the calculation of vertical
velocities, radial velocities for ~70% of the stars were taken from Barbier, Brossat & Figon (2000).

Hipparcos new reduction: http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=+I%2F311%2Fhip2&-from=nav&-nav=cat%3AI%2F311%26tab%3A%7BI%2F311%2Fhip2%7D%26key%3Asource%3DI%2F311%2Fhip2%26HTTPPRM%3A%26-out.max%3D50%26-out.form%3DHTML+Table%26-out.add%3D_r%26-out.add%3D_RAJ%2C_DEJ%26-sort%3D_r%26-oc.form%3Dsexa%26
Radial velocity catalog: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=III%2F213
"""

import math, numpy as np
from scipy import interpolate
from sympy.mpmath import *
from astropy import constants as const
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Distance
import gala.coordinates as gc
np.seterr(divide='ignore', invalid='ignore')

#----------------constants & paths-------------------------

w_sun = 7.25 #* (u.km/u.s)
u_sun = 11.1 #* (u.km/u.s)
v_sun = 12.24 #* (u.km/u.s)
rc_A_stars = 0.170 * u.kpc
hc_A_stars = 0.170 * u.kpc
rc_F_stars = 0.050 * u.kpc
hc_F_stars_min =  - 0.109 * u.kpc
hc_F_stars_max = 0.04 * u.kpc  
infilepath = '/Users/Jatan/Desktop/DDDM/Data/'
outfilepath = '/Users/Jatan/Desktop/DDDM/Results/'

#--------------------input---------------------------------

data_hip = infilepath + 'hipparcos_data.tsv'   
hd, RA, DE, plx, err_plx, pm_RA, pm_DE, Hpmag, BV = np.loadtxt(data_hip, delimiter= ";", skiprows=3, usecols=(2, 6, 7, 8, 9, 10, 11, 12, 13), unpack=True)
data_rad = infilepath + 'radial_velocity.tsv'	#len = 36145
hd_rad, rad_vel, err_rad_vel = np.genfromtxt(data_rad, delimiter= ";", filling_values=0, skip_header=3, usecols=(3, 8, 9), unpack=True)  
hip_id = [int(x) for x in hd]
hip_id_rad = [int(x) for x in hd_rad]

#--------------------main()---------------------------------

temp_coord = SkyCoord(RA, DE, frame='icrs', unit='deg')
Gal_long = temp_coord.galactic.l.degree
Gal_lat = temp_coord.galactic.b.degree

distance = abs(1/plx) #Unit -> kpc
z_coord = (SkyCoord(l = Gal_long * u.degree, b = Gal_lat * u.degree, distance = distance*u.kpc, frame='galactic')).cartesian.z
y_coord = (SkyCoord(l = Gal_long * u.degree, b = Gal_lat * u.degree, distance = distance*u.kpc, frame='galactic')).cartesian.y
x_coord = (SkyCoord(l = Gal_long * u.degree, b = Gal_lat * u.degree, distance = distance*u.kpc, frame='galactic')).cartesian.x
r_coord = np.sqrt(pow(x_coord, 2) + pow(y_coord, 2))         #positivity of r follows from radial symmetry
Mv = Hpmag + 5 * np.log10(abs(plx)/100) #to obtain absolute magnitude; formula from https://www.cosmos.esa.int/web/hipparcos/sample-tables
pm_long, pm_lat = gc.pm_icrs_to_gal(temp_coord, [pm_RA, pm_DE] * u.mas/u.yr).value
velz_all = ((4.74047 * pm_lat * np.cos(Gal_lat*u.deg))/plx) + (w_sun * pow(np.cos(Gal_lat*u.deg), 2)) - (u_sun * np.cos(Gal_long*u.deg) * np.cos(Gal_lat*u.deg) * np.sin(Gal_lat*u.deg)) - (v_sun * np.sin(Gal_long*u.deg) * np.cos(Gal_lat*u.deg) * np.sin(Gal_lat*u.deg)) #* (u.km/u.s)

										
										#cuts to obtain A star population and velocities for low-lat stars

indexes_A_stars = []
indexes_vel_A_stars = []
for i in range(len(hd)):
	if ((0.0 < Mv[i] < 1.0) & (-0.2 < BV[i] < 0.6)):
		if((r_coord[i] < rc_A_stars) & ((- hc_A_stars) < z_coord[i] < (hc_A_stars))):	#should zsun be included? No, because the reference point for all distances is the sun! zsun is included while constructing densities.
			indexes_A_stars.append(i)

for i in indexes_A_stars:
 	if (abs(Gal_lat[i]) < 12):
 		indexes_vel_A_stars.append(i)

hip_vel_A = [hip_id[i] for i in indexes_vel_A_stars]		  #determine hipparcos id for low lat A-stars
temp_list_A = list(set(hip_vel_A).intersection(hip_id_rad))   #find elements common to both (hipparcos & radial velocity) catalogs
vel_indices_A = [[hip_vel_A.index(i), hip_id_rad.index(i)] for i in temp_list_A]   #[element in indexes_vel_A_stars, element in radial_vel]

flag_A = np.zeros(len(indexes_vel_A_stars))
rad_list_A = [vel_indices_A[i][0] for i in range(len(vel_indices_A))] 
for i in rad_list_A: 
	flag_A[i]+=1

vel_A = []
for i in indexes_vel_A_stars:
	j = indexes_vel_A_stars.index(i)
	if flag_A[j] == 0:
		vel_A.append(((4.74047 * pm_lat[i] * np.cos(Gal_lat[i]*u.deg))/plx[i]) + (w_sun * pow(np.cos(Gal_lat[i]*u.deg), 2)) - (u_sun * np.cos(Gal_long[i]*u.deg) * np.cos(Gal_lat[i]*u.deg) * np.sin(Gal_lat[i]*u.deg)) - (v_sun * np.sin(Gal_long[i]*u.deg) * np.cos(Gal_lat[i]*u.deg) * np.sin(Gal_lat[i]*u.deg))) #* (u.km/u.s)
	else:
		temp_index = vel_indices_A[temp_list_A.index(hip_vel_A[j])][1]
		#print "calculating using radial velocity!!"
		vel_A.append((w_sun + (4.74047 * pm_lat[i] * np.cos(Gal_lat[i]*u.deg))/plx[i]) + rad_vel[temp_index] * np.sin(Gal_lat[i]*u.deg))

velz_A_lowlat = [vel_A[i].value for i in range(len(vel_A))]

									#cuts to obtain F star population and velocities for low-lat stars
indexes_F_stars = []
indexes_vel_F_stars = []
for i in range(len(hd)):
	if ((1.0 < Mv[i] < 2.5) & (-0.2 < BV[i] < 0.6)):
		if((r_coord[i] < rc_F_stars) & ((hc_F_stars_min) < z_coord[i] < (hc_F_stars_max))):
			indexes_F_stars.append(i)

for i in indexes_F_stars:
	if (abs(Gal_lat[i]) < 12):
		indexes_vel_F_stars.append(i)

hip_vel_F = [hip_id[i] for i in indexes_vel_F_stars]
temp_list_F = list(set(hip_vel_F).intersection(hip_id_rad))    
vel_indices_F = [[hip_vel_F.index(i), hip_id_rad.index(i)] for i in temp_list_F]   #[element in indexes_vel_F_stars, element in radial_vel]

flag_F = np.zeros(len(indexes_vel_F_stars))
rad_list_F = [vel_indices_F[i][0] for i in range(len(vel_indices_F))] 
for i in rad_list_F: 
	flag_F[i]+=1

vel_F = []
for i in indexes_vel_F_stars:
	j = indexes_vel_F_stars.index(i)
	if flag_F[j] == 0:
		vel_F.append(((4.74047 * pm_lat[i] * np.cos(Gal_lat[i]*u.deg))/plx[i]) + (w_sun * pow(np.cos(Gal_lat[i]*u.deg), 2)) - (u_sun * np.cos(Gal_long[i]*u.deg) * np.cos(Gal_lat[i]*u.deg) * np.sin(Gal_lat[i]*u.deg)) - (v_sun * np.sin(Gal_long[i]*u.deg) * np.cos(Gal_lat[i]*u.deg) * np.sin(Gal_lat[i]*u.deg))) #* (u.km/u.s)
	else:
		temp_index = vel_indices_F[temp_list_F.index(hip_vel_F[j])][1]
		#print "calculating using radial velocity!!"
		vel_F.append((w_sun + (4.74047 * pm_lat[i] * np.cos(Gal_lat[i]*u.deg))/plx[i]) + rad_vel[temp_index] * np.sin(Gal_lat[i]*u.deg))

velz_F_lowlat = [vel_F[i].value for i in range(len(vel_F))]


 								#saving into new files for A and F stars separately

# Gal_long_A = [Gal_long[i] for i in indexes_A_stars]; Gal_long_F = [Gal_long[i] for i in indexes_F_stars]
# Gal_lat_A = [Gal_lat[i] for i in indexes_A_stars]; Gal_lat_F = [Gal_lat[i] for i in indexes_F_stars]
# plx_A = [plx[i] for i in indexes_A_stars]; plx_F = [plx[i] for i in indexes_F_stars]
# err_plx_A = [err_plx[i] for i in indexes_A_stars]; err_plx_F = [err_plx[i] for i in indexes_F_stars]
# Mag_A = [Mv[i] for i in indexes_A_stars]; Mag_F = [Mv[i] for i in indexes_F_stars]
# BV_A = [BV[i] for i in indexes_A_stars]; BV_F = [BV[i] for i in indexes_F_stars]
# z_coord_A = [z_coord.value[i] for i in indexes_A_stars]; z_coord_F = [z_coord.value[i] for i in indexes_F_stars]; 
vel_z_A = [velz_all[i] for i in indexes_vel_A_stars]; vel_z_F = [velz_all[i] for i in indexes_vel_F_stars]; 

# newfile1 = outfilepath + 'hipparcos_data_Astars.txt'   #distance.txt
# np.savetxt(newfile1, np.c_[Gal_long_A, Gal_lat_A, plx_A, err_plx_A, z_coord_A, vel_z_A, Mag_A, BV_A], fmt='%e', header='created in python with cuts following 1604.01407 \nLongitude      Latitude       plx      err_plx      z_coord    Velocity_z      Hpmag      (B-V)')
# newfile2 = outfilepath + 'hipparcos_data_Fstars.txt'  
# np.savetxt(newfile2, np.c_[Gal_long_F, Gal_lat_F, plx_F, err_plx_F, z_coord_F, vel_z_F, Mag_F, BV_F], fmt='%e', header='created in python with cuts following 1604.01407 \nLongitude      Latitude       plx      err_plx      z_coord    Velocity_z      Hpmag      (B-V)')
newfile3 = outfilepath + 'vel_Astars_uncorr.txt' 
np.savetxt(newfile3, np.c_[vel_z_A, velz_A_lowlat], fmt='%e', header='created in python with cuts following 1604.01407 \nvel_z   vel_z (w/rad vel)')
newfile4 = outfilepath + 'vel_Fstars_uncorr.txt' 
np.savetxt(newfile4, np.c_[vel_z_F, velz_F_lowlat], fmt='%e', header='created in python with cuts following 1604.01407 \nvel_z   vel_z (w/rad vel)')
print "Your new files are %s" %newfile3 +" and %s" %newfile4 #%s" %newfile1 +", %s" %newfile2 +", 


#print len(indexes_A_stars)
# print len(indexes_vel_A_stars)
