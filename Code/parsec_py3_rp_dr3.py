import emcee
import pickle
from pylab import *
import glob
import scipy
from scipy import interpolate
from astropy import constants as cons
import corner
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u
import sys

# Path to isochrone files
default_isopath = '/Users/rbrahm/data/PARSEC-2/'
default_isopath = 'PARSEC_2/'

def mag_to_flux(mag,name):
	dct = {'Johnson_B':6.491e-9,
		   'Johnson_V':3.734e-9,
		   'SDSS_g':5.056e-9,
		   'SDSS_r':2.904e-9,
		   'SDSS_i':1.967e-9,
		   'Gmag':18851516970.802105,
		   'BPmag':13821504093.71292,
		   'RPmag':8030969471.986631,
		   '2mass_J':3.129e-10,
		   '2mass_H':1.133e-10,
		   '2mass_Ks':4.283e-11,
		   'W1':8.081e-12,
		   'W2':2.397e-12,
		   'W3':7.112e-14,
		   'W4':4.507e-15,
	}

	flux  = dct[name] * 10**(-mag/2.5)
	return flux

def flux_to_mag(flux,name):
	"""
	Convert flux to magnitude for a variety of filters.

	Parameters:
	flux (float): stellar flux
	name (str): filter name

	Returns:
	mag (float): stellar magnitude in the filter
	"""

	# Define conversion factors for different filters
	dct = {'Johnson_B':6.491e-9,
		   'Johnson_V':3.734e-9,
		   'SDSS_g':5.056e-9,
		   'SDSS_r':2.904e-9,
		   'SDSS_i':1.967e-9,
		   'Gmag':18851516970.802105,
		   'BPmag':13821504093.71292,
		   'RPmag':8030969471.986631,
		   '2mass_J':3.129e-10,
		   '2mass_H':1.133e-10,
		   '2mass_Ks':4.283e-11,
		   'W1':8.081e-12,
		   'W2':2.397e-12,
		   'W3':7.112e-14,
		   'W4':4.507e-15,
	}
	# Compute magnitude from flux
	mag = -2.5*np.log10(flux/dct[name])
	# Return the magnitude
	return mag


def rhostar(M,R):
	"""
	Compute stellar density from mass and radius.

	Parameters:
	M (float): stellar mass in Msun units
	R (float): stellar radius in Rsun units

	Returns:
	rho (float): stellar density in cgs units
	"""
	rho = M*cons.M_sun.cgs.value/(4.*np.pi*(R*cons.R_sun.cgs.value)**3/3.)
	return rho

def loggtoRs(logg,Ms):
	"""
	Compute stellar radius from surface gravity and mass.

	Parameters:
	logg (float): log of the stellar surface gravity
	Ms (float): stellar mass in Msun units

	Returns:
	Rs (float): stellar radius in Rsun units
	"""
	# compute surface gravity from logg
	g = 10**logg
	# compute stellar radius from Newton's law in cgs
	Rs = np.sqrt(cons.G.cgs.value*Ms*cons.M_sun.cgs.value/g)
	# convert to Rsun units
	Rs = Rs / cons.R_sun.cgs.value
	# return stellar radius in Rsun
	return Rs

def FEHtoZ(FEH,Zsun=0.0152):
	"""
	Convert Fe/H to Z (metal mass fraction).

	Parameters:
	FEH (float): stellar metallicity
	Zsun (float): solar Z

	Returns:
	(float): stellar Z
	"""
	return Zsun*10**(FEH)

def ZtoFEH(Z,Zsun=0.0152):
	return np.log10(Z/Zsun)

def get_Zs(isopath=default_isopath):
	fls = glob.glob(isopath+'parsec_av0.0*.dat')
	fehs = []
	for fl in fls:
		z = fl.split('_z')[-1][:-4]
		fehs.append(float(z))
	return np.sort(np.array(fehs))

def stringed(Z):
	"""
	Convert Z isochrone value back into key string
	Necessary because prior float conversion strips right zeroes.

	Parameters:
	Z (float): Z isochrone key value previously converted to float

	Returns
	sZ (string): isochrone key corresponding to the Z value
	"""

	# round Z to 4 decimal places and convert value to str
	sZ = str(np.around(Z,4))
	# get difference of length from 5
	diff = 5 - len(sZ)
	i=0
	# while length is less than 5, add right zeroes
	while i < diff:
		sZ = sZ + '0'
		i+=1
	# return string value with right zeroes added if necessary
	return sZ

def make_dct(isopath=default_isopath):
	Zs = get_Zs(isopath=isopath)

	#avs = ['0.2','0.4','0.6','0.8','1.0']
	avs = ['1.0']
	dct = {}
	for Z in Zs:
		sZ = stringed(Z)
		dct[sZ] = {}

		path = isopath + 'parsec_av0.0_z' + sZ + '.dat'
		#print(path, dsa)

		tage =-999
		f = open(path,'r')
		lines = f.readlines()
		f.close()


		olines = {}
		for av in avs:
			opath = isopath + 'parsec_av'+av+'_z' + sZ + '.dat'
			of = open(opath,'r')
			olines[av] =  of.readlines()
			of.close()

		count = 2
		for iline in np.arange(len(lines)):
			line = lines[iline]
			print(line)
			if line[0]!='#':
				if float(line.split()[1]) != tage:
					if tage != -999:
						dct[sZ][str(tage)]['imass'] = np.array(imass)
						dct[sZ][str(tage)]['mass'] = np.array(mass)
						dct[sZ][str(tage)]['logL'] = np.array(logL)
						dct[sZ][str(tage)]['logT'] = np.array(logT)
						dct[sZ][str(tage)]['logg'] = np.array(logg)
						dct[sZ][str(tage)]['Gmag_0.0'] = np.array(Gmag)
						dct[sZ][str(tage)]['BPmag_0.0'] = np.array(BPmag)
						dct[sZ][str(tage)]['RPmag_0.0'] = np.array(RPmag)
						dct[sZ][str(tage)]['Jmag_0.0'] = np.array(Jmag)
						dct[sZ][str(tage)]['Hmag_0.0'] = np.array(Hmag)
						dct[sZ][str(tage)]['Kmag_0.0'] = np.array(Kmag)

						count += 1

						for av in avs:
							dct[sZ][str(tage)]['Gmag_'+av] = np.array(Gmagd[av])
							dct[sZ][str(tage)]['BPmag_'+av] = np.array(BPmagd[av])
							dct[sZ][str(tage)]['RPmag_'+av] = np.array(RPmagd[av])
							dct[sZ][str(tage)]['Jmag_'+av] = np.array(Jmagd[av])
							dct[sZ][str(tage)]['Hmag_'+av] = np.array(Hmagd[av])
							dct[sZ][str(tage)]['Kmag_'+av] = np.array(Kmagd[av])


					tage = float(line.split()[1])
					dct[sZ][str(tage)] = {}

					imass = [float(line.split()[2])]
					mass = [float(line.split()[3])]
					logL = [float(line.split()[4])]
					logT = [float(line.split()[5])]
					logg = [float(line.split()[6])]
					Gmag = [float(line.split()[8])]
					BPmag = [float(line.split()[9])]
					RPmag = [float(line.split()[10])]
					Jmag = [float(line.split()[13])]
					Hmag = [float(line.split()[14])]
					Kmag = [float(line.split()[15])]

					Gmagd, BPmagd, RPmagd, Jmagd, Hmagd, Kmagd = {},{},{},{},{},{}
					for av in avs:
						Gmagd[av] = [float(olines[av][iline+count].split()[8])]
						BPmagd[av] = [float(olines[av][iline+count].split()[9])]
						RPmagd[av] = [float(olines[av][iline+count].split()[10])]
						Jmagd[av] = [float(olines[av][iline+count].split()[13])]
						Hmagd[av] = [float(olines[av][iline+count].split()[14])]
						Kmagd[av] = [float(olines[av][iline+count].split()[15])]

				else:
					if not float(line.split()[2]) in imass:
						imass.append(float(line.split()[2]))
						mass.append(float(line.split()[3]))
						logL.append(float(line.split()[4]))
						logT.append(float(line.split()[5]))
						logg.append(float(line.split()[6]))
						Gmag.append(float(line.split()[8]))
						BPmag.append(float(line.split()[9]))
						RPmag.append(float(line.split()[10]))
						Jmag.append(float(line.split()[13]))
						Hmag.append(float(line.split()[14]))
						Kmag.append(float(line.split()[15]))

						for av in avs:
							Gmagd[av].append(float(olines[av][iline+count].split()[8]))
							BPmagd[av].append(float(olines[av][iline+count].split()[9]))
							RPmagd[av].append(float(olines[av][iline+count].split()[10]))
							Jmagd[av].append(float(olines[av][iline+count].split()[13]))
							Hmagd[av].append(float(olines[av][iline+count].split()[14]))
							Kmagd[av].append(float(olines[av][iline+count].split()[15]))

		dct[sZ][str(tage)]['imass'] = np.array(imass)
		dct[sZ][str(tage)]['mass'] = np.array(mass)
		dct[sZ][str(tage)]['logL'] = np.array(logL)
		dct[sZ][str(tage)]['logT'] = np.array(logT)
		dct[sZ][str(tage)]['logg'] = np.array(logg)
		dct[sZ][str(tage)]['Gmag_0.0'] = np.array(Gmag)
		dct[sZ][str(tage)]['BPmag_0.0'] = np.array(BPmag)
		dct[sZ][str(tage)]['RPmag_0.0'] = np.array(RPmag)
		dct[sZ][str(tage)]['Jmag_0.0'] = np.array(Jmag)
		dct[sZ][str(tage)]['Hmag_0.0'] = np.array(Hmag)
		dct[sZ][str(tage)]['Kmag_0.0'] = np.array(Kmag)

		count += 1

		for av in avs:
			dct[sZ][str(tage)]['Gmag_'+av] = np.array(Gmagd[av])
			dct[sZ][str(tage)]['BPmag_'+av] = np.array(BPmagd[av])
			dct[sZ][str(tage)]['RPmag_'+av] = np.array(RPmagd[av])
			dct[sZ][str(tage)]['Jmag_'+av] = np.array(Jmagd[av])
			dct[sZ][str(tage)]['Hmag_'+av] = np.array(Hmagd[av])
			dct[sZ][str(tage)]['Kmag_'+av] = np.array(Kmagd[av])


	print(isopath+'parsec.pkl')
	pickle.dump(dct,open(isopath+'parsec.pkl','w'))

def presus(x,vec):
	"""
	Find closest elements of a sorted array to a value.

	Parameters:
	x (float): value to test
	vec (array): sorted array in which to find elements

	Returns:
	vec[I2] (float): closest element of vec that is smaller than x
	vec[I1] (float): closest element of vec that is larger than x
	"""
	# compute difference from x for each element in vec
	res = vec - x
	# check if x is contained in range spanned by vec
	if res[-1] >= 0 and res[0]<=0:
		# get index of closest element of vec that is larger than x
		I1 = np.where(res>=0)[0][0]
		# get index of closest element of vec that is smaller than x
		I2 = np.where(res<=0)[0][-1]
	# if x is larger than vec range, select largest two elements of vec
	elif res[-1]<0:
		I1 = -1
		I2 = -2
	# if x is smaller than vec range, select smallest two elements of vec
	elif res[0]>0:
		I1 = 1
		I2 = 0
	# return the two closest elements of vec to x
	return vec[I2],vec[I1]


def trilinear(x,y,z,x0,x1,y0,y1,z0,z1,c000,c100,c010,c001,c110,c011,c101,c111):
	if x1 != x0:
		xd = (x-x0) / (x1 - x0)
	else:
		xd = 0
	if y1 != y0:
		yd = (y-y0) / (y1 - y0)
	else:
		yd = 0
	if z1 != z0:
		zd = (z-z0) / (z1 - z0)
	else:
		zd = 0

	c00 = c000*(1-xd) + c100*xd
	c01 = c001*(1-xd) + c101*xd
	c10 = c010*(1-xd) + c110*xd
	c11 = c011*(1-xd) + c111*xd
	c0 = c00*(1-yd) + c10*yd
	c1 = c01*(1-yd) + c11*yd
	c = c0 * (1-zd) + c1*zd
	return c

def bilinear(x,y,x0,x1,y0,y1,q11,q21,q12,q22):
	"""
	Perform a bilinear interpolation of a function f(x,y),
	given four coordinates and corresponding f values.

	Parameters:
	x (float): x coordinate at which to interpolate the function
	y(float): y coordinate at which to interpolate the function
	x0 (float): x coordinate of the first and third data points
	x1 (float): x coordinate of the second and fourth data points
	y0 (float): y coordinate of the first and second data points
	y1 (float): y coordinate of the third and fourth data points
	q11 (float): value of f at (x0, y0)
	q21 (float): value of f at (x1, y0)
	q12 (float): value of f at (x0, y1)
	q22 (float): value of f at (x1, y1)

	Returns:
	float: the value of the interpolated function at (x,y)
	"""
	# check the input coordinates are different
	if x1 != x0 and y1 != y0:
		# do linear interpolation in x
		cxa = (x1-x) / (x1 - x0)
		cxb = (x-x0) / (x1 - x0)

		fxy1 = cxa*q11 + cxb*q21
		fxy2 = cxa*q12 + cxb*q22

		# do linear interpolation in y
		cya = (y1-y) / (y1 - y0)
		cyb = (y-y0) / (y1 - y0)

		fxy = cya*fxy1 + cyb*fxy2

		# return the interpolated value
		return fxy

	# if all input coordinates are the same
	elif x1 == x0 and y1 == y0:
		# return the value at the point
		return q11

	# if x coords are the same but y coords differ
	elif x1 == x0 and y1 != y0:
		# do linear interpolation in y
		m = (q11 - q12) / (y1 - y0)
		n = q11 - m*y0
		# return the interpolated value
		return m*y + n

	# if x coords differ and y coords are the same
	elif x1 != x0 and y1 == y0:
		# do linear interpolation in x
		m = (q11 - q21) / (x1 - x0)
		n = q11 - m*x0
		# return the interpolated value
		return m*x + n



def get_vals(Ms, Age, FeH, isopath=default_isopath,keys =['mass','logL','logT','logg','Gmag_0.0','BPmag_0.0','RPmag_0.0','Jmag_0.0','Hmag_0.0','Kmag_0.0']):
	"""
	Interpolate values for a set of parameters based on the isochrones.

	Parameters:
	Ms (float): initial mass
	Age (float): stellar age
	FeH (float): metallicity
	isopath (str): file path to isochrones
	keys (list of str): parameters for which to interpolate values

	Returns:
	out (dict): values for each key
	"""
	# convert stellar metallicity to stellar Z
	Z = FEHtoZ(FeH)
	#print(Ms, Age,Z)
	# list the Z values for each isochrone (main dct keys)
	zetas = list(dct.keys())
	#list the times for the first isochrone (keys to the dict for each Z)
	times = list(dct[zetas[0]].keys())
	# sort the Z and time arrays
	zetas = np.sort(np.array(zetas).astype('float'))
	times = np.sort(np.array(times).astype('float'))
	#print(Age, times)
	# initialize output dictionary
	out = {}
	# If Z or age are outside isochrone values return -99999 for all keys
	if Z < zetas.min() or Z > zetas.max() or Age < times.min() or Age > times.max():
		for key in keys:
			out[key] = -999999999999999
		return out

	# If Z and age are within isochrone values continue
	# get the closest two elements of zetas to stellar Z
	Z1,Z2 = presus(Z,zetas)
	# get the closest two elements of times to stellar age
	T1,T2 = presus(Age,times)
	# get imass values for all combinations of Z1, Z2, T1, T2
	M11 = dct[stringed(Z1)][str(T1)]['imass']
	M12 = dct[stringed(Z1)][str(T2)]['imass']
	M21 = dct[stringed(Z2)][str(T1)]['imass']
	M22 = dct[stringed(Z2)][str(T2)]['imass']

	# get the maximum of the minimum imass values
	m1 = np.max([np.min(M11),np.min(M12),np.min(M21),np.min(M22)])
	# get the minimum of the maximum imass values
	m2 = np.min([np.max(M11),np.max(M12),np.max(M21),np.max(M22)])

	# get the sorted index values for each imass array
	I11 = np.argsort(M11)
	I12 = np.argsort(M12)
	I21 = np.argsort(M21)
	I22 = np.argsort(M22)

	# combine and sort all unique imass values from the four arrays into one array
	mtot = np.sort(np.unique(np.hstack((M11,M12,M21,M22))))
	# select only values within the m1, m2 range as reference imass values
	I = np.where((mtot>=m1) & (mtot<=m2))[0]
	mref = mtot[I]

	# initialize output dictionary
	out = {}
	# If stellar mass is outside isochrone values return -99999 for all keys
	if Ms < mref.min() or Ms > mref.max():
		for key in keys:
			out[key] = -999999999999999
		return out

	#If stellar mass is within isochrone values continue
	# Ensure m1, m2 are in mref
	if mref[0]>m1:
		mref = np.hstack((np.array([m1]),mref))
	if mref[-1]<m2:
		mref = np.hstack((mref,np.array([m2])))

	#M1,M2 = presus(Ms,mref)

	# initialize output dictionary
	out = {}
	if True:
		# Loop through the keys
		for key in keys:
			# get key isochrone values for all combinations of Z1, Z2, T1, T2
			V11 = dct[stringed(Z1)][str(T1)][key]
			V12 = dct[stringed(Z1)][str(T2)][key]
			V21 = dct[stringed(Z2)][str(T1)][key]
			V22 = dct[stringed(Z2)][str(T2)][key]
			# interpolate key values by spline over corresponding isochrone imass values
			# then evaluate at the stellar initial mass value
			tck = interpolate.splrep(M11[I11],V11[I11],k=1)
			logT11 = interpolate.splev(Ms,tck)
			#logT112 = interpolate.splev(M2,tck)
			tck = interpolate.splrep(M12[I12],V12[I12],k=1)
			logT12 = interpolate.splev(Ms,tck)
			#logT122 = interpolate.splev(M2,tck)
			tck = interpolate.splrep(M21[I21],V21[I21],k=1)
			logT21 = interpolate.splev(Ms,tck)
			#logT212 = interpolate.splev(M2,tck)
			tck = interpolate.splrep(M22[I22],V22[I22],k=1)
			logT22 = interpolate.splev(Ms,tck)
			#logT222 = interpolate.splev(M2,tck)

			#logT = trilinear(Ms,Age,Z,M1,M2,T1,T2,Z1,Z2,logT111,logT211,logT121,logT112,logT221,logT122,logT212,logT222)
			# do bilinear interpolation over all combinations
			# of Z1, Z2, T1, T2, and the spline-interpolated values
			# to get the value of the key at (Z, Age)
			logT = bilinear(Z,Age,Z1,Z2,T1,T2,logT11,logT21,logT12,logT22)
			# save the value to output dictionary
			out[key] = logT
	# return the output dictionary with the values of all the keys
	return out


def get_model(Ms, Age, FeH, Av=0., isopath=default_isopath):
	dct = pickle.load(open(isopath+'parsec.pkl','rb'), encoding = 'latin1')

	aage = np.log10(Age*1e9)
	#avs = ['0.0','0.2','0.4','0.6','0.8','1.0']
	avs = ['0.0','1.0']

	vavs = np.array(avs).astype('float')
	av1,av2 = presus(Av,vavs)

	out = get_vals(Ms, aage, FeH, isopath=isopath, keys =['mass','logL','logT','logg','Gmag_'+str(av1),'BPmag_'+str(av1),'RPmag_'+str(av1),'Jmag_'+str(av1),'Hmag_'+str(av1),'Kmag_'+str(av1),\
						'Gmag_'+str(av2),'BPmag_'+str(av2),'RPmag_'+str(av2),'Jmag_'+str(av2),'Hmag_'+str(av2),'Kmag_'+str(av2)])
	Rs = loggtoRs(out['logg'],out['mass'])
	Teff = 10**out['logT']
	print(Rs, Teff)


def lnprior(theta):
	"""
	Define parameter space priors for the emcee posterior probability function.

	Parameters:
	theta (array): age, mass, Av, Fe/H if fitting metallicity too

	Returns:
	(float): 0 if parameters are within constraints, -inf if not
	"""

	# check if metallicity is being fitted
	if feh_free:
		# split age, mass, Av, FE/H from theta
		AGE, MASS, Av, feh = theta
		# check if parameters are within defined constraints
		if 0.05 < AGE < 20 and 0.4 < MASS < 4.5  and 0.0 < Av < 2. and -1 < feh < 0.5:
			# if yes return zero
			return 0.0
	# if metallicity is fixed
	else:
		# split age, mass, Av from theta
		AGE, MASS, Av = theta
		# check if parameters are within defined constraints
		if 0.05 < AGE < 20 and 0.4 < MASS < 4.5  and 0.0 < Av < 2.:
			# if yes return zero
			return 0.0

	# if the parameters are not within the constraints return -infinity
	return -np.inf

def lnlike(theta, y, yerr,isopath=default_isopath):
	"""
	Compute the log likelihood

	Parameters:
	theta (array): age, mass, Av, Fe/H if fitting metallicity too
	y: test values
	yerr: test value errors
	isopath (str): filepath to isochrones

	Returns:
	ret(float): log likelihood
	"""

	#print(theta)
	# check if metallicity is being fitted
	if feh_free:
		# split age, mass, Av, FE/H from theta
		AGE, MASS, Av, feh2 = theta
	# if metallicity is fixed
	else:
		# split age, mass, Av from theta
		AGE, MASS, Av = theta
		# set metallicity to fixed value
		feh2 = feh

	# get the log10 of the age
	aage = np.log10(AGE*1e9)

	# initialize Av limit values
	avs = ['0.0','1.0']
	# convert to array
	vavs = np.array(avs).astype('float')
	# get closest elements of vavs to Av
	# TODO these will always be 0 and 1 as it only has those elements?
	av1,av2 = presus(Av,vavs)
	# initialize list of keys
	keys = ['mass','logL','logT','logg']
	# append Av limit values to keys
	for av in avs:
		for band in obands:
			keys.append(band+'_'+av)

	# interpolate values from the isochrones for the list of keys
	out = get_vals(MASS, aage, feh2, isopath=isopath, keys =keys)
	# save output logT values to model array
	model = np.array(10**out['logT'])
	# loop through the bands
	for band in obands:
		# calculate values TODO what are these?
		m = (out[band+'_'+str(av2)] - out[band+'_'+str(av1)]) / float(av2-av1)
		n = out[band+'_'+str(av2)] - m * av2
		ys = 10**(-0.2*(m*Av + n))
		# add values to model
		model = np.hstack((model,np.array([ys])))
	# compute inverse of the error
	inv_sigma2 = 1.0/(yerr**2)
	#print(y, model)
	#ret = -np.log(2*np.pi) + np.log(np.sum(np.exp(-0.5*((y-model)/yerr)**2)/yerr))
	# compute the log likelihood of values y given the model
	ret = -0.5*(np.sum(inv_sigma2*(y-model)**2 - np.log(inv_sigma2)))

	#print(MASS,AGE,ret)
	#print(out)

	#print(theta, y, model, ret)
	# if the log likelihood is NaN or the output mass is bad
	if np.isnan(ret) or out['mass'] == -999999999999999:
		# return -infinity
		return -np.inf
	# otherwise return the log likelihood
	else:
		return ret

def lnprob(theta, y, yerr):
	"""
	Define the posterior probability for the emcee sampler.
	(emcee doc) A function that takes a vector in the parameter space as input
	and returns the natural logarithm of the posterior probability (up to an
	additive constant) for that position.

	Parameters:
	theta (array): age, mass, Av, Fe/H if fitting metallicity too
	y: test values
	yerr: test value errors

	Returns:
	(float): posterior probability
	"""

	# get the prior probabilities
	lp = lnprior(theta)
	# check if prior is infinite (i.e. parameters outside constraints in lnprior)
	if not np.isfinite(lp):
		# if it is return - infinity
		return -np.inf
	# if not return the posterior probability
	return lp + lnlike(theta, y, yerr)

def get_sum(vec):
	"""
	Get the median and the 1-sigma error bars.

	Parameters
	vec (array): a numpy array

	Returns:
	fval (float): the median of the array
	vali (float): the lower 1-sigma error bar
	valf (float): the upper 1-sigma error bar
	"""
	# sort the input array
	fvec   = np.sort(vec)
	# get the median
	fval  = np.median(fvec)
	# get the 1sigma array index
	nn = int(np.around(len(fvec)*0.15865))
	# get the 1-sigma error bars
	vali,valf = fval - fvec[nn],fvec[-nn] - fval
	# return the values
	return fval,vali,valf

def get_mags_2mass(RA,DEC,band):
	"""
	Get 2MASS magnitudes in a band from Vizier using input coordinates

	Parameters:
	RA (float): the right ascension
	DEC (float): the declination
	band (str): the band to query

	Returns:
	float(dat[band]) (float): the magnitude in the band
	float(dat['e_'+band]) (float): the error of the magnitude in the band
	"""

	# Query 2MASS catalogue on Vizier in a 10s radius region around the coordinates
	result = Vizier.query_region(coord.SkyCoord(ra=RA,dec=DEC,unit=(u.deg, u.deg),\
		frame="icrs"),width='10s',catalog='II/246')

	dat = result[0]
	# Compute the angular distance of the results to the input coordinates
	# dist = np.sqrt((dat['RAJ2000'] - RA)**2 + (dat['DEJ2000'] - DEC)**2)
	dist1 = np.sin(np.deg2rad(dat['DEJ2000']))*np.sin(np.deg2rad(DEC))
	dist2 = np.cos(np.deg2rad(dat['DEJ2000']))*np.cos(np.deg2rad(DEC))
	dist3 = np.cos(np.deg2rad(dat['RAJ2000'])-np.deg2rad(RA))
	dist = np.arccos(dist1 + dist2 * dist3)
	# Get the closest result to the input coordinates
	imin = np.argmin(dist)
	dat = dat[imin]

	# Return magnitude and error in the band for closest result to input coords
	return float(dat[band]),float(dat['e_'+band])

def get_mass(Ms, P, e, K, i):
	"""
	Compute planet mass from stellar mass and orbital parameters.

	Parameters:
	Ms (array): set of stellar mass posteriors
	P (float): orbital period
	e (float): orbital eccentricity
	K (float): RV amplitude
	i (float): orbital inclination

	Returns:
	(float): planet mass
	"""
	# if the amplitude is negative set to 1
	if K<0:
		K=1.

	# define constants
	G = cons.G.cgs.value
	Mjup = cons.M_jup.cgs.value
	Msun = cons.M_sun.cgs.value

	# create initial array of masses from 0.0001 to 100 MJ
	masses = np.arange(0.0001,100.,0.1)*Mjup

	# define function of mass, P, Ms, i, K from definition of K
	# (function will be equal to 1 for Mp)
	C1 = (2.*np.pi*G/P)**(1./3.)
	C2 = masses * np.sin(i) / ( (Ms * Msun + masses)**(2./3.) )
	C3 = 1./np.sqrt(1.-e**2)
	cte = C1 * C2 * C3 / K
	# interpolate the function over initial array of masses
	tck = interpolate.splrep(cte,masses,k=3)
	# return function evaluated at 1
	return interpolate.splev(1.,tck) / Mjup

def get_teq(L,e,a,A=0.,beta=0.5,eps=1.):
	"""
	Compute planet equilibrium temperature
	Equation from Mendez and Rivera-Valentin 2017 for elliptical orbits

	Parameters:
	L (float): stellar luminosity
	e (float): orbital eccentricity
	a (float): orbital semi-major axis
	A (float): planet bond albedo
	beta (float): heat distribution
	eps (float): emissivity

	Returns:
	(float): equilibrium temperature
	"""
	# set T0 = equilibrium temperature of Earth for zero albedo
	T0= 278.5
	# cumpute auxiliary parameter eff
	eff = (1.-A) / (beta*eps)
	# comupte and return equilibrium temperature
	return (T0*(eff * L / (a*a))**(1./4.)) * (1. - e*e*1./16. - e*e*e*e*15./1024.)

def runit(obname,isopath=default_isopath, fehfree=False, pars_path='pars/', stepnumber = 10000):
	"""
	Main function to compute stellar parameters

	Parameters:
	obname (str): name of stellar target
	isopath (str): filepath to isochrones
	fehfree (boolean): whether to recompute Fe/H from isochrones
	par_path (str): path to stellar files
	stepnumber: Numbers of steps in the MCMC

	Returns:
	(to screen and pkl file): imass, age, mass, Teff, logg, radius, rho, Av,
								  Fe/H median values and 1-sigma error bars
	"""

	# setup

	# set global parameters
	global dct,obands,feh, feh_free

	# set whether to recompute Fe/H from isochrones
	feh_free=fehfree

	# set whether to apply gaia correction
	gaia_corr = True
	# set number of resamples from normal distributions
	nsim = 1000
	# set photometric bands to query
	obands = ['Gmag','BPmag','RPmag','Jmag','Hmag','Kmag']
	# initialize magnitude and error dictionaries
	mags = {}
	emags = {}

	# read stellar Teff, Fe/H, GAIA source id from file
	fl = open(pars_path+obname+'.txt','r').readlines()
	for line in fl:
		cos = line.split()
		if cos[0] == 'teff':
			teff = float(cos[1])
			eteff = float(cos[2])
		if cos[0] == 'feh':
			feh = float(cos[1])
			efeh = float(cos[2])
		if cos[0] == 'source_id':
			sid = cos[1]

	# Query GAIA DR3 for coords, parallax, magnitudes using source ID
	job = Gaia.launch_job("select top 500 source_id,ra,dec,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_flux,phot_g_mean_flux_error,phot_g_mean_mag,phot_bp_mean_flux,phot_bp_mean_flux_error,phot_bp_mean_mag,phot_rp_mean_flux,phot_rp_mean_flux_error,phot_rp_mean_mag,radial_velocity,radial_velocity_error,teff_gspphot,ag_gspphot FROM gaiadr3.gaia_source  where (source_id="+sid+")")
	r = job.get_results()
	ra,dec = r['ra'], r['dec']
	parallax,eparallax = r['parallax'], r['parallax_error']
	pmra,epmra = r['pmra'], r['pmra_error']
	pmdec,epmdec = r['pmdec'], r['pmdec_error']
	ra2000,dec2000 = ra - pmra*15.5/(1000.*3600.), r['dec'] - pmdec*15.5/(1000.*3600.)

	# loop through the bands
	for band in obands:
		key = None
		if band == 'Gmag':
			key = 'phot_g_mean'
			systematic = 0.002
		elif band == 'BPmag':
			key = 'phot_bp_mean'
			systematic = 0.005

		elif band == 'RPmag':
			key = 'phot_rp_mean'
			systematic = 0.003
		# For GAIA magnitudes compute the error from the flux
		if key != None:
			mags[band] = float(r[key+'_mag'])
			e1 = flux_to_mag(float(r[key+'_flux'])+float(r[key+'_flux_error']),band)
			e2 = flux_to_mag(float(r[key+'_flux'])-float(r[key+'_flux_error']),band)
			err = 0.5*(e2 - e1)
			emags[band] = np.sqrt(err**2+systematic**2)
		# For 2MASS band query Vizier
		if band in ['Jmag','Hmag','Kmag']:
			mag,emag = get_mags_2mass(ra2000[0],dec2000[0],band)
			mags[band] = mag
			emags[band] = emag

	# Print name, Teff, Fe/H, magnitudes to terminal
	print('\n\t'+obname+'\n')
	print('\tTeff [K]: '+str(teff)+' +/- '+str(eteff))
	print('\t[Fe/H]: '+str(feh)+' +/- '+str(efeh))
	print('\n')
	for band in obands:
		print('\t'+band+': '+str(mags[band])+' +/- '+str(emags[band]))

	# Correct GAIA parallax
	if gaia_corr:
		parallax += 0.082
		eparallax = eparallax * 10. / 7.

	# Get distance from parallax (using random normal sampling for error)
	diss = 1. / np.random.normal(parallax/1000.,eparallax/1000.,nsim)
	print('\n\t D [pc] = ',np.mean(diss),' +/- ', np.sqrt(np.var(diss)))

	# Create data and error arrays with teff, Fe/H, vec
	# TODO what is vec???
	data = [teff]
	error = [eteff]
	for band in obands:
		ms  = np.random.normal(mags[band],emags[band],nsim)
		vec = diss*10**(-0.2*(ms+5))
		data.append(np.mean(vec))
		error.append(np.sqrt(np.var((vec))))
	data = np.array(data)
	error = np.array(error)

	# load isochrone data
	dct = pickle.load(open(isopath+'parsec.pkl','rb'), encoding = 'latin1')

	# set number of walkers for mcmc
	nwalkers =10

	# initialize values
	pos = []
	# initial values for age, mass, Av
	guess = [2.,1.2,0.1]
	# If metallicity is free initialize that too
	if feh_free:
		guess.append(feh)
	ndim = len(guess)

	# loop over the number of walkers
	while len(pos) < nwalkers:
		# set initial guesses for each walker
		# by adding random samples from gaussian
		vala = guess[0] + 0.1*np.random.randn()
		valm = guess[1] + 0.1*np.random.randn()
		valv = guess[2] + 0.1*np.random.randn()
		if feh_free:
			valf = guess[3] + 0.1*np.random.randn()
		if feh_free:
			pos.append(np.array([vala,valm,valv,valf]))
		else:
			pos.append(np.array([vala,valm,valv]))

	# Instantiate the emcee sampler
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(data, error))
	# iterate the mcmc over 10000 iterations
	sampler.run_mcmc(pos, stepnumber, progress = True)
	
	
	fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
	samples = sampler.get_chain()
	labels = ["Age", "Mass", "AV"]
	for i in range(ndim):
    		ax = axes[i]
    		ax.plot(samples[:, :, i], "k", alpha=0.3)
    		ax.set_xlim(0, len(samples))
    		ax.set_ylabel(labels[i])
    		ax.yaxis.set_label_coords(-0.1, 0.5)

	axes[-1].set_xlabel("step number");
	fig.savefig(pars_path+obname+"_walker_out.png")
	plt.show()
	
	burnIn = input("Please Enter BurnIn phase (previous steps in chain get discarded): ") 
	burnIn = int(burnIn) 
	
	# get the results, discarding first 50 chains (burn-in?)
	samples = sampler.chain[:, burnIn:, :].reshape((-1, ndim))
	print(np.shape(samples))
	
	

	# create and save corner plots of the mcmc results
	if feh_free:
		fig = corner.corner(samples, labels=["$AGE$", "$MASS$", "$Av$","[Fe/H]"])
	else:
		fig = corner.corner(samples, labels=["$AGE$", "$MASS$", "$Av$"])
	fig.savefig(pars_path+obname+"_out.png")

	# created sorted arrays for age, initial mass, Av, Fe/H if applicable
	fages   = np.sort(samples[:,0])
	fimasses = np.sort(samples[:,1])
	fav = np.sort(samples[:,2])
	if feh_free:
		fehs = np.sort(samples[:,3])

	# initialize mass, Teff, logg, Rad, Lum, Rho, lists
	fmasses,ftefs, floggs, frads, flums, rhos = [],[],[],[],[],[]
	# loop over the mcmc chains
	for i in range(len(samples)):
		# use the init mass, log10(age*10**8), Fe/H from each chain
		# to get mass, logL, logT, logg by interpolating the isochrones
		results = get_vals(samples[i,1],np.log10(samples[i,0]*1e9),feh,keys =['mass','logL','logT','logg'])
		# save interpolated values
		floggs.append(results['logg'])
		ftefs.append(10**results['logT'])
		flums.append(10**results['logL'])
		fmasses.append(results['mass'])
		# compute stellar radius from mass and logg
		Rs = loggtoRs(results['logg'],results['mass'])
		# save logg again - TODO check why!!
		floggs.append(results['logg'])
		# save radius
		frads.append(Rs)
		# comupte and save Rho
		rhos.append(rhostar(results['mass'],Rs))

	# convert results lists to arrays
	fmasses,ftefs, floggs, frads, flums, rhos = np.array(fmasses), np.array(ftefs), np.array(floggs), np.array(frads), np.array(flums), np.array(rhos)
	# get the median and 1-sigma error bars for all the parameters
	fage,age1,age2 = get_sum(fages)
	fimass,imass1,imass2 = get_sum(fimasses)
	fmass,mass1,mass2 = get_sum(fmasses)
	ftef,tef1,tef2 = get_sum(ftefs)
	flogg, logg1, logg2 = get_sum(floggs)
	frad,rad1,rad2 = get_sum(frads)
	flum,lum1,lum2 = get_sum(flums)
	rho,rho1,rho2 = get_sum(rhos)
	av,av1,av2 = get_sum(fav)
	if feh_free:
		fehv,feh1,feh2 = get_sum(fehs)

	# print all parameters and intervals to screen
	print( '\n')
	print('IMASS =', fimass, '(',imass1,imass2,')')
	print('AGE =', fage, '(',age1, age2,')')
	print('MASS =', fmass, '(',mass1,mass2,')')
	print('Teff =', ftef, '(',tef1, tef2,')')
	print('log(g) =', flogg, '(',logg1, logg2,')')
	print('Rs =', frad, '(',rad1,rad2,')')
	print('L =', flum, '(',lum1, lum2,')')
	print('Rho =', rho, '(',rho1,rho2,')')
	print('Av =', av, '(',av1,av2,')')
	if feh_free:
		print('[Fe/H] =', fehv, '(',feh1,feh2,')')

	# creat save dictionary
	dct = {}
	# add all the parameters to dictionary
	dct['imass'] = fimasses
	dct['age'] = fages
	dct['mass'] = fmasses
	dct['teff'] = ftefs
	dct['logg'] = floggs
	dct['rstar'] = frads
	dct['lstar'] = flums
	dct['rhostar'] = rhos
	dct['av'] = fav
	if feh_free:
		dct['feh'] = fehs
	# save the dictionary to file
	pickle.dump(dct,open(pars_path+'/'+obname+"_samples.pkl",'wb'))

def trunc(val,force=-999):
	"""
	Round decimals with automatic selection of decimals to round to

	Parameters:
	val (float): parameter to be rounded
	force (float): optional, decimals to round to

	Returns:
	ret: rounded value
	ret2: decimal position rounded to (not returned if force given)
	"""
	# if the position to round to is not forced
	if force == -999:
		# convert input value to string
		st = str(val)
		j = 0
		# check if input value is smaller than one
		if st[0] == '0':
			for i in range(len(st)):
				# find first significant decimal
				if st[i] != '0' and st[i]!='.':
					j = i
					break
			# round to second significant decimal
			ret = np.around(val,j)
			# save position rounded to
			ret2 = j
		# if inupt value is larger than one
		else:
			for i in range(len(st)):
				# find decimal point
				if st[i] == '.':
					j=i
					break
			# round to 2 - decimal point index
			ret =  np.around(val,2-j)
			# save position rounded to
			ret2  = 2 - j
		# check if decimal part is zero
		if str(ret)[-2:] == '.0':
			# if yes convert to int, return rounded value and pos
			return int(ret), ret2
		else:
			# return rounded value and position rounded to
			return ret,ret2

	# if the decimal position is forced
	else:
		# convert input value to string
		st = str(val)
		# round it to forced position
		ret = np.around(val,force)
		#print(ret)
		# check if decimal part is zero
		if str(ret)[-2:] == '.0':
			# if yes convert to int
			return int(ret)
		else:
			# return rounded value
			return ret



def planet_props(star_sample, juliet_samples, connect_rho=False,mps=False,smf=1.4):
	"""
	Compute planet parameters from juliet planet fit posteriors and
	parsec stellar parameters posteriors.

	Parameters:
	star_sample (str): path to parsec star parameter posteriors pkl file
	juliet_samples (str): path to juliet fit posterios pkl file
	connect_rho (boolean): use parsec posterior on rho to constrain the juliet one
	mps (boolean): K in m/s
	smf (float): factor by which to broaden the random noise added to mass, rho

	Returns:
	(to screen): Mp, Rp, a, Teq values
	"""
	#dp = pickle.load(open(juliet_samples))
	# load the juliet planet fit posteriors
	with open(juliet_samples, 'rb') as f:
		dp = pickle.load(f, encoding = 'latin1')

	# Check if the fit used the Espinoza 18 r1,r2 parametrization
	if 'pu' in dp.keys():
		# get the maximum planet-to-star radius ratio
		pu = dp['pu']
		# get the smallest planet-to-star radius ratio
		pl = dp['pl']
		Ar = (pu - pl)/(2. + pl + pu)
		# set a reference key
		if 'r1_p1' in dp['posterior_samples'].keys():
			trans = True
		else:
			trans = False
	else:
		trans = False

	# keep only the posterior samples from the juliet fit
	dp = dp['posterior_samples']

	# Check if the fit used transits or not
	if 'a_p' in dp.keys() or 'rho' in dp.keys():
		rv_only = False
	else:
		rv_only = True

	# load the parsec stellar parameters posteriors
	ds = pickle.load(open(star_sample, 'rb'), encoding = 'latin1')

	# put stellar parameters into arrays
	mstar   = ds['mass']
	rstar   = ds['rstar']
	rhostar = ds['rhostar']
	lstar   = ds['lstar']
	age	 = ds['age']
	# create indices array
	I=np.arange(len(mstar))
	# select indices
	#I = np.where((age<5.3)&(mstar>1.2))[0]
	# cut parameter arrays to I indices
	# (currently no effect as no filter on I)
	mstar,rstar,rhostar,lstar,age = mstar[I],rstar[I],rhostar[I],lstar[I],age[I]
	# get standard deviation of stellar masses
	devmstar = np.sqrt(np.var(mstar))
	# add random gaussian noise to masses using std-dev and smf factor
	mstar = mstar + np.random.normal(0,smf*devmstar,len(mstar))
	# get standard deviation of stellar densities
	devrhostar = np.sqrt(np.var(rhostar))
	# add random gaussian noise to densities using std-dev and smf factor
	rhostar = rhostar + np.random.normal(0,smf*devrhostar,len(rhostar))
	# get the median and 1-sigma error bars for stellar mass
	fms,ms1,ms2 = get_sum(mstar)
	# get rounded values and decimal positions rounded to for error bars
	ms1,t1 = trunc(ms1)
	ms2,t2 = trunc(ms2)
	# round median mass to maximum number of decimals kept for error bars
	fms = trunc(fms,force=max(t1,t2))
	# print rounded mass and error bars
	print('\nMs =', fms, '(',ms1, ms2,')')
	# get the median and 1-sigma error bars for stellar radius
	frs,rs1,rs2 = get_sum(rstar)
	# get rounded values and decimal positions rounded to for error bars
	rs1,t1 = trunc(rs1)
	rs2,t2 = trunc(rs2)
	# round median radius to maximum number of decimals kept for error bars
	frs = trunc(frs,force=max(t1,t2))
	# print rounded radius and error bars
	print('Rs =', frs, '(',rs1, rs2,')')
	# get the median and 1-sigma error bars for stellar density
	fds,ds1,ds2 = get_sum(rhostar)
	# get rounded values and decimal positions rounded to for error bars
	ds1,t1 = trunc(ds1)
	ds2,t2 = trunc(ds2)
	# round median density to maximum number of decimals kept for error bars
	fds = trunc(fds,force=max(t1,t2))
	# print rounded density and error bars
	print('RHOs =', fds, '(',ds1, ds2,')')
	# get the median and 1-sigma error bars for stellar luminosity
	fls,ls1,ls2 = get_sum(lstar)
	# get rounded values and decimal positions rounded to for error bars
	ls1,t1 = trunc(ls1)
	ls2,t2 = trunc(ls2)
	# round median luminosity to maximum number of decimals kept for error bars
	fls = trunc(fls,force=max(t1,t2))
	# print rounded luminosity and error bars
	print('Ls =', fls, '(',ls1, ls2,')')
	# get the median and 1-sigma error bars for stellar age
	fas,as1,as2 = get_sum(age)
	# get rounded values and decimal positions rounded to for error bars
	as1,t1 = trunc(as1)
	as2,t2 = trunc(as2)
	# round median age to maximum number of decimals kept for error bars
	fas = trunc(fas,force=max(t1,t2))
	# print rounded age and error bars
	print('AGEs =', fas, '(',as1, as2,')\n')

	# get the keys for the planet fit posteriors
	keys = dp.keys()
	# intialize number of planets npl
	npl = 0
	# loop through the keys
	for key in keys:
		# loop through possible number of planets (1-7)
		for i in [1,2,3,4,5,6,7]:
			si = str(i)
			# if planet i is in this key and npl is smaller than i
			if '_p'+si in key and npl < i:
				# set number of planets to i
				npl = i
	# print number of planets
	print('Number of planets:', npl)

	# create array of planet indices
	pls = np.arange(npl).astype('int') + 1

	# loop through the planets
	for plt in pls:
		# print planet being processed
		print('Planet', plt)
		# convert period to seconds
		P = dp['P_p'+str(plt)] * 24. * 3600.
		# convert K to centimetres assuming km/s
		K = dp['K_p'+str(plt)] * 100000.
		# if K was in m/s fix conversion
		if mps:
			K /= 1000.

		# check if the fit included light curves
		if not rv_only:
			# check if the scaled semi-major axis a/Rs was fitted
			if 'a_p'+str(plt) in dp.keys():
				# if yes obtain it from dict
				a = dp['a_p'+str(plt)]
			# if not estimate it now
			else:
				# get the stellar density prior given to juliet
				rho = dp['rho']
				# compute a from Kepler's third law
				a = rho**(1./3.) * (cons.G.value * P * P /(3.*np.pi))**(1./3.)
		# if it didn't sample from rhostar
		else:
			# change to cgs
			rhostar_cgs = rhostar = rhostar*1000
			# sample len(P) values from rhostar
			rho = np.random.choice(rhostar_cgs, size = len(P), replace = False)
			# compute a from Kepler's third law
			a = rho**(1./3.) * (cons.G.value * P * P /(3.*np.pi))**(1./3.)

		# check if the eccentricity was fitted
		if 'ecc_p'+str(plt) in dp.keys():
			# if yes obtain it from dict
			e = dp['ecc_p'+str(plt)]
		# if not set it to zero
		else:
			e = np.zeros(len(P))

		# If light curves were used, compute the inclination
		if not rv_only:
			# Check if the fit used the Espinoza 18 r1,r2 parametrization
			if trans:
				# get r1, r2
				r1 = dp['r1_p'+str(plt)]
				r2 = dp['r2_p'+str(plt)]
				# initialize b impact parameter, p Rp/Rs ratio
				b,p = np.zeros(len(r1)),np.zeros(len(r1))
				# compute b, p from r1, r2
				for i in range(len(r1)):
					if r1[i] > Ar:
						b[i] = (1.+pl)*(1. + (r1[i]-1.)/(1.-Ar))
						p[i] = (1.-r2[i])*pl + r2[i]*pu
					else:
						b[i],p[i] = (1. + pl) + np.sqrt(r1[i]/Ar)*r2[i]*(pu-pl), pu + (pl-pu)*np.sqrt(r1[i]/Ar)*(1.-r2[i])
			# if not get b and p from dict
			else:
				p = dp['p_p'+str(plt)]
				b = dp['b_p'+str(plt)]

			# compute the orbital inclination
			inc = np.arccos(b/a)

		# else set it to 90
		inc = np.ones(np.shape(a))*np.pi/2.

		# initialize a, Mp, Rp, Teq
		As,MPs,RPs,TEQs = [],[],[],[]

		# compute stellar densities in cgs
		# TODO why recompute??
		densities = a**3 * 3*np.pi/(cons.G.cgs.value*P*P)
		#print(densities)

		# loop over the length of P (ie number of juliet samples)
		for i in range(len(P)):
			# if connecting the densities
			# this needs to be false if rv only
			if connect_rho and not rv_only:
				# get difference between rho from parsec and densities
				res = np.absolute(rhostar - densities[i])
				# get indices of 20 closest values
				I = np.argsort(res)[:20]
				# shuffle index values randomly
				i2 = I[np.random.randint(len(I))]
			# if not connecting the densities
			else:
				# shuffle all parsec indices
				i2 = np.random.randint(len(rstar))

			# use the i-th planet posteriors and the random i2 star posteriors
			# get the planet mass for this (i,i2) set
			MP = get_mass(mstar[i2],P[i],e[i],K[i],inc[i])#*cons.M_jup.cgs.value
			if not rv_only:
				# get the semi-major axis for this (i,i2) set
				A  = a[i]*rstar[i2] * cons.R_sun.value / cons.au.value
				# get the planet radius for this (i,i2) set if not rv-only
				RP = rstar[i2] * p[i] * cons.R_sun.value / cons.R_jup.value
			else:
				# get the semi-major axis for this (i, i2) set if RV-only
				# get cgs constants
				Msun = cons.M_sun.cgs.value
				G = cons.G.cgs.value
				Mjup = cons.M_jup.cgs.value
				# compute from K
				A = (G/(1-e[i]**2))*((MP*Mjup)**2/((mstar[i2]*Msun+MP*Mjup)*K[i]**2))
				# convert to AU
				A = A/cons.au.cgs.value
			# save a, Rp, Mp to lists
			As.append(A)
			if not rv_only:
				RPs.append(RP)
			MPs.append(MP)
			# get and save Teq for this (i,i2) set
			TEQs.append(get_teq(lstar[i2],e[i],A,A=0.,beta=0.5,eps=1.))

		# convert parameter lists to arrays
		As,MPs,TEQs = np.array(As), np.array(MPs), np.array(TEQs)
		if not rv_only:
			RPs = np.array(RPs)
		#print(MPs)
		# get median values and 1-sigma error bars for Mp, Rp, a, Teq
		fmp,mp1,mp2 = get_sum(MPs)
		if not rv_only:
				frp,rp1,rp2 = get_sum(RPs)
		fa,a1,a2	= get_sum(As)
		ftq,tq1,tq2    = get_sum(TEQs)

		if not rv_only:
			# Create a histogram of the radii
			hist(RPs)

		# get rounded values and decimal positions rounded to for Mp error bars
		mp1,t1 = trunc(mp1)
		mp2,t2 = trunc(mp2)
		# round median Mp to maximum number of decimals kept for error bars
		fmp = trunc(fmp,force=max(t1,t2))

		if not rv_only:
			# get rounded values and decimal positions rounded to for Rp error bars
			rp1,t1 = trunc(rp1)
			rp2,t2 = trunc(rp2)
			# round median Rp to maximum number of decimals kept for error bars
			frp = trunc(frp,force=max(t1,t2))

		# get rounded values and decimal positions rounded to for a error bars
		a1,t1 = trunc(a1)
		a2,t2 = trunc(a2)
		# round median a to maximum number of decimals kept for error bars
		fa = trunc(fa,force=max(t1,t2))

		# get rounded values and decimal positions rounded to for Teq error bars
		tq1,t1 = trunc(tq1)
		tq2,t2 = trunc(tq2)
		# round median Teq to maximum number of decimals kept for error bars
		ftq = trunc(ftq,force=max(t1,t2))

		# print rounded median values and 1-sigma error bars
		print('\n')
		print('a =', fa, '(',a1,a2,')')
		print('Mp =', fmp, '(',mp1, mp2,')')
		if not rv_only:
			print('Rp =', frp, '(',rp1,rp2,')')
		print('Teq =', ftq, '(',tq1,tq2,')')
		#show()

def planet_props_lc(star_sample, juliet_samples, connect_rho=False,mps=False,smf=1.4):
	"""
	Compute planet parameters from juliet planet fit posteriors and
	parsec stellar parameters posteriors - lc only.

	Parameters:
	star_sample (str): path to parsec star parameter posteriors pkl file
	juliet_samples (str): path to juliet fit posterios pkl file
	connect_rho (boolean): use parsec posterior on rho to constrain the juliet one
	mps (boolean): K in m/s
	smf (float): factor by which to broaden the random noise added to mass, rho

	Returns:
	(to screen): Mp, Rp, a, Teq values
	"""
	#dp = pickle.load(open(juliet_samples))
	# load the juliet planet fit posteriors
	with open(juliet_samples, 'rb') as f:
		dp = pickle.load(f, encoding = 'latin1')

	# Check if the fit used the Espinoza 18 r1,r2 parametrization
	if 'pu' in dp.keys():
		# get the maximum planet-to-star radius ratio
		pu = dp['pu']
		# get the smallest planet-to-star radius ratio
		pl = dp['pl']
		Ar = (pu - pl)/(2. + pl + pu)
		# set a reference key
		if 'r1_p1' in dp['posterior_samples'].keys():
			trans = True
		else:
			trans = False
	else:
		trans = False

	# keep only the posterior samples from the juliet fit
	dp = dp['posterior_samples']

	# load the parsec stellar parameters posteriors
	ds = pickle.load(open(star_sample, 'rb'), encoding = 'latin1')

	# put stellar parameters into arrays
	mstar   = ds['mass']
	rstar   = ds['rstar']
	rhostar = ds['rhostar']
	lstar   = ds['lstar']
	age	 = ds['age']
	# create indices array
	I=np.arange(len(mstar))
	# select indices
	#I = np.where((age<5.3)&(mstar>1.2))[0]
	# cut parameter arrays to I indices
	# (currently no effect as no filter on I)
	mstar,rstar,rhostar,lstar,age = mstar[I],rstar[I],rhostar[I],lstar[I],age[I]
	# get standard deviation of stellar masses
	devmstar = np.sqrt(np.var(mstar))
	# add random gaussian noise to masses using std-dev and smf factor
	mstar = mstar + np.random.normal(0,smf*devmstar,len(mstar))
	# get standard deviation of stellar densities
	devrhostar = np.sqrt(np.var(rhostar))
	# add random gaussian noise to densities using std-dev and smf factor
	rhostar = rhostar + np.random.normal(0,smf*devrhostar,len(rhostar))
	# get the median and 1-sigma error bars for stellar mass
	fms,ms1,ms2 = get_sum(mstar)
	# get rounded values and decimal positions rounded to for error bars
	ms1,t1 = trunc(ms1)
	ms2,t2 = trunc(ms2)
	# round median mass to maximum number of decimals kept for error bars
	fms = trunc(fms,force=max(t1,t2))
	# print rounded mass and error bars
	print('\nMs =', fms, '(',ms1, ms2,')')
	# get the median and 1-sigma error bars for stellar radius
	frs,rs1,rs2 = get_sum(rstar)
	# get rounded values and decimal positions rounded to for error bars
	rs1,t1 = trunc(rs1)
	rs2,t2 = trunc(rs2)
	# round median radius to maximum number of decimals kept for error bars
	frs = trunc(frs,force=max(t1,t2))
	# print rounded radius and error bars
	print('Rs =', frs, '(',rs1, rs2,')')
	# get the median and 1-sigma error bars for stellar density
	fds,ds1,ds2 = get_sum(rhostar)
	# get rounded values and decimal positions rounded to for error bars
	ds1,t1 = trunc(ds1)
	ds2,t2 = trunc(ds2)
	# round median density to maximum number of decimals kept for error bars
	fds = trunc(fds,force=max(t1,t2))
	# print rounded density and error bars
	print('RHOs =', fds, '(',ds1, ds2,')')
	# get the median and 1-sigma error bars for stellar luminosity
	fls,ls1,ls2 = get_sum(lstar)
	# get rounded values and decimal positions rounded to for error bars
	ls1,t1 = trunc(ls1)
	ls2,t2 = trunc(ls2)
	# round median luminosity to maximum number of decimals kept for error bars
	fls = trunc(fls,force=max(t1,t2))
	# print rounded luminosity and error bars
	print('Ls =', fls, '(',ls1, ls2,')')
	# get the median and 1-sigma error bars for stellar age
	fas,as1,as2 = get_sum(age)
	# get rounded values and decimal positions rounded to for error bars
	as1,t1 = trunc(as1)
	as2,t2 = trunc(as2)
	# round median age to maximum number of decimals kept for error bars
	fas = trunc(fas,force=max(t1,t2))
	# print rounded age and error bars
	print('AGEs =', fas, '(',as1, as2,')\n')

	# get the keys for the planet fit posteriors
	keys = dp.keys()
	# intialize number of planets npl
	npl = 0
	# loop through the keys
	for key in keys:
		# loop through possible number of planets (1-7)
		for i in [1,2,3,4,5,6,7]:
			si = str(i)
			# if planet i is in this key and npl is smaller than i
			if '_p'+si in key and npl < i:
				# set number of planets to i
				npl = i
	# print number of planets
	print('Number of planets:', npl)

	# create array of planet indices
	pls = np.arange(npl).astype('int') + 1

	# loop through the planets
	for plt in pls:
		# print planet being processed
		print('Planet', plt)
		# convert period to seconds
		P = dp['P_p'+str(plt)] * 24. * 3600.

		# check if the scaled semi-major axis a/Rs was fitted
		if 'a_p'+str(plt) in dp.keys():
			# if yes obtain it from dict
			a = dp['a_p'+str(plt)]
		# if not estimate it now
		else:
			# get the stellar density prior given to juliet
			rho = dp['rho']
			# compute a from Kepler's third law
			a = rho**(1./3.) * (cons.G.value * P * P /(3.*np.pi))**(1./3.)

		# check if the eccentricity was fitted
		if 'ecc_p'+str(plt) in dp.keys():
			# if yes obtain it from dict
			e = dp['ecc_p'+str(plt)]
		# if not set it to zero
		else:
			e = np.zeros(len(P))

		# If light curves were used, compute the inclination
		# Check if the fit used the Espinoza 18 r1,r2 parametrization
		if trans:
			# get r1, r2
			r1 = dp['r1_p'+str(plt)]
			r2 = dp['r2_p'+str(plt)]
			# initialize b impact parameter, p Rp/Rs ratio
			b,p = np.zeros(len(r1)),np.zeros(len(r1))
			# compute b, p from r1, r2
			for i in range(len(r1)):
				if r1[i] > Ar:
					b[i] = (1.+pl)*(1. + (r1[i]-1.)/(1.-Ar))
					p[i] = (1.-r2[i])*pl + r2[i]*pu
				else:
					b[i],p[i] = (1. + pl) + np.sqrt(r1[i]/Ar)*r2[i]*(pu-pl), pu + (pl-pu)*np.sqrt(r1[i]/Ar)*(1.-r2[i])
		# if not get b and p from dict
		else:
			p = dp['p_p'+str(plt)]
			b = dp['b_p'+str(plt)]

		# compute the orbital inclination
		inc = np.arccos(b/a)


		# initialize a, Mp, Rp, Teq
		As,RPs,TEQs = [],[],[]

		# compute stellar densities in cgs
		# TODO why recompute??
		densities = a**3 * 3*np.pi/(cons.G.cgs.value*P*P)
		#print(densities)

		# loop over the length of P (ie number of juliet samples)
		for i in range(len(P)):
			# if connecting the densities
			# this needs to be false if rv only
			if connect_rho:
				# get difference between rho from parsec and densities
				res = np.absolute(rhostar - densities[i])
				# get indices of 20 closest values
				I = np.argsort(res)[:20]
				# shuffle index values randomly
				i2 = I[np.random.randint(len(I))]
			# if not connecting the densities
			else:
				# shuffle all parsec indices
				i2 = np.random.randint(len(rstar))

			# use the i-th planet posteriors and the random i2 star posteriors
			# get the semi-major axis for this (i,i2) set
			A  = a[i]*rstar[i2] * cons.R_sun.value / cons.au.value
			# get the planet radius for this (i,i2) set if not rv-only
			RP = rstar[i2] * p[i] * cons.R_sun.value / cons.R_jup.value
			# save a, Rp, Mp to lists
			As.append(A)
			RPs.append(RP)
			# get and save Teq for this (i,i2) set
			TEQs.append(get_teq(lstar[i2],e[i],A,A=0.,beta=0.5,eps=1.))

		# convert parameter lists to arrays
		As,RPs,TEQs = np.array(As), np.array(RPs), np.array(TEQs)
		frp,rp1,rp2 = get_sum(RPs)
		fa,a1,a2	= get_sum(As)
		ftq,tq1,tq2    = get_sum(TEQs)

		# get rounded values and decimal positions rounded to for Rp error bars
		rp1,t1 = trunc(rp1)
		rp2,t2 = trunc(rp2)
		# round median Rp to maximum number of decimals kept for error bars
		frp = trunc(frp,force=max(t1,t2))

		# get rounded values and decimal positions rounded to for a error bars
		a1,t1 = trunc(a1)
		a2,t2 = trunc(a2)
		# round median a to maximum number of decimals kept for error bars
		fa = trunc(fa,force=max(t1,t2))

		# get rounded values and decimal positions rounded to for Teq error bars
		tq1,t1 = trunc(tq1)
		tq2,t2 = trunc(tq2)
		# round median Teq to maximum number of decimals kept for error bars
		ftq = trunc(ftq,force=max(t1,t2))

		# print rounded median values and 1-sigma error bars
		print('\n')
		print('a =', fa, '(',a1,a2,')')
		print('Rp =', frp, '(',rp1,rp2,')')
		print('Teq =', ftq, '(',tq1,tq2,')')
		#show()

def make_plot(syst,masses = [1.01,1.14,1.3],isopath=default_isopath,pars_path='pars/'):
	"""
	Make and save plot

	Parameters:
	syst (str): system to plot
	masses (array): test masses
	isopath (string): path to isochrones
	pars_path (string): path to folder with stellar parameter files

	Returns
	(to pdf): plot of R vs Teff
	"""

	# Set up plot parameters
	rc('font', **{'family': 'Helvetica'})
	rc('text', usetex=True)
	matplotlib.rcParams.update({'font.size':22})
	rc('legend', **{'fontsize':7})

	# make the dct keyword global
	global dct

	# decide whether to use GAIA correction or not
	gaia_corr = True
	# set number of resamples from normal distributions
	nsim = 1000
	# set photometric bands to query
	obands = ['Gmag','BPmag','RPmag','Jmag','Hmag','Kmag']
	# initialize magnitude and error dictionaries
	mags = {}
	emags = {}

	# load the isochrones
	dct = pickle.load(open(isopath+'parsec.pkl','rb'), encoding = 'latin1')
	# create array of ages with 0.01 step
	ages1 = np.arange(1.01,9.9,0.01)
	# create array of ages with 1 step
	ages2 = np.arange(1.01,9.9,1.)
	# open txt file with Teff, Fe/H, GAIA source ID for the system
	fl = open(pars_path+syst+'.txt','r').readlines()
	# read Teff, Fe/H, GAIA source ID
	for line in fl:
			cos = line.split()
			if cos[0] == 'teff':
				teff = float(cos[1])
				eteff = float(cos[2])
			if cos[0] == 'feh':
				feh = float(cos[1])
				efeh = float(cos[2])
			if cos[0] == 'source_id':
				sid = cos[1]

	# Query GAIA DR3 for coords, parallax, magnitudes using source ID
	job = Gaia.launch_job("select top 500 source_id,ra,dec,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_flux,phot_g_mean_flux_error,phot_g_mean_mag,phot_bp_mean_flux,phot_bp_mean_flux_error,phot_bp_mean_mag,phot_rp_mean_flux,phot_rp_mean_flux_error,phot_rp_mean_mag,radial_velocity,radial_velocity_error,teff_gspphot,ag_gspphot FROM gaiadr3.gaia_source  where (source_id="+sid+")")
	r = job.get_results()
	ra,dec = r['ra'], r['dec']
	parallax,eparallax = r['parallax'], r['parallax_error']
	pmra,epmra = r['pmra'], r['pmra_error']
	pmdec,epmdec = r['pmdec'], r['pmdec_error']
	ra2000,dec2000 = ra - pmra*15.5/(1000.*3600.), r['dec'] - pmdec*15.5/(1000.*3600.)

	# loop through the bands
	for band in obands:
		key = None
		if band == 'Gmag':
			key = 'phot_g_mean'
			systematic = 0.002
		elif band == 'BPmag':
			key = 'phot_bp_mean'
			systematic = 0.005

		elif band == 'RPmag':
			key = 'phot_rp_mean'
			systematic = 0.003

		# For GAIA magnitudes compute the error from the flux
		if key != None:
			mags[band] = float(r[key+'_mag'])
			e1 = flux_to_mag(float(r[key+'_flux'])+float(r[key+'_flux_error']),band)
			e2 = flux_to_mag(float(r[key+'_flux'])-float(r[key+'_flux_error']),band)
			err = 0.5*(e2 - e1)
			emags[band] = np.sqrt(err**2+systematic**2)

		# For 2MASS band query Vizier
		if band in ['Jmag','Hmag','Kmag']:
			mag,emag = get_mags_2mass(ra2000[0],dec2000[0],band)
			mags[band] = mag
			emags[band] = emag

	# load stellar posterior samples
	samples = pickle.load(open(pars_path+syst+'_samples.pkl','rb'), encoding = 'latin1')
	# interpolate values for log T and magnitudes from stellar parameters and isochrones
	out = get_vals(np.mean(samples['mass']), np.log10(np.mean(samples['age'])*1e9), feh, isopath=default_isopath,keys=['logT','Gmag_0.0','BPmag_0.0','RPmag_0.0','Jmag_0.0','Hmag_0.0','Kmag_0.0','Gmag_1.0','BPmag_1.0','RPmag_1.0','Jmag_1.0','Hmag_1.0','Kmag_1.0'])
	# put logT into array
	model = np.array([10**out['logT']])
	av1,av2 = 0.0,1.0
	# loop over the bands
	for band in obands:
		# correct for extinction? (TODO check)
		m = (out[band+'_'+str(av2)] - out[band+'_'+str(av1)]) / float(av2-av1)
		n = out[band+'_'+str(av2)] - m * av2
		ys = 10**(-0.2*(m*np.mean(samples['av']) + n))

		model = np.hstack((model,np.array([ys])))

	nsim=1000

	# correct GAIA parallax
	if gaia_corr:
		parallax += 0.082
		eparallax = eparallax * 10. / 7.

	# Get distance from parallax (using random normal sampling for error)
	diss = 1. / np.random.normal(parallax/1000.,eparallax/1000.,nsim)
	print('\n\t D [pc] = ',np.mean(diss),' +/- ', np.sqrt(np.var(diss)))

	# Create data and error arrays with teff, Fe/H, TODO
	# TODO what is vec???
	data = [teff]
	error = [eteff]
	for band in obands:
		ms  = np.random.normal(mags[band],emags[band],nsim)
		vec = diss*10**(-0.2*(ms+5))
		data.append(np.mean(vec))
		error.append(np.sqrt(np.var((vec))))
	data = np.array(data)
	error = np.array(error)
	# print values
	print(data, error,model,)
	print(data - model)
	print((data - model)/error)
	print(fds)
	clf()

	# get mass, logT, logg for 1.14 mass, np.log10(6.7*1e9) age, 0.26 metallicity
	# TODO why these values?
	out =  get_vals(1.14, np.log10(6.7*1e9), 0.26, isopath=default_isopath,keys=['mass','logT','logg'])
	#print(out['mass'],10**out['logT'],out['logg'],loggtoRs(out['logg'],out['mass']))
	#print(gfd)

	# loop over test masses
	for mass in masses:
		# print the mass
		print(mass)
		# set up lists for R and Teff
		rads,teffs =[],[]
		rads2,teffs2 = [],[]
		# loop over fine array of ages
		for age in ages1:
			# get mass, logT, logg from isochrones for these values
			out = get_vals(mass, np.log10(age*1e9), feh, isopath=default_isopath,keys=['mass','logT','logg'])
			try:
				# compute and save R, Teff
				Rs = loggtoRs(out['logg'],out['mass'])
				teff = 10**out['logT']
				teffs.append(teff)
				rads.append(Rs)
			except:
				# if they fail do nothing
				'nothing'
		# loop over broad array of ages
		for age in ages2:
			# get mass, logT, logg from isochrones for these values
			out = get_vals(mass, np.log10(age*1e9), feh, isopath=default_isopath,keys=['mass','logT','logg'])
			try:
				# compute and save R, Teff
				Rs = loggtoRs(out['logg'],out['mass'])
				teff = 10**out['logT']
				teffs2.append(teff)
				rads2.append(Rs)
			except:
				# if they fail do nothing
				'nothing'
		rads,teffs = np.array(rads),np.array(teffs)
		# plot the R and Teff computed with the fine array of ages
		plot(teffs,rads)
		# overplot (as scatter) the R and Teff computed with the broad array of ages
		scatter(teffs2,rads2)

	# get the radius and Teff from the stellar posteriors
	st,sr = samples['teff'],samples['rstar']
	# plot the posteriors (save first 5000 points)
	plot(st[-5000:],sr[-5000:],'k.',alpha=0.01)
	# set axis and limits
	gca().invert_xaxis()
	ylim([0.8,2.5])
	# set labels
	xlabel(r'T$_{eff}$ [K]')
	ylabel(r'R$_{\star}$ [R$_{\odot}$]')
	# save the plot
	savefig('test.pdf',dpi=300,bbox_inches='tight')
