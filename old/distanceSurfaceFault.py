import math
import numpy as np
#import scipy as scp
#from scipy import optimize
#import matplotlib.pyplot as plt
#from matplotlib.mlab import griddata

def MinDistToFault(lat0, lon0, H, strike0, dip, Mw, lonE, latE):
    '''
    Funcion mediante la cual se calculara la minima distancia Rrup a un plano de ruptura dado
    '''
    L = 10**(-2.9 + 0.63*Mw) # Calculamos el Largo [Km] teorico del plano de falla a partir de esta relacion de escala
    W = 10**(-0.86 + 0.35*Mw) # Calculamos el Ancho [Km] teorico del plano de falla a partir de esta relacion de escala
    
    # Con el fin de conseguir el menor error posible se consideraran subfallas de ~1[Km], por lo que se utilizara la siguiente estamcion de subfallas en el alrgo y ancho por terremoto.
    nx = int(L) # Se considera el entero del largo y ancho en kilometros, para asi tener subfallas de ~1[Km]
    ny = int(W)
    
    L = L * 1000 # Pasamos las unidades a metros
    W = W * 1000
    H = H * 1000
    # Procedemos a calcular las coordenadas geograficas de los 4 extremos que daran lugar a nuestro plano de falla
    latt1, lont1, alpha21 = vinc_pt(lat0,lon0, (strike0 + 90), (W/2)*math.cos(math.radians(dip)))
    latt2, lont2, alpha21 = vinc_pt(lat0,lon0, (strike0 - 90), (W/2)*math.cos(math.radians(dip)))
    # Los 4 puntos que definen la falla grande correponden a los siguientes:
    lat1, lon1, alpha21 = vinc_pt(latt1,lont1, strike0, L/2. )
    lat2, lon2, alpha21 = vinc_pt(latt1,lont1, (strike0 - 180.), L/2. )
    lat3, lon3, alpha21 = vinc_pt(latt2,lont2, strike0, L/2. )
    lat4, lon4, alpha21 = vinc_pt(latt2,lont2, (strike0 - 180.), L/2. )
    # Pasamos a calcular las coordenadas de cada una de las subfallas
    l = L/nx
    w = W/ny
    Hmax = H + (W/2)*math.sin(math.radians(dip))
    DistSubFallas = np.zeros(nx*ny).reshape(nx,ny)
    epsilon = 1e-7
    for ix in range(nx):
        lla1, llo1, alpha21 = vinc_pt(lat1, lon1, (strike0 - 180), l*ix+epsilon)
        for iy in range(ny):
            lla2, llo2, alpha21 = vinc_pt(lla1, llo1, (strike0 - 90), w*iy+epsilon)
            dist, alpha12,  alpha21 = vinc_dist(lla2, llo2, latE, lonE)
            Hi = Hmax - w*math.sin(math.radians(dip))*iy
            DistSubFallas[ix, iy] = math.sqrt(dist**2 + Hi**2)
            
    disref, alpha12,  alpha21 = vinc_dist(lat0, lon0, latE, lonE)
    Rref = math.sqrt(disref**2 + H**2)/1000
	#print(DistSubFallas, Rref)
    Rrup = np.amin(DistSubFallas)/1000
#    print(Rrup)
    if Rref < Rrup:
        print("Revisa la cochinada de codigo")
    return Rrup



# ----------------------------------------------------------------------
# | Algrothims from Geocentric Datum of Australia Technical Manual	    |
# | 								                                    |
# | http://www.anzlic.org.au/icsm/gdatum/chapter4.html	        		|
# | 								                                    |
# | This page last updated 11 May 1999 	                				|
# | 								                                    |
# | Computations on the Ellipsoid	                    				|
# | 								                                    |
# | There are a number of formulae that are available           		|
# | to calculate accurate geodetic positions, 		            		|
# | azimuths and distances on the ellipsoid.			                |
# | 								                                    |
# | Vincenty's formulae (Vincenty, 1975) may be used 		            |
# | for lines ranging from a few cm to nearly 20,000 km, 	            |
# | with millimetre accuracy. 					                        |
# | The formulae have been extensively tested 		                    |
# | for the Australian region, by comparison with results       		|
# | from other formulae (Rainsford, 1955 & Sodano, 1965). 	            |
# |								                                        |
# | * Inverse problem: azimuth and distance from known 	        		|
# |			latitudes and longitudes 			                        |
# | * Direct problem: Latitude and longitude from known 	            |
# |			position, azimuth and distance. 		                    |
# | * Sample data 						                                |
# | * Excel spreadsheet 			                            		|
# | 								                                    |
# | Vincenty's Inverse formulae				                    		|
# | Given: latitude and longitude of two points                 		|
# |			(phi1, lembda1 and phi2, lembda2), 	                        |
# | Calculate: the ellipsoidal distance (s) and 	            		|
# | forward and reverse azimuths between the points (alpha12, alpha21).	|
# |									                                    |
# ---------------------------------------------------------------------- 

def vinc_dist(  phi1,  lembda1,  phi2,  lembda2 ) :
        """ 
        Returns the distance between two geographic points on the ellipsoid
        and the forward and reverse azimuths between these points.
        lats, longs and azimuths are in decimal degrees, distance in metres 

        Returns ( s, alpha12,  alpha21 ) as a tuple
        """

        f = 1.0 / 298.257223563		# WGS84
        a = 6378137.0 			# metres
        if (abs( phi2 - phi1 ) < 1e-8) and ( abs( lembda2 - lembda1) < 1e-8 ) :
                return 0.0, 0.0, 0.0

        piD4   = math.atan( 1.0 )
        two_pi = piD4 * 8.0

        phi1    = phi1 * piD4 / 45.0
        lembda1 = lembda1 * piD4 / 45.0		# unfortunately lambda is a key word!
        phi2    = phi2 * piD4 / 45.0
        lembda2 = lembda2 * piD4 / 45.0

        b = a * (1.0 - f)

        TanU1 = (1-f) * math.tan( phi1 )
        TanU2 = (1-f) * math.tan( phi2 )

        U1 = math.atan(TanU1)
        U2 = math.atan(TanU2)

        lembda = lembda2 - lembda1
        last_lembda = -4000000.0		# an impossibe value
        omega = lembda

        # Iterate the following equations, 
        #  until there is no significant change in lembda 

        while ( last_lembda < -3000000.0 or lembda != 0 and abs( (last_lembda - lembda)/lembda) > 1.0e-9 ) :

                sqr_sin_sigma = pow( math.cos(U2) * math.sin(lembda), 2) + \
                        pow( (math.cos(U1) * math.sin(U2) - \
                        math.sin(U1) *  math.cos(U2) * math.cos(lembda) ), 2 )

                Sin_sigma = math.sqrt( sqr_sin_sigma )

                Cos_sigma = math.sin(U1) * math.sin(U2) + math.cos(U1) * math.cos(U2) * math.cos(lembda)
        
                sigma = math.atan2( Sin_sigma, Cos_sigma )

                Sin_alpha = math.cos(U1) * math.cos(U2) * math.sin(lembda) / math.sin(sigma)
                alpha = math.asin( Sin_alpha )

                Cos2sigma_m = math.cos(sigma) - (2 * math.sin(U1) * math.sin(U2) / pow(math.cos(alpha), 2) )

                C = (f/16) * pow(math.cos(alpha), 2) * (4 + f * (4 - 3 * pow(math.cos(alpha), 2)))

                last_lembda = lembda

                lembda = omega + (1-C) * f * math.sin(alpha) * (sigma + C * math.sin(sigma) * \
                        (Cos2sigma_m + C * math.cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2) )))

        u2 = pow(math.cos(alpha),2) * (a*a-b*b) / (b*b)

        A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))

        B = (u2/1024) * (256 + u2 * (-128+ u2 * (74 - 47 * u2)))

        delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B/4) * \
                (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2) ) - \
                (B/6) * Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) * \
                (-3 + 4 * pow(Cos2sigma_m,2 ) )))

        s = b * A * (sigma - delta_sigma)

        alpha12 = math.atan2( (math.cos(U2) * math.sin(lembda)), \
                (math.cos(U1) * math.sin(U2) - math.sin(U1) * math.cos(U2) * math.cos(lembda)))

        alpha21 = math.atan2( (math.cos(U1) * math.sin(lembda)), \
                (-math.sin(U1) * math.cos(U2) + math.cos(U1) * math.sin(U2) * math.cos(lembda)))

        if ( alpha12 < 0.0 ) : 
                alpha12 =  alpha12 + two_pi
        if ( alpha12 > two_pi ) : 
                alpha12 = alpha12 - two_pi

        alpha21 = alpha21 + two_pi / 2.0
        if ( alpha21 < 0.0 ) : 
                alpha21 = alpha21 + two_pi
        if ( alpha21 > two_pi ) : 
                alpha21 = alpha21 - two_pi

        alpha12    = alpha12    * 45.0 / piD4
        alpha21    = alpha21    * 45.0 / piD4
        return s, alpha12,  alpha21 

   # END of Vincenty's Inverse formulae 


#-------------------------------------------------------------------------------
# Vincenty's Direct formulae							|
# Given: latitude and longitude of a point (phi1, lembda1) and 			|
# the geodetic azimuth (alpha12) 						|
# and ellipsoidal distance in metres (s) to a second point,			|
# 										|
# Calculate: the latitude and longitude of the second point (phi2, lembda2) 	|
# and the reverse azimuth (alpha21).						|
# 										|
#-------------------------------------------------------------------------------

def  vinc_pt( phi1, lembda1, alpha12, s ) :
        """

        Returns the lat and long of projected point and reverse azimuth
        given a reference point and a distance and azimuth to project.
        lats, longs and azimuths are passed in decimal degrees

        Returns ( phi2,  lambda2,  alpha21 ) as a tuple 

        """
 
        f = 1.0 / 298.257223563		# WGS84
        a = 6378137.0 			# metres
        piD4 = math.atan( 1.0 )
        two_pi = piD4 * 8.0

        phi1    = phi1    * piD4 / 45.0
        lembda1 = lembda1 * piD4 / 45.0
        alpha12 = alpha12 * piD4 / 45.0
        if ( alpha12 < 0.0 ) : 
                alpha12 = alpha12 + two_pi
        if ( alpha12 > two_pi ) : 
                alpha12 = alpha12 - two_pi

        b = a * (1.0 - f)

        TanU1 = (1-f) * math.tan(phi1)
        U1 = math.atan( TanU1 )
        sigma1 = math.atan2( TanU1, math.cos(alpha12) )
        Sinalpha = math.cos(U1) * math.sin(alpha12)
        cosalpha_sq = 1.0 - Sinalpha * Sinalpha

        u2 = cosalpha_sq * (a * a - b * b ) / (b * b)
        A = 1.0 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * \
                (320 - 175 * u2) ) )
        B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2) ) )

        # Starting with the approximation
        sigma = (s / (b * A))

        last_sigma = 2.0 * sigma + 2.0	# something impossible

        # Iterate the following three equations 
        #  until there is no significant change in sigma 

        # two_sigma_m , delta_sigma
        while ( abs( (last_sigma - sigma) / sigma) > 1.0e-9 ) :
                two_sigma_m = 2 * sigma1 + sigma

                delta_sigma = B * math.sin(sigma) * ( math.cos(two_sigma_m) \
                        + (B/4) * (math.cos(sigma) * \
                        (-1 + 2 * pow( math.cos(two_sigma_m), 2 ) -  \
                        (B/6) * math.cos(two_sigma_m) * \
                        (-3 + 4 * pow(math.sin(sigma), 2 )) *  \
                        (-3 + 4 * pow( math.cos (two_sigma_m), 2 ))))) \

                last_sigma = sigma
                sigma = (s / (b * A)) + delta_sigma

        phi2 = math.atan2 ( (math.sin(U1) * math.cos(sigma) + math.cos(U1) * math.sin(sigma) * math.cos(alpha12) ), \
                ((1-f) * math.sqrt( pow(Sinalpha, 2) +  \
                pow(math.sin(U1) * math.sin(sigma) - math.cos(U1) * math.cos(sigma) * math.cos(alpha12), 2))))

        lembda = math.atan2( (math.sin(sigma) * math.sin(alpha12 )), (math.cos(U1) * math.cos(sigma) -  \
                math.sin(U1) *  math.sin(sigma) * math.cos(alpha12)))

        C = (f/16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq ))

        omega = lembda - (1-C) * f * Sinalpha *  \
                (sigma + C * math.sin(sigma) * (math.cos(two_sigma_m) + \
                C * math.cos(sigma) * (-1 + 2 * pow(math.cos(two_sigma_m),2) )))

        lembda2 = lembda1 + omega

        alpha21 = math.atan2 ( Sinalpha, (-math.sin(U1) * math.sin(sigma) +  \
                math.cos(U1) * math.cos(sigma) * math.cos(alpha12)))

        alpha21 = alpha21 + two_pi / 2.0
        if ( alpha21 < 0.0 ) :
                alpha21 = alpha21 + two_pi
        if ( alpha21 > two_pi ) :
                alpha21 = alpha21 - two_pi

        phi2       = phi2       * 45.0 / piD4
        lembda2    = lembda2    * 45.0 / piD4
        alpha21    = alpha21    * 45.0 / piD4

        return phi2,  lembda2,  alpha21 

  # END of Vincenty's Direct formulae

#--------------------------------------------------------------------------
# Notes: 
# 
# * "The inverse formulae may give no solution over a line 
# 	between two nearly antipodal points. This will occur when 
# 	lembda ... is greater than pi in absolute value". (Vincenty, 1975)
#  
# * In Vincenty (1975) L is used for the difference in longitude, 
# 	however for consistency with other formulae in this Manual, 
# 	omega is used here. 
# 
# * Variables specific to Vincenty's formulae are shown below, 
# 	others common throughout the manual are shown in the Glossary. 
# 
# 
# alpha = Azimuth of the geodesic at the equator
# U = Reduced latitude
# lembda = Difference in longitude on an auxiliary sphere (lembda1 & lembda2 
# 		are the geodetic longitudes of points 1 & 2)
# sigma = Angular distance on a sphere, from point 1 to point 2
# sigma1 = Angular distance on a sphere, from the equator to point 1
# sigma2 = Angular distance on a sphere, from the equator to point 2
# sigma_m = Angular distance on a sphere, from the equator to the 
# 		midpoint of the line from point 1 to point 2
# u, A, B, C = Internal variables
# 
# 

#*******************************************************************

def est_dist( phi1,  lembda1,  phi2,  lembda2 ) :
        """ 

        Returns an estimate of the distance between two geographic points
        This is a quick and dirty vinc_dist 
        which will generally estimate the distance to within 1%
        Returns distance in metres

        """

        f = 1.0 / 298.257223563		# WGS84
        a = 6378137.0 			# metres
        piD4   = 0.785398163397 

        phi1    = phi1 * piD4 / 45.0
        lembda1 = lembda1 * piD4 / 45.0
        phi2    = phi2 * piD4 / 45.0
        lembda2 = lembda2 * piD4 / 45.0

        c = math.cos((phi2+phi1)/2.0)

        return math.sqrt( pow(math.fabs(phi2-phi1), 2) + \
                pow(math.fabs(lembda2-lembda1)*c, 2) ) * a * ( 1.0 - f + f * c )

   # END of rough estimate of the distance.

