#import sys
#import os
import numpy as np
#import matplotlib.pyplot as plt
import argparse
#from tqdm import tqdm
import time
import pysiaf
from pysiaf.utils.rotations import attitude
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits

#######################################
# Create command line argument parser
#######################################

def create_parser():

	# Handle user input with argparse
    parser = argparse.ArgumentParser(
        description="Flags and options from user.")

    parser.add_argument('-i', '--input',
        dest='input',
        default='input.pointing',
        metavar='input',
        type=str,
        help='APT pointing file to process.')

    parser.add_argument('--nirspec',
        dest='nirspec',
        default='nirspec.txt',
        metavar='nirspec',
        type=str,
        help='Output file for nirspec pointings.')

    parser.add_argument('--nircam',
        dest='nircam',
        default='nircam.txt',
        metavar='nircam',
        type=str,
        help='Output file for nircam pointings.')

    parser.add_argument('--nc_color',
        dest='nc_color',
        default='#C0C0C0',
        metavar='nc_color',
        type=str,
        help='Default color for nircam overlay.')

    parser.add_argument('--nc_fill',
        dest='nc_fill',
        default=0.01,
        metavar='nc_fill',
        type=float,
        help='Default fill for nircam overlay.')

    parser.add_argument('--ns_color',
        dest='ns_color',
        default='#FFD700',
        metavar='ns_color',
        type=str,
        help='Default color for nirspec overlay.')

    parser.add_argument('--ns_fill',
        dest='ns_fill',
        default=0.01,
        metavar='ns_fill',
        type=float,
        help='Default fill for nirspec overlay.')

    parser.add_argument('--prime',
        dest='prime',
        default='NS',
        metavar='prime',
        type=str,
        help='Prime instrument defining the aperture PA (default NS, opts NC, MIRI).')

    parser.add_argument('--fits',
        dest='fits',
        default=None,
        metavar='fits',
        type=str,
        help='Input fits file with WCS to convert to pixel positions.')

    parser.add_argument('-v', '--verbose',
        dest='verbose',
        action='store_true',
        help='Print helpful information to the screen? (default: False)',
        default=False)

    return parser


#######################################
# get_aperture_boxes() function
#######################################
def get_aperture_boxes(lin,instrument='NC',prime='NS',wcs=None):

    #warning -- target name is split always in two

    l = lin.split()
    print(l)
    boxes = []
    pixel_aper = []

    i_exp = 18
    print(l[i_exp])

    flag_exp = False
    if(l[i_exp]=='SCIENCE'):
        pointing = {'Tar':l[0], 'Tile':l[1], 'Exp':l[2], 
                'Dith':l[3], 'Aperture':l[4], 'Target':l[6],
                'RA':float(l[7]), 'DEC':float(l[8]),
                'BaseX':float(l[9]), 'BaseY':float(l[10]),
                'DithX':float(l[11]), 'DithY':float(l[12]),
                'V2':float(l[13]),'V3':float(l[14]),
                'IdlX':float(l[15]),'IdlY':float(l[16]),
                'dDist':float(l[20]),'PA':float(l[22])}
        flag_exp = True
    elif(l[i_exp]=='PARALLEL'):
        pointing = {'Tar':l[0], 'Tile':l[1], 'Exp':l[2], 
                'Dith':l[3], 'Aperture':l[4], 'Target':l[6],
                'RA':float(l[7]), 'DEC':float(l[8]),
                'BaseX':float(l[9]), 'BaseY':float(l[10]),
                'DithX':float(l[11]), 'DithY':float(l[12]),
                'V2':float(l[13]),'V3':float(l[14]),
                'IdlX':float(l[15]),'IdlY':float(l[16]),
                'dDist':float(0),'PA':float(l[21])}
        flag_exp = True

    #print(pointing)

    if(prime=='NS'):
        siaf_prime = pysiaf.Siaf('NIRSpec')
        aper_prime = siaf_prime['NRS_FULL_MSA']
    elif(prime=='NC'):
        siaf_prime = pysiaf.Siaf('NIRCam') 
        aper_prime = siaf_prime['NRCALL_FULL']

    if(flag_exp):


        if(instrument=='NC'):
            siaf = pysiaf.Siaf('NIRCam')
            nc_list = ['NRCA5_FULL','NRCB5_FULL']
            apers = [siaf[nc] for nc in nc_list] # NRCA5_FULL is the SIAF aperture name for NIRCam Module A LW channel
        elif(instrument=='NS'):
            siaf = pysiaf.Siaf('NIRSpec')
            msa_list = ["NRS_FULL_MSA1","NRS_FULL_MSA2",
                        "NRS_FULL_MSA3", "NRS_FULL_MSA4"]
            apers = [siaf[msa] for msa in msa_list] # NRCA5_FULL is the SIAF aperture name for NIRCam Module A LW channel
           # msa_dict = {"NRS_FULL_MSA1": "NRS2_FULL", "NRS_FULL_MSA2": "NRS2_FULL",
           #             "NRS_FULL_MSA3": "NRS1_FULL", "NRS_FULL_MSA4": "NRS1_FULL",
           #             "NRS_S1600A1_SLIT": "NRS1_FULL"}
            #apers = [siaf[nc] for nc in nc_list] # NRCA5_FULL is the SIAF aperture name for NIRCam Module A LW channel

        #v2_1, v3_1 = aper.idl_to_tel(ap.XIdlVert1, ap.YIdlVert1)
        #v2_2, v3_2 = aper.idl_to_tel(ap.XIdlVert2, ap.YIdlVert2)
        #v2_3, v3_3 = aper.idl_to_tel(ap.XIdlVert3, ap.YIdlVert3)
        #v2_4, v3_4 = aper.idl_to_tel(ap.XIdlVert4, ap.YIdlVert4)

        ra = pointing['RA']
        dec = pointing['DEC']
        pa = pointing['PA']
        V2 = pointing['V2']
        V3 = pointing['V3']
        IdlX = pointing['IdlX']
        IdlY = pointing['IdlY']

        print(pointing)
#        print(f'IdlX {IdlX}')
#        print(f'IdlY {IdlY}')
#        print(f'PA = {pa}')


        for aper in apers:
            pixels = []

            v23 = [aper.idl_to_tel(aper.XIdlVert1, aper.YIdlVert1), 
                      aper.idl_to_tel(aper.XIdlVert2, aper.YIdlVert2),
                      aper.idl_to_tel(aper.XIdlVert3, aper.YIdlVert3),
                      aper.idl_to_tel(aper.XIdlVert4, aper.YIdlVert4)]

#            v23 = [[aper.XIdlVert1, aper.YIdlVert1],[aper.XIdlVert2, aper.YIdlVert2],
#                      [aper.XIdlVert3, aper.YIdlVert3],
#                      [aper.XIdlVert4, aper.YIdlVert4]]
#            v23 = [aper.idl_to_tel(aper.XIdlVert1+IdlX, aper.YIdlVert1+IdlY), 
#                      aper.idl_to_tel(aper.XIdlVert2+IdlX, aper.YIdlVert2+IdlY),
#                      aper.idl_to_tel(aper.XIdlVert3+IdlX, aper.YIdlVert3+IdlY),
#                      aper.idl_to_tel(aper.XIdlVert4+IdlX, aper.YIdlVert4+IdlY)]
#53097, -27.8883
        #    matt = [attitude(V2,V3,ra,dec,pa) for v23i in v23]

# target is in lower right MSA quadarant


            matt = attitude(V2,V3,ra,dec,pa-aper_prime.V3IdlYAngle)
#            matt = attitude(V2,V3,ra,dec,0) # V3 PA range 0.5254 when PA = 139.1

            aper.set_attitude_matrix(matt)

            for i in range(len(v23)):
                #aper.set_attitude_matrix(matt)
                v23i = v23[i]
                print(f'{i} Vertex in v23 {v23[i]}')
                print(f'V2 {V2} V3 {V3} V2Ref {aper.V2Ref} V3Ref {aper.V3Ref} V3IdlYAngle {aper.V3IdlYAngle}')
#                matt = attitude(v23i[0],v23i[1],ra,dec,pa)
#                aper.set_attitude_matrix(matt)
                sc = SkyCoord( *aper.tel_to_sky(v23i[0], v23i[1]), unit='deg' )
#                sc = SkyCoord( *aper.idl_to_sky(v23i[0], v23i[1]), unit='deg' )

                #sc.to_string('decimal')
                rao, deco = sc.ra.value, sc.dec.value
                boxes.append([rao,deco])
                if(wcs is not None):
                    psc = sc.to_pixel(wcs)
                    pixels.append([float(psc[1]),float(psc[0])])
                #print(v23i[0],v23i[1],rao,deco)#.data.value)

            #for p in pixels:
            #print(pixels)
            pixel_aper.append(pixels)


    if(flag_exp==False):
        pointing = None
    return boxes,pixel_aper,pointing

#######################################
# main() function
#######################################
def main():

    #begin timer
    time_global_start = time.time()

    #create the command line argument parser
    parser = create_parser()

    #store the command line arguments
    args   = parser.parse_args()

    # open files
    fp = open(args.input,'r')

    n_ns = 0
    n_nc = 0
    n_v  = 0
    n_o  = 0

    fl = fp.readlines()
    fp.close()

    # initial scan
    for l in fl:

        try:
            if(l.split()[0]=='*'):
                n_o += 1
        except:
            pass

        try:
            if(l.split()[0]=='**'):
                n_v += 1
        except:
            pass

        try:
            if(l.split()[4]=='NRS_FULL_MSA'):
                print(l.split())
                n_ns += 1
        except:
            pass

        try:
            if(l.split()[4]=='NRCALL_FULL'):
                print(l.split())
                n_nc += 1
        except:
            pass

    print(f'n_o  {n_o}')
    print(f'n_v  {n_v}')
    print(f'n_ns {n_ns}')
    print(f'n_nc {n_nc}')

    wcs = None
    if(args.fits is not None):
        header = fits.getheader(args.fits,'SCI')
        wcs = WCS(header)


    obs = ''
    visit = ''

    nc_pixels = []
    nc_poly = []

    ns_pixels = []
    ns_poly = []

    nc_names = []
    ns_names = []

    # reparse
    k = 0
    for l in fl:

        print(k,len(fl))

        flag = False
        flag_ns = False
        flag_nc = False

        #get observation
        try:
            if(l.split()[0]=='*'):
                obs = l.split()[2]
                #flag = True
        except:
            pass

        #get visit
        try:
            if(l.split()[0]=='**'):
                visit = l.split()[2].replace(':','_')
                flag = True
        except:
            pass

        #get NS
        try:
            if(l.split()[4]=='NRS_FULL_MSA'):
                flag = True
                flag_ns = True
        except:
            pass

        #get NC
        try:
            if(l.split()[4]=='NRCALL_FULL'):
                flag = True
                flag_nc = True
        except:
            pass


        if(flag):
            try:
                fpns.close()
            except:
                pass
            try:
                fpnc.close()
            except:
                pass

            #print(obs,visit)



            # OK we have some visits coming
            # https://jwst-docs.stsci.edu/jwst-observatory-hardware/jwst-telescope/jwst-focal-plane
            if(flag_ns):
                print(l.split())
                #ns_poly = get_nirspec_boxes(l)
                ns_poly_i, ns_pixels_i, ns_pointing_i = get_aperture_boxes(l,instrument='NS',prime=args.prime,wcs=wcs)

                if(ns_pointing_i is not None):
                    ns_pixels.append(ns_pixels_i)
                    ns_poly.append(ns_poly_i)
                    fpns = open(f'nirspec_{obs}_{visit}.txt','a')
                    ns_pixels.append(ns_pixels_i)
                    ns_poly.append(ns_poly_i)
                    Exp  = ns_pointing_i["Exp"]
                    Dith = ns_pointing_i['Dith']
                    Targ = ns_pointing_i['Target'].strip('.')
                    Targ = Targ.replace('-','')
                    name = f'{Targ}_ns_{obs}_{visit}_{Exp}_{Dith}_regions'
                    s = f'{name}={ns_pixels_i}\n'
                    fpns.write(s)
                    poly_name = f'{name}_poly'
                    ns_names.append(poly_name)
                    s = f'const {poly_name}=L.polygon({name},{{color:"{args.ns_color}","fillOpacity":{args.ns_fill}}},).bindPopup("{Targ} MSA Pointing {visit}.{Exp}.{Dith}");\n'
                    fpns.write(s)
                    fpns.close()

            if(flag_nc):
                nc_poly_i, nc_pixels_i, nc_pointing_i = get_aperture_boxes(l,instrument='NC',prime=args.prime,wcs=wcs)
                if(nc_pointing_i is not None):
                    fpnc = open(f'nircam_{obs}_{visit}.txt','a')
                    nc_pixels.append(nc_pixels_i)
                    nc_poly.append(nc_poly_i)
                    Exp = nc_pointing_i["Exp"]
                    Dith = nc_pointing_i['Dith']
                    Targ = nc_pointing_i['Target'].strip('.')
                    Targ = Targ.replace('-','')
                    name = f'{Targ}_nc_{obs}_{visit}_{Exp}_{Dith}_regions'
                    s = f'{name}={nc_pixels_i}\n'
                    fpnc.write(s)
                    poly_name = f'{name}_poly'
                    nc_names.append(poly_name)
                    s = f'const {poly_name}=L.polygon({name},{{color:"{args.nc_color}","fillOpacity":{args.nc_fill}}},).bindPopup("{Targ} NC Pointing {visit}.{Exp}.{Dith}");\n'
                    fpnc.write(s)
                    fpnc.close()




        k+=1


    # if we have nirspec regions, create the
    # final layer group for fitsmap
    if(len(ns_names)>0):
        fpns = open('ns_layer_group.txt','w')
        s = f'const {Targ}_ns_regions=L.layerGroup({ns_names})\n'
        s = s.replace("'",'')
        fpns.write(s)
        s = f'"{Targ} NS MSA Regions":{Targ}_ns_regions,'
        fpns.write(s)
        fpns.close()


    # if we have nircam regions, create the 
    # final layer group for fitsmap
    if(len(nc_names)>0):
        fpnc = open('nc_layer_group.txt','w')
        s = f'const {Targ}_nc_regions=L.layerGroup({nc_names})\n'
        s = s.replace("'",'')
        fpnc.write(s)
        s = f'"{Targ} NC Regions":{Targ}_nc_regions,'
        fpnc.write(s)
        fpnc.close()




# example JS for fitsmap overlay
#const deep_hst_gs_triple_1_pointing_7_ds9_slit = L.polygon(deep_hst_gs_triple_1_pointing_7_ds9_slit_regions_layer,{color:"#CACFD2","fillOpacity":0},).bindPopup("Deep HST Pointing 7");
#deep_hst_gs_triple_1_pointing_13_ds9_slit_regions_layer=[[[],[]]]
#const deep_hst_gs_triple_1_pointing_7_ds9_slit = L.polygon(deep_hst_gs_triple_1_pointing_7_ds9_slit_regions_layer,{color:"#CACFD2","fillOpacity":0},).bindPopup("Deep HST Pointing 7");
#const deep_hst_gs_slits = L.layerGroup([deep_hst_gs_triple_1_pointing_7_ds9_slit,deep_hst_gs_triple_1_pointing_8_ds9_slit,deep_hst_gs_triple_1_pointing_13_ds9_slit])
#     "JADES NIRSpec Deep HST":deep_hst_gs_slits,

    #end timer
    time_global_end = time.time()
    if(args.verbose):
    	print(f"Time to execute program: {time_global_end-time_global_start}s.")

#######################################
# Run the program
#######################################
if __name__=="__main__":
	main()