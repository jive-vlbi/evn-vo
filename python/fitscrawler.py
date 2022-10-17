import csv
import os
import urllib
import sys
import time as tm

from astropy.io import fits
from astropy import constants as const
from astropy import units as u
from astropy.time import Time
import numpy as np

import vlbi

diameters = {
    # EVN and associates
    'AR': 305,
    'BD': 32,
    'CM': 32,
    'DA': 25,
    'DE': 25,
    'DW': 25,
    'EB': 100,
    'EF': 100,
    'EV': 70,		# Evpatoria
    'HH': 26,
    'IB': 16,
    'IR': 32,
    'JB': 76,
    'JL': 76,
    'JV': 30.82,	# Geometric mean for Joddrell Mark2
    'KM': 40,
    'KN': 25,
    'KS': 34,		# Kashima
    'KT': 22,
    'KU': 22,
    'KY': 22,
    'MA': 20,		# Matera
    'MC': 32,
    'MH': 14,
    'NK': 32,
    'NT': 32,
    'O6': 20,
    'O8': 25,
    'ON': 25,
    'PI': 25,
    'RO': 70,
    'SH': 25,
    'SR': 64,
    'SV': 32,
    'T6': 65,
    'TA': 25,
    'TR': 32,
    'UR': 25,
    'VM': 20,		# VERA Mizusawa
    'VS': 20,		# VERA Ishigakijima
    'WB': 25,
    'WN': 13,
    'WZ': 20,
    'YM': 32,		# Yamaguchi
    'YS': 40,
    'ZC': 32,

    # VLBA
    'BR': 25,
    'FD': 25,
    'GB': 104.88,	# Geometric mean for Green Bank
    'HN': 25,
    'KP': 25,
    'LA': 25,
    'MK': 25,
    'NL': 25,
    'OV': 25,
    'PT': 25,
    'SC': 25,
    'YY': 0,

    # LBA
    'AT': 0,
    'CD': 30,
    'HO': 26,
    'KE': 12,
    'MP': 22,
    'PA': 64,
    'TD': 34,
    'TI': 70,
    'WA': 30,
    'YG': 12,

    # Test stations
    'ED': 0,
    'EX': 0,
    'HD': 0,
    'HX': 0,
    'JD': 0,
    'KD': 0,
    'MD': 0,
    'ND': 0,
    'OD': 0,
    'OX': 0,
    'SD': 0,
    'WD': 0,
    'YD': 0,
    }

class IdiHDU(fits.PrimaryHDU):
    @classmethod
    def match_header(cls, header):
        try:
            keyword = header.cards[0].keyword
        except:
            keyword = header.ascard[0].key
            pass
        return (keyword == 'SIMPLE' and 'GROUPS' in header and
                header['GROUPS'] == True and 'NAXIS' in header and
                header['NAXIS'] == 0)

fits.register_hdu(IdiHDU)

os.environ['TZ'] = 'UTC'
tm.tzset()

def parse_fitsidi(exp_id, product_id, idifiles, csvfile):
    nparts = 0
    for file in idifiles:
        part = int(os.path.splitext(file)[1][4:])
        if not part == nparts + 1:
            raise RuntimeError("non-sequenctial FITS-IDI files")
        nparts += 1
        continue

    hdulist = fits.open(idifiles[0])
    #hdulist.info()

    hdu = hdulist['FREQUENCY']
    #print(hdu.columns)

    f_min = 1e100
    f_max = 0
    f_resolution_min = 1e100
    f_resolution_max = 0

    ref_freq = hdu.header['REF_FREQ']
    for row in hdu.data:
        if hdu.header['NO_BAND'] == 1:
            zipped = zip([row['BANDFREQ']], [row['SIDEBAND']],
                         [row['TOTAL_BANDWIDTH']])
        else:
            zipped = zip(row['BANDFREQ'], row['SIDEBAND'],
                         row['TOTAL_BANDWIDTH'])
            pass
        for band in zipped:
            if band[1]:
                f_min = min(f_min, ref_freq + band[0])
                f_max = max(f_max, ref_freq + band[0] + band[2])
            else:
                f_min = min(f_min, ref_freq + band[0] - band[2])
                f_max = max(f_max, ref_freq + band[0])
                pass
            continue
        f_resolution_min = min(f_resolution_min, row['CH_WIDTH'].min())
        f_resolution_max = max(f_resolution_max, row['CH_WIDTH'].max())
        continue

    f_min = f_min * u.Hz
    f_max = f_max * u.Hz
    f_resolution_min = f_resolution_min * u.Hz
    f_resolution_max = f_resolution_max * u.Hz
    f_resolution = (f_resolution_min + f_resolution_max) / 2

    em_min = f_max.to(u.meter, equivalencies=u.spectral())
    em_max = f_min.to(u.meter, equivalencies=u.spectral())
    em_res_power = ref_freq * u.Hz / f_resolution
    em_res_power_min = f_min / f_resolution_max
    em_res_power_max = f_max / f_resolution_min
    em_resolution = f_resolution.to(u.meter, equivalencies=u.spectral())

    B = 0
    D = 25
    stabxyz = []
    hdu = hdulist['ARRAY_GEOMETRY']
    for row in hdu.data:
        anname = row['ANNAME'].upper()
        if anname in diameters:
            D = max(D, diameters[anname])
        else:
            print('%s: Unknown diameter' % anname)
            pass
        stabxyz.append(row['STABXYZ'])
        continue
    for i in range(len(stabxyz)):
        for j in range(i):
            baseline = stabxyz[i] - stabxyz[j]
            B = max(B, np.linalg.norm(baseline))
            continue
        continue

    hdu = hdulist['SOURCE']
    #print(hdu.columns)

    max_source_id = 0
    for row in hdu.data:
        if row['SOURCE_ID'] > max_source_id:
            max_source_id = row['SOURCE_ID']
            pass
        continue

    target_name = (max_source_id + 1) * [""] 
    s_ra = (max_source_id + 1) * [0.0]
    s_dec = (max_source_id + 1) * [0.0]
    s_fov = (max_source_id + 1) * [0.0]
    s_region = (max_source_id + 1) * [""]
    s_resolution = (max_source_id + 1) * [0.0]
    for row in hdu.data:
        source_id = row['SOURCE_ID']
        target_name[source_id] = row['SOURCE']
        # XXX Are these indeed ICRS coordinates?
        s_ra[source_id] = row['RAEPO'] * u.deg
        s_dec[source_id] = row['DECEPO'] * u.deg
        continue
    s_xel1 = -1
    s_xel2 = -1

    hdu = hdulist['UV_DATA']
    #print(hdu.columns)

    obs_id = hdu.header['OBSCODE']
    em_xel = hdu.header['NO_BAND'] * hdu.header['NO_CHAN']
    pol_xel = hdu.header['NO_STKD']

    prev_source_id = -1
    prev_jd = 0.0
    jd_max = (max_source_id + 1) * [0.0]
    jd_min = (max_source_id + 1) * [1e100]
    nvis = (max_source_id + 1) * [0]
    inttim = (max_source_id + 1) * [0.0]
    inttim_min = (max_source_id + 1) * [1e100]
    inttim_max = (max_source_id + 1) * [0.0]
    access_estsize = 0
    for file in idifiles:
        hdulist = fits.open(file)
        hdu = hdulist['UV_DATA']
        for row in hdu.data:
            # Skip bogus rows
            if row['DATE'] == 0.0 and row['TIME'] == 0.0:
                continue
            if row['SOURCE_ID'] < 0 or row['SOURCE_ID'] > max_source_id:
                continue
            jd = row['DATE'] + row['TIME']
            jd_min[source_id] = min(jd_min[source_id], jd)
            jd_max[source_id] = max(jd_max[source_id], jd)
            source_id = row['SOURCE_ID']
            if source_id != prev_source_id or jd != prev_jd:
                nvis[source_id] += 1
                inttim[source_id] += row['INTTIM']
                prev_source_id = source_id
                prev_jd = jd
                pass
            inttim_min[source_id] = min(inttim_min[source_id], row['INTTIM'])
            inttim_max[source_id] = max(inttim_max[source_id], row['INTTIM'])
            continue
        access_estsize += (os.path.getsize(file) + 999) // 1000
        continue

    t_xel = nvis
    t_min = Time(jd_min, format='jd')
    t_min.format = "mjd"
    t_max = Time(jd_max, format='jd')
    t_max.format = "mjd"
    t_exptime = inttim * u.s
    t_resolution = (np.add(inttim_min, inttim_max) / 2) * u.s

    dataproduct_type = "visibility"
    calib_level = 1

    obs_collection = "EVN"
    obs_publisher_did = "ivo://jive.eu/~?" + exp_id + product_id

    access_url = "http://archive.jive.nl/exp/" + exp_id + "/fits"
    access_format = "application/x-fits-idi"

    o_ucd = "stat.uncalib"

    instrument_name = "EVN"

    record = {}
    record['dataproduct_type'] = dataproduct_type
    record['calib_level'] = calib_level
    record['target_name'] = None
    record['obs_id'] = obs_id
    record['obs_collection'] = obs_collection
    record['obs_publisher_did'] = None
    record['access_url'] = access_url
    record['access_format'] = access_format
    record['access_estsize'] = access_estsize
    record['s_ra'] = None
    record['s_dec'] = None
    record['s_fov'] = None
    record['s_region'] = None
    record['s_resolution'] = None
    record['s_xel1'] = s_xel1
    record['s_xel2'] = s_xel2
    record['t_xel'] = None
    record['t_min'] = None
    record['t_max'] = None
    record['t_exptime'] = None
    record['t_resolution'] = None
    record['em_xel'] = em_xel
    record['em_min'] = em_min.to_value(u.m)
    record['em_max'] = em_max.to_value(u.m)
    record['em_res_power'] = em_res_power
    record['o_ucd'] = o_ucd
    record['pol_xel'] = pol_xel
    record['instrument_name'] = instrument_name
    record['_nparts'] = nparts
    record['_product_id'] = exp_id + product_id

    fieldnames = record.keys()
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    if csvfile.tell() == 0:
        writer.writeheader()
        pass

    for source_id in range(len(t_xel)):
        if t_xel[source_id] == 0:
            continue
        nu = ((f_min + f_max) / 2).to_value(u.Hz)
        delta_nu = ((f_resolution_min + f_resolution_max) / 2).to_value(u.Hz)
        tau = t_resolution[source_id].to_value(u.s)
        s_fov[source_id] = vlbi.fov(nu, delta_nu, D, B, tau) * u.rad
        ra = s_ra[source_id].to_value(u.deg)
        dec = s_dec[source_id].to_value(u.deg)
        fov = s_fov[source_id].to_value(u.deg)
        s_region[source_id] = "Circle J2000 %f %f %f" % (ra, dec, fov)
        s_resolution[source_id] = vlbi.resolution(nu, B) * u.rad
        continue

    for source_id in range(len(t_xel)):
        if t_xel[source_id] == 0:
            continue
        record['target_name'] = target_name[source_id]
        record['s_ra'] = s_ra[source_id].to_value(u.deg)
        record['s_dec'] = s_dec[source_id].to_value(u.deg)
        record['s_fov'] = s_fov[source_id].to_value(u.deg)
        record['s_region'] = s_region[source_id]
        record['s_resolution'] = s_resolution[source_id].to_value(u.arcsec)
        record['t_xel'] = t_xel[source_id]
        record['t_min'] = t_min[source_id].value
        record['t_max'] = t_max[source_id].value
        record['t_exptime'] = t_exptime[source_id].to_value(u.s)
        record['t_resolution'] = t_resolution[source_id].to_value(u.s)
        target = urllib.parse.quote(target_name[source_id])
        freq = "%.2fMHz" % f_min.to_value(u.MHz)
        record['obs_publisher_did'] = obs_publisher_did + '_' + target + '_' + freq
        writer.writerow(record)
        continue
    return

idifiles = sys.argv[1:]
idifiles = sorted(idifiles, key=lambda s: int(s[s.rfind('IDI') + 3:]))
exp_id = os.path.split(os.path.split(os.path.split(idifiles[0])[0])[0])[1]
product_id = os.path.splitext(os.path.basename(idifiles[0]))[0]
product_id = product_id[product_id.find('_'):]
csvfile = open('records.csv', 'a')
parse_fitsidi(exp_id, product_id, idifiles, csvfile)
