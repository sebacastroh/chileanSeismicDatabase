"""import urllib
from stations_csn import stations_csn

names = {}

for station in stations_csn:
    fid = urllib.urlopen('http://evtdb.csn.uchile.cl/station/' + station)
    webpage = fid.read().decode('utf-8')
    name = ''
    for line in webpage.splitlines():
        if station + ' (' in line:
            i = line.find('(') + 1
            j = line.find(')')
            name = line[i:j]
            break
    stations_csn[station] = [name] + stations_csn[station]

"""
for station in sorted(stations_csn):
    if len(stations_csn[station]) == 5:
        print("'" + station + "': " + str( stations_csn[station]) + ',')
    else:
        print("'" + station + "': " + str( [stations_csn[station][0]]  + stations_csn[station][1:]) + ',')

