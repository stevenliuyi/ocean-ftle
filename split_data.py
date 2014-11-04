def split(filename, s):
    f = open(filename)
    lon = -89.875
    lat = 25.125
    day = 1
    fw = open(s + '_' + str(day).zfill(3) + '.ascii', 'w')
    for line in f:
        if (lon == 9.875):
            if (lat == 74.875):
                lon = -89.875
                lat = 25.125
                fw.write(line)
                fw.close()
                day += 1
                fw = open(s + '_' + str(day).zfill(3) + '.ascii', 'w')
                continue
            lat += .25
            lon = -89.875
        else:
            lon += .25
        fw.write(line)
    f.close()

ufilename = "u_NATL_2005.ascii"
vfilename = "v_NATL_2005.ascii"

split(ufilename, 'u')
split(vfilename, 'v')
