import glob, numpy, math, datetime, time, pickle, os, h5py

debug = False

def new_gdf(dirn,dtype="<i2",itemsize=4):
    result = {}
    if os.path.isfile("%s/sampler_config.hdf5"%(dirn)):
        if debug:
            print("Version >1.0")
        f = h5py.File("%s/sampler_config.hdf5"%(dirn),"r")

        for k in f.keys():
            result[k] = f[k].value

        if debug:
            print(result["version"])
            print(result["dtype"])
            print(result["itemsize"])
            print(result["n_channels"])
        result["data_dirs"] = []
        for i in range(result["n_channels"]):
            result["data_dirs"].append("%s/%03d"%(dirn,i))
        f.close()
    else:
        if debug:
            print("Version 1.0")
        result["version"] = "1.0"
        result["itemsize"] = itemsize
        result["dtype"] = dtype
        result["data_dirs"] = [dirn]
        result["n_channels"] = 1

    files = []
    for i in range(len(result["data_dirs"])):
        ddir = "%s/*.gdf"%(result["data_dirs"][i])
        if debug:
            print(ddir)
        fl = glob.glob(ddir)
        fl.sort()
        files.append(fl)
    
    result['file_size'] = os.path.getsize(files[0][0])/result["itemsize"]
    result['file_list'] = files
    result['max_n'] = result['file_size']*len(files[0])
    result['cache'] = []
    result['cache_idx'] = []
    for i in range(result["n_channels"]):
        result["cache"].append(numpy.zeros(result['file_size'],dtype=numpy.complex64))
        result['cache_idx'].append(-1)

    result['scale'] = 1.0
    result['dirn'] = dirn
    result['re_idx'] = numpy.arange(result['file_size'],dtype=numpy.int64)*2
    result['im_idx'] = numpy.arange(result['file_size'],dtype=numpy.int64)*2+1

    if result["version"] == "1.0":
        tstamp_file = open("%s/timestamps.log"%(dirn),'r')
        result['t0'] = float(tstamp_file.readline().split(" ")[1])
        tstamp_file.close()

    if debug:
        print("Read gdf dir. Version %s Nchannels %d Nfiles=%d Fsize=%d t0 %1.2f"%(result["version"],len(result['file_list']),len(result["file_list"][0]),result['file_size'],result['t0']))
    return(result)

def read_vec(gdf, idx, length, chan=0):
    idx = int(idx)
    length = int(length)
    files = gdf['file_list'][chan]
    f0_idx = int(math.floor(idx/gdf['file_size']))
    fn_idx = int(math.floor((idx+length-1)/gdf['file_size']))+1
    res_vec = numpy.zeros([length],dtype=numpy.complex64)
    if debug:
        print("f0 %d %d"%(f0_idx,fn_idx))
    c0_idx = idx % gdf['file_size']
    c1_idx = (idx+length-1) % gdf['file_size']+1
    if debug:
        print("%d %d"%(c0_idx,c1_idx))
    n_read = 0
    for f_idx in range(f0_idx,fn_idx):
        c0 = 0
        c1 = gdf['file_size']
        if f_idx == f0_idx:
            c0 = c0_idx
        if f_idx+1 == fn_idx:
            c1 = c1_idx
        if debug:
            print("open file %d %s"%(f_idx,files[f_idx]))

        if gdf['cache_idx'][chan] != f_idx:
            a = numpy.fromfile(files[f_idx],dtype=gdf['dtype'])
            gdf['cache'][chan] = numpy.array(a[gdf['re_idx']]*gdf['scale'] + 1j*a[gdf['im_idx']]*gdf['scale'],dtype=numpy.complex64)
            gdf['cache_idx'][chan] = f_idx
        if debug:
            print("%d %d %d"%(len(res_vec),c0,c1))
        res_vec[numpy.arange(c1-c0)+n_read] = gdf['cache'][chan][c0:c1]
        n_read = n_read + (c1-c0)
        if debug:
            print("indices %d %d %d"%(c0,c1,n_read))
    return(res_vec)

