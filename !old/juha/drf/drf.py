import glob, numpy, math, datetime, time, pickle, os, h5py, re

debug = True

def new_drf(dirn):
    result = {}

    f = h5py.File("%s/sampler_config.hdf5"%(dirn),"r")
    for k in f.keys():
        result[k] = f[k].value
    f.close()

    files = glob.glob("%s/*/*.drf"%(dirn))
    files.sort()
    
    result["file_timestamps"]=[]
    
    for f in files:
        result["file_timestamps"].append(numpy.uint64(float(re.search("chunck-(.*).drf",f).group(1))))

    result["file_timestamps"] = numpy.array(result["file_timestamps"],dtype=numpy.uint64)*numpy.uint64(result["sample_rate"])

    # global sample index.
    result["t0"] = result["file_timestamps"][0]
    result["t1"] = result["file_timestamps"][len(files)-2]

    result['file_size'] = numpy.uint64(os.path.getsize(files[0]))/numpy.uint64(result["itemsize"])
    result['file_list'] = files
    
    # create a continous list of files. missing files indicated with ""
    n_virtual_files = numpy.uint64((result["t1"]-result["t0"])/result["file_size"])
    formal_list = [""]*n_virtual_files
    fi = 0
    for i in range(len(files)-2):
        f_idx = (result["file_timestamps"][i]-result["t0"])/result["file_size"]
        formal_list[f_idx] = files[i]

    result["formal_file_list"]=formal_list
    nsecs_per_file = numpy.uint64(result["file_size"]/result["sample_rate"])

    result["cache"] = numpy.zeros(result['file_size'],dtype=numpy.complex64)
    result['cache_idx'] = -1

    result['scale'] = 1.0
    result['dirn'] = dirn
    result['re_idx'] = numpy.arange(result['file_size'],dtype=numpy.int64)*2
    result['im_idx'] = numpy.arange(result['file_size'],dtype=numpy.int64)*2+1

    if debug:
        print("Read digital_rf dir. Version %s n_files=%d file_size=%d t0 %d t1 %d"%(result["version"],len(result["file_list"]),result['file_size'],result['t0'],result['t0']))
    return(result)

#
# using global sample index, read vector. 
# throw exception if file missing.
# 
def read_vec(drf, idx, length):
    idx = numpy.uint64(idx)-drf["t0"]
    length = numpy.uint64(length)
    files = drf['formal_file_list']    
    f0_idx = int(math.floor(idx/drf['file_size']))
    fn_idx = int(math.floor((idx+length-1)/drf['file_size']))+1
    res_vec = numpy.zeros([length],dtype=numpy.complex64)
    if debug:
        print "f0 ",f0_idx," ",fn_idx
    c0_idx = numpy.mod(idx,drf['file_size'])
    c1_idx = numpy.mod((idx+length-numpy.uint64(1)),drf['file_size'])+numpy.uint64(1)
    if debug:
        print c0_idx," ",c1_idx
    n_read = numpy.uint64(0)
    for f_idx in range(f0_idx,fn_idx):
        c0 = numpy.uint64(0)
        c1 = drf['file_size']
        if f_idx == f0_idx:
            c0 = c0_idx
        if f_idx+1 == fn_idx:
            c1 = c1_idx
        if debug:
            print "open file ",f_idx," ",files[f_idx],"\n"

        if drf['cache_idx'] != f_idx:
            drf["cache"] = None
            if files[f_idx] == "":
                print("file missing!")
                drf["cache"] = numpy.zeros([drf["file_size"]],dtype=numpy.complex64)
            else:
                a = numpy.fromfile(files[f_idx],dtype=drf['dtype'])
                drf['cache'] = numpy.array(a[drf['re_idx']]*drf['scale'] + 1j*a[drf['im_idx']]*drf['scale'],dtype=numpy.complex64)
            drf['cache_idx'] = f_idx
        if debug:
            print len(res_vec)," c0 ",c0," c1 ",c1
        res_vec[numpy.arange(c1-c0,dtype=numpy.uint64)+n_read] = drf['cache'][c0:c1]
        n_read = n_read + (c1-c0)
        if debug:
            print "indices ",c0,"-",c1," n_read ",n_read
    return(res_vec)
