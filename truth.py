#-------------Wrap verosimilitud-------------_#

def pyfunc(np.ndarray[np.int32_t, ndim=1] in_array):    
    cdef np.ndarray out_array = np.zeros((512,), dtype = np.int32)
    n = len(in_array)
    mymodule.c_func(<int *> in_array.data, n, <int *> out_array.data)
    return out_array
