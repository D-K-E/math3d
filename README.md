# math3d
yet another 3d math library

This one doesn't throw any exceptions. It doesn't allocate memory on the heap.
It doesn't contain any virtual function, nor a recursive function.
It also doesn't use any function pointers. No concurrency and threading
libraries are used either.

Hence it is highly compatible with OpenCL-C++ and CUDA at the same time.
For CUDA, one would need to add `__host__`, `__device__` attributes in front
of the methods. 
For OpenCL-C++ you probably don't need anything, just copy-paste the stuff
into your kernel and it should be good to go.

# Usage

Here are some usage examples to get an idea:
