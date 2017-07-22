# linear-optics
A fast algorithm for simulating linear optical quantum circuits written in C++.

Compile with the Makefile.

For the GCC compiler, use command line flags

  -Ofast -funroll-loops

For the intel compiler, use the command line flags

  -O3 -xavx -inline-forceinline -funroll-loops

The linear algebra, Eigen, is required:

  http://eigen.tuxfamily.org/index.php?title=Main_Page

In main.cpp, we demonstrate how to use the LinearOpticalTransform object. We simulate the action of a Mach-Zehnder interferometer acting on a simple optical quantum state. Follow the comments for details.

Any requests, suggestions, , questions, comments etc, please email

  jsmith74@tulane.edu
