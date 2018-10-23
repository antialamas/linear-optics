This is a fork from https://github.com/jsmith74/linear-optics.git

# linear-optics
A fast, hybrid algorithm for simulating linear optical quantum circuits written in C++.

For more details and author contact information, refer to

https://arxiv.org/abs/1711.01319

or

https://journals.aps.org/pra/abstract/10.1103/PhysRevA.97.012320 

It is preferable to use Ryser's algorithm over the methods described in the links above for evaluating certain elements of the matrix A(U). This code automatically determines which method to use for each element of A(U) and will construct the full matrix using an optimal hybrid of methods.

Compile with the included Makefile.

The linear algebra, Eigen, is required:

  http://eigen.tuxfamily.org/index.php?title=Main_Page

In main.cpp, we demonstrate how to use the LinearOpticalTransform object. We simulate the action of a Mach-Zehnder interferometer acting on a simple optical quantum state. Follow the comments in the code for details.

Feel free to fork, change, modify etc.


