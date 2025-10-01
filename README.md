# Ray Tracer

I followed the Ray Tracing in one Weekend tutorial at <https://raytracing.github.io/>.

Following along with the explanations, I implemented my own version of the code.

`gpng_src.cpp` and `gpng.h` contain my own code to output an image into a file that follows the PNG standard. Currently there is no compression.

`gmath_src.cpp` and `gmath.h` contain mathematical tools useful for ray tracing, including random number generation and support for 3D vectors and their operations.

`main.cpp` contains the ray tracing code.

To compile, run `g++ -Wall -o out *.cpp` in the project's root directory.

### Example Image
![alt text](https://github.com/suspicious-salmon/Ray-Tracing-in-One-Weekend/blob/master/cover_image.png?raw=true)
