docker run -it --rm -v $(pwd):/sosi  conanio/clang14-ubuntu16.04:latest  bash -c 'cd /sosi; clang++  affine.cpp'
