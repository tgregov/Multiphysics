# Compile on Mac
cd ../..
#. ./envs/linux-macos.sh
rm -rf ./build
mkdir build
cd build
cmake ..
make
cd ..
cd ./run/mac/