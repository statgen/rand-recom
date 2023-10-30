#rand-recom


## Build
```
cmake -P cmake/get-dependencies.cmake cget
mkdir build; cd build
cmake -DCMAKE_PREFIX_PATH="../cget" -DCMAKE_CXX_FLAGS="-I../cget/include" -DCMAKE_BUILD_TYPE=Release ..
make
```
