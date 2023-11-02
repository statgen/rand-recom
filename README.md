# rand-recom
A tool for randomly applying artificial recombination.

## Build and Install
```
cmake -P cmake/get-dependencies.cmake cget
mkdir build; cd build
cmake -DCMAKE_PREFIX_PATH="../cget" -DCMAKE_CXX_FLAGS="-I../cget/include" -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

## Usage
```
TARGET_SEGMENT_LENGTH=25000000
RANDOM_SEED=12345

# Recombination will occur once every $TARGET_SEGMENT_LENGTH on average
rand-recom input.bcf --seed $RANDOM_SEED --target-length $TARGET_SEGMENT_LENGTH -O bcf -o output.bcf
```
