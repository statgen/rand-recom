# rand-recom
A tool for randomly applying artificial recombination.

## Build
```
cmake -P cmake/get-dependencies.cmake cget
mkdir build; cd build
cmake -DCMAKE_PREFIX_PATH="../cget" -DCMAKE_CXX_FLAGS="-I../cget/include" -DCMAKE_BUILD_TYPE=Release ..
make
```

## Usage
```
# (target_seg_legnth * records_in_chr20) / chr20_length
RECOM_PARAM=$(( (25000000 * 9291832) / 64444167 ))

SEED=12345

# Currently defaults to uncompressed BCF output to stdout. Will add CLI parameters for output later.
rand-recom $SEED $RECOM_PARAM input.bcf | bcftools view -Ob -o out_file.bcf
```
