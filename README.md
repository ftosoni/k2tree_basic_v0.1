# k2tree_basic_v0.1
k²-tree implementation by the University of A Coruna

## Install
```bash
mkdir build
cd build
cmake ..
make -j 12
```

## Create plain format
In the `gen_graph` folder
```bash
python3 edgedescr2plain.py  ../data/example14.graph-txt ../data/example14.kt-plain
```

## Compress
In the `build` folder
```bash
./build_tree ../data/example14.kt-plain ../data/example14
```

## Multiply
In the `build` folder
```bash
./multiply ../data/example14.kt ../data/invec14
```
