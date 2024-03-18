"""
This script converts a binary matrix in this format:

5
0 1 2 3 4
1 2 4
3 4
0 1 2
1 2 3 4

Into a plain representation for the k²-tree
"""

import struct, sys

def read_graph(infilepath):
    infile = open(infilepath, 'r')
    line = infile.readline()
    assert(line and line[-1]=='\n')

    num_nodes = int(line[:-1])
    yield num_nodes
    line = infile.readline()

    while line :
        assert(line[-1]=='\n')
        yield list(map(int, line[:-1].split()))
        line = infile.readline()

    infile.close()
    return

def num_nodes_edge(infilepath) :
    nnodes,nedges = 0,0
    for i,e in enumerate(read_graph(infilepath)) :
        if i==0 :
            nnodes = e
            continue
        nedges += len(e)
    return nnodes,nedges

def write_binary_graph(infilepath, outfilepath) :
    outfile = open(outfilepath, 'wb')
    nnodes,nedges = num_nodes_edge(infilepath)
    outfile.write(struct.pack('i', nnodes))
    outfile.write(struct.pack('q', nedges))

    for i,adjl in enumerate(read_graph(infilepath)) :
        if i==0 :
            continue

        outfile.write(struct.pack('i', -1))  # Mark the start of a new adjacency list
        for neighbour in adjl:
            outfile.write(struct.pack('i', neighbour+1)) #1-based!
    outfile.close()
    return

def main():
    if len(sys.argv) != 2+1 :
        print('Usage is:', sys.argv[0], '<infilepath> <outfilepath>')
        exit(1)
    
    #args
    infilepath = sys.argv[1]
    outfilepath = sys.argv[2]

    write_binary_graph(infilepath, outfilepath)

    print(f"Binary graph file '{outfilepath}' has been created.")

if __name__ == '__main__':
    main()