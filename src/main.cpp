#include <iostream>
#include "BCHCodes.h"
#include <fstream>

void usage()
{
    std::cout<<"Usage: bch_matrices n d filename\nn - length of the code;\nd - minimal distance of the code;\nfilename - optional, if given - outputs matrices in the file, else - outputs in console"<<std::endl;
}

int main(int argc,  char** argv) {
    if(!(argc == 3 || argc == 4))
    {
        usage();
    }
    else
    {
        unsigned n, d;
        n = atoi(argv[1]);
        d = atoi(argv[2]);
        std::streambuf* buf;
        std::ofstream fout;
        if(argc == 3)
        {
            buf = std::cout.rdbuf();
        }
        else
        {
            fout.open(argv[3]);
            buf = fout.rdbuf();
        }
        std::ostream output(buf);
        build_bch_matrices(n, d, output);
        build_bch_matrices_bit_order(n, d, output);
        fout.close();
    }
    return 0;
}