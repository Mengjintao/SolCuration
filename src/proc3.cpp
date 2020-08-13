#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
using namespace std;

#define CHUNKSIZE 1
#define ThreadNum 50

int main()
{
        vector<string> smiles;
        string line;

        ifstream myfile("extend.smi");
        while(getline(myfile, line))
                smiles.push_back(line);
        myfile.close();

        int i, tid;
        int chunk = CHUNKSIZE;
        omp_set_num_threads(ThreadNum);

        #pragma omp parallel for schedule(dynamic, chunk) \
           shared(smiles) private(i, tid)
        for(i=0;i<smiles.size();i++)
        {
                tid = omp_get_thread_num();
                string line = smiles[i];
                char str[1000];
                char idx[100];
                sscanf(line.c_str(), "%s %s", str, &idx);
                printf("Thread %d: %s %s\n", tid, str, idx);
                string cmd = string("echo ") + string("\"") + line + string("\"") +  string("| xedex -m 10 -i l -o s >") + string(idx) + string(".sdf 2> ") + string(idx) + string(".log");
                printf("%s\n", cmd.c_str());
                system(cmd.c_str());
        }
        return 1;
}
~
