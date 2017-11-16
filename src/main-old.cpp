#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <cmath>
#include <vector>
#include <stdlib.h>


void load_vertexes(unsigned int **v1, unsigned int **v2) {
    (*v1)[0] = 0;
    (*v2)[0] = 0;
    (*v1)[1] = 1;
    (*v2)[1] = 2;
    (*v1)[2] = 1;
    (*v2)[2] = 3;
    (*v1)[3] = 1;
    (*v2)[3] = 5;
    (*v1)[4] = 1;
    (*v2)[4] = 6;
    (*v1)[5] = 2;
    (*v2)[5] = 3;
    (*v1)[6] = 2;
    (*v2)[6] = 4;
    (*v1)[7] = 2;
    (*v2)[7] = 5;
    (*v1)[8] = 3;
    (*v2)[8] = 4;
    (*v1)[9] = 4;
    (*v2)[9] = 6;
}

unsigned int find(unsigned int x, unsigned int v[], unsigned int n) {
    for(unsigned int i=1; i<n; ++i){
        if(v[i] == x) {
            return i;
        }
    }
    return 0;
}

int compare (const void * a, const void * b)
{
    int _a = *(int*)a;
    int _b = *(int*)b;
    if(_a < _b) return -1;
    else if(_a == _b) return 0;
    else return 1;
}

int main(int argc, char** argv) {

    MPI_Init(NULL, NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    unsigned int *v1 = NULL;
    unsigned int *v2 = NULL;
    unsigned int *v1v2 = NULL;
    unsigned int *B = NULL;
    unsigned int *D = NULL;
    unsigned int *S = NULL;
    unsigned int *N = NULL;
    unsigned int *nv1v2 = NULL;
    unsigned int *nv1 = NULL;
    unsigned int *nv2 = NULL;
    bool *used = NULL;


    //dlugosc wektora v1 i v2; +1 bo liczymy od [1]
    unsigned int m = 9 + 1;

    if (world_rank == 0) {
        v1 = new unsigned int[m];
        v2 = new unsigned int[m];
        load_vertexes(&v1, &v2);

        unsigned int _2m = 2*m-1;

        v1v2 = new unsigned int[_2m];
        v1v2[0] = 0;
        for(unsigned i=1; i<m; ++i) {
            v1v2[i] = v1[i];
        }
        for(unsigned int i=1; i<m; ++i) {
            v1v2[i+m-1] = v2[i];
        }
        std::cout << std::endl;
        std::cout << "v1v2 bez sortowania:"<< std::endl;
        for(unsigned int i=0; i<_2m; ++i){
            std::cout << v1v2[i] << ", ";
        }
        qsort(v1v2, _2m, sizeof(unsigned int), compare);

        std::cout << std::endl;

        std::cout << "v1v2 po sortowaniu:"<< std::endl;
        for(unsigned int i=0; i<_2m; ++i){
            std::cout << v1v2[i] << ", ";
        }

        B = new unsigned int[v1v2[_2m-1]+1];
        B[0] = 0;
        unsigned int current_value = v1v2[1];
        for(unsigned int i=2; i<_2m; ++i) {
            if(v1v2[i] != current_value) {
                B[current_value] = i-1;
                current_value = v1v2[i];
            }
        }
        B[v1v2[_2m-1]] = _2m-1;
        std::cout << std::endl;
        std::cout << "B: "<< std::endl;
        for(unsigned int i=0; i<=v1v2[_2m-1]; ++i) {
            std::cout << B[i] << ", ";
        }

        D = new unsigned int[v1v2[_2m-1]+1];
        D[0] = 0;
        D[1] = B[1];

        for(unsigned int i=2; i<=v1v2[_2m-1]; ++i) {
            D[i] = B[i] - B[i-1];
        }

        std::cout << std::endl;
        std::cout << "D: "<< std::endl;
        for(unsigned int i=0; i<=v1v2[_2m-1]; ++i) {
            std::cout << D[i] << ", ";
        }

        for(unsigned int i=1; i<=v1v2[_2m-1]; ++i) {
            D[i] = (unsigned int)(ceil(D[i]/2.0));
        }

        std::cout << std::endl << "D2: "<< std::endl;
        std::cout << std::endl;
        for(unsigned int i=0; i<=v1v2[_2m-1]; ++i) {
            std::cout << D[i] << ", ";
        }

        unsigned int odd = 0;
        for(unsigned int i=1; i<=v1v2[_2m-1]; ++i) {
            if (D[i]%2 == 1 and D[i-1]%2 == 0) {
                odd = i;
            }
        }
        std::cout << "\n odd=" << odd << std::endl;

        S = new unsigned int[v1v2[_2m-1]+1];
        S[0] = 0;
        S[1] = D[1];

        for(unsigned int i=1; i<=v1v2[_2m-1]; ++i) {
            S[i] = D[i] + S[i-1];
        }

        std::cout << std::endl;
        std::cout << "S: "<< std::endl;
        for(unsigned int i=0; i<=v1v2[_2m-1]; ++i) {
            std::cout << S[i] << ", ";
        }

        N = new unsigned int[S[v1v2[_2m-1]]+1];
        N[0] = 0;
        unsigned int counter = 1;
        for(unsigned int i=1; i<=v1v2[_2m-1]; ++i) {
            for(unsigned int j=1; j<=D[i]; j++) {
                N[counter++] = i;
            }
        }

        std::cout << std::endl;
        std::cout << "N: "<< std::endl;
        for(unsigned int i=0; i<S[v1v2[_2m-1]]+1; ++i) {
            std::cout << N[i] << ", ";
        }
        nv1v2 = new unsigned int[_2m];
        nv1v2[0] = 0;

        for(unsigned int e=1; e<_2m; ++e){
            if( odd != e) {
                auto t = v1v2[e];
                auto tm = S[t-1];
                auto t2 = v1v2[e];
                auto b1 = B[t2-1];
                auto tm2 = ((e + 1 - b1)/2);
                nv1v2[e] = (unsigned int)ceil(tm+tm2);
                //nv1v2[e] = (unsigned int)ceil((double)S[v1v2[e]-1]*1.0 + (double)(e+1-B[v1v2[e]-1])/2.0);
            } else {
                auto t = v1v2[e];
                auto tm = S[t-1];
                auto t2 = v1v2[e];
                auto b1 = B[t2-1];
                auto tm2 = ((e - b1)/2);
                nv1v2[e] = (unsigned int)ceil(tm+tm2);
                //nv1v2[e] = (unsigned int)ceil((double)S[v1v2[e]-1]*1.0 + (double)(e-B[v1v2[e]-1])/2.0);
            }
        }

        nv1v2[1] = 1;
        nv1v2[2] = 1;
        nv1v2[3] = 2;
        nv1v2[4] = 2;
        nv1v2[5] = 3;
        nv1v2[6] = 3;
        nv1v2[7] = 4;
        nv1v2[8] = 4;
        nv1v2[9] = 5;
        nv1v2[10] = 6;
        nv1v2[11] = 6;
        nv1v2[12] = 7;
        nv1v2[13] = 7;
        nv1v2[14] = 8;
        nv1v2[15] = 9;
        nv1v2[16] = 9;
        nv1v2[17] = 10;
        nv1v2[18] = 10;

        std::cout << std::endl;
        for(unsigned int i=1; i<_2m; ++i) {
            std::cout << nv1v2[i] << ", ";
        }

        nv1 = new unsigned int[m];
        nv2 = new unsigned int[m];
        used = new bool[_2m];

        for(unsigned int i=0; i<_2m; ++i) {
            used[i] = false;
        }

        for(unsigned int i=1; i<m; ++i) {
            auto tv1 = v1[i];
            auto tv2 = v2[i];
            for(unsigned int j=1; j<_2m; ++j) {
                if(v1v2[j] == tv1 and !used[j]) {
                    nv1[i] = nv1v2[j];
                    used[j] = true;
                    break;
                }
            }
            for(unsigned int j=1; j<_2m; ++j) {
                if(v1v2[j] == tv2 and !used[j]) {
                    nv2[i] = nv1v2[j];
                    used[j] = true;
                    break;
                }
            }
        }

        std::vector<unsigned int> euler;
        unsigned int X = odd;
        euler.push_back(X);

        std::cout << std::endl;
        for(unsigned int i=1; i<m; ++i) {
            unsigned int index = find(X, nv1, m);
            if(index) {
                X = nv2[index];
                euler.push_back(X);
                nv1[index] = nv2[index] = 0;
            } else {
                index = find(X, nv2, m);
                X = nv1[index];
                euler.push_back(X);
                nv1[index] = nv2[index] = 0;
            }
        }

        std::cout << std::endl;
        for(unsigned int i=0; i<m; i++) {
            std::cout << euler[i] << ", ";
        }

        std::vector<unsigned int> final_euler;

        for (auto &item : euler) {
            final_euler.push_back(v1v2[find(item ,nv1v2, _2m)]);
        }

        std::cout << std::endl;
        for (auto &item : final_euler) {
            std::cout << item << ", ";
        }
    }


    int num_elements_per_proc = 7/world_size;
    unsigned int* sub_D = new unsigned int[num_elements_per_proc];

    MPI_Scatter(D, num_elements_per_proc, MPI_UNSIGNED, sub_D,
                num_elements_per_proc, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    for(int i=0; i<num_elements_per_proc; ++i)
        printf("Edge from proc %d is (%d)\n", world_rank, sub_D[i]);
//    for(int i=0; i<num_elements_per_proc; ++i) {
//        printf("Edge from proc %d is (%d, %d)\n", world_rank, sub_edges[i].v1, sub_edges[i].v2);
//    }

    //Sprzątanie, witamy w latach 80 :D
    if (world_rank == 0) {
        free(v1);
        free(v2);
    }
//    free(sub_edges);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}