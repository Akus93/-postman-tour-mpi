#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <fstream>
#include <algorithm>

bool Printing = 0;


//LOAD VERTICES

unsigned int get_file_length(const char* file_name) {
    unsigned int count = 0;
    char c;
    FILE *fp = fopen(file_name, "r");
    if (fp == NULL) {
        printf("Could not open file %s", file_name);
        return 0;
    }
    for (c = getc(fp); c != EOF; c = getc(fp))
        if (c == '\n')
            count = count + 1;
    fclose(fp);
    return count;
}

void load_vertices(const char* file_name, unsigned int **v1, unsigned int **v2, unsigned int &size) {
    size = get_file_length(file_name);
    *v1 = new unsigned int[size];
    *v2 = new unsigned int[size];
    (*v1)[0] = 0;
    (*v2)[0] = 0;
    std::ifstream file(file_name);
    unsigned int x1, x2;
    int i = 1;
    while (file >> x1 >> x2) {
        (*v1)[i] = x1;
        (*v2)[i] = x2;
        i++;
    }
}


// SORT

unsigned int * merge(unsigned int *A, unsigned int asize, unsigned int *B, unsigned int bsize) {
    unsigned int ai, bi, ci, i;
    unsigned int *C;
    unsigned int csize = asize+bsize;

    ai = 0;
    bi = 0;
    ci = 0;

    C = (unsigned int *)malloc(csize*sizeof(unsigned int));
    while ((ai < asize) && (bi < bsize)) {
        if (A[ai] <= B[bi]) {
            C[ci] = A[ai];
            ci++; ai++;
        } else {
            C[ci] = B[bi];
            ci++; bi++;
        }
    }

    if (ai >= asize)
        for (i = ci; i < csize; i++, bi++)
            C[i] = B[bi];
    else if (bi >= bsize)
        for (i = ci; i < csize; i++, ai++)
            C[i] = A[ai];

    for (i = 0; i < asize; i++)
        A[i] = C[i];
    for (i = 0; i < bsize; i++)
        B[i] = C[asize+i];

    return C;
}

void swap(unsigned int *v, unsigned int i, unsigned int j) {
    unsigned int t;
    t = v[i];
    v[i] = v[j];
    v[j] = t;
}

void m_sort(unsigned int *A, unsigned int min, unsigned int max)
{
    unsigned int *C;
    unsigned int mid = (min+max)/2;
    unsigned int lowerCount = mid - min + 1;
    unsigned int upperCount = max - mid;

    if (max == min) {
        return;
    } else {
        m_sort(A, min, mid);
        m_sort(A, mid+1, max);
        C = merge(A + min, lowerCount, A + mid + 1, upperCount);
    }
}

//END SORT

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

int compare (const void * a, const void * b)
{
    int _a = *(int*)a;
    int _b = *(int*)b;
    if(_a < _b) return -1;
    else if(_a == _b) return 0;
    else return 1;
}

unsigned int* generatev1v2(int world_rank, int world_size, unsigned int *v1, unsigned int *v2, unsigned int msize, unsigned int m2size){
    if(world_rank == 0){
        unsigned int *v1v2 = new unsigned int[m2size];
        v1v2[0] = 0;
        // code from pdf. not working :(
        // for (unsigned int e = 1; e < m2size; e++){
        //  v1v2[2*(e-1)+1] = v1[e];
        //  v1v2[2*e] = v2[e];
        // }
        for(unsigned i=1; i<msize; ++i) {
            v1v2[i] = v1[i];
        }
        for(unsigned int i=1; i<msize; ++i) {
            v1v2[i+msize-1] = v2[i];
        }
        //qsort(v1v2, m2size, sizeof(unsigned int), compare);
        if (Printing){
            std::cout << "\nv1v2:\n\t\t";
            for(unsigned int i=0; i<m2size; ++i){
                std::cout << v1v2[i] << ", ";
            }
        }
        return v1v2;
    }
    else {
        return (unsigned int*) malloc(sizeof(unsigned int) * m2size);
    }
} 

unsigned int* generateB(int world_rank, int world_size, unsigned int *v1v2, unsigned int msize, unsigned int m2size){
    if(world_rank == 0){
        unsigned int* B = new unsigned int[v1v2[m2size-1]+1];
        B[0] = 0;
        unsigned int current_value = v1v2[1];
        for(unsigned int i=2; i<m2size; ++i) {
            if(v1v2[i] != current_value) {
                B[current_value] = i-1;
                current_value = v1v2[i];
            }
        }
        B[v1v2[m2size-1]] = m2size-1;
        if(Printing){
            std::cout << std::endl;
            std::cout << "B:\n\t\t ";
            for(unsigned int i=0; i<=v1v2[m2size-1]; ++i) {
                std::cout << B[i] << ", ";
            }
        }
        return B;
    }
    else {
        return NULL;
    }
}
unsigned int* generateBparallel(int world_rank, int world_size, unsigned int *v1v2, unsigned int msize, unsigned int m2size){
    MPI_Bcast(v1v2, m2size, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    int Bsize = v1v2[m2size-1]+1;
    unsigned int* Bbuffer = new unsigned int[Bsize];
    Bbuffer[0] = 0;
    for(int i=1; i<Bsize; i++){ Bbuffer[i] = 0; }
    unsigned int current_value;
    for(unsigned int i=2+world_rank; i<m2size; i+=world_size) {
        current_value = v1v2[i-1];
        if(v1v2[i] != current_value) {
            Bbuffer[current_value] = i-1;
        }
    };

    Bbuffer[v1v2[m2size-1]] = m2size-1;
    MPI_Barrier(MPI_COMM_WORLD);
    /* Merge Bbuffers */
    unsigned int* Brecvbuffer = new unsigned int[Bsize];

    if (world_rank)
        MPI_Send(Bbuffer, Bsize, MPI_UNSIGNED, 0, 50, MPI_COMM_WORLD);
//    if(!world_rank)
//        MPI_Send(Bbuffer, Bsize, MPI_UNSIGNED, 0, 50, MPI_COMM_WORLD);

    if(world_rank == 0){
        for(int i=1; i<world_size; i++){
            MPI_Status status;
            MPI_Recv(Brecvbuffer, Bsize, MPI_UNSIGNED, MPI_ANY_SOURCE, 50, MPI_COMM_WORLD, &status);
            for(int j=0; j<Bsize; j+=1){
                if(Brecvbuffer[j] > Bbuffer[j]){
                    Bbuffer[j] = Brecvbuffer[j];
                }
            }
        }
        /* print calculated B*/
        if(Printing){
            std::cout << "\nParallel B: \n\t\t";
            for(int i = 0; i<Bsize; i++){
                std::cout << Bbuffer[i] << ", ";
            }
        }

        return Bbuffer;
    }
}

unsigned int* generateD(int world_rank, int world_size, unsigned int *v1v2, unsigned int *B, unsigned int msize, unsigned int m2size){
    if(world_rank == 0){
        unsigned int Dsize = v1v2[m2size-1]+1;
        unsigned int* D = new unsigned int[Dsize];
        D[0] = 0;
        D[1] = B[1];

        for(unsigned int i=2; i<Dsize; i++) {
            D[i] = B[i] - B[i-1];
        }

        std::cout << std::endl;
        if(Printing){
            std::cout << "D: ";
            for(unsigned int i=0; i<Dsize; ++i) {
                std::cout << D[i] << ", ";
            }
        }
        return D;
    }
    else {
        return NULL;
    }
}
unsigned int* generateDparallel(int world_rank, int world_size, unsigned int *v1v2, unsigned int *B, unsigned int msize, unsigned int m2size, unsigned int Dsize){
    if(world_rank != 0){
        B = (unsigned int*) malloc(sizeof(unsigned int)*Dsize);
    }
    MPI_Bcast(B, Dsize, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    unsigned int* D = new unsigned int[Dsize];
    D[0] = 0;
    D[1] = B[1];
    for(int i=2; i<Dsize; i++){D[i] = 0;}
    for(int i=2+world_rank; i<Dsize; i+=world_size){
        D[i] = B[i] - B[i-1];
    }

    /* Merge all D */
    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank)
        MPI_Send(D, Dsize, MPI_UNSIGNED, 0, 51, MPI_COMM_WORLD);
//    MPI_Barrier(MPI_COMM_WORLD);
    if(world_rank == 0){
        unsigned int* Drecvbuffer = new unsigned int[Dsize];
        for (int i=1; i<world_size; i++){
            MPI_Status status;
            MPI_Recv(Drecvbuffer, Dsize, MPI_UNSIGNED, MPI_ANY_SOURCE, 51, MPI_COMM_WORLD, &status);
            for(unsigned j=0; j<Dsize; j++){
                if(Drecvbuffer[j] != 0){
                    D[j] = Drecvbuffer[j];
                }
            }
        }
        if(Printing){
            std::cout << std::endl;
            std::cout << "D parallel: ";
            for(unsigned int i=0; i<Dsize; ++i) {
                std::cout << D[i] << ", ";
            }
        }
        return D;
    }
}

unsigned int* generateS(int world_rank, int world_size, unsigned int *v1v2, unsigned int *D, unsigned int m2size){
    if(world_rank == 0){

        unsigned int Ssize = v1v2[m2size-1]+1;
        unsigned int* S = new unsigned int[Ssize];

        /* prepare DD table from D */
        unsigned int* DD = new unsigned int[Ssize];
        if(Printing){std::cout << "\nDD: ";}
        for(unsigned int i=0; i<Ssize; ++i) {
            D[i] = (unsigned int)(ceil(D[i]/2.0));
            if(Printing){std::cout << D[i] << ", ";}
        }

        S[0] = 0;
        S[1] = D[1];

        for(unsigned int i=1; i<Ssize; ++i) {
            S[i] = D[i] + S[i-1];
        }

        if(Printing){
            std::cout << std::endl;
            std::cout << "S: ";
            for(unsigned int i=0; i<Ssize; ++i) {
                std::cout << S[i] << ", ";
            }
        }
        return S;
    }
    else {
        return NULL;
    }
}

unsigned int getOdd(int world_rank, int world_size, unsigned int* v1v2, unsigned int *D, unsigned int m2size){
    if(world_rank == 0){
        unsigned int odd = 0;
        unsigned int Dsize = v1v2[m2size-1]+1;
        for(unsigned int i=1; i<Dsize; ++i) {
            if (D[i]%2 == 1 and D[i-1]%2 == 0) {
                odd = i;
            }
        }
        if(Printing){std::cout << "\nodd=" << odd;}
        return odd;
    }
    else {
        return 0;
    }
}

unsigned int* generateN(int world_rank, int world_size, unsigned int* v1v2, unsigned int *D, unsigned int *S, unsigned int m2size){
    if(world_rank == 0){
        unsigned int Nsize = S[v1v2[m2size-1]]+1;
        unsigned int* N = new unsigned int[Nsize];
        N[0] = 0;
        unsigned int counter = 1;
        for(unsigned int i=1; i<=v1v2[m2size-1]; ++i) {
            for(unsigned int j=1; j<=D[i]; j++) {
                N[counter++] = i;
            }
        }

        if(Printing){
            std::cout << std::endl;
            std::cout << "N: ";
            for(unsigned int i=0; i<S[v1v2[m2size-1]]+1; ++i) {
                std::cout << N[i] << ", ";
            }
        }
        return N;
    }
    else {
        return NULL;
    }
}

unsigned int* generateNv1v2(int world_rank, int world_size, unsigned int* v1v2, unsigned int* S, unsigned int* B, unsigned int odd, unsigned int m2size){
    if(world_rank == 0){
        unsigned int* nv1v2 = new unsigned int[m2size];
        nv1v2[0] = 0;

        for(unsigned int e=1; e<m2size; ++e){
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
        if(Printing){
            std::cout << "\nnv1v2: ";
            for(unsigned int i=1; i<m2size; ++i) {
                std::cout << nv1v2[i] << ", ";
            }
        }
        return nv1v2;
    }
    else {
        return NULL;
    }
}

unsigned int* mockupNv1v2(int world_rank){
    if(world_rank == 0){
        unsigned int* nv1v2 = new unsigned int[19];
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
        std::cout << "\nnv1v2: ";
        for(unsigned int i=1; i<19; ++i) {
            std::cout << nv1v2[i] << ", ";
        }
        std::cout << " (mocked)";
        return nv1v2;
    } else {
        return NULL;
    }
}

unsigned int find(unsigned int x, unsigned int v[], unsigned int n) {
    for(unsigned int i=1; i<n; ++i){
        if(v[i] == x) {
            return i;
        }
    }
    return 0;
}

void generateNvs(unsigned int* nv1, unsigned int* nv2, int world_rank, int world_size, unsigned int* v1, unsigned int* v2, unsigned int* v1v2, unsigned int* nv1v2, unsigned int m2size, unsigned int m){
    if(world_rank == 0){
        bool* used = new bool[m2size];
        for(unsigned int i=0; i<m2size; ++i) {
            used[i] = false;
        }        
        for(unsigned int i=1; i<=m; ++i) {
            auto tv1 = v1[i];
            auto tv2 = v2[i];
            for(unsigned int j=1; j<m2size; ++j) {
                if(v1v2[j] == tv1 and !used[j]) {
                    nv1[i] = nv1v2[j];
                    used[j] = true;
                    break;
                }
            }
            for(unsigned int j=1; j<m2size; ++j) {
                if(v1v2[j] == tv2 and !used[j]) {
                    nv2[i] = nv1v2[j];
                    used[j] = true;
                    break;
                }
            }
        }
    }
}

void getEuler(int world_rank, int world_size, unsigned int odd, unsigned int m, unsigned int* nv1, unsigned int* nv2, unsigned int* nv1v2, unsigned int* v1v2, unsigned int m2size, unsigned int msize){
    if(world_rank ==0){
        std::vector<unsigned int> euler;
        unsigned int X = odd;
        euler.push_back(X);

        for(unsigned int i=1; i<m; ++i) {
            unsigned int index = find(X, nv1, msize);
            if(index) {
                X = nv2[index];
                euler.push_back(X);
                nv1[index] = nv2[index] = 0;
            } else {
                index = find(X, nv2, msize);
                X = nv1[index];
                euler.push_back(X);
                nv1[index] = nv2[index] = 0;
            }
        }
        if(Printing){
            std::cout << "\nEulerian path in N graph:\n";
            for(unsigned int i=1; i<m; i++) {
                std::cout << euler[i] << "-";
            }
        }
        std::vector<unsigned int> final_euler;
        for (auto &item : euler) {
            final_euler.push_back(v1v2[find(item ,nv1v2, m2size)]);
        }
        if(Printing){
            std::cout << "\nEulerian path in main graph:\n";
            for (auto &item : final_euler) {
                std::cout << item << "-";
            }
            std::cout << "\n";
        }
    }
}

int main(int argc, char** argv) {

    MPI_Init(NULL, NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    /* vectors */
    unsigned int *v1 = NULL;
    unsigned int *v2 = NULL;
    unsigned int *v1v2;

    unsigned int *B = NULL;
    unsigned int *D = NULL;
    unsigned int *S = NULL;
    unsigned int *N = NULL;
    unsigned int *nv1v2 = NULL;
    unsigned int *nv1 = NULL;
    unsigned int *nv2 = NULL;
    unsigned int odd = 0;

    /* sizes... */
    unsigned int m;             // number of edges
    unsigned int msize;     // size of v1 or v2
    unsigned int m2;           // length of v1v2, nv1v2
    unsigned int m2size;   // size of v1v2, nv1v2
    unsigned int sort_chunk_size = 0;

    /* sort */
    unsigned int *sort_v1v2 = NULL;
    unsigned int * chunk;
    unsigned int * other;
    int step;
    MPI_Status status;
    unsigned int M = m2;

    /* time measure */
    double starttime, endtime, sort_start_time, sort_end_time,
            b_start_time, b_end_time, d_start_time, d_end_time,
            odd_start_time, odd_end_time, s_start_time, s_end_time,
            n_start_time, n_end_time, nv1v2_start_time, nv1v2_end_time;

    if(world_rank == 0) {
        starttime = MPI_Wtime();
    }

    /* load graph data. current: mockup */
    if(world_rank == 0) {
        const char* file_name = "output.txt";
        load_vertices(file_name, &v1, &v2, m);
        std::cout<< "Zaladowano plik " << file_name << " z " << m << " krawedziami." << std::endl;
        msize = m+1;
        m2 = 2*m;
        m2size = m2 + 1;
        MPI_Bcast(&m,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        MPI_Bcast(&msize,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        MPI_Bcast(&m2,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        MPI_Bcast(&m2size,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    } else {
        MPI_Bcast(&m,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        MPI_Bcast(&msize,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        MPI_Bcast(&m2,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        MPI_Bcast(&m2size,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
    }

    /* just in case... */
    MPI_Barrier(MPI_COMM_WORLD);

    
    v1v2 = generatev1v2(world_rank, world_size, v1, v2, msize, m2size);
    

    //SORT

    if(world_rank == 0) {
        if (m2%world_size) {
            std::cout << std::endl << "Liczba elementow ( " << m2 << " )nie jest podzielna przez liczbe procesow. Zmien liczbe procesow i sproboj ponownie." << std::endl;
            exit(0);
        }
        std::cout << "Rozpoczecie sortowania v1v2... \n";
        sort_start_time   = MPI_Wtime();
        sort_chunk_size = m2/world_size;
        MPI_Bcast(&sort_chunk_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        chunk = (unsigned int *)malloc(sort_chunk_size*sizeof(unsigned int));
        sort_v1v2 = (unsigned int *) malloc((m2) * sizeof(unsigned int));
        memcpy(sort_v1v2, v1v2+1, m2 * sizeof(unsigned int));
        MPI_Scatter(sort_v1v2,sort_chunk_size,MPI_UNSIGNED,chunk,sort_chunk_size,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        m_sort(chunk, 0, sort_chunk_size-1);
    } else {
        MPI_Bcast(&sort_chunk_size,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        chunk = (unsigned int *)malloc(sort_chunk_size*sizeof(unsigned int));
        MPI_Scatter(sort_v1v2,sort_chunk_size,MPI_UNSIGNED,chunk,sort_chunk_size,MPI_UNSIGNED,0,MPI_COMM_WORLD);
        m_sort(chunk, 0, sort_chunk_size-1);
    }


    step = 1;
    while(step<world_size)
    {
        if(world_rank%(2*step)==0)
        {
            if(world_rank+step<world_size)
            {
                MPI_Recv(&M,1,MPI_UNSIGNED,world_rank+step,0,MPI_COMM_WORLD,&status);
                other = (unsigned int *)malloc(M*sizeof(unsigned int));
                MPI_Recv(other,M,MPI_UNSIGNED,world_rank+step,0,MPI_COMM_WORLD,&status);
                chunk = merge(chunk,sort_chunk_size,other,M);
                sort_chunk_size = sort_chunk_size+M;
            }
        }
        else
        {
            int near = world_rank-step;
            MPI_Send(&sort_chunk_size,1,MPI_UNSIGNED,near,0,MPI_COMM_WORLD);
            MPI_Send(chunk,sort_chunk_size,MPI_UNSIGNED,near,0,MPI_COMM_WORLD);
            break;
        }
        step = step*2;
    }


    if (world_rank==0) {
        memcpy(v1v2+1, chunk, m2 * sizeof(unsigned int));
        sort_end_time  = MPI_Wtime();
        std::cout << "Posortowano tablice v1v2 w " << sort_end_time-sort_start_time << "\n";

        if(Printing){
            std::cout << "\nSorted v1v2:\n\t";
            for(int i=0;i<m2size;i++){
                std::cout << v1v2[i] << ", ";
            }
        }
    }
    //SORT END


    if (!world_rank) {
        std::cout << "Rozpoczęcie generowania tablicy B... \n";
        b_start_time = MPI_Wtime();
    }
    B = generateBparallel(world_rank, world_size, v1v2, msize, m2size);
    //    B = generateB(world_rank, world_size, v1v2, msize, m2size);
    if (!world_rank) {
        b_end_time = MPI_Wtime();
        std::cout << "Wygenerowano tablice B w "<< b_end_time-b_start_time << "\n";
    }


    unsigned int Dsize;
    if(world_rank == 0){
        Dsize = v1v2[m2size-1]+1;
    }

    MPI_Bcast(&Dsize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (!world_rank){
        std::cout << "Rozpoczęcie generowania tablicy D... \n";
        d_start_time = MPI_Wtime();
    }

    D = generateDparallel(world_rank, world_size, v1v2, B, msize, m2size, Dsize);
    //    D = generateD(world_rank, world_size, v1v2, B, msize, m2size);
    if (!world_rank) {
        d_end_time = MPI_Wtime();
        std::cout << "Wygenerowano tablice D w " << d_end_time-d_start_time << "\n";
    }


    if (!world_rank) {
        std::cout << "Rozpoczęcie generowania odd... \n";
        odd_start_time = MPI_Wtime();
    }
    odd = getOdd(world_rank, world_size, v1v2, D, m2size);
    if (!world_rank) {
        odd_end_time = MPI_Wtime();
        std::cout << "Wygenerowano odd w " << odd_end_time-odd_start_time << "\n";
    }

    if (!world_rank) {
        std::cout << "Rozpoczęcie generowania tablicy S... \n";
        s_start_time = MPI_Wtime();
    }
    S = generateS(world_rank, world_size, v1v2, D, m2size);
    if (!world_rank) {
        s_end_time = MPI_Wtime();
        std::cout << "Wygenerowano tablice S w " << s_end_time-s_start_time << "\n";
    }

    if (!world_rank) {
        std::cout << "Rozpoczęcie generowania tablicy N... \n";
        n_start_time = MPI_Wtime();
    }
    N = generateN(world_rank, world_size, v1v2, D, S, m2size);
    if (!world_rank) {
        n_end_time = MPI_Wtime();
        std::cout << "Wygenerowano tablice N w " << n_end_time-n_start_time << "\n";
    }

    if (!world_rank) {
        std::cout << "Rozpoczęcie generowania tablicy nv1v2... \n";
        nv1v2_start_time = MPI_Wtime();
    }
    nv1v2 = generateNv1v2(world_rank, world_size, v1v2, S, B, odd, m2size);
    if (!world_rank) {
        nv1v2_end_time = MPI_Wtime();
        std::cout << "Wygenerowano tablice nv1v2 w " << nv1v2_end_time-nv1v2_start_time << "\n";
    }
    // nv1v2 = mockupNv1v2(world_rank);

    if(world_rank == 0){
        nv1 = new unsigned int[msize];
        nv2 = new unsigned int[msize];
    }
    if (!world_rank) std::cout << "Rozpoczęcie generowania tablic nv1 i nv2... \n";
    generateNvs(nv1, nv2, world_rank, world_size, v1, v2, v1v2, nv1v2, m2size, m);
    if (!world_rank) std::cout << "Wyenerowano tablie nv1 i nv2. \n";
    if(world_rank == 0 and Printing){
        std::cout << "\nnv1: ";
        for (unsigned int i=1; i<=m; i++){
            std::cout << nv1[i] << ", ";
        }
        std::cout << "\nnv2: ";
        for (unsigned int i=1; i<=m; i++){
            std::cout << nv2[i] << ", ";
        }
    }
     getEuler(world_rank, world_size, odd, m, nv1, nv2, nv1v2, v1v2, m2size, msize);
    //getEuler(world_rank, world_size, 5, m, nv1, nv2, nv1v2, v1v2, m2size, msize);

    if(world_rank == 0){
        endtime   = MPI_Wtime();
        std::cout << "m=" << m;
        std::cout << "\nTime: " << endtime-starttime;
        std::cout << "\n\nGoodbye!\n\n";
//        free(v1);
//        free(v2);
//        free(v1v2);
//        free(sort_v1v2);
//        free(B);
//        free(D);
//        free(nv1v2);
//        free(nv1);
//        free(nv2);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}