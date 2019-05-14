#ifndef __Common_Utilities_h__
#define __Common_Utilities_h__

#include <cstdlib>


template <class T, size_t Nx, size_t Ny, size_t Nz>
void Populate_Stiffness_Matrix(T (&K_matrix)[Nx][Ny][Nz][3][3][3][3][3], const bool use_random_number) {
        constexpr int d = 3;
        //std::srand(time(NULL));
         // Set up stiffness matrix
        //#pragma omp parallel for
        for(int i=0;i<Nx;i++)
        for(int j=0;j<Ny;j++)
        for(int k=0;k<Nz;k++)
            for(int di=-1;di<=1;di++)
            for(int dj=-1;dj<=1;dj++)
            for(int dk=-1;dk<=1;dk++)
                for (int v=0; v<d; v++)
                for (int w=0; w<d; w++)
                    if( (i+di < 0) || (i+di >= Nx) || (j+dj < 0) || (j+dj >= Ny) || (k+dk < 0) || (k+dk >= Nz) )
                        K_matrix[i][j][k][di+1][dj+1][dk+1][v][w] = T();
                    else if ( (di==0) && (dj == 0) && (dk == 0) && (v==w))
                        if (!use_random_number)
                            K_matrix[i][j][k][di+1][dj+1][dk+1][v][w] = T(80)*T(2);
                        else
                            K_matrix[i][j][k][di+1][dj+1][dk+1][v][w] = T(82*2)+(T(std::rand())/T(RAND_MAX))*T(10);
                    else
                        if (!use_random_number)
                            K_matrix[i][j][k][di+1][dj+1][dk+1][v][w]= T(-1);
                        else
                            K_matrix[i][j][k][di+1][dj+1][dk+1][v][w] = -(T(std::rand())/T(RAND_MAX));

        // Symmetric K
        for(int i=0;i<Nx;i++)
        for(int j=0;j<Ny;j++)
        for(int k=0;k<Nz;k++)
            for(int di=-1;di<=1;di++)
            for(int dj=-1;dj<=1;dj++)
            for(int dk=-1;dk<=1;dk++)
                for (int v=0;v<d;v++)
                for (int w=0;w<d;w++)
                    if ((i+di>=0)&&(i+di<Nx)&&(j+dj>=0)&&(j+dj<Ny)&&(k+dk>=0)&&(k+dk<Nz)) {
                        K_matrix[i][j][k][di+1][dj+1][dk+1][v][w] +=  K_matrix[i+di][j+dj][k+dk][1-di][1-dj][1-dk][w][v];
                        K_matrix[i+di][j+dj][k+dk][1-di][1-dj][1-dk][w][v] = K_matrix[i][j][k][di+1][dj+1][dk+1][v][w];}
    }

template <size_t Nx, size_t Ny, size_t Nz>
void PrintArray(int (&arr)[Nx][Ny][Nz]) {
  const int offset = 5;

  for(size_t i=0;i<Nx;i++)
    for(size_t j=0;j<Ny;j++) {
      for (size_t v=0;v<(offset-1)*j;v++)
        std::cout<<" ";
      for(size_t k=0;k<Nz;k++) {
        auto data = arr[i][j][k];
        std::cout<<std::setw(offset)<<data;
      }
      std::cout<<std::endl;
    }
}

template <class T, int d, size_t Nx, size_t Ny, size_t Nz>
    void Fill_Variable(T (&x)[Nx][Ny][Nz][d], const T fill_value, const bool use_random_number) {
    //std::srand(time(NULL));
    for(int i=0;i<Nx;i++)
        for(int j=0;j<Ny;j++)
            for(int k=0;k<Nz;k++)
                for(int v=0;v<d;v++){
                    if (use_random_number) {
                        x[i][j][k][v]=(T(std::rand())/T(RAND_MAX));
                    }
                    else {
                        x[i][j][k][v]=fill_value;}}
}

template <class T, int d, size_t Nx, size_t Ny, size_t Nz>
    void Multiply(const T (&K_matrix)[Nx][Ny][Nz][3][3][3][d][d], const T (&x)[Nx][Ny][Nz][d], T (&y)[Nx][Ny][Nz][d]) {
    using variable_type = double (&)[Nx][Ny][Nz][d];
    double * temp_raw = new double[Nx*Ny*Nz*d]();
    variable_type temp = reinterpret_cast<variable_type>(*temp_raw);
    // y=K*x
#pragma omp parallel for collapse(4)
    for(int i=0;i<Nx;i++)
    for(int j=0;j<Ny;j++)
    for(int k=0;k<Nz;k++)
        for(int v=0;v<d;v++) {
//            y[i][j][k][v] = T();
            for(int di=-1;di<=1;di++)
            for(int dj=-1;dj<=1;dj++)
            for(int dk=-1;dk<=1;dk++)
                if( (i+di >= 0) && (i+di < Nx) && (j+dj >= 0) && (j+dj < Ny) && (k+dk >= 0) && (k+dk < Nz) )
                    for(int w=0;w<d;w++)
                        temp[i][j][k][v] += (double)K_matrix[i][j][k][di+1][dj+1][dk+1][v][w] * x[i+di][j+dj][k+dk][w];}
#pragma omp parallel for collapse(4)
    for(int i=0;i<Nx;i++)
    for(int j=0;j<Ny;j++)
    for(int k=0;k<Nz;k++)
        for(int v=0;v<d;v++)
            y[i][j][k][v] = (T)temp[i][j][k][v];
}

#endif // __Common_Utilities_h__
