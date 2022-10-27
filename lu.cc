#include <chrono>

#include <cstdio>

#include "matrix.h"

//#include "linearalgebra.h"

#include <iostream>

#include <vector>

using namespace std;


/* Global variables holding the matrix data. To complete this assignment
 * you are requested to only use arrays and access these arrays with
 * subscripts. Do not use pointers.
 */

const int max_n_elements = 131072;
const int max_n_rows = 16384;
const int test_vector_count = 5;

const int n = 765;

static double values[max_n_elements];

static int col_ind[max_n_elements];
static int row_ptr_begin[max_n_rows];
static int row_ptr_end[max_n_rows];




template<typename T, size_t n>
void print_array(T const(& arr)[n])
{
    for (size_t i = 0; i < n; i++) {
        std::cout << arr[i] << ' ';
    }
}


int main(int argc, char **argv)
{
    void partialPivoting(float [][n],float [][5], int e[n]);
    void lu(float[][n], float[][n], float[][n], int n);
    void output(float[][n], int);

    cout<< "Start"<< endl;
  if (argc != 2)
    {
      fprintf(stderr, "usage: %s <filename>\n", argv[0]);
      return -1;
    }

  int nnz, n_rows, n_cols;
  bool ok(false);

  ok = load_matrix_market(argv[1], max_n_elements, max_n_rows,
                          nnz, n_rows, n_cols,
                          values, col_ind, row_ptr_begin, row_ptr_end);
  if (!ok)
    {
      fprintf(stderr, "failed to load matrix.\n");
      return -1;
    }

  // For debugging, can be removed when implementation is finished.
 // dump_nonzeros(n_rows, values, col_ind, row_ptr_begin, row_ptr_end);




    //////////////////////////////////////////////

  auto factorization_start_time = std::chrono::high_resolution_clock::now();

  // Perform LU factorization here

    float l[n][n], u[n][n];

    int i = 0, j = 0, k = 0;

    // Printing the matrix before LU Factorization after partial pivoting
    float x[n][5]= {};
    int e[n] = {};
    //vector<float> x(n, 0.1);
    //Making alternating vector
    // Initialising 5 vectors for further implementation because x values also need to shift during pivoting

    for (i = 0; i < n; i++)
    {
        x[i][0] = 1;
        x[i][1] = 0.1;

        // e matrix is for specifying index of b values
        e[i]= i+1;
        if(i%2 == 0){
            x[i][2] = 1;
            x[i][3] = 5;
            x[i][4] = 100;
        }else{
            x[i][2] = -1;
            x[i][3] = -5;
            x[i][4] = -100;
        }
    }

    cout<<"x:  "<<endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < 5; j++)
        {
            printf("%f ", x[i][j]);
        }
        cout << "\n";
    }

    float a[n][n] = {};

    for( k = 0; k< n;k++){
        for (int idx = row_ptr_begin[k]; idx <= row_ptr_end[k]; ++idx){
            if(col_ind[idx]>=0 && (values[idx]<=-0.005||values[idx]>0.005)){
                //cout<<values[idx];
                a[k][col_ind[idx]] = values[idx];
            }else{
                a[k][col_ind[idx]] = 0;
            }
        }
    }
    cout<< endl;


    cout<< "Matrix A"<< endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%f ", a[i][j]);
        }
        cout << "\n";
    }
    partialPivoting(a, x, e);
    cout<< "A after partial pivoting"<<endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%f ", a[i][j]);
        }
        cout << "\n";
    }
    cout<<"x after partial pivoting:  "<<endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < 5; j++)
        {
            printf("%f ", x[i][j]);
        }
        cout << "\n";
    }

    /*for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%f ", b[i][j]);
        }
        cout << "\n";
    }*/
    //

    lu(a, l, u, n);
    cout << "\nL Decomposition\n\n";
    output(l, n);
    cout << "\nU Decomposition\n\n";
    output(u, n);

    auto factorization_end_time = std::chrono::high_resolution_clock::now();

    auto solve_start_time = std::chrono::high_resolution_clock::now();

    // Compute all 5 solution vectors here
    // Ax = L(U * x) = b
    // U * x = c  matrix
    float c[n][5] = {};

    for (int m = 0; m < 5 ; ++m) {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                c[i][m]= c[i][m] + u[i][j] * x[j][m];
            }

        }
    }


    // L * c = b
    float b[n][5]= {};
    for (int m = 0; m < 5 ; ++m) {
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            b[i][m]= b[i][m] + l[i][j] * c[j][m];
        }
    }
    }

    cout<< "Result: " <<endl;
    for (i = 0; i < n; i++) {
        for (j = 0; j < 5; j++){

            cout << "b"<<e[i]<<": "<<b[i][j]<<", ";
    }
        cout<<endl;


    }


  auto solve_end_time = std::chrono::high_resolution_clock::now();
  
  
  double relative_errors[test_vector_count] = {0};
  
  // Compute relative errors here
  std::chrono::duration<double> factorization_elapsed_time = factorization_end_time - factorization_start_time;
  std::chrono::duration<double> solve_elapsed_time = solve_end_time - solve_start_time;
  
  
  // Print results
  cout<< "\nMeasurement Results 1: " << endl;
  fprintf(stdout, "%.20f\n", factorization_elapsed_time.count());
  cout<< "Measurement 2"<< endl;
  fprintf(stdout, "%.20f\n", solve_elapsed_time.count());
  cout<< "Measurement 3"<< endl;
  for (size_t vector_idx = 0; vector_idx < test_vector_count; ++vector_idx)
    {
      fprintf(stdout, "%.20f\n", relative_errors[vector_idx]);
    }
  return 0;
}


void partialPivoting(float a[][n], float x[][5],int e[n]){
    float c,d;
    int f;
    for(int i=n-1;i>0;i--) {
        if (a[i - 1][0] < a[i][0]){

            for (int j = 0; j < n; j++) {

                //Pivoting
                c = a[i][j];
                a[i][j] = a[i - 1][j];
                a[i - 1][j] = c;


            }
            for (int j = 0; j < 5; j++) {
            // We also need to shift multiplier positions on partial pivoting
            d = x[i][j];
            x[i][j] = x[i - 1][j];
            x[i - 1][j] = d;

            // Also swapping e matrix values for printing the b indexes.
            f = e[i];
            e[i] =  e[i-1];
            e[i-1] = f;

            }
        }
    }

}
void lu(float a[][n], float l[][n], float u[][n], int n)
{
    int i = 0, j = 0, k = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (j < i)
                l[j][i] = 0;
            else
            {
                l[j][i] = a[j][i];
                for (k = 0; k < i; k++)
                {
                    l[j][i] = l[j][i] - l[j][k] * u[k][i];
                }
            }
        }
        for (j = 0; j < n; j++)
        {
            if (j < i)
                u[i][j] = 0;
            else if (j == i)
                u[i][j] = 1;
            else
            {
                u[i][j] = a[i][j] / l[i][i];
                for (k = 0; k < i; k++)
                {
                    u[i][j] = u[i][j] - ((l[i][k] * u[k][j]) / l[i][i]);
                }
            }
        }
    }
}
void output(float x[][n], int n)
{
    int i = 0, j = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%f ", x[i][j]);
        }
        cout << "\n";
    }
}