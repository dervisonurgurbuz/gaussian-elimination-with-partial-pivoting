/* ----------------------- norm ----------------------- */
/*  Given an array and its length, this function 
    computes the 2-norm of the array.
    
    Input variables:
        x     : pointer to array for which the 2-norm should
                 be computed.
        length: number of entries in x.

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.      */

double norm (double * x, int length) {
    int i, length5;
    double a, sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += x[i] * x[i];
    }
    for(; i < length; i += 5) {
        sum += x[i] * x[i] + x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2]
                           + x[i + 3] * x[i + 3] + x[i + 4] * x[i + 4];
    }

    return sqrt(sum);
}


/* ------------------ subinfnorm_index ----------------- */
/*  Given an array, an index, and the length of the 
    subarray, this function returns the index of the 
    subarray which attains the value of the infity norm
    over that subarray.
    
    Input variables:
        x     : pointer to array for which index should
                 be returned.
        sub   : index of x at which subarray begins.
        length: number of entries in x.

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory. If
    more than one value attains the infinity norm then
    the smallest index for which the infinity norm was
    attained is returned.                                */

int subinfnorm_index (double * x, int sub, int length) {
    int i, index;
    double t, inf;

    inf = abs(x[sub]);
    index = sub;

    for(i = sub + 1; i < length + sub; i++) {
        t = abs(x[i]);
        if(t > inf) {
            inf = t;
            index = i;
        }
    }

    return index;
}


/* ----------------------- col_swap ----------------------- */
/*  Given a matrix, and two indices this function exchanges
    the columns of the matrix with the given indices.
    
    Input variables:
        x: pointer to array whose entries are the 
            columns of matrix A.
        m: dimension of matrix A (m by m).
        i: first index.
        j: second index.

    Features: This implementation has time complexity 
    O(1) and requires O(1) additional memory.               */

void col_swap (double ** x, int i, int j) {
    double * t;
    
    t = x[i];
    x[i] = x[j];
    x[j] = t;
}

/* ----------------------- row_swap ----------------------- */
/*  Given a matrix, its dimension, and two indices this 
    function exchanges the rows of the matrix with the 
    given indices.
    
    Input variables:
        x: pointer to array whose entries are the 
            columns of matrix A.
        m: dimension of matrix A (m by m).
        i: first index.
        j: second index.

    Features: This implementation has time complexity 
    O(m) and requires O(1) additional memory.               */

void row_swap (double ** x, int m, int i, int j) {
    int k;
    double t;
    
    for(k = 0; k < m; k++) {
        t = x[k][i];
        x[k][i] = x[k][j];
        x[k][j] = t;
    }
}

/* ----------------------- vec_copy ----------------------- */
/*  Given two arrays of the same length and their length, 
    this function stores the values from the first array
    in the second array.
    
    Input variables:
        x     : pointer to array whose entries are to be
                 copied.
        y     : pointer to array in which the components
                 of x are to be stored.
        length: number of entries in x and in y.

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.          */

void vec_copy (double * x, double * y, int length) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i];
    }
    for(; i < length; i += 5) {
        y[i] = x[i];
        y[i + 1] = x[i + 1];
        y[i + 2] = x[i + 2];
        y[i + 3] = x[i + 3];
        y[i + 4] = x[i + 4];
    }
}


/* ------------------- partialvec_copy ------------------- */
/*  Given two arrays, the length of the second array, and
    an index this function stores the values from the
    subarray x[index : index + length] in the array
    y[0 : length].
    
    Input variables:
        x     : pointer to array whose entries are to be
                 copied.
        y     : pointer to array in which the components
                 of x are to be stored.
        length: number of entries in y.
        index : starting index of subarray of x to be
                copied to y.

    Example: Suppose x is a pointer to the array 
    {1, 2, 3, 4, 5}, y is a pointer to the array {0, 0, 0}, 
    length = 3, and index = 2. Then after executing
    partialvec_copy(x, y, 3, 2), the array pointed to by 
    y is now {3, 4, 5}.                         

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.         */

void partialvec_copy (double * x, double * y, int length, int index) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i + index];
    }
    for(; i < length; i += 5) {
        y[i] = x[i + index];
        y[i + 1] = x[i + index + 1];
        y[i + 2] = x[i + index + 2];
        y[i + 3] = x[i + index + 3];
        y[i + 4] = x[i + index + 4];
    }
}


/* ----------------------- scalar_div ----------------------- */
/*  Given two arrays of the same length, their length, and a
    scalar value this function divides the values from the 
    first array by the scalar value and stores the computed
    number in the second array.
    
    Input variables:
        x     : pointer to array whose components are to be
                 divided by r and stored in second array, y.
        r     : scalar used in division.
        length: number of entries in x and in y.
        y     : pointer to array in which the components
                 of x are to be stored.


    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.            */

void scalar_div (double * x, double r, int length, double * y) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i]/r;
    }
    for(; i < length; i += 5) {
        y[i] = x[i]/r;
        y[i + 1] = x[i + 1]/r;
        y[i + 2] = x[i + 2]/r;
        y[i + 3] = x[i + 3]/r;
        y[i + 4] = x[i + 4]/r;
    }
}


/* ----------------------- scalar_sub ----------------------- */
/*  Given two arrays of the same length, their length, and a
    scalar value this function multiplies the values from the 
    first array by the scalar value and then subtracts the 
    computed components from the components the second array.
    
    Input variables:
        x     : pointer to array whose components are to be
                 multiplied by r then subtracted from the
                 components of the second array, y.
        r     : scalar used in multiplication.
        length: number of entries in x and in y.
        y     : pointer to array in which the components
                 of x are to be stored.


    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.            */

void scalar_sub (double * x, double r, int length, double * y) {
    int i, length5;

    length5 = length % 5;

    if(r == 1) {
        for(i = 0; i < length5; i++) {
            y[i] -= x[i];
        }
        for(; i < length; i += 5) {
            y[i] -= x[i];
            y[i + 1] -= x[i + 1];
            y[i + 2] -= x[i + 2];
            y[i + 3] -= x[i + 3];
            y[i + 4] -= x[i + 4];
        }
    }

    else{
        for(i = 0; i < length5; i++) {
            y[i] -= r * x[i];
        }
        for(; i < length; i += 5) {
            y[i] -= r * x[i];
            y[i + 1] -= r * x[i + 1];
            y[i + 2] -= r * x[i + 2];
            y[i + 3] -= r * x[i + 3];
            y[i + 4] -= r * x[i + 4];
        }
    }
}


/* --------------------- partialscalar_sub --------------------- */
/*  Given two arrays, the length of the second array, a scalar 
    value, and an index, this function multiplies the values 
    starting at the given index from the first array by the 
    scalar value and then subtracts the computed components from 
    the components the second array.
    
    Input variables:
        x     : pointer to array whose components are to be
                 multiplied by r then subtracted from the
                 components of the second array, y.
        r     : scalar used in multiplication.
        length: number of entries in y.
        index : 
        y     : pointer to array in which the components
                 of x are to be stored.

    Example: Suppose x is a pointer to the array 
    {1, 2, 3, 4, 5}, y is a pointer to the array {0, 0, 0}, 
    length = 3, r = -1, and index = 2. Then after executing
    partialscalar_sub(x, -1, 3, 2, y), the array pointed to 
    by y is now {-3, -4, -5}. 

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.               */

void partialscalar_sub (double * x, double r, int length, 
                                              int index, double * y) 
{
    int i, length5;

    length5 = length % 5;

    if(r == 1) {
        for(i = 0; i < length5; i++) {
            y[i + index] -= x[i];
        }
        for(; i < length; i += 5) {
            y[i + index] -= x[i];
            y[i + index + 1] -= x[i + 1];
            y[i + index + 2] -= x[i + 2];
            y[i + index + 3] -= x[i + 3];
            y[i + index + 4] -= x[i + 4];
        }
    }
    else if(r == -1) {
        for(i = 0; i < length5; i++) {
            y[i + index] += x[i];
        }
        for(; i < length; i += 5) {
            y[i + index] += x[i];
            y[i + index + 1] += x[i + 1];
            y[i + index + 2] += x[i + 2];
            y[i + index + 3] += x[i + 3];
            y[i + index + 4] += x[i + 4];
        }
    }
    else{
        for(i = 0; i < length5; i++) {
            y[i + index] -= r * x[i];
        }
        for(; i < length; i += 5) {
            y[i + index] -= r * x[i];
            y[i + index + 1] -= r * x[i + 1];
            y[i + index + 2] -= r * x[i + 2];
            y[i + index + 3] -= r * x[i + 3];
            y[i + index + 4] -= r * x[i + 4];
        }
    }
}


/* ------------------------- rowsubrow ------------------------- */
/*  Given two arrays, the length of the second array, a scalar 
    value, a row number, and an index, this function multiplies 
    the values starting at the given index from the first array
    by the scalar value and then subtracts the computed 
    components from the components the second array.
    
    Input variables:
        x     : pointer to array whose components are to be
                 multiplied by r then subtracted from the
                 components of the second array, y. The 
                 elements of the array should be the columns
                 of the matrix A.
        r     : scalar used in multiplication.
        length: number of entries in y.
        index : index at which to start subtraction
        row1  : index of first row 
        row2  : index of second row at which to perform 
                 subtraction

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.               */

void rowsubrow (double ** x, double r, int length, int index, 
                                       int row1, int row2) 
{
    int i, length5;

    length5 = length % 5;

    if(r == 1) {
        for(i = 0; i < length5; i++) {
            x[i + index][row1] -= x[i + index][row2];
        }
        for(; i < length; i += 5) {
            x[i + index][row1] -= x[i + index][row2];
            x[i + index + 1][row1] -= x[i + index + 1][row2];
            x[i + index + 2][row1] -= x[i + index + 2][row2];
            x[i + index + 3][row1] -= x[i + index + 3][row2];
            x[i + index + 4][row1] -= x[i + index + 4][row2];
        }
    }

    else{
        for(i = 0; i < length5; i++) {
            x[i + index][row1] -= r * x[i + index][row2];
        }
        for(; i < length; i += 5) {
            x[i + index][row1] -= r * x[i + index][row2];
            x[i + index + 1][row1] -= r * x[i + index + 1][row2];
            x[i + index + 2][row1] -= r * x[i + index + 2][row2];
            x[i + index + 3][row1] -= r * x[i + index + 3][row2];
            x[i + index + 4][row1] -= r * x[i + index + 4][row2];
        }
    }
}


/* ----------------------- matrixrow_sub ----------------------- */
/*  Given two arrays, the length of the second array, a scalar 
    value, a row number, and an index, this function multiplies 
    the values starting at the given index from the first array
    by the scalar value and then subtracts the computed 
    components from the components the second array.
    
    Input variables:
        x     : pointer to array whose components are to be
                 multiplied by r then subtracted from the
                 components of the second array, y. The 
                 elements of the array should be the columns
                 of the matrix A.
        r     : scalar used in multiplication.
        length: number of entries in y.
        index : index at which to start subtraction
        row   : index of row at which to perform subtraction
        y     : pointer to array in which the components
                 of x are to be stored.

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.               */

void matrixrow_sub (double * x, double r, int length, 
                               int index, int row, double ** y) 
{
    int i, length5;

    length5 = length % 5;

    if(r == 1) {
        for(i = 0; i < length5; i++) {
            y[i + index][row] -= x[i];
        }
        for(; i < length; i += 5) {
            y[i + index][row] -= x[i];
            y[i + index + 1][row] -= x[i + 1];
            y[i + index + 2][row] -= x[i + 2];
            y[i + index + 3][row] -= x[i + 3];
            y[i + index + 4][row] -= x[i + 4];
        }
    }

    else{
        for(i = 0; i < length5; i++) {
            y[i + index][row] -= r * x[i];
        }
        for(; i < length; i += 5) {
            y[i + index][row] -= r * x[i];
            y[i + index + 1][row] -= r * x[i + 1];
            y[i + index + 2][row] -= r * x[i + 2];
            y[i + index + 3][row] -= r * x[i + 3];
            y[i + index + 4][row] -= r * x[i + 4];
        }
    }
}


/* --------------------- dot_product --------------------- */
/*  Given two arrays of the same length and their length, 
    this function returns the dot product of the two 
    arrays.
    
    Input variables:
        x     : pointer to first array.
        y     : pointer to second array.
        length: number of entries in x and in y.

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.         */

double dot_product (double * x, double * y, int length) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += x[i] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
                           + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
    }

    return sum;
}


/* ------------------ partialdot_product ------------------ */
/*  Given two arrays of the same length, their length, and
    an index this function returns the dot product of the 
    two subarrays x[index : length] and y[index : length].
    
    Input variables:
        x     : pointer to first array.
        y     : pointer to second array.
        length: number of entries in x and in y.
        index : starting index for subarrays.

    Example: Suppose x is a pointer to the array 
    {1, 2, 3, 4}, y is a pointer to the array {5, 6, 7, 8}, 
    length = 4, and index = 2. Then the value returned by
    executing partialdot_product(x, y, 4, 2) is 53, which
    is computed by
        x[2] * y[2] + x[3] * y[3] = 3 * 7 + 4 * 8
                                  = 21 + 32
                                  = 53.

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.          */

double partialdot_product (double * x, double * y, int length, int index) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = index; i < length5; i++) {
        sum += x[i] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
                           + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
    }

    return sum;
}


/* -------------------- subdot_product -------------------- */
/*  Given two arrays, the length of the second array, and
    an index this function returns the dot product of the 
    two subarrays x[index : index + length] and 
    y[0 : length]. It is necessary that index + length is
    at most the length of the first array.
    
    Input variables:
        x     : pointer to first array.
        y     : pointer to second array.
        length: number of entries in y.
        index : starting index for subarray of x.

    Example: Suppose x is a pointer to the array 
    {1, 2, 3, 4, 5}, y is a pointer to the array 
    {-1, -2, -3}, length = 3, and index = 2. Then the value 
    returned by executing subdot_product(x, y, 3, 2) is 53, 
    which is computed by
            x[2] * y[0] + x[3] * y[1] + x[4] * y[2] 

          =  3   *  -1  +  4   *  -2  +  5   *  -3

          = -    3      -      8      -      15

          = -26.

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.          */

double subdot_product (double * x, double * y, int length, int index) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += x[i + index] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i + index] * y[i] + x[i + index + 1] * y[i + 1] 
                                   + x[i + index + 2] * y[i + 2]
                                   + x[i + index + 3] * y[i + 3]
                                   + x[i + index + 4] * y[i + 4];
    }

    return sum;
}


/* -------------------- submatrow_product -------------------- */
/*  Given two arrays, the length of the second array, an index
    and a row number, this function returns the dot product of
    the two subarrays x[index : index + length][row] and 
    y[0 : length]. It is necessary that index + length is
    at most the length of the first array.
    
    Input variables:
        x     : pointer to first array whose entries should be
                the columns of the matrix A.
        y     : pointer to second array.
        length: number of entries in y.
        index : starting index for subarray of x.
        row   : row index of matrix A

    Features: This implementation has time complexity 
    O(length) and requires O(1) additional memory.          */

double submatrow_product (double ** x, double * y, int length, 
                                       int index, int row) 
{
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += x[i + index][row] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i + index][row] * y[i] + x[i + index + 1][row] * y[i + 1] 
                                        + x[i + index + 2][row] * y[i + 2]
                                        + x[i + index + 3][row] * y[i + 3]
                                        + x[i + index + 4][row] * y[i + 4];
    }

    return sum;
}
