//Linear Algebra by Andrew Fowler

/** Various functions dealing with vectors and matrices. */
class LinearAlgebra {

    /**
     * Returns the magnitude of the vector v (which may be of any length).
     * This is found by adding up the squares of all of the elements of v
     * and taking the square root of the total.
     */

    static double magnitude(double[] v) { //finished

        double sum = 0.0;
        for (int i = 0; i < v.length; i++) {
            sum += v[i]*v[i];
        }

        return Math.sqrt(sum);
    }

    /**
     * Returns the sum of vectors v and w. This is a vector of the same
     * length, each of whose elements is the sum of the corresponding
     * elements in v and w.
     */

    static double[] sum(double[] v, double[] w) { //finished
        int n = w.length;
        double [] add = new double[n];
        for (int i = 0; i < w.length; i++) {
            add[i] = v[i] + w[i];
        }

        return add;
    }

    /**
     * Returns the difference between vectors v and w. This is a vector
     * of the same length, each of whose elements is the difference
     * between the corresponding elements in v and w.
     */
    static double[] difference(double[] v, double[] w) { //finished
        int n = w.length;
        double [] add = new double[n];
        for (int i = 0; i < w.length; i++) {
            add[i] = v[i] - w[i];
        }

        return add;
    }

    /**
     * Returns the element-wise between vectors v and w. This is a vector
     * of the same length, each of whose elements is the product of the
     * corresponding elements in v and w.
     */
    static double[] elementwiseProduct(double[] v, double[] w) { //finished
        int n = w.length;
        double [] add = new double[n];
        for (int i = 0; i < w.length; i++) {
            add[i] = v[i] * w[i];
        }

        return add;
    }

    /**
     * Returns the dot product of vectors v and w. This is the sum of
     * the products of the corresponding elements.
     */
    static double dotProduct(double[] v, double[] w) {
        double sum = 0.0;
        for (int i = 0; i < w.length; i++) {
            sum += v[i] * w[i];
        }
        return sum;
    }

    /**
     * Returns, as an array of two elements, the dimensions of matrix m.
     */
    static int[] dimensions(double[][] m) {
        int[] dim = new int [2];
        dim[0] = m.length;
        dim[1] = m[0].length;
        return dim;
    }

    /**
     * Returns the element-wise sum of matrices m and n.
     */
    static double[][] sum(double[][] m, double[][] n) {
        int[] mdim = new int[2];
        int[] ndim = new int[2];
        int[] matrixSize = new int[2];

        mdim[0] = m.length;
        mdim[1] = m[0].length;
        ndim[0] = n.length;
        ndim[1] = n[0].length;

        for (int i = 0; i < 2; i++) {
            if (mdim[i] >= ndim[i]) {
                matrixSize[i] = mdim[i];
            } else {
                matrixSize[i] = ndim[i];
            }
        }

        double[][] sum = new double[matrixSize[0]][matrixSize[1]];

        for (int i = 0; i < matrixSize[0]; i++) {
            for (int j = 0; j < matrixSize[1]; j++) {
                sum[i][j] = m[i][j] + n[i][j];
            }
        }
        return sum;
    }

    /**
     * Returns the element-wise product of matrices m and n.
     */
        static double[][] elementwiseProduct(double[][] m, double[][] n) {
            int[] mdim = new int[2];
            int[] ndim = new int[2];
            int[] matrixSize = new int[2];

            mdim[0] = m.length;
            mdim[1] = m[0].length;
            ndim[0] = n.length;
            ndim[1] = n[0].length;

            for (int i = 0; i < 2; i++) {
                if (mdim[i] >= ndim[i]) {
                    matrixSize[i] = mdim[i];
                } else {
                    matrixSize[i] = ndim[i];
                }
            }

            double [][] product = new double[matrixSize[0]][matrixSize[1]];

            for (int i =0; i < matrixSize[0]; i++){
                for(int j =0; j < matrixSize[1]; j++){
                    product[i][j] = m[i][j] * n[i][j];
                }
            }
            return product;
        }

    /**
     * Returns the transpose of m, that is, a matrix where element
     * i, j is element j, i from m.
     */
    static double[][] transpose(double[][] m) {
        int[] dim = new int[2];
        dim[0] = m.length;
        dim[1] = m[0].length;
        double [][] transpose = new double[dim[1]][dim[0]];

        for (int i = 0; i < dim[1]; i ++) {
            for (int j = 0; j < dim[0]; j++) {
                transpose[i][j] = m[j][i];
            }
        }
        return transpose;
    }

    /**
     * Returns the matrix product of m and n. (Search the web for a
     * definition.)
     */
    static double[][] product(double[][] m, double[][] n) {
        int[] mdim = new int[2];
        int[] ndim = new int[2];
        int[] matrixSize = new int[2];

        mdim[0] = m[0].length;
        mdim[1] = m.length;
        ndim[0] = n[0].length;
        ndim[1] = n.length;
        matrixSize[0] = ndim[0];
        matrixSize[1] = mdim[1];

        double[][] product = new double[matrixSize[1]][matrixSize[0]];

        //initialization
        for (int i = 0; i < matrixSize[1]; i++) {
            for (int j = 0; j < matrixSize[0]; j++) {
                product[i][j] = 0;
            }
        }

        for (int j = 0; j < matrixSize[1]; j++) {
            for (int i = 0; i < matrixSize[0]; i++) {
                double sum = 0;
                for (int k =0; k < mdim[0]; k++) {
                    sum = sum + (m[j][k] * n[k][i]);
                    product[j][i] = sum;
                }
            }
        }
        return product;
    }

}
