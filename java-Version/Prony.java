package com.chinomars.prony;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import android.util.Log;
import android.widget.CompoundButton;

/**
 * Created by Chino on 4/11/16.
 */
public class Prony {

    /**
     * Vector dot times each other
     */
    public double DotTimes(double[] var1, double[] var2) {
        if (var1.length != var2.length) return 0;

        int len = var1.length;
        double result = 0;
        for (int i = 0; i < len; ++i) {
            result += var1[i] * var2[i];    
        }

        return result;
    }

    /**
     * Calculate the C matrix in Marple Algorithm of AR analysis module
     * i = 0:p; j = 0:p;
     */
    public double[][] CalcCMat(double[] x, int p) {
        double[][] cMat = new double[p][p];
        int N = x.length;
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < p; ++j) {
                double[] vecTmp1 = new double[N-p];
                System.arraycopy(x, p-i-1, vecTmp1, 0, N-p);
                double[] vecTmp2 = new double[N-p];
                System.arraycopy(x, p-j-1, vecTmp2, 0, N-p);
                cMat[i][j] = DotTimes(vecTmp1, vecTmp2);

                System.arraycopy(x, i, vecTmp1, 0, N-p);
                System.arraycopy(x, j, vecTmp2, 0, N-p);

                cMat[i][j] += DotTimes(vecTmp1, vecTmp2);
            }
        }

        return cMat;
    }
    
    /**
     * Marple algorithm to calculate a_i
     */
    public Matrix Marple(double[] x, int p) {
        double[][] cMat = CalcCMat(x, p);
        Matrix C = new Matrix(cMat);

        double[][] cP = new double[p][1];
        int N = x.length;
        for (int i = 0; i < p; ++i) {
            double[] vecTmp1 = new double[N-p];
            System.arraycopy(x, p-i-1, vecTmp1, 0, N-p);
            double[] vecTmp2 = new double[N-p];
            System.arraycopy(x, p, vecTmp2, 0, N-p);

            cP[i][0] = DotTimes(vecTmp1, vecTmp2);

            System.arraycopy(x, i, vecTmp1, 0, N-p);
            System.arraycopy(x, 0, vecTmp2, 0, N-p);

            cP[i][0] += DotTimes(vecTmp1, vecTmp2);

            cP[i][0] *= -1;
        }
        Matrix CP = new Matrix(cP);

        return C.solve(CP); // return a_i Matrix

//        double[][] a = C.solve(CP).getArrayCopy();
//        double[][] res = new double[p+1][1];
//        res[0][0] = 1;
//        try{
//            System.arraycopy(a, 0, res, 1, p);
//        } catch (Exception e) {
//            Log.e("error","dimension error");
//        }
//
//        return res;

    }

    /**
     * solve the poly
     * @param a
     * @return
     */
    public EigenvalueDecomposition NRoots(double[] a) {
        int len = a.length - 1;
        double[][] genA = new double[len][len];
        for (int i = 0; i < len; ++i) {
            genA[0][i] = a[i+1] / a[0];
        }

        for (int i = 0;  i < len; ++i) {
            genA[len-1][i] = 0;
        }

        for (int i = 1; i < len-1; ++i) {
            genA[i][i] = 1;
        }

        Matrix AMat = new Matrix(genA);
        AMat.print(0,0);

        return AMat.eig();

    }

    /**
     * complex pow
     * @param c
     * @param i
     * @return c^i
     */
    public Complex mPow(Complex c, int i) {
        Complex result = null;
        for (int j = 0; j < i; ++j) {
            result = c.times(c);
        }

        return result;
    }









}













