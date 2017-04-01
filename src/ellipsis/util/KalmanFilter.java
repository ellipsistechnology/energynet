package ellipsis.util;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;

public class KalmanFilter
{

    private RealVector[] a;
    private RealMatrix H;
    private double forgettingFactor = 1;

    public double[] kalmanFilter(double[] y, ArrayRealVector x)
    {
        // Initialise:
        if(H == null)
        {
            // Setup H matrix:
            int featureCount = x.getDimension();
            double[][] array = new double[featureCount][featureCount];
            H = new Array2DRowRealMatrix(array, false);
            
            // Setup coefficient vector:
            int resultCount = y.length;
            a = new ArrayRealVector[resultCount];
            for (int i = 0; i < resultCount; i++)
            {
                a[i] = new ArrayRealVector(new double[featureCount]);
            };
        }
        
        // Update coefficients:
        H = H.scalarMultiply(forgettingFactor).add(x.outerProduct(x));
        RealMatrix H_inv = new SingularValueDecomposition(H).getSolver().getInverse();
        
        double approx[] = new double[a.length];
        for (int j = 0; j < approx.length; j++)
        {
            approx[j] = evaluateAndUpdate(y[j], x, H_inv, a[j]);
        }
        
        // Adjust forgetting factor:
//        double error = 0;
//        for (int i = 0; i < y.length; i++)
//        {
//            double e = Math.abs(y[i] - approx[i]);
//            error += e/y[i];
//        }
//        error /= y.length;
//        forgettingFactor = (forgettingFactor + Math.exp(-error))/2;
        
        return approx;
    }

    private double evaluateAndUpdate(double y, ArrayRealVector x, RealMatrix H_inv, RealVector a)
    {
        double approx = evaluate(x, a);
        double error = y - approx;
        if(error == 0)
            return approx;
        
        RealVector num = x.mapMultiply(error);
        RealVector new_a = a.add(H_inv.operate(num));

//for (int i = 0; i < x.getDimension(); i++)
//{
//	System.out.print(x.getEntry(i));
//	System.out.print(",");
//}
//System.out.println();
        for (int i = 0; i < a.getDimension(); i++)
        {
            a.setEntry(i, new_a.getEntry(i));
//System.out.print(new_a.getEntry(i));
//System.out.print(",");
        }
//System.out.println();

        return approx;
    }
    
    public double[] evaluate(ArrayRealVector x)
    {
        double approx[] = new double[a.length];
        for (int j = 0; j < approx.length; j++)
        {
            approx[j] = evaluate(x, a[j]);
        }
        return approx;
    }

    private double evaluate(ArrayRealVector x, RealVector a)
    {
        return a.dotProduct(x);
    }
}
