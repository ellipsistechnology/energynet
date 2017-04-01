package ellipsis.energy.calculation;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class FastDecoupledSolver implements PowerFlowSolver
{
    private boolean[] isSlackSource;
    private Complex[] power;
    private double error;
    private RealMatrix B_inverse;
    private AdmittanceMatrix admittances;
    int slackCount = 0;
    
    // Variables to be stored between iterations:
    private double[] deltas;
    
    public FastDecoupledSolver(Complex[] power, boolean[] isSlackSource, AdmittanceMatrix admittances)
    {
        this.power = power;
        this.isSlackSource = isSlackSource;
        this.admittances = admittances;
        
        // Remove slack busses:
        int rows = admittances.getRowCount();
        int cols = admittances.getColumnCount();
        for (int i = 0; i < rows; i++)
        {
            if(isSlackSource[i])
                ++slackCount;
        }
        
        // Build and invert susceptance matrix:
        double[][] data = new double[rows-slackCount][cols-slackCount];
System.out.println("Susceptance:");
        int iOffset = 0;
        for (int i = 0; i < rows; i++)
        {
            if(isSlackSource[i])
            {
                ++iOffset;
                continue;
            }
            
            int jOffset = 0;
            for (int j = 0; j < cols; j++)
            {
                if(isSlackSource[j])
                {
                    ++jOffset;
                    continue;
                }
                
                data[i-iOffset][j-jOffset] = admittances.get(i, j).getImaginary();
System.out.print(data[i-iOffset][j-jOffset]+"\t");
            }
System.out.println();
        }
        RealMatrix B = new Array2DRowRealMatrix(data, false);
        B_inverse = new QRDecomposition(B).getSolver().getInverse();
        
System.out.println("\nB:");
for (int i = 0; i < data.length; i++)
{
    for (int j = 0; j < data[i].length; j++)
    {
        System.out.print(B_inverse.getEntry(i, j)+"\t");
    }
    System.out.println();
}
    }
    
    @Override
    public Complex[] solve(Complex[] Vk)
    {
        if(deltas == null)
            initDeltas(Vk);
        
        error = 0;
        
        // Build P, Q, Delta P, and Delta Q vectors:
//        RealVector P = new ArrayRealVector();
//        RealVector Q = new ArrayRealVector();
        RealVector DeltaP = new ArrayRealVector();
        RealVector DeltaQ = new ArrayRealVector();
        for (int i = 0; i < Vk.length; i++)
        {
            if(isSlackSource[i])
            {
                continue;
            }
            
            double P_i = 0;
            double Q_i = 0;

            for (int j = 0; j < Vk.length; j++)
            {
                Complex s = Vk[i].multiply(Vk[j]).multiply(admittances.get(i, j));
                double m = s.abs(); // magnitude
                double a = s.getArgument(); // angle
                // TODO The above is wrong: theta_ij is from Y_ij only
                
                // (6.52):
                // P_i = sum[j=1..n](|V_i||V_j||Y_ij|cos(theta_ij-delta_i+delta_j)
                P_i += m * Math.cos(a-deltas[i]+deltas[j]);
                
                // (6.53):
                // Q_i = -sum[j=1..n](|V_i||V_j||Y_ij|sin(theta_ij-delta_i+delta_j)
                Q_i += -m * Math.sin(a-deltas[i]+deltas[j]);
            }
            
            // (6.63):
            // Delta P_i = P_i(sched) - P_i
            double DeltaP_i = power[i].getReal() - P_i;
            
            // (6.64):
            // Delta Q_i = Q_i(sched) - Q_i
            double DeltaQ_i = power[i].getImaginary() - Q_i;
            
            DeltaP.append(DeltaP_i);
            DeltaQ.append(DeltaQ_i);
            
System.out.println("P,Q="+DeltaP_i+","+DeltaQ_i);
        }
        
        // (6.77):
        // Delta delta = -B^-1 Delta P/|V|
        
        // (6.78):
        // Delta |V| = -B^-1 Delta Q/|V|
        
        // Get the maximum error:
        double newError = 0;
        if(newError > error)
            error = newError;
        
        return Vk;
    }

    public void initDeltas(Complex[] Vk)
    {
        deltas = new double[Vk.length];
    }

    @Override
    public double getError()
    {
        return error;
    }
}
