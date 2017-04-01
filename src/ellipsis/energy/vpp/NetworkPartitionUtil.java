package ellipsis.energy.vpp;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import ellipsis.energy.calculation.AdmittanceMatrix;
import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;

public class NetworkPartitionUtil
{
    private double[][] alphas;
    private RealMatrix jacobianVQ;
    private boolean approximate = false;
    private Grid grid;
    private AnalysisResults results;
    private Map<String, Integer> busNumbers;
    private AdmittanceMatrix admittanceMatrix;
    
    /**
     * 
     * @param grid The power grid containing all busses of interest.
     * @param results May be null if approximate is true. Otherwise must contain results of call to {@link LoadFlowAnalyser#analyse()}.
     * @param approximate If true, inverse of admittance matrix's imaginary components are used instead of the jacobian matrix.
     * @param baseImpedance Only required if results are null (used for the construction of the bus admittance matrix).
     */
    public NetworkPartitionUtil(Grid grid, AnalysisResults results, boolean approximate, double baseImpedance)
    {
        this.grid = grid;
        
        this.results = results;
        if(results != null)
        {
            busNumbers = results.getBusNumbers();
            admittanceMatrix = results.getAdmittanceMatrix();
        }
        else
        {
            busNumbers = createBusNumbers(grid);
            admittanceMatrix = new AdmittanceMatrix(grid.getBusses(), busNumbers, grid.getLines(), baseImpedance);
        }
        
        this.approximate = approximate;
        int n = grid.getBusses().size();
        alphas = new double[n][n];
        setupJacobianInverse();
    }
    
    private static Map<String, Integer> createBusNumbers(Grid grid)
    {
        Map<String, Integer> busNumbers = new HashMap<String, Integer>();
        int i = 0;
        for (Bus bus : grid.getBusses())
        {
            busNumbers.put(bus.getName(), i);
            ++i;
        }
        return busNumbers;
    }

    public double distance(Bus bus_i, Bus bus_j)
    {
        return -Math.log10(alpha(bus_i, bus_j)*alpha(bus_j, bus_i));
    }

    public double distance(Bus bus_i, Collection<Bus> busses_j)
    {
        double min = Double.MAX_VALUE;
        double max = 0;
        
        if(busses_j.size() == 1 && busses_j.contains(bus_i))
            return 0;
        
        for (Bus bus_j : busses_j)
        {
            if(bus_i == bus_j)
                continue;
            
            double d = distance(bus_i, bus_j);
            if(d > max)
                max = d;
            if(d < min)
                min = d;
        }
        
        return (min+max)/2;
    }
    
    public double distance(VPP vpp_i, VPP vpp_j) 
    {
        Collection<Bus> busses_i = vpp_i.getBusses();
        Collection<Bus> busses_j = vpp_j.getBusses();
        
        return distance(busses_i, busses_j);
    }

    public double distance(Collection<Bus> busses_i, Collection<Bus> busses_j)
    {
        double min = Double.MAX_VALUE;
        double max = 0;
        
        for (Bus bus_i : busses_i)
        {
            for (Bus bus_j : busses_j)
            {
                double d = distance(bus_i, bus_j);
                if(d > max)
                    max = d;
                if(d < min)
                    min = d;
            }
        }
        
        return (min+max)/2;
    }

    public double radius(Collection<Bus> busses)
    {
        switch (busses.size())
        {
        case 0:
            throw new RuntimeException("VPP of 0 busses cannot have it's radius calculated.");
        case 1:
            return 0;
        case 2:
            Iterator<Bus> iterator = busses.iterator();
            return 0.5*distance(iterator.next(), iterator.next());
        }
        
        // For 3 or more busses find the maximum radius:
        double radius = 0;
        for (Bus bus : busses)
        {
            radius = Math.max(radius, distance(bus, busses));
        }
        return radius;
    }
    
    private double[][] jacobianMatrixQV()
    {
        Collection<Bus> busses = grid.getBusses();
        int n = busses.size();
        double[][] QV = new double[n][n];
        
        // Fill in diagonal elements (del Q_i/del V_i):
        for (Bus bus_i : busses)
        {
            String name_i = bus_i.getName();
            int i = busNumbers.get(name_i);
            Complex Y_ii = admittanceMatrix.get(i, i);
            
            if(approximate)
            {
                QV[i][i] = -/*V_i.abs()**/Y_ii.getImaginary();
            }
            else
            {
                Complex V_i = results.getBusVoltage(name_i);
                double sum_j = 0;
                for (Bus bus_j : busses)
                {
                    String name_j = bus_j.getName();
                    int j = busNumbers.get(name_j);
                    if(j == i)
                        continue;
                    
                    Complex V_j = results.getBusVoltage(name_j);
                    Complex Y_ij = admittanceMatrix.get(i, j);
                    
                    sum_j += V_j.abs()*Y_ij.abs()*Math.sin(Y_ij.getArgument() - V_i.getArgument() + V_j.getArgument());
                }
                
                QV[i][i] = -2*V_i.abs()*Y_ii.abs()*Math.sin(Y_ii.getArgument()) - sum_j;
            }
        }
        
        // Fill in off-diagonal elements:
        for (Bus bus_i : busses)
        {
            String name_i = bus_i.getName();
            int i = busNumbers.get(name_i);
            
            for (Bus bus_j : busses)
            {
                String name_j = bus_j.getName();
                int j = busNumbers.get(name_j);
                
                if(i == j)
                    continue;
                
                Complex Y_ij = admittanceMatrix.get(i, j);
                
                if(approximate)
                {
                    QV[i][j] = -/*V_i.abs()**/Y_ij.getImaginary();
                }
                else
                {
                    Complex V_i = results.getBusVoltage(name_i);
                    Complex V_j = results.getBusVoltage(name_j);
                    QV[i][j] = -V_i.abs()*Y_ij.abs()*Math.sin(Y_ij.getArgument() - V_i.getArgument() + V_j.getArgument());
                }
            }
        }
        
        // Construct and return matrix:
        return QV;
    }

    private void setupJacobianInverse()
    {
        // Create dQ/dV matrix:
        double[][] QV = jacobianMatrixQV();
        
        // Remove slack bus:
        Collection<Bus> busses = grid.getBusses();
        int slackIndex = -1;
        for (Bus bus : busses)
        {
            if(!bus.getSlackVoltage().equals(Complex.ZERO))
            {
                slackIndex = busNumbers.get(bus.getName());
                break;
            }
        }
        if(slackIndex == -1)
            throw new RuntimeException("Grid does not have a slack bus.");
        double[][] QV2 = new double[QV.length-1][QV.length-1];
        int iOffset = 0;
        for (int i = 0; i < QV.length; i++)
        {
            if(i == slackIndex)
            {
                iOffset = -1;
            }
            else
            {
                int jOffset = 0;
                for (int j = 0; j < QV[i].length; j++)
                {
                    if(j == slackIndex)
                        jOffset = -1;
                    else
                        QV2[i+iOffset][j+jOffset] = QV[i][j];
                }
            }
        }
        Array2DRowRealMatrix jacobianMatrixQV = new Array2DRowRealMatrix(QV2, false);
        
        // Invert to get dV/dQ matrix:
        // WARNING: Use of SVD here can lead to issues if the matrix is singular.
        //          It will produce a pseudo-inverse which is quite different to the actual inverse.
        this.jacobianVQ = new LUDecomposition(jacobianMatrixQV).getSolver().getInverse();
        
        // Insert row and column for the slack bus:
        double[][] VQ2 = jacobianVQ.getData();
        double[][] VQ = new double[QV.length][QV.length];
        
        iOffset = 0;
        for (int i = 0; i < VQ.length; i++)
        {
            if(i == slackIndex)
                iOffset = -1;

            int jOffset = 0;
            for (int j = 0; j < VQ[i].length; j++)
            {
                if(i == slackIndex || j == slackIndex)
                {
                    VQ[i][j] = 0;
                    if(j == slackIndex)
                        jOffset = -1;
                }
                else
                {
                    VQ[i][j] = VQ2[i+iOffset][j+jOffset];
                }
            }
        }
        
        this.jacobianVQ = new Array2DRowRealMatrix(VQ, false);
    }
    
    private double alpha(Bus bus_i, Bus bus_j)
    {
        int i = busNumbers.get(bus_i.getName());
        int j = busNumbers.get(bus_j.getName());
        
        // Check cache:
        if(alphas[i][j] > 0)
            return alphas[i][j];

        // Not in cache so calculate and store before returning:
        alphas[i][j] = jacobianVQ.getEntry(i, j) / jacobianVQ.getEntry(j, j);
        
        return alphas[i][j];
    }
}
