package ellipsis.energy.calculation;

import java.util.Collection;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.FieldMatrix;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Line;
import ellipsis.energy.grid.ShuntAdmittance;

public class AdmittanceMatrix
{
    public FieldMatrix<Complex> Y;
    
    AdmittanceMatrix(FieldMatrix<Complex> yMatrix)
    {
        Y = yMatrix;
    }
    
    public AdmittanceMatrix(Collection<Bus> busses, Map<String, Integer> busNumbers, Collection<Line> lines, double baseImpedance)
    {
        int nBusses = busses.size();

        // Get admittances:
        Complex[] lineAdmittances = new Complex[lines.size()];
        int i = 0;
        for (Line line : lines)
        {
            double length = line.getLength();
            double resistance = line.getResistancePerMetre()*length;
            double inductance = line.getInductancePerMetre()*length;
            lineAdmittances[i] = new Complex(resistance/baseImpedance, inductance/baseImpedance).reciprocal();
            
            ++i;
        }
        
        // Initialise Y matrix:
        Complex[][] Y_array = new Complex[nBusses][nBusses];
        for (i = 0; i < Y_array.length; i++)
        {
            for (int j = 0; j < Y_array[i].length; j++)
            {
                Y_array[i][j] = Complex.ZERO;
            }
        }
        Y = new Array2DRowFieldMatrix<Complex>(Y_array, /*copyArray=*/false); // TODO delay creation of matrix and use array for longer - more efficient
     
        // Calculate Y matrix off-diagonal elements from line data:
        i = 0;
        for (Line line : lines)
        {
            int fromNumber = busNumbers.get(line.getFromBus().getName());
            int toNumber = busNumbers.get(line.getToBus().getName());
            
            Complex value = Y.getEntry(fromNumber, toNumber).subtract(lineAdmittances[i]);
            Y.setEntry(fromNumber, toNumber, value);
            Y.setEntry(toNumber, fromNumber, value);
            
            ++i;
        }
        
        // Calculate diagonal elements of Y:
        for (Bus bus : busses)
        {
            i = busNumbers.get(bus.getName());

            // Consider all lines attached to this bus:
            int j = 0;
            for (Line line : lines)
            {
                int fromNumber = busNumbers.get(line.getFromBus().getName());
                int toNumber = busNumbers.get(line.getToBus().getName());
                if(fromNumber == i || toNumber == i)
                {
                    Complex value = Y.getEntry(i, i).add(lineAdmittances[j]);
                    Y.setEntry(i, i, value);
                }
                
                ++j;
            }
            
            // Consider all shunt admittances:
            for (ShuntAdmittance shuntAdmittance : bus.getShuntAdmittances())
            {
                Complex value = Y.getEntry(i, i).add(shuntAdmittance.getAdmittance().multiply(baseImpedance));
                Y.setEntry(i, i, value);
            }
        }
    }
    
    public Complex get(int i, int j)
    {
        return Y.getEntry(i, j);
    }

    public FieldMatrix<Complex> getAdmittances()
    {
        return Y;
    }

    public int getRowCount()
    {
        return Y.getRowDimension();
    }
    
    public int getColumnCount()
    {
        return Y.getColumnDimension();
    }
}
