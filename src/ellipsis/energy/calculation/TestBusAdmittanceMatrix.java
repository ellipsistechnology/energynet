package ellipsis.energy.calculation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.ArrayFieldVector;
import org.apache.commons.math3.linear.FieldDecompositionSolver;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldVector;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Line;
import ellipsis.energy.grid.SlackSource;

/**
 * Ref. Hadi Saadat, Power System Analysis Third Edition, Chapter 6, PSA Publishing, 2010.
 * @author bmillar
 *
 */
@Deprecated
public class TestBusAdmittanceMatrix
{
    public static void main(String[] args)
    {
        Grid grid = Grid.grid().

            Bus("Bus 1").
                SlackSource("HV1", 1.1, 0, 0, 1).
                Line("Line 1-2", 1, 0, 0.4).
                    Bus("Bus 2").
                        SlackSource("HV2", 1, 0, 0, 0.8).
                        Line("Line 2-3", "Bus 3", 1, 0, 0.2).
                    terminate(). // Bus 2
                terminate(). // Line 1-2
                Line("Line 1-3", 1, 0, 0.2).
                    Bus("Bus 3").
                        Line("Line 3-4", 1, 0, 0.08).
                            Bus("Bus 4").terminate().
                        terminate(). // Line 3-4
                    terminate(). // Bus 3
                terminate(). // Line 1-3
            terminate(). // Bus 1
                
        grid();
        
        // Initialise impedance and admittance vectors:
        List<Line> lines = new ArrayList<Line>(grid.getLines());
        List<Bus> busses = new ArrayList<Bus>(grid.getBusses());
        List<SlackSource> hvSources = new ArrayList<SlackSource>(grid.get(SlackSource.class));
        int nBusses = busses.size();
        int nLines = lines.size();
        //Complex[] Z_array = new Complex[lines.size()+hvSources.size()];
        List<Complex> y = new ArrayList<Complex>(lines.size()+hvSources.size());
        for (int i = 0; i < nLines; ++i)
        {
            Line line = lines.get(i);
            double length = line.getLength();
            Complex z = new Complex(line.getResistancePerMetre()*length, line.getInductancePerMetre()*length);
            //Z_array[i] = z;
            y.add(z.reciprocal());
        }
        for(int i = 0; i < hvSources.size(); ++i)
        {
            SlackSource hvSource = hvSources.get(i);
            Complex z = new Complex(hvSource.getResistance(), hvSource.getInductance());
            //Z_array[nLines+i] = z;
            y.add(z.reciprocal());
        }
        
        // Assign numbers to busses:
        Map<String, Integer> busNumbers = new HashMap<String, Integer>();
        for (int i = 0; i < nBusses; ++i)
        {
            Bus bus = busses.get(i);
            busNumbers.put(bus.getName(), i);
        }
        
        // Initialise Y matrix:
        Complex[][] Y_array = new Complex[nBusses][nBusses];
        for (int i = 0; i < Y_array.length; i++)
        {
            for (int j = 0; j < Y_array[i].length; j++)
            {
                Y_array[i][j] = Complex.ZERO;
            }
        }
        Array2DRowFieldMatrix<Complex> Y = new Array2DRowFieldMatrix<Complex>(Y_array, /*copyArray=*/false);
        
        // Calculate Y matrix off-diagonal elements from line data:
        for (int i = 0; i < lines.size(); i++)
        {
            Line line = lines.get(i);
            int fromNumber = busNumbers.get(line.getFromBus().getName());
            int toNumber = busNumbers.get(line.getToBus().getName());
            
            Complex value = Y.getEntry(fromNumber, toNumber).subtract(y.get(i));
            Y.setEntry(fromNumber, toNumber, value);
            Y.setEntry(toNumber, fromNumber, value);
        }
        
        // Calculate diagonal elements of Y:
        for (int i = 0; i < busses.size(); i++)
        {
            // Consider all lines attached to this bus:
            for (int j = 0; j < lines.size(); j++)
            {
                Line line = lines.get(j);
                int fromNumber = busNumbers.get(line.getFromBus().getName());
                int toNumber = busNumbers.get(line.getToBus().getName());
                if(fromNumber == i || toNumber == i)
                {
                    Complex value = Y.getEntry(i, i).add(y.get(j));
                    Y.setEntry(i, i, value);
                }
            }
            
            // Consider all HV Sources attached to this bus:
            for (int j = 0; j < hvSources.size(); j++)
            {
                SlackSource hvSource = hvSources.get(j);
                String parentBusName = hvSource.getParent().getName();
                if(i == busNumbers.get(parentBusName))
                {
                    Complex value = Y.getEntry(i, i).add(y.get(nLines+j));
                    Y.setEntry(i, i, value);
                }
            }
        }
        
        // Make current vector:
        Complex[] I_array = new Complex[nBusses];
        for(int i = 0; i < nBusses; ++i)
        {
            I_array[i] = Complex.ZERO;
        }
        for(int i = 0; i < hvSources.size(); ++i)
        {
            SlackSource hvSource = hvSources.get(i);
            int busNumber = busNumbers.get(hvSource.getParent().getName());
            Complex v = hvSource.getVoltage();
            Complex z = new Complex(hvSource.getResistance(), hvSource.getInductance());
            I_array[busNumber] = v.divide(z);
        }
        ArrayFieldVector<Complex> I = new ArrayFieldVector<Complex>(I_array, /*copyArray=*/false);
        
        // Calculate voltages:
        //FieldMatrix<Complex> Z = new FieldLUDecomposition<Complex>(Y).getSolver().getInverse();
        //FieldVector<Complex> V = Z.operate(I);
        FieldDecompositionSolver<Complex> solver = new FieldLUDecomposition<Complex>(Y).getSolver();
        FieldVector<Complex> V = solver.solve(I);
        
        // Print results:
        System.out.println("\nNumber|Name|Current|Voltage");
        for (int i = 0; i < nBusses; ++i)
        {
            Bus bus = busses.get(i);
            String name = bus.getName();
            System.out.println(i+"|"+name+"|"+I.getEntry(i)+"|"+V.getEntry(i));
        }
        
        System.out.println("\nY:");
        for(int i = 0; i < nBusses; ++i)
        {
            for(int j = 0; j < nBusses; ++j)
            {
                System.out.print(Y.getEntry(i, j)+"\t");
            }
            System.out.println();
        }
        
        
        /* Calculate voltages allowing for DGs and loads: */
        
        // Set up loads & DGs:
        Complex loads[] = new Complex[nBusses];
        Complex outputs[] = new Complex[nBusses];
System.out.println("\nLoads and Sources:");
        for(int i = 0; i < nBusses; ++i)
        {
            Bus bus = busses.get(i);
            loads[i] = bus.getLoadPower();
            outputs[i] = bus.getGeneratedPower();
System.out.println(bus.getName()+": load="+loads[i]+", output="+outputs[i]);
        }
        
        // Find voltages:
System.out.println("\nCalculating voltages:");
        int iterations = 20;
        double targetError = 0.00001; // NOTE base error on real and imaginary power mismatch
        FieldVector<Complex> Vk = new ArrayFieldVector<Complex>(nBusses, Complex.ONE);
        boolean hasError = true;
        for(int k = 0; k < iterations && hasError; ++k)
        {
            hasError = false;
            Complex vSum = Complex.ZERO;
            for (int i = 0; i < nBusses; ++i) // NOTE this is incorrect as y are the admittances of ALL lines, not just for that particular bus; see LoadFlowAnalyser
            {
                vSum = vSum.add(y.get(i).multiply(Vk.getEntry(i)));
            }
            
            FieldVector<Complex> newV = new ArrayFieldVector<Complex>(nBusses, Complex.ONE);
            for (int i = 0; i < nBusses; ++i)
            {
                // V_i(k+1) = {S_i/V_i*(k) + sum(y_ij*Vj(k), i != j)}/sum(y_ij)
                Complex v_i = Vk.getEntry(i);
                Complex vk = outputs[i].divide(v_i.conjugate()).add(vSum).subtract(v_i).divide(Y.getEntry(i, i));
                newV.setEntry(i, vk);
                double error = Math.abs(v_i.subtract(vk).abs());
                if(error > targetError)
                    hasError = true;
System.out.println(busses.get(i).getName()+"'s error = "+error);
            }
            Vk = newV;
            
System.out.println(k+": hasError="+hasError+","+
		" Bus 1:"+Vk.getEntry(busNumbers.get("Bus 1"))+
        " Bus 2:"+Vk.getEntry(busNumbers.get("Bus 2"))+
        " Bus 3:"+Vk.getEntry(busNumbers.get("Bus 3"))+
        " Bus 4:"+Vk.getEntry(busNumbers.get("Bus 4"))
        );

        }
        
        // Print results:
        System.out.println("\nVoltage:");
        if(hasError)
            System.out.println("Did not converge!");
        for (int i = 0; i < nBusses; ++i)
        {
            Bus bus = busses.get(i);
            System.out.println(bus.getName()+": "+Vk.getEntry(i));
        }
    }
}
