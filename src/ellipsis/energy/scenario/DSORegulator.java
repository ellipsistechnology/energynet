package ellipsis.energy.scenario;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Capacitor;
import ellipsis.energy.grid.Line;
import ellipsis.energy.grid.Transformer;

public class DSORegulator implements Regulator
{
    private Set<Transformer> oltcs = new HashSet<Transformer>();
    private Set<Capacitor> capacitors = new HashSet<Capacitor>();
    private double basePower;
    
    @Override
    public void regulate(AnalysisResults results)
    {
        if(basePower <= 0)
            throw new RuntimeException("Invalid base power "+basePower);
        
        for (Transformer txr : oltcs)
        {
            regulate(txr, results);
        }
        
        for (Capacitor c : capacitors)
        {
            regulate(c, results);
        }
    }

    /**
     * Currently, this regulates by approximating the best transformer
     * coil/voltage ratio. However, in a scenario with a shorter time 
     * between iterations it would be better to more realistically 
     * simulate the transformer's behaviour by adjusting step by step
     * when the voltage is too high or too low. However this approach
     * is too slow when time increments are an hour.
     * @param txr
     * @param results
     */
    private void regulate(Transformer txr, AnalysisResults results)
    {
        // Check valid transformer:
        if(txr.getBaseVoltageRatio() <= 0)
            throw new RuntimeException("Invalid base voltage ratio for transformer "+txr.getName()+" of "+txr.getBaseVoltageRatio());
        if(txr.getImpedance().abs() <= 0)
            throw new RuntimeException("Invalid impedance for transformer "+txr.getName()+" of "+txr.getImpedance());
        if(txr.getMinimumT() <= 0 || txr.getMaximumT() <= txr.getMinimumT())
            throw new RuntimeException("Invalid min/max T for transformer "+txr.getName()+" of "+txr.getMinimumT()+" and "+txr.getMaximumT());
        
        // Get transformer parameters:
        double tStep = txr.getTStepSize();
        double t = txr.getT();
        double maximumT = txr.getMaximumT();
        double minimumT = txr.getMinimumT();
        double deadBandStepPercent = txr.getDeadBandStepPercent();
        Bus toBus = txr.getToBus();
        double nextT;
        
        if(txr.isIncrementalTapMode())
        {
            // Calculate transformer dead band:
            double deadBandWidth = deadBandStepPercent*tStep;
            double lowerDeadBandLimit = 1.0 - deadBandWidth;
            double upperDeadBandLimit = 1.0 + deadBandWidth;
            
            // Get the bus voltage:
            String toBusName = toBus.getName();
            Complex toBusVoltage = results.getBusVoltage(toBusName);
            double voltage = toBusVoltage.getReal(); // Regulate active only
            
            // Calculate next t:
            nextT = 
                    voltage < lowerDeadBandLimit ? t+tStep : // Voltage too low
                    voltage > upperDeadBandLimit ? t-tStep : // Voltage too high
         /*lower limit <= voltage <= upper limit*/ t;        // Within dead band
        }
        else
        {
            // Get the bus voltage:
            String toBusName = toBus.getName();
            
            Complex impedance = txr.getImpedance();
            Complex toBusPower = results.getBusPower(toBusName);
            Complex fromBusVoltage = results.getBusVoltage(txr.getFromBus().getName());

            double a = 
                    toBusPower.conjugate()
                    .add(impedance.reciprocal())
                    .divide(fromBusVoltage.divide(impedance))
                    .getReal();
            nextT = a/txr.getBaseVoltageRatio();
            
            // Correct for discrete steps:
            int n = (int) Math.round((nextT - minimumT)/tStep);
            nextT = minimumT + n*tStep;
        }
        
   /*System.out.println(
           "t = "+t+
           ", U1 = "+lowerDeadBandLimit+
           ", U2 = "+upperDeadBandLimit+
           ", voltage(real) = "+voltage+
           ", nextT = "+nextT);*/
                
        // Constrain t to min and max:
        if(nextT > maximumT)
            nextT = maximumT;
        else if(nextT < minimumT)
                nextT = minimumT;
        
        // Update t:
        txr.setT(nextT);
    }

    // NOTE ALSO, should probably include the load so as to correct for that too
    private void regulate(Capacitor c, AnalysisResults results)
    {
        // Calculate inflowing power to parent bus:
        Bus bus = c.getParent();
        Complex powerIn = Complex.ZERO;
        for (Line line : bus.getLines())
        {
            String lineName = line.getName();
            if(line.getFromBus() == bus)
            {
                Complex s = results.getLinePowerReverse(lineName);
                if(s.getReal() > 0)
                    powerIn = powerIn.add(s);
            }
            else
            {
                Complex s = results.getLinePower(lineName);
                if(s.getReal() > 0)
                    powerIn = powerIn.add(s);
            }
        }
        
        // Check power factor:
        double Q = powerIn.getImaginary();
        double P = powerIn.getReal();
        double powerFactor = P/powerIn.abs();
        double targetPowerFactor = c.getTargetPowerFactor();
        if(powerFactor >= targetPowerFactor) // Already good enough
            return;

        // Calculate required capacitor power value to obtain target power factor:
        double targetLagAngle = Math.acos(targetPowerFactor);
        double targetQ = Math.abs(P*Math.tan(targetLagAngle));
        
        double requiredAdditional = basePower*(Q-targetQ);
        double cVar = requiredAdditional+c.getVAR();

        if(cVar > c.getMax())
        	c.setVAR(c.getMax());
        else if(cVar > 0)
            c.setVAR(cVar);
        else
            c.setVAR(0);
    }
    
    
    /* Accessors */

    public double getBasePower()
    {
        return basePower;
    }

    public void setBasePower(double basePower)
    {
        this.basePower = basePower;
    }

    public boolean add(Transformer t)
    {
        return oltcs.add(t);
    }

    public boolean add(Capacitor c)
    {
        return capacitors.add(c);
    }

    public Collection<Transformer> getTransformers()
    {
        return oltcs;
    }

    public Collection<Capacitor> getCapacitors()
    {
        return capacitors;
    }
}
