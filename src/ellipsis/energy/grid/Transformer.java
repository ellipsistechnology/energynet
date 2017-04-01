package ellipsis.energy.grid;

import org.apache.commons.math3.complex.Complex;

/**
 * @author bmillar
 *
 */
public class Transformer extends Line
{
    // Transformer parameters:
    private double baseVoltageRatio;
    private Complex impedance;
    
    // Equivalent circuit with shunt and series impedances:
    private ShuntAdmittance fromShunt = new ShuntAdmittance() { public Unit copy(Grid grid) { return null; }; };
    private ShuntAdmittance toShunt = new ShuntAdmittance() { public Unit copy(Grid grid) { return null; }; };
    
    // Regulation parameters:
    private double t = 1;
    private double tStepSize;
    private double minimumT;
    private double maximumT;
    private double deadBandStepPercent = 0.75;
    private boolean incrementalTapMode = false;
    
    @Override
    public void wasRemovedFromBus(Bus parent)
    {
        super.wasRemovedFromBus(parent);
        Bus fromBus = getFromBus();
        if(parent.equals(fromBus))
        {
            fromBus.removeChild(fromShunt);
        }
        else
        {
            Bus toBus = getToBus();
            if(parent.equals(toBus))
                toBus.removeChild(toShunt);
//            else
//                throw new RuntimeException("Attempt was made to remove a Transformer from a bus it was not connected to: "+parent);
        }
    }
    
    @Override
    public void setFromBus(Bus fromBus)
    {
        super.setFromBus(fromBus); // note that this will trigger wasRemovedFromBus
        
        if(fromBus != null)
            fromBus.addChild(fromShunt);
    }
    
    @Override
    public void setToBus(Bus toBus)
    {
        super.setToBus(toBus); // note that this will trigger wasRemovedFromBus
        
        if(toBus != null)
            toBus.addChild(toShunt);
    }
    
    @Override
    public void setName(String name)
    {
        super.setName(name);
        fromShunt.setName(name+":FromShunt");
        toShunt.setName(name+":ToShunt");
    }

    /**
     * Equivalent circuit is a shunt admittance on the from and to busses and 
     * a series admittance on the line.
     * 
     *  Yi = Yt(a-1)/a
     *  Yj = Yt(1-a)/a^2
     * Yij = Yt/a
     * 
     * Where Yi, Yj and Yij are the from and to shunt admitances and line admittance respectively,
     * and 1:a is the transformer turn ratio.
     */
    private void updateEquivalentUnits()
    {
        if(impedance == null || baseVoltageRatio == 0)
            return;
        
        double voltageRatio = voltageRatio();
        
        Complex admittance = impedance.reciprocal();
        
        fromShunt.setAdmittance(admittance.multiply((voltageRatio-1)/voltageRatio));
        toShunt.setAdmittance(admittance.multiply((1-voltageRatio)/voltageRatio/voltageRatio));
        
        Complex lineImpedance = impedance.multiply(voltageRatio);
        setResistancePerMetre(lineImpedance.getReal());
        setInductancePerMetre(lineImpedance.getImaginary());
        setLength(1);
    }
    
    public double voltageRatio()
    {
        return baseVoltageRatio*t;
    }

    public void swapEndpoints()
    {
        Bus from = getFromBus();
        Bus to = getToBus();
        
        setFromBus(null);
        setToBus(from);
        setFromBus(to);
    }
    
    
    /* ACCESSORS */

    /**
     * Note that the actual voltage ratio will be the base voltage
     * ratio multiplied by t.
     * @return The ratio of output to input voltage for t = 1.
     * @see {@link #getT()}
     * @see {@link #voltageRatio()}
     */
    public double getBaseVoltageRatio()
    {
        return baseVoltageRatio;
    }

    public void setBaseVoltageRatio(double voltageRatio)
    {
        this.baseVoltageRatio = voltageRatio;
        updateEquivalentUnits();
    }

    public Complex getImpedance()
    {
        return impedance;
    }

    public void setImpedance(Complex impedance)
    {
        this.impedance = impedance;
        updateEquivalentUnits();
    }

    public ShuntAdmittance getFromShunt()
    {
        return fromShunt;
    }

    public ShuntAdmittance getToShunt()
    {
        return toShunt;
    }

    /**
     * The per unit adjustment for the turn ratio.
     * The turn ratio is t*baseTurnRatio.
     * @return The per unit turn ratio adjustment.
     */
    public double getT()
    {
        return t;
    }

    public void setT(double t)
    {
        if(t < minimumT)
            throw new RuntimeException("Invalid t: "+t+", "+minimumT+ " < t < "+maximumT);
        if(t > maximumT)
            throw new RuntimeException("Invalid t: "+t+", "+minimumT+ " < t < "+maximumT);
        this.t = t;
        updateEquivalentUnits();
    }

    public double getTStepSize()
    {
        return tStepSize;
    }

    public void setTStepSize(double tStep)
    {
        this.tStepSize = tStep;
    }

    public double getMinimumT()
    {
        return minimumT;
    }

    public void setMinimumT(double minimumT)
    {
        this.minimumT = minimumT;
    }

    public double getMaximumT()
    {
        return maximumT;
    }

    public void setMaximumT(double maximumT)
    {
        this.maximumT = maximumT;
    }

    public double getDeadBandStepPercent()
    {
        return deadBandStepPercent;
    }

    public void setDeadBandStepPercent(double deadBandStepPercent)
    {
        this.deadBandStepPercent = deadBandStepPercent;
    }

    public boolean isIncrementalTapMode()
    {
        return incrementalTapMode;
    }

    public void setIncrementalTapMode(boolean incrementalTapMode)
    {
        this.incrementalTapMode = incrementalTapMode;
    }
    
    @Override
    public Unit copy(Grid grid) 
    {
    	Transformer copy = new Transformer();
    	copyInto(copy);
    	copy.setBaseVoltageRatio(getBaseVoltageRatio());
    	copy.setImpedance(getImpedance());

    	copy.setDeadBandStepPercent(getDeadBandStepPercent());
    	copy.setMinimumT(getMinimumT());
    	copy.setMaximumT(getMaximumT());
    	copy.setT(getT());
    	copy.setTStepSize(getTStepSize());
    	copy.setIncrementalTapMode(isIncrementalTapMode());
    	
    	return copy;
    }
}
