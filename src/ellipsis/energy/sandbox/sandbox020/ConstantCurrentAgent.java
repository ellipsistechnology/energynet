package ellipsis.energy.sandbox.sandbox020;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import static ellipsis.util.VectorHelper.*;
import static ellipsis.util.MatrixHelper.*;
import static ellipsis.util.Sum.*;
import static ellipsis.util.ArrayHelper.*;

public class ConstantCurrentAgent extends Agent
{
    private double current;

    public double getCurrent()
    {
        return current;
    }

    public void setCurrent(double current)
    {
        this.current = current;
    }
    
    public ConstantCurrentAgent(double current)
    {
        this.current = current;
    }
    
    @Override
    public RealVector g()
    {
        return vector(
            sum(neighbour -> (vplus-neighbour.vplus)*conductance(neighbour), neighbours()) - current,
            sum(neighbour -> -(vminus-neighbour.vminus)*conductance(neighbour), neighbours()) - current
        );
    }

    @Override
    public RealMatrix gGrad(Agent wrt)
    {
        if(wrt == this)
            return matrix(
                array(conductanceSum(), 0.0              ), // dg/dv+
                array(0.0,              -conductanceSum()) // dg/dv-
            );
        else // neighbour
            return matrix(
                array(-conductance(wrt), 0.0             ), // dg/dv+
                wrt.isGrounded() ? array(0.0, 0.0) : array(0.0, conductance(wrt)) // dg/dv-
            );
    }

    @Override
    public double inequalityDelta(Agent wrt, RealVector step)
    {
        return 0;
    }
    
    @Override
    public double h()
    {
        return 0;
    }

    @Override
    public RealVector hGrad(Agent wrt)
    {
        return new ArrayRealVector(STATE_DIMENSION);
    }

    @Override
    public RealVector costGrad(Agent wrt)
    {
        return vector(STATE_DIMENSION, 0.0);
    }
    
    @Override
    public double costDelta(Agent wrt, RealVector step)
    {
        return 0;
    }
    
    @Override
    public RealVector constraintDelta(Agent wrt, RealVector step)
    {
        double vplusStep = step.getEntry(0);
        double vminusStep = step.getEntry(1);
        if(wrt == this)
            return
                vector(
                    vplusStep*conductanceSum(),
                    -vminusStep*conductanceSum()
                );
        else
            return
                vector(
                    -vplusStep*conductance(wrt),
                    vminusStep*conductance(wrt)
                );
    }
    
    @Override
    public String primaryVariable()
    {
        return "current";
    }
    
    @Override
    public double primaryValue()
    {
        return getCurrent();
    }
    
    @Override
    public double resistance()
    {
        return voltage()/getCurrent();
    }
}