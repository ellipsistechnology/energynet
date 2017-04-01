package ellipsis.energy.sandbox.sandbox020;

import static ellipsis.util.ArrayHelper.array;
import static ellipsis.util.MatrixHelper.matrix;
import static ellipsis.util.Sum.sum;
import static ellipsis.util.VectorHelper.vector;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class ConstantPowerAgent extends Agent
{
    private double power;
    
    public double getPower()
    {
        return power;
    }

    public void setPower(double power)
    {
        this.power = power;
    }
    
    public ConstantPowerAgent(double power)
    {
        this.power = power;
    }

    @Override
    public RealVector g()
    {
        double vDiff = vplus - vminus;
        return vector(
            sum(neighbour -> vDiff*(vplus-neighbour.vplus)*conductance(neighbour), neighbours()) - power,
            sum(neighbour -> -vDiff*(vminus-neighbour.vminus)*conductance(neighbour), neighbours()) - power
        );
    }

    @Override
    public RealMatrix gGrad(Agent wrt)
    {
        if(wrt == this)
            return matrix(
                array(dgdvp() ), // dg/dv+
                array(dgdvm() ) // dg/dv-
            );
        else // neighbour
            return matrix(
                array(-(vplus-vminus)*conductance(wrt), 0.0                            ), // dg/dv+
                wrt.isGrounded() ? array(0.0, 0.0) : array(0.0, (vplus-vminus)*conductance(wrt)) // dg/dv-
            );
    }

    private double[] dgdvm()
    {
        return array(
                -sum(neighbour -> (vplus-neighbour.vplus)*conductance(neighbour), neighbours()),
                (2*vminus-vplus)*conductanceSum() - sum(neighbour -> neighbour.vminus*conductance(neighbour), neighbours())
            );
    }

    protected double[] dgdvp()
    {
        return array(
            (2*vplus-vminus)*conductanceSum() - sum(neighbour -> neighbour.vplus*conductance(neighbour), neighbours()),
            -sum(neighbour -> (vminus-neighbour.vminus)*conductance(neighbour), neighbours())
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
                    sum(neighbour -> (vplus+vplusStep-neighbour.vplus)*conductance(neighbour), neighbours()),
                    -sum(neighbour -> (vminus+vminusStep-neighbour.vminus)*conductance(neighbour), neighbours())
                ).mapMultiply((vplusStep-vminusStep))
                .add(vector(
                    vplusStep,
                    -vminusStep
                ).mapMultiply(conductanceSum()*(vplus-vminus)));
        else
            return
                vector(
                    -vplusStep,
                    vminusStep
                ).mapMultiply(conductance(wrt)*(vplus-vminus));
    }
    
    @Override
    public String primaryVariable()
    {
        return "power";
    }
    
    @Override
    public double primaryValue()
    {
        return getPower();
    }
    
    @Override
    public double resistance()
    {

        double v = voltage();
        return v*v/getPower();
    }
}