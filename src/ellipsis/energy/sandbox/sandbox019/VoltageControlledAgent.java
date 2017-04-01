package ellipsis.energy.sandbox.sandbox019;

import static ellipsis.util.ArrayHelper.array;
import static ellipsis.util.MatrixHelper.matrix;
import static ellipsis.util.Sum.sum;
import static ellipsis.util.VectorHelper.vector;

import java.util.Random;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class VoltageControlledAgent extends Agent
{
    private static final Random rand = new Random(0);
    
    private double a = 0.0025; // This gives 80% efficiency at 100W DC output
    private double b = 1;
    private double iMax;
    
    
    public VoltageControlledAgent()
    {
        // FIXME do not commit
        z = rand.nextDouble();
    }
    
    @Override
    public RealVector g()
    {
        return vector(
            currentPlus() + currentMinus(),
            0.0
        );
    }

//@Override
//public double h()
//{
//    return 0.0;
//}
//
//@Override
//public RealVector hzGrad(Agent wrt)
//{
//    return vector(0.0, 0.0, 0.0);
//}
//
//@Override
//protected double inequalityDelta(Agent wrt, RealVector step)
//{
//    return 0.0;
//}
    
    @Override
    protected double inequalityDelta(Agent wrt, RealVector step)
    {
        if(wrt == this)
        {
            double stepVPlus = step.getEntry(0);
            double stepZ = step.getEntry(2);
            return stepVPlus*conductanceSum() + stepZ*stepZ + 2*z*stepZ;
        }
        else
        {
            double stepVPlus = step.getEntry(0);
            return -stepVPlus*conductance(wrt);
        }
    }
    
    @Override
    public double h()
    {
        return currentPlus() - iMax;
    }
    
    @Override
    public RealVector hzGrad(Agent wrt)
    {
        if(wrt == this)
            return vector(conductanceSum(), 0.0, 2*z);
        else
            return vector(-conductance(wrt), 0.0, 0.0);
    }

    @Override
    public RealMatrix gGrad(Agent wrt)
    {
        if(wrt == this)
            return matrix(
                array(conductanceSum(),  0.0), // dg/dv+
                isGrounded() ? array(0.0, 0.0) : array(conductanceSum(), 0.0), // dg/dv-
                array(0.0, 0.0)                // dg/dz
            );
        else // neighbour
            return matrix(
                array(-conductance(wrt), 0.0), // dg/dv+
                wrt.isGrounded() ? array(0.0, 0.0) : array(-conductance(wrt), 0.0), // dg/dv-
                array(0.0, 0.0)                // dg/dz
            );
    }
    
    /**
     * The power transfer function gives the input power required to provide the given power output.
     * It is an increasing convex function.
     * @return The gradient of the power transfer function.
     */
    public double powerTransferGrad()
    {
        double power = power();
        validate(power, "power");
        return 2*a*power + b;
    }
    
    protected double powerTransfer(double vplus, double vminus)
    {
        double p = power(vplus, vminus);
        return a*p*p + b*p;
    }

    protected double power(double vplus, double vminus)
    {
        return (vplus-vminus)*sum(neighbour -> (vplus-neighbour.vplus)*conductance(neighbour), neighbours());
    }
    
    protected double power()
    {
        return power(vplus, vminus);
    }
    
    @Override
    public double costLocal()
    {
//        double power = power();
//        return (a*power + b)*power;
        return powerTransfer(vplus, vminus);
    }

    @Override
    public RealVector costGrad(Agent wrt)
    {
        double powerTransferGrad = powerTransferGrad();
        validate(powerTransferGrad, "powerTransferGrad");
        if(wrt == this)
        {
            double dcdvplus = (2*vplus-vminus)*conductanceSum() - sum(neighbour -> neighbour.getVplus()*conductance(neighbour), neighbours());
            validate(dcdvplus, "c grad wrt vplus (self)");
            double dcdvminus;
            if(isGrounded())
                dcdvminus = 0.0;
            else
                dcdvminus = -sum(neighbour -> (vplus - neighbour.getVplus())*conductance(neighbour), neighbours());
            validate(dcdvminus, "c grad wrt vminus (self)");
            return 
                vector(
                    dcdvplus,
                    dcdvminus,
                    0.0
                ).mapMultiply(powerTransferGrad);
        }
        else
        {
            double dcdvplus = -voltage()*conductance(wrt);
            validate(dcdvplus, "c grad wrt vplus (neighbour "+wrt.getName()+")");
            return 
                vector(
                    dcdvplus,
                    0.0,
                    0.0
                ).mapMultiply(powerTransferGrad);
        }
    }
    
    @Override
    public double costDelta(Agent wrt, RealVector step)
    {
        double vplusStep = step.getEntry(0);
        double vminusStep = step.getEntry(1);
        if(wrt == this)
        {
            return powerTransfer(vplus+vplusStep, vminus+vminusStep) - powerTransfer(vplus, vminus);
        }
        else
        {
            double yij = conductance(wrt);
//            return voltage()*vplusStep*yij*(a*voltage()*vplusStep*yij - b - 2*a*power());
            double delta_jP_i = -voltage()*vplusStep*yij;
            return (a*(delta_jP_i + 2*power()) + b)*delta_jP_i;
        }
    }
    
    @Override
    public RealVector constraintDelta(Agent wrt, RealVector step)
    {
        double vplusStep = step.getEntry(0);
        double vminusStep = step.getEntry(1);
        if(wrt == this)
            return 
                vector(
                    (vplusStep+vminusStep)*conductanceSum(),
                    0.0
                );
        else
            return 
                vector(
                    -(vplusStep+vminusStep)*conductance(wrt),
                    0.0
                );
    }
    
    @Override
    public String primaryVariable()
    {
        return "voltage";
    }
    
    @Override
    public double primaryValue()
    {
        return voltage();
    }
    
    @Override
    public double resistance()
    {
        double v = voltage();
        return v*v/power();
    }

    public void setIMax(double maxCurrent)
    {
        this.iMax = maxCurrent;
    }
}