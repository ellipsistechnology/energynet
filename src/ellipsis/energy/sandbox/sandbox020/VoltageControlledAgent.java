package ellipsis.energy.sandbox.sandbox020;

import static ellipsis.util.ArrayHelper.array;
import static ellipsis.util.MatrixHelper.matrix;
import static ellipsis.util.Sum.sum;
import static ellipsis.util.VectorHelper.vector;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

public class VoltageControlledAgent extends Agent
{
    private double a = 0;//0.0025; // This gives 80% efficiency at 100W DC output
    private double b = 1;
    private double iMax;
    
    public VoltageControlledAgent()
    {
    }
    
    @Override
    public RealVector g()
    {
        return vector(
            currentPlus() + currentMinus(),
            0.0
        );
    }

    @Override
    public RealMatrix gGrad(Agent wrt)
    {
        if(wrt == this)
            return matrix(
                array(conductanceSum(),  0.0), // dg/dv+
                isGrounded() ? array(0.0, 0.0) : array(conductanceSum(), 0.0) // dg/dv-
            );
        else // neighbour
            return matrix(
                array(-conductance(wrt), 0.0), // dg/dv+
                wrt.isGrounded() ? array(0.0, 0.0) : array(-conductance(wrt), 0.0) // dg/dv-
            );
    }
    
    public double h()
    {
//        if(!enableInequalities)
//            return 0;
        
        return sum(n -> (vplus-n.vplus)*conductance(n), neighbours()) - iMax;
    }

    public RealVector hGrad(Agent wrt)
    {
        if(!enableInequalities)
            return new ArrayRealVector(STATE_DIMENSION);
        
        if(wrt == this)
            return vector(conductanceSum(), 0);
        else
            return vector(-conductance(wrt), 0);
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
                    dcdvminus
                ).mapMultiply(powerTransferGrad);
        }
        else
        {
            double dcdvplus = -voltage()*conductance(wrt);
            validate(dcdvplus, "c grad wrt vplus (neighbour "+wrt.getName()+")");
            return 
                vector(
                    dcdvplus,
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
    
    public double hDelta(Agent wrt, RealVector step)
    {
        if(!enableInequalities)
            return 0;
        
        double deltaVPlus = step.getEntry(0);
        if(wrt == this)
        {
            return deltaVPlus*conductanceSum();
        }
        else
        {
            return -deltaVPlus*conductance(wrt);
        }
    }

    @Override
    public double inequalityDelta(Agent wrt, RealVector step)
    {
        if(!enableInequalities)
            return 0;
        
        double hDelta = hDelta(wrt, step);
        double h = h();
        if((mu + alpha*h > 0) && (mu + alpha*(h + hDelta) > 0))
        {
            return (alpha/2.0)*hDelta*(hDelta + 2*h) + mu*hDelta;
        }
        else if((mu + alpha*h > 0) && (mu + alpha*(h + hDelta) <= 0))
        {
            return (alpha/2)*h*h + mu*h + mu*mu/(2*alpha);
        }
        else if((mu + alpha*h <= 0) && (mu + alpha*(h + hDelta) > 0))
        {
            return -(alpha/2)*h*h - mu*h - mu*mu/(2*alpha);
        }
        else
        {
            return 0;
        }
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