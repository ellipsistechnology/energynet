package ellipsis.energy.sandbox.sandbox019;
import static ellipsis.util.VectorHelper.vector;
import static ellipsis.util.ApproximateComparator.approxEquals;
import static ellipsis.util.Sum.sum;
import static ellipsis.util.Sum.sumV;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.Sum;

public abstract class Agent
{
    //// Members ////
    
    public static final int EQUALITY_CONSTRAINT_COUNT = 2;
    public static final int STATE_DIMENSION = 3; // [v+, v-, z]
    
    protected double alphaMultiplier = 1.0001;
    protected double epsilon = 100;
    protected double epsilonMultiplier = 1;//0.9998;
    protected double epsilonMin = 1e-3;
    
    protected double vminus, vplus;  // Primal state
    protected double z; // Linearisation variable for inequality constraints
    protected double vmin, vmax; // Admissible voltage range
    protected RealVector lambda; // Dual state for equality constraints
    protected double mu; // Dual variable for inequality constraints
    protected double alpha;
    private Map<Agent, Double> conductances = new HashMap<>(); // Constants
    private double ySum = 0; // Sum of neighbouring conductances
    private boolean grounded = false;
    private String name;
    
    protected double lambdaMax = 1e6;
    protected double muMax = 1e6;
    
    public String getName()
    {
        return name;
    }

    public void setName(String name)
    {
        this.name = name;
    }

    public boolean isGrounded()
    {
        return grounded;
    }
    
    public void setGrounded(boolean grounded)
    {
        this.grounded = grounded;
    }
    
    public void addNeighbour(Agent neighbour, double conductance)
    {
        ySum = 0;
        conductances.put(neighbour, conductance);
    }
    
    public double conductance(Agent neighbour)
    {
        return conductances.get(neighbour);
    }
    
    public Iterable<Agent> neighbours()
    {
        return conductances.keySet();
    }
    
    public int neighbourCount()
    {
        return conductances.size();
    }

    public double getVplus()
    {
        return vplus;
    }

    public void setVplus(double vplus)
    {
        this.vplus = vplus;
    }

    public double getVminus()
    {
        return vminus;
    }

    public void setVminus(double vminus)
    {
        this.vminus = vminus;
    }

    public double getZ()
    {
        return z;
    }
    
    public void setZ(double z)
    {
        this.z = z;
    }
    
    public double getVmin()
    {
        return vmin;
    }

    public void setVmin(double vmin)
    {
        this.vmin = vmin;
    }

    public double getVmax()
    {
        return vmax;
    }

    public void setVmax(double vmax)
    {
        this.vmax = vmax;
    }
    
    public RealVector getLambda()
    {
        return lambda;
    }
    
    public double getMu()
    {
        return mu;
    }

    public void setLambda(RealVector lambda)
    {
        this.lambda = lambda;
    }

    public double getAlpha()
    {
        return alpha;
    }

    public void setAlpha(double alpha)
    {
        this.alpha = alpha;
    }
    
    public double voltage()
    {
        return getVplus() - getVminus();
    }

    @Override
    public String toString()
    {
        return name;
    }

    protected double conductanceSum()
    {
        if(ySum == 0)
        {
            ySum = sum(neighbour -> conductance(neighbour), neighbours());
        }
        return ySum;
    }
    
    public double getEpsilon()
    {
        return epsilon;
    }
    
    
    //// Initialisation ////
    
    public Agent()
    {
        lambda = new ArrayRealVector(EQUALITY_CONSTRAINT_COUNT);
        alpha = 1.0;
    }
    
    
    //// Abstract Functions ////
    
    /**
     * Equality constraints for this agent.
     * @return A 2 dimensional vector.
     */
    public abstract RealVector g();
    
    /**
     * Inequality constraints for this agent.
     * @return 
     */
    public abstract double h();
    
    /**
     * Equality constraint gradient: For this agent's constraints 
     * with respect to the given agent's state.
     * @param wrt The agent the gradient is to be with respect to. This is assumed to be a neighbour.
     * @return A 3x2 matrix.
     */
    public abstract RealMatrix gGrad(Agent wrt);
    
    /**
     * Inequality constraint (h(v)+z^2) gradient: For this agent's constraints 
     * with respect to the given agent's state [v+, v-, z].
     * @param wrt The agent the gradient is to be with respect to. This is assumed to be a neighbour.
     * @return A 3 dimensional vector.
     */
    public abstract RealVector hzGrad(Agent wrt);
    
    /**
     * The gradient of the cost function with respect to the state: [v+ v- v].
     * @param wrt The agent whose state the gradient is with respect to.
     * @return The gradient vector.
     */
    public abstract RealVector costGrad(Agent wrt);
    
    /**
     * Calculates the change in this agent's cost due to the change in the given agent's state 
     * specified by step.
     * If the given agent is not equal to this agent then neighbour agents must not be used in
     * the calculation.
     * @param wrt The agent whose state is changed. This is assumed to be a neighbour.
     * @param step The state change.
     * @return The change in cost.
     */
    public abstract double costDelta(Agent wrt, RealVector step);
    
    /**
     * Calculates the change in this agent's lambda/constraint dot product due to the change in 
     * the given agent's state specified by step.
     * If the given agent is not equal to this agent then neighbour agents must not be used in
     * the calculation.
     * @param wrt The agent whose state is changed. This is assumed to be a neighbour.
     * @param step The state change.
     * @return The change in cost.
     */
    public abstract RealVector constraintDelta(Agent wrt, RealVector step);

    /**
     * Calculates the change in this agent's mu*h(v) due to the change in the given agent's
     * state specified by step.
     * If the given agent is not equal to this agent then neighbour agents must not be used in
     * the calculation.
     * @param step
     * @return
     */
    protected abstract double inequalityDelta(Agent wrt, RealVector step);

    /**
     * The primary value controller or specified by this agent.
     * @return
     */
    public abstract double primaryValue();

    /**
     * The name of the primary value controller or specified by this agent.
     * @return
     */
    public abstract String primaryVariable();

    public abstract double resistance();
    
    
    //// Optimisation ////
    
    public void step()
    {
        // Step primal state [v+, v-, z]:
        RealVector grad = vector(Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE);
int i = 0;
double oldVPlus = vplus;
double oldVMinus = vminus;
double oldZ = z;
        while(grad.getNorm() > epsilon)
        {
//            grad = minimizeApproximate(1, 1, 1);
            grad.setEntry(0, minimizeApproximate(1, 0, 0).getEntry(0)); // step v+
            grad.setEntry(1, minimizeApproximate(0, 1, 0).getEntry(1)); // step v-
            grad.setEntry(2, minimizeApproximate(0, 0, 1).getEntry(2)); // step z
//if(i > 1000 && i%1000 == 0)
//    System.out.println(i+",\t"+grad);
++i;
        }

        validate(vplus, "v+");
        validate(vminus, "v-");
        validate(z, "z");
        
        // FIXME uncomment
        // Step dual variabls:
//        RealVector g = g();
//        lambda = constrain(lambda.add(g.mapMultiply(alpha)), lambdaMax);
//        validate(lambda, "lambda");
//        
//        double h = h();
//        mu += constrain(alpha*(h + z*z), muMax);
//        validate(mu, "mu");
        
        // Increase alpha:
        alpha *= alphaMultiplier;
//if(alpha > 40)
//    alpha = 40;
        validate(alpha, "alpha");
        
        // Decrease epsilon:
        epsilon *= epsilonMultiplier;
        validate(epsilon, "epsilon");
        if(epsilon < epsilonMin)
            epsilon = epsilonMin;
    }

    protected double constrain(double d, double max)
    {
        if(d > max)
            return max;
        if(d < -max)
            return -max;
        return d;
    }

    protected RealVector constrain(RealVector v, double max)
    {
        int dimension = v.getDimension();
        for(int i = 0; i < dimension; ++i)
            v.setEntry(i, constrain(v.getEntry(i), max));
        return v;
    }

    private RealVector minimizeApproximate(double vplusFlag, double vminusFlag, double zFlag)
    {
        RealVector grad = lagrangeGrad();
        RealVector gradFull = grad.copy();
        grad.setEntry(0, grad.getEntry(0)*vplusFlag);
        grad.setEntry(1, grad.getEntry(1)*vminusFlag);
        grad.setEntry(2, grad.getEntry(2)*zFlag);
        validate(grad, "Gradient");
        
        double stepSize = backtrack(grad);
        validate(stepSize, "step size");
//System.out.println("\ngrad_v+:\t"+grad.getEntry(0));
        
        RealVector step = grad.mapMultiply(-stepSize);
        if(step.getDimension() != STATE_DIMENSION)
            throw new RuntimeException("wrong dimension for step");
        double lagrangeBefore = lagrange();
        stepState(step);
        
        projectVoltages(gradFull);
        
        // FIXME uncomment
//        double lagrangeAfter = lagrange();
//        if(approxLessThan(lagrangeBefore, lagrangeAfter, -10e-3))
//            throw new RuntimeException(
//                    "Lagrange increased after step: From "+lagrangeBefore+" to "+lagrangeAfter
//                    +", with step "+step);
        
        return gradFull;
    }

    /**
     * Add the given values to the agent's state.
     * @param step [vStepPlus, vStepMinus, zStep]
     */
    public void stepState(RealVector step)
    {
        vplus += step.getEntry(0);  // delta v+
        vminus += step.getEntry(1); // delta v-
        z += step.getEntry(2);      // delta z
    }

    /**
     * If the resultant voltage (v+ - v-) is not
     * within the admissible set defined by vmin and vmax then its value will be projected.
     * @param gradToMask Elements will be set to zero if projection occurs in that dimension.
     */
    protected void projectVoltages(RealVector gradToMask)
    {
        if(isGrounded())
        {
            if(vplus > vmax)
            {
                vplus = vmax;
                gradToMask.setEntry(0, 0.0);
            }
            else if(vplus < vmin)
            {
                vplus = vmin;
                gradToMask.setEntry(0, 0.0);
            }
        }
        else
        {
            double v = voltage();
            double vlimit = 0;
            if(v > vmax)
                vlimit = vmax;
            else if(v < vmin)
                vlimit = vmin;
            if(vlimit != 0)
            {
                double delta = (vminus - vplus + vlimit)/2;
                vplus += delta;
                vminus -= delta;
                gradToMask.setEntry(0, 0.0);
                gradToMask.setEntry(1, 0.0);
            }
        }
    }
    
    protected void validateEquals(double a, double b, double precision, String vars)
    {
        if(!approxEquals(a, b, precision))
            throw new RuntimeException(getName()+": "+vars+" were not equal");
    }
    
    protected void validate(double d, String var)
    {
        if(Double.isNaN(d))
            throw new RuntimeException(getName()+": "+var+" was NaN");
        if(Double.isInfinite(d))
            throw new RuntimeException(getName()+": "+var+" was Infinite");
    }

    protected void validate(RealVector v, String var)
    {
        if(v.isNaN())
            throw new RuntimeException(getName()+": "+var+" was NaN");
        if(v.isInfinite())
            throw new RuntimeException(getName()+": "+var+" was Infinite");
    }

    protected void validate(RealMatrix m, String var)
    {
        int dimension = m.getRowDimension();
        for (int i = 0; i < dimension; ++i)
        {
            RealVector v = m.getRowVector(i);
            if(v.isNaN())
                throw new RuntimeException(getName()+": "+var+" was NaN");
            if(v.isInfinite())
                throw new RuntimeException(getName()+": "+var+" was Infinite");
        }
    }

    /**
     * Lagrange function gradient with respect to this agent's state.
     * @return The gradient vector.
     */
    public RealVector lagrangeGrad()
    {
        // Cost gradient:
        RealVector costGrad = costGrad();
        validate(costGrad, "cost gradient");
        
        // Equality constraint gradients:
        RealVector constraintGrad = constraintGrad();
        
        // Inequality constraints:
        RealVector inequalityGrad = inequalityGrad();

        return costGrad.add(constraintGrad).add(inequalityGrad);
    }

    private RealVector inequalityGrad()
    {
        RealVector hGradSelf = hzGrad(this).mapMultiply(mu + alpha*(h() + z*z));
        RealVector hGradNeighbours = 
            sumV(neighbour ->
                neighbour.hzGrad(this).mapMultiply(neighbour.mu + neighbour.alpha*(neighbour.h() + neighbour.z*neighbour.z)),
                neighbours(), 
                STATE_DIMENSION);

        RealVector inequalityGrad = hGradSelf.add(hGradNeighbours);
        return inequalityGrad;
    }

    private RealVector constraintGrad()
    {
        RealMatrix gGradSelf = gGrad(this);
        validate(gGradSelf, "g gradient (self)");
        RealVector constraintGradSelf = gGradSelf.operate(lambda.add(g().mapMultiply(alpha)));
        validate(constraintGradSelf, "constraint grad (self)");
        RealVector constraintGradNeighbours = 
            sumV(neighbour -> 
                neighbour.gGrad(this).operate(neighbour.lambda.add(neighbour.g().mapMultiply(neighbour.alpha))),
                neighbours(), 
                STATE_DIMENSION
            );
        validate(constraintGradNeighbours, "constraint grad (neighbours)");
        
        RealVector constraintGrad = constraintGradSelf.add(constraintGradNeighbours);
        return constraintGrad;
    }

    protected RealVector costGrad()
    {
        RealVector costGradSelf = costGrad(this);
        validate(costGradSelf, "cost grad (self)");
        RealVector costGradNeighbours = sumV(neighbour -> neighbour.costGrad(this), neighbours(), STATE_DIMENSION);
        validate(costGradNeighbours, "cost grad (neighbours)");
        return costGradSelf.add(costGradNeighbours);
    }

    private static final double DELTA = 0.5;
    private static final double BETA = 0.5;
    /**
     * Approximately finds the best step size for the given gradient assuming
     * gradient decent will be used to minimise the Lagrange function.
     * @param grad
     * @return
     */
    protected double backtrack(RealVector grad)
    {
        double stepSize = 1;
        double d = -DELTA*grad.dotProduct(grad);
        double lagDelta = Double.MAX_VALUE;
        while(lagDelta > stepSize*d || Double.isNaN(lagDelta))
        {
            lagDelta = lagrangeDelta(stepSize, grad);
            stepSize *= BETA;
        }
        return stepSize;
    }

    protected double lagrangeDelta(double stepSize, RealVector grad)
    {
        RealVector step = grad.mapMultiply(-stepSize);
        double costDelta = costDelta(step);
        double constraintDelta = constraintDelta(step);
        double constraintSquaredDelta = constraintSquaredDelta(step);
        double inequalityDelta = inequalityDelta(step);
        double inequalitySquaredDelta = inequalitySquaredDelta(step);
        return costDelta + constraintDelta + constraintSquaredDelta + inequalityDelta  + inequalitySquaredDelta;
    }

    /**
     * Change in cost due to a change in this agent's state.
     * @param step The change in state.
     * @return The change in cost.
     */
    protected double costDelta(RealVector step)
    {
        double costDelta = costDelta(this, step);
        double neighbourCostDelta = sum(neighbour -> neighbour.costDelta(this, step), neighbours());
        return
            costDelta +
            neighbourCostDelta;
    }

    /**
     * Change in the lambda/constraint dot product due to a change in 
     * this agent's state.
     * @param step The change in state.
     * @return The change in the lambda/constraint dot product.
     */
    protected double constraintDelta(RealVector step)
    {
        return
            lambda.dotProduct(constraintDelta(this, step)) +
            sum(neighbour -> neighbour.lambda.dotProduct(neighbour.constraintDelta(this, step)), neighbours());
    }

    /**
     * 
     * @param step
     * @return
     */
    protected double inequalityDelta(RealVector step)
    {
        return 
            mu*inequalityDelta(this, step) +
            sum(neighbour -> neighbour.mu*neighbour.inequalityDelta(this, step), neighbours());
    }

    /**
     * Change in the augmentation term due to a change in 
     * this agent's state.
     * @param step The change in state.
     * @return The change in the augmentation term.
     */
    protected double constraintSquaredDelta(RealVector step)
    {
        return 
            constraintSquaredDelta(this, step) +
            sum(neighbour -> neighbour.constraintSquaredDelta(this, step), neighbours());
    }

    /**
     * Change in this agent's augmentation term due to a change in
     * the given agent's state.
     * @param wrt The agent whose state is changed.
     * @param step The change in state.
     * @return The change in the augmentation term.
     */
    protected double constraintSquaredDelta(Agent wrt, RealVector step)
    {
        RealVector dg = constraintDelta(wrt, step);
        return (alpha/2)*dg.dotProduct(dg.add(g().mapMultiply(2.0)));
    }
    
    /**
     * Change in the inequality term due to a change in
     * this agent's state.
     * @param step The change in state.
     * @return The change in the inequality term.
     */
    protected double inequalitySquaredDelta(RealVector step)
    {
        return 
            inequalitySquaredDelta(this, step) +
            sum(neighbour -> neighbour.inequalitySquaredDelta(this, step), neighbours());
    }

    /**
     * Change in this agent's inequality term due to a change in
     * the given agent's state.
     * @param wrt The agent whose state is changed.
     * @param step The change in state.
     * @return The change in the inequality term.
     */
    protected double inequalitySquaredDelta(Agent wrt, RealVector step)
    {
        double deltaHZ = inequalityDelta(wrt, step);
        return (alpha/2)*deltaHZ*(deltaHZ + 2*(h()+z*z));
    }
    
    
    //// Information ////

    public double currentPlus()
    {
        return sum(n -> (vplus - n.vplus)*conductance(n), neighbours());
    }
    
    public double currentMinus()
    {
        return sum(n -> (vminus - n.vminus)*conductance(n), neighbours());
    }

    public double lagrange()
    {
        return lagrangeLocal() + sum(Agent::lagrangeLocal, neighbours());
    }
    
    public double lagrangeLocal()
    {
        return costLocal() + lagrangeLocal_constraints() + lagrangeLocal_inequalities();
    }
    
    private double lagrangeLocal_inequalities()
    {
        double hz = h() + z*z;
        return mu*hz + (alpha/2)*hz*hz;
    }

    private double lagrangeLocal_constraints()
    {
        RealVector g = g();
        return lambda.dotProduct(g) + (alpha/2)*g.dotProduct(g);
    }

    /**
     * Overwrite in subclass.
     * @return 0 by default
     */
    public double costLocal()
    {
        return 0;
    }
    
    
    //// Debug ////

    @SuppressWarnings("unused")
    private Object zz__showGrad_gSelf = new Object()
    {
        public String toString() 
        {
            return ""+gGrad(Agent.this);
        }
    };
    
    @SuppressWarnings("unused")
    private Object zz__showConstraintGrad = new Object()
    {
        public String toString() 
        {
            return ""+gGrad(Agent.this).operate(lambda.add(g().mapMultiply(alpha)))
                    .add(Sum.sumV(v -> v.gGrad(Agent.this).operate(v.lambda.add(v.g().mapMultiply(v.alpha))), neighbours(), STATE_DIMENSION));
        }
    };

    @SuppressWarnings("unused")
    private Object zz__showGrad_costNeighbours = new Object()
    {
        public String toString() 
        {
            StringBuffer sb = new StringBuffer();
            for (Agent n : neighbours())
            {
                sb.append(n.getName());
                sb.append(":\n");
                sb.append(n.costGrad(Agent.this));
                sb.append('\n');
            }
            return sb.toString();
        }
    };
    
    @SuppressWarnings("unused")
    private Object zz__showGrad_costSelf = new Object()
    {
        public String toString() 
        {
            return ""+costGrad(Agent.this);
        }
    };

    @SuppressWarnings("unused")
    private Object zz__showGrad_gNeighbours = new Object()
    {
        public String toString() 
        {
            StringBuffer sb = new StringBuffer();
            for (Agent n : neighbours())
            {
                sb.append(n.getName());
                sb.append(":\n");
                sb.append(n.gGrad(Agent.this));
                sb.append('\n');
            }
            return sb.toString();
        }
    };

    private String zz__gradEstimate_g(Agent wrt)
    {
        double oldVPlus = wrt.vplus;
        double oldVMinus = wrt.vminus;
        double oldZ = wrt.z;
        double step = 0.0001;
        double[][] grad = new double[3][];
        
        // vplus:
        {
            wrt.vplus -= step;
            RealVector left = this.g();
            wrt.vplus = oldVPlus+step;
            RealVector right = this.g();
            wrt.vplus = oldVPlus;
            grad[0] = right.subtract(left).mapDivide(2*step).toArray();
        }

        // vminus:
        if(!wrt.grounded)
        {
            wrt.vminus -= step;
            RealVector left = this.g();
            wrt.vminus = oldVMinus+step;
            RealVector right = this.g();
            wrt.vminus = oldVMinus;
            grad[1] = right.subtract(left).mapDivide(2*step).toArray();
        }
        else
        {
            grad[1] = new double[2];
        }
        
        // z:
        {
            wrt.z -= step;
            RealVector left = this.g();
            wrt.z = oldZ+step;
            RealVector right = this.g();
            wrt.z = oldZ;
            grad[2] = right.subtract(left).mapDivide(2*step).toArray();
        }
        
        return ""+new Array2DRowRealMatrix(grad);
    };
    
    @SuppressWarnings("unused")
    private Object zz__showGradEstimate_g = new Object()
    {
        public String toString() 
        {
            Agent wrt = Agent.this;
            
            return zz__gradEstimate_g(wrt);
        }
    };
    
    @SuppressWarnings("unused")
    private Object zz__showGradEstimate_gNeighbours = new Object()
    {
        public String toString() 
        {
            StringBuffer sb = new StringBuffer();
            for (Agent n : neighbours())
            {
                sb.append(n.getName());
                sb.append(":\n");
                
                sb.append(n.zz__gradEstimate_g(Agent.this));
                
                sb.append('\n');
            }
            return sb.toString();
        }
    };
    
    @SuppressWarnings("unused")
    private Object zz__showGradEstimate_constraint = new Object()
    {
        public String toString() 
        {
            Agent wrt = Agent.this;
            
            return zz__gradEstimate_constraint(wrt);
        }
    };
    
    private String zz__gradEstimate_constraint(Agent wrt)
    {
        double oldVPlus = wrt.vplus;
        double oldVMinus = wrt.vminus;
        double oldZ = wrt.z;
        double step = 0.0001;
        double[] grad = new double[3];
        
        // vplus:
        {
            wrt.vplus -= step;
            double left = zz__lagrange_constraints();
            wrt.vplus = oldVPlus+step;
            double right = zz__lagrange_constraints();
            wrt.vplus = oldVPlus;
            grad[0] = (right-left)/(2*step);
        }

        // vminus:
        if(!wrt.grounded)
        {
            wrt.vminus -= step;
            double left = zz__lagrange_constraints();
            wrt.vminus = oldVMinus+step;
            double right = zz__lagrange_constraints();
            wrt.vminus = oldVMinus;
            grad[1] = (right-left)/(2*step);
        }
        
        // z:
        {
            wrt.z -= step;
            double left = zz__lagrange_constraints();
            wrt.z = oldVMinus+step;
            double right = zz__lagrange_constraints();
            wrt.z = oldZ;
            grad[2] = (right-left)/(2*step);
        }

        return ""+new ArrayRealVector(grad);
    }
    
    @SuppressWarnings("unused")
    private Object zz__showGradEstimate_inequalities = new Object()
    {
        public String toString() 
        {
            Agent wrt = Agent.this;
            
            return zz__gradEstimate_inequalities(wrt);
        }
    };
    
    private String zz__gradEstimate_inequalities(Agent wrt)
    {
        double oldVPlus = wrt.vplus;
        double oldVMinus = wrt.vminus;
        double oldZ = wrt.z;
        double step = 0.0001;
        double[] grad = new double[3];
        
        // vplus:
        {
            wrt.vplus -= step;
            double left = zz__lagrange_inequalities();
            wrt.vplus = oldVPlus+step;
            double right = zz__lagrange_inequalities();
            wrt.vplus = oldVPlus;
            grad[0] = (right-left)/(2*step);
        }

        // vminus:
        if(!wrt.grounded)
        {
            wrt.vminus -= step;
            double left = zz__lagrange_inequalities();
            wrt.vminus = oldVMinus+step;
            double right = zz__lagrange_inequalities();
            wrt.vminus = oldVMinus;
            grad[1] = (right-left)/(2*step);
        }
        
        // z:
        {
            wrt.z -= step;
            double left = zz__lagrange_inequalities();
            wrt.z = oldZ+step;
            double right = zz__lagrange_inequalities();
            wrt.z = oldZ;
            grad[2] = (right-left)/(2*step);
        }

        return ""+new ArrayRealVector(grad);
    }
    
    @SuppressWarnings("unused")
    private Object zz__showLagrange_inequalities_z = new Object()
    {
        public String toString() 
        {
            double zOld = z;
            StringBuffer sb = new StringBuffer();
            sb.append("z,mu(h()+z^2)+(a/2)(h()+z^2)^2\n");
            for(double x = 0.5; x < 1; x += 0.0025)
            {
                sb.append(x);
                sb.append(',');
                z = x;
                sb.append(lagrangeLocal_inequalities());
                z = zOld;
                sb.append('\n');
            }
            return sb.toString();
        };
    };

    private double zz__lagrange_constraints()
    {
        return lagrangeLocal_constraints() + sum(Agent::lagrangeLocal_constraints, neighbours());
    }

    private double zz__lagrange_inequalities()
    {
        return lagrangeLocal_inequalities() + sum(Agent::lagrangeLocal_inequalities, neighbours());
    }
    
    @SuppressWarnings("unused")
    private Object zz__showGradEstimate_cost = new Object()
    {
        public String toString() 
        {
            Agent wrt = Agent.this;
            
            return zz__gradEstimate_cost(wrt);
        }
    };
    
    private String zz__gradEstimate_cost(Agent wrt)
    {
        double oldVPlus = wrt.vplus;
        double oldVMinus = wrt.vminus;
        double oldZ = wrt.z;
        double step = 0.0001;
        double[] grad = new double[3];
        
        // vplus:
        {
            wrt.vplus -= step;
            double left = this.costLocal();
            wrt.vplus = oldVPlus+step;
            double right = this.costLocal();
            wrt.vplus = oldVPlus;
            grad[0] = (right-left)/(2*step);
        }

        // vminus:
        if(!wrt.grounded)
        {
            wrt.vminus -= step;
            double left = this.costLocal();
            wrt.vminus = oldVMinus+step;
            double right = this.costLocal();
            wrt.vminus = oldVMinus;
            grad[1] = (right-left)/(2*step);
        }
        
        // z:
        {
            wrt.z -= step;
            double left = this.costLocal();
            wrt.z = oldZ+step;
            double right = this.costLocal();
            wrt.z = oldZ;
            grad[2] = (right-left)/(2*step);
        }

        return ""+new ArrayRealVector(grad);
    }
    
    @SuppressWarnings("unused")
    private Object zz__showGradEstimate_costNeighbours = new Object()
    {
        public String toString() 
        {
            StringBuffer sb = new StringBuffer();
            for (Agent n : neighbours())
            {
                sb.append(n.getName());
                sb.append(":\n");
                
                sb.append(n.zz__gradEstimate_cost(Agent.this));
                
                sb.append('\n');
            }
            return sb.toString();
        }
    };
    
    @SuppressWarnings("unused")
    private Object zz__showBacktrack = new Object()
    {
        public String toString() 
        {
            StringBuffer sb = new StringBuffer();
            sb.append("t,delta L(), true delta L(),target,lower,delta part(), true delta part()\n");
            
            double vplusOld = vplus;
            double vminusOld = vminus;
            double zOld = z;
            
            RealVector grad = lagrangeGrad();
            double d = -DELTA*grad.dotProduct(grad);
            double before = zz__LagrangeLocal();
            double partBefore = zz__LagrangeLocal_inequalitySquared(); // WARNING: Don't forget to change all three
            
            double tmin = -2e-7;
            double tmax = 2e-6;
            double tinc = (tmax-tmin)/50.0; // Note that more than this truncates results
            for(double t = tmin; t < tmax; t+=tinc)
            {
                RealVector step = grad.mapMultiply(-t);
                
                vplus = vplusOld + step.getEntry(0);
                vminus = vminusOld + step.getEntry(1);
                z = zOld + step.getEntry(2);
                
                double after = zz__LagrangeLocal();
                double lagrangeDelta = after - before;
                double truePartDelta = zz__LagrangeLocal_inequalitySquared() - partBefore;
            
                // Reset:
                vplus = vplusOld;
                vminus = vminusOld;
                z = zOld;
                
                double partDelta = inequalitySquaredDelta(step);
                
                sb.append(t+","+lagrangeDelta(t, grad)+","+lagrangeDelta+","+(t*d)+","+(t*d/DELTA)+","+partDelta+","+truePartDelta);
                sb.append('\n');
            }
            return sb.toString();
        };
    };
    
    /**
     * Suitable only for calculating differences due to local changes.
     * @return Lagrange cost only including local and neighbour values.
     */
    private double zz__LagrangeLocal()
    {
        double cost = zz__costFull();
        
        double lambdaG = zz__LagrangeLocal_constraint();
        
        double alphaGSquared = zz__LagrangeLocal_constraintSquared();
        
        double muHZ = zz__LagrangeLocal_inequality();
        
        double alphaHZSquared = zz__LagrangeLocal_inequalitySquared();
        
        return cost + lambdaG + alphaGSquared + muHZ + alphaHZSquared;
    }

    private double zz__LagrangeLocal_inequalitySquared()
    {
        double alphaHZSquared = (alpha/2)*(h()+z*z)*(h()+z*z);
        alphaHZSquared += sum(n -> (n.alpha/2)*(n.h()+n.z*n.z)*(n.h()+n.z*n.z), neighbours());
        return alphaHZSquared;
    }

    private double zz__LagrangeLocal_constraintSquared()
    {
        double alphaGSquared = (alpha/2)*g().dotProduct(g());
        alphaGSquared += sum(n -> (n.alpha/2)*n.g().dotProduct(n.g()), neighbours());
        return alphaGSquared;
    }

    private double zz__LagrangeLocal_inequality()
    {
        double muHZ = mu*(h()+z*z);
        muHZ += sum(n -> n.mu*(n.h()+n.z*n.z), neighbours());
        return muHZ;
    }

    private double zz__LagrangeLocal_constraint()
    {
        double lambdaG = lambda.dotProduct(g());
        lambdaG += sum(n -> n.lambda.dotProduct(n.g()), neighbours());
        return lambdaG;
    }
    
    @SuppressWarnings("unused")
    private Object zz__showLagrange_vplus = new Object()
    {
        public String toString() 
        {
            double vOld = vplus;
            StringBuffer sb = new StringBuffer();
            sb.append("v+,L(v+)\n");
            
            double low = 12.06001;
            double high = 12.06002;
            double step = (high-low)/100;
            for(double v = low; v < high; v += step)
            {
                sb.append(v);
                sb.append(',');
                vplus = v;
                sb.append(lagrange());
                vplus = vOld;
                sb.append('\n');
            }
            return sb.toString();
        };
    };
    
    @SuppressWarnings("unused")
    private Object zz__showLagrange_vminus = new Object()
    {
        public String toString() 
        {
            double vOld = vminus;
            StringBuffer sb = new StringBuffer();
            sb.append("v-,L(v-)\n");
            for(double v = -1; v < 1; v += 0.1)
            {
                sb.append(v);
                sb.append(',');
                vminus = v;
                sb.append(lagrange());
                vminus = vOld;
                sb.append('\n');
            }
            return sb.toString();
        };
    };
    
    @SuppressWarnings("unused")
    private Object zz__showLagrange_z = new Object()
    {
        public String toString() 
        {
            double vOld = z;
            StringBuffer sb = new StringBuffer();
            sb.append("z,L(z)\n");
            for(double x = -5; x < 5; x += 0.1)
            {
                sb.append(x);
                sb.append(',');
                z = x;
                sb.append(lagrange());
                z = vOld;
                sb.append('\n');
            }
            return sb.toString();
        };
    };
    
    @SuppressWarnings("unused")
    private Object zz__showCost_vplus = new Object()
    {
        public String toString() 
        {
            double vOld = vplus;
            StringBuffer sb = new StringBuffer();
            sb.append("v+,c_i(v+)\n");
            double low = 11.7;
            double high = 12.7;
            double step = (high - low)/50;
            for(double x = low; x < high; x += step)
            {
                sb.append(x);
                sb.append(',');
                vplus = x;
                sb.append(costLocal());
                vplus = vOld;
                sb.append('\n');
            }
            return sb.toString();
        };
    };
    
    @SuppressWarnings("unused")
    private Object zz__showCostFull_vplus = new Object()
    {
        public String toString() 
        {
            double vOld = vplus;
            StringBuffer sb = new StringBuffer();
            sb.append("v+,c(v+)\n");
            double low = 11.7;
            double high = 12.7;
            double step = (high - low)/50;
            for(double x = low; x < high; x += step)
            {
                sb.append(x);
                sb.append(',');
                vplus = x;
                sb.append(zz__costFull());
                vplus = vOld;
                sb.append('\n');
            }
            return sb.toString();
        }
    };

    private double zz__costFull()
    {
        return costLocal()+sum(Agent::costLocal, neighbours());
    };
    
    @SuppressWarnings("unused")
    private Object zz__showCost_vminus = new Object()
    {
        public String toString() 
        {
            double vOld = vminus;
            StringBuffer sb = new StringBuffer();
            sb.append("v-,c(v+)\n");
            for(double x = -2; x < 2; x += 0.1)
            {
                sb.append(x);
                sb.append(',');
                vminus = x;
                sb.append(costLocal());
                vminus = vOld;
                sb.append('\n');
            }
            return sb.toString();
        };
    };
    
    @SuppressWarnings("unused")
    private Object zz__showCost_z = new Object()
    {
        public String toString() 
        {
            double vOld = z;
            StringBuffer sb = new StringBuffer();
            sb.append("z,c(v+)\n");
            for(double x = -5; x < 5; x += 0.1)
            {
                sb.append(x);
                sb.append(',');
                z = x;
                sb.append(costLocal());
                z = vOld;
                sb.append('\n');
            }
            return sb.toString();
        };
    };
}