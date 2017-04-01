package ellipsis.energy.sandbox.sandbox020;
import static ellipsis.util.ApproximateComparator.approxEquals;
import static ellipsis.util.Sum.sum;
import static ellipsis.util.Sum.sumV;
import static ellipsis.util.VectorHelper.vector;

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
    
    public static final boolean enableInequalities = true;
    
    public static final int EQUALITY_CONSTRAINT_COUNT = 2;
    public static final int STATE_DIMENSION = 2; // [v+, v-]
    
//    public double alphaInit = 0.5e-2;
    public double alphaMultiplier = 1.0001; // 2 node (1.004)
    public double alphaMax = 1e1;
    
    public static final double epsilonInit = 10.0;
    public static final double epsilonMultiplier = 1.0;//0.998;//0.998;
    public static final double epsilonMin = 1e-3;

    public double lambdaMax = 1e6;
    public double lambdaMultiplier = 1.0;

    public double muMax = 1e6;
    public double muMultiplier = 1.0;
    
    protected double vminus, vplus;  // Primal state
    protected double vmin, vmax; // Admissible voltage range
    protected RealVector lambda; // Dual variable for equality constraints
    protected double mu; // Dual variable for inequality constraints
    protected double alpha;
    protected double epsilon;
    private Map<Agent, Double> conductances = new HashMap<>(); // Constants
    private double ySum = 0; // Sum of neighbouring conductances
    private boolean grounded = false;
    private String name;
    
    public void init()
    {
        // Disable constraints:
        double oldAlpha = alpha;
        RealVector oldLambda = lambda;
        double oldMu = mu;
        alpha = 0;
        lambda = vector(STATE_DIMENSION, 0.0);
        mu = 0;
        
        // Step primal state [v+, v-]:
//        RealVector grad = vector(STATE_DIMENSION, Double.MAX_VALUE);

//        while(grad.getNorm() > 0.1)
//        {
//            grad = minimizeApproximate(1, 1);
//        }
        
        // Reset constraints:
        alpha = oldAlpha;
        lambda = oldLambda;
        mu = oldMu;
    }
    
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
        if(grounded)
        	vminus = 0.0;
    }
    
    public void addNeighbour(Agent neighbour, double conductance)
    {
        ySum = 0;
        conductances.put(neighbour, conductance);
    }
    
    public double conductance(Agent neighbour)
    {
        Double y_ij = conductances.get(neighbour);
		return y_ij == null ? 0.0 : y_ij;
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
    	if(!grounded)
    		this.vminus = vminus;
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

    public void setLambda(RealVector lambda)
    {
        this.lambda = lambda;
    }
    
    public double getMu()
    {
        return mu;
    }
    
    public void setMu(double mu)
    {
        this.mu = mu;
    }

    public double getAlpha()
    {
        return alpha;
    }

    public void setAlpha(double alpha)
    {
        this.alpha = alpha;
    }
    
    public void setEpsilon(double epsilon)
    {
        this.epsilon = epsilon;
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

    public double conductanceSum()
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
    
    public void setAlphaMultiplier(double alphaMultiplier)
    {
        this.alphaMultiplier = alphaMultiplier;
    }
    
    public void setAlphaMax(double alphaMax)
    {
        this.alphaMax = alphaMax;
    }
    
    public double getAlphaMultiplier()
    {
        return alphaMultiplier;
    }
    
    public double getAlphaMax()
    {
        return alphaMax;
    }
    
    public void setLambdaMax(double lambdaMax)
    {
        this.lambdaMax = lambdaMax;
    }
    
    public double getLambdaMax()
    {
        return lambdaMax;
    }
    
    public void setLambdaMultiplier(double lambdaMultiplier)
    {
        this.lambdaMultiplier = lambdaMultiplier;
    }
    
    public double getLambdaMultiplier()
    {
        return lambdaMultiplier;
    }
    
    public double getMuMax()
    {
        return muMax;
    }
    
    public void setMuMax(double muMax)
    {
        this.muMax = muMax;
    }
    
    public double getMuMultiplier()
    {
        return muMultiplier;
    }
    
    public void setMuMultiplier(double muMultiplier)
    {
        this.muMultiplier = muMultiplier;
    }
    
    
    //// Initialisation ////
    
    public Agent()
    {
        lambda = new ArrayRealVector(EQUALITY_CONSTRAINT_COUNT);
//        alpha = alphaInit;
        epsilon = epsilonInit;
    }
    
    
    //// Abstract Functions ////
    
    /**
     * Equality constraints for this agent.
     * @return A 2 dimensional vector.
     */
    public abstract RealVector g();
    
    /**
     * Equality constraint gradient: For this agent's constraints 
     * with respect to the given agent's state.
     * @param wrt The agent the gradient is to be with respect to. This is assumed to be a neighbour.
     * @return A 3x2 matrix.
     */
    public abstract RealMatrix gGrad(Agent wrt);
    
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
    
    public abstract double inequalityDelta(Agent wrt, RealVector step);

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
    
    public abstract double h();
    
    public abstract RealVector hGrad(Agent wrt); 
    
    
    //// Optimisation ////
    
    public void step()
    {
        // Step primal state [v+, v-]:
        RealVector grad = vector(STATE_DIMENSION, Double.MAX_VALUE);

        while(grad.getNorm() > epsilon)
        {
            grad = minimizeApproximate(1, 1);
//            double g0 = minimizeApproximate(1, 0).getEntry(0);
//            grad.setEntry(0, g0); // step v+
//            double g1 = minimizeApproximate(0, 1).getEntry(1);
//            grad.setEntry(1, g1); // step v-
        }

        validate(vplus, "v+");
        validate(vminus, "v-");

        // Step dual variables:
        RealVector g = g();
        lambda = constrain(lambda.add(g.mapMultiply(alpha*lambdaMultiplier)), lambdaMax);
        validate(lambda, "lambda");
        
        double h = h();
        mu = constrain(Math.max(0, mu + h*alpha*muMultiplier), muMax);
        validate(mu, "mu");
        
        // Increase alpha:
        alpha *= alphaMultiplier;
        if(alpha > alphaMax)
            alpha = alphaMax;
        validate(alpha, "alpha");
        
        // Decrease epsilon:
        epsilon *= epsilonMultiplier;
        validate(epsilon, "epsilon");
        if(epsilon < epsilonMin)
            epsilon = epsilonMin;
    }

    public static double constrain(double d, double max)
    {
        if(d > max)
            return max;
        if(d < -max)
            return -max;
        return d;
    }

    public static RealVector constrain(RealVector v, double max)
    {
        int dimension = v.getDimension();
        for(int i = 0; i < dimension; ++i)
            v.setEntry(i, constrain(v.getEntry(i), max));
        return v;
    }

    @SuppressWarnings("deprecation")
    private RealVector minimizeApproximate(double vplusFlag, double vminusFlag)
    {
        RealVector grad = lagrangeGrad();
        grad.setEntry(0, grad.getEntry(0)*vplusFlag);
        grad.setEntry(1, grad.getEntry(1)*vminusFlag);
        validate(grad, "Gradient");
        
        double stepSize = backtrack(grad);
        validate(stepSize, "step size");
        
        RealVector step = grad.mapMultiply(-stepSize);
        if(step.getDimension() != STATE_DIMENSION)
            throw new RuntimeException("wrong dimension for step");
//        double lagrangeBefore = lagrange();
        stepState(step);
        
        RealVector projected = projectVoltages();
        grad = lagrangeGrad();
        grad = grad.ebeMultiply(projected);
        
        // FIXME uncomment
//        double lagrangeAfter = lagrange();
//        if(lagrangeBefore < lagrangeAfter)
//            throw new RuntimeException(
//                    "Lagrange increased after step: From "+lagrangeBefore+" to "+lagrangeAfter
//                    +", with step "+step+" for agent "+getName());
        
        return grad;
    }

    /**
     * Add the given values to the agent's state.
     * @param step [vStepPlus, vStepMinus]
     */
    public void stepState(RealVector step)
    {
        vplus += step.getEntry(0);  // delta v+
        vminus += step.getEntry(1); // delta v-
    }

    /**
     * If the resultant voltage (v+ - v-) is not
     * within the admissible set defined by vmin and vmax then its value will be projected.
     * @return Vector with elements equal to 0.0 if the index was projected, and equal to 1.0 otherwise.
     */
    public RealVector projectVoltages()
    {
        RealVector projected = new ArrayRealVector(STATE_DIMENSION, 1.0);
        if(isGrounded())
        {
            if(vplus > vmax)
            {
                vplus = vmax;
                projected.setEntry(0, 0.0);
            }
            else if(vplus < vmin)
            {
                vplus = vmin;
                projected.setEntry(0, 0.0);
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
                projected.setEntry(0, 0.0);
                projected.setEntry(1, 0.0);
            }
        }
        
        return projected;
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
        
        // Constraint gradients:
        RealVector constraintGrad = constraintGrad();

        return costGrad.add(constraintGrad);
    }

    private RealVector constraintGrad()
    {
        RealVector eqConGrad, inConGrad;
        
        // sum{ grad g(lambda + alpha*g()) }
        eqConGrad = equalityGrad();
        
        // Inequalities:
        inConGrad = inequalityGrad();
        
        return eqConGrad.add(inConGrad);
    }

    private RealVector equalityGrad()
    {
        RealVector eqConGrad;
        {
            RealMatrix gGradSelf = gGrad(this);
            validate(gGradSelf, "g gradient (self)");
            RealVector eqConGradSelf = gGradSelf.operate(lambda.add(g().mapMultiply(alpha)));
            validate(eqConGradSelf, "constraint grad (self)");
            RealVector eqConGradNeighbours = 
                sumV(neighbour -> 
                    neighbour.gGrad(this).operate(neighbour.lambda.add(neighbour.g().mapMultiply(neighbour.alpha))),
                    neighbours(), 
                    STATE_DIMENSION
                );
            validate(eqConGradNeighbours, "constraint grad (neighbours)");
            
            eqConGrad = eqConGradSelf.add(eqConGradNeighbours);
        }
        return eqConGrad;
    }

    private RealVector inequalityGrad()
    {
        if(!enableInequalities)
            return new ArrayRealVector(STATE_DIMENSION);

        RealVector hGradSelf = hGrad(this);
        validate(hGradSelf, "h gradient (self)");
        RealVector inConGradSelf = (mu + alpha*h() > 0) ? hGradSelf.mapMultiply(mu + alpha*h()) : vector(0.0, 0.0);
        validate(inConGradSelf, "inequality cons grad (self)");
        RealVector inConGradNeighbours = 
            sumV(neighbour -> 
                (neighbour.mu + neighbour.alpha*neighbour.h() > 0) ?
                        neighbour.hGrad(this).mapMultiply(neighbour.mu + neighbour.h()*neighbour.alpha) :
                        vector(0.0, 0.0),
                neighbours(), 
                STATE_DIMENSION
            );
        validate(inConGradNeighbours, "inequality cons grad (neighbours)");
        
        return inConGradSelf.add(inConGradNeighbours);
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
        return costDelta + constraintDelta + constraintSquaredDelta + inequalityDelta;
    }

    public double inequalityDelta(RealVector step)
    {
        if(!enableInequalities)
            return 0;
        
        double deltaSelf = inequalityDelta(this, step);
        double deltaNeighbours = sum(n -> n.inequalityDelta(this, step), neighbours());
        
        return deltaSelf + deltaNeighbours;
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
        return costLocal() + lagrangeLocal_equalities() + lagrangeLocal_inequalities();
    }

    public double lagrangeLocal_inequalities()
    {
        if(!enableInequalities)
            return 0;
        
        double h = h();
        if(mu + alpha*h > 0)
        {
            return (alpha/2.0)*h*h + mu*h;
        }
        else
        {
            return -mu*mu/(2*alpha);
        }
    }

    private double lagrangeLocal_equalities()
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
    
    public double lagrangeCheck()
    {
        return cost() + equalities() + inequalities();
    }
    
    public double cost()
    {
        return costLocal() + sum(Agent::costLocal, neighbours());
    }
    
    public double equalities()
    {
        return lagrangeLocal_equalities() + sum(Agent::lagrangeLocal_equalities, neighbours());
    }
    
    public double inequalities()
    {
        return lagrangeLocal_inequalities() + sum(Agent::lagrangeLocal_inequalities, neighbours());
    }
    
    public RealVector gradEstimate()
    {
        double oldVPlus = vplus;
        double oldVMinus = vminus;
        double step = 0.0001;
        double[] grad = new double[2];
        
        // vplus:
        {
            vplus -= step;
            double left = this.lagrange();
            vplus = oldVPlus+step;
            double right = this.lagrange();
            vplus = oldVPlus;
            grad[0] = (right - left)/(2*step);
        }

        // vminus:
        if(!grounded)
        {
            vminus -= step;
            double left = this.lagrange();
            vminus = oldVMinus+step;
            double right = this.lagrange();
            vminus = oldVMinus;
            grad[1] = (right - left)/(2*step);
        }
        else
        {
            grad[1] = 0;
        }
        
        return new ArrayRealVector(grad);
    }
    
    public RealVector _gradEstimate()
    {
        return gradEstimate(this::lagrange);
    }
    
    public RealVector gradEstimate_cost()
    {
        return gradEstimate(this::cost);
    }
    
    public RealVector gradEstimate_equalities()
    {
        return gradEstimate(this::equalities);
    }
    
    public RealVector gradEstimate_inequalities()
    {
        return gradEstimate(this::inequalities);
    }
    
    public static interface ValueFunction
    {
        double value();
    }
    public RealVector gradEstimate(ValueFunction f)
    {
        double oldVPlus = vplus;
        double oldVMinus = vminus;
        double step = 0.0001;
        double[] grad = new double[2];
        
        // vplus:
        {
            vplus -= step;
            double left = f.value();
            vplus = oldVPlus+step;
            double right = f.value();
            vplus = oldVPlus;
            grad[0] = (right - left)/(2*step);
        }

        // vminus:
        if(!grounded)
        {
            vminus -= step;
            double left = f.value();
            vminus = oldVMinus+step;
            double right = f.value();
            vminus = oldVMinus;
            grad[1] = (right - left)/(2*step);
        }
        else
        {
            grad[1] = 0;
        }
        
        return new ArrayRealVector(grad);
    }

    
    //// Debug ////
    
    @SuppressWarnings("unused")
    private Object zz__showGrad = new Object()
    {
        public String toString() 
        {
            return lagrangeGrad().toString();
        };
    };

    @SuppressWarnings("unused")
    private Object zz__showGrad_gSelf = new Object()
    {
        public String toString() 
        {
            return ""+gGrad(Agent.this);
        }
    };
    
    @SuppressWarnings("unused")
    private Object zz__showGrad_Constraint = new Object()
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
        double step = 0.0001;
        double[][] grad = new double[2][];
        
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

    private String zz__gradEstimate(Agent wrt)
    {
        double oldVPlus = wrt.vplus;
        double oldVMinus = wrt.vminus;
        double step = 0.0001;
        double[] grad = new double[2];
        
        // vplus:
        {
            wrt.vplus -= step;
            double left = this.lagrange();
            wrt.vplus = oldVPlus+step;
            double right = this.lagrange();
            wrt.vplus = oldVPlus;
            grad[0] = (right - left)/(2*step);
        }

        // vminus:
        if(!wrt.grounded)
        {
            wrt.vminus -= step;
            double left = this.lagrange();
            wrt.vminus = oldVMinus+step;
            double right = this.lagrange();
            wrt.vminus = oldVMinus;
            grad[1] = (right - left)/(2*step);
        }
        else
        {
            grad[1] = 0;
        }
        
        return ""+new ArrayRealVector(grad);
    };
    
    @SuppressWarnings("unused")
    private Object zz__showGradEstiamte = new Object()
    {
       public String toString() 
       {
           StringBuffer sb = new StringBuffer();
           sb.append(zz__gradEstimate(Agent.this));
           sb.append('\n');
           for (Agent a : neighbours())
           {
               sb.append(zz__gradEstimate(a));
               sb.append('\n');
           }
           
           return sb.toString();
       }; 
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
        double step = 0.0001;
        double[] grad = new double[2];
        
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

        return ""+new ArrayRealVector(grad);
    }

    private double zz__lagrange_constraints()
    {
        return lagrangeLocal_equalities() + sum(Agent::lagrangeLocal_equalities, neighbours());
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
        double step = 0.0001;
        double[] grad = new double[2];
        
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
            
            RealVector grad = lagrangeGrad();
            double d = -DELTA*grad.dotProduct(grad);
            double before = lagrange();
            double partBefore = 0;//zz__LagrangeLocal_constraintSquared(); // WARNING: Don't forget to change all three
            
            double tmin = -2e-5;
            double tmax = -10*tmin;//2e-7;
            double tinc = (tmax-tmin)/50.0; // Note that more than this truncates results
            for(double t = tmin; t < tmax; t+=tinc)
            {
                RealVector step = grad.mapMultiply(-t);
                
                vplus = vplusOld + step.getEntry(0);
                vminus = vminusOld + step.getEntry(1);
                
                double after = lagrange();
                double lagrangeDelta = after - before;
                double truePartDelta = 0;//zz__LagrangeLocal_constraintSquared() - partBefore;
            
                // Reset:
                vplus = vplusOld;
                vminus = vminusOld;
                
                double partDelta = 0;//constraintSquaredDelta(step);
                
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
    @SuppressWarnings("unused")
    private double zz__LagrangeLocal()
    {
        double cost = zz__costFull();
        
        double lambdaG = zz__LagrangeLocal_constraint();
        
        double alphaGSquared = zz__LagrangeLocal_constraintSquared();
        
        return cost + lambdaG + alphaGSquared + lagrangeLocal_inequalities();
    }

    private double zz__LagrangeLocal_constraintSquared()
    {
        double alphaGSquared = (alpha/2)*g().dotProduct(g());
        alphaGSquared += sum(n -> (n.alpha/2)*n.g().dotProduct(n.g()), neighbours());
        return alphaGSquared;
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
            
            double low = 11;
            double high = 12.5;
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
}