package ellipsis.energy.sandbox;

import static ellipsis.util.Sum.sum;
import static ellipsis.util.VectorHelper.vector;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Util;

import ellipsis.util.TeeOutputStream;

public class Sandbox018C_Centralized extends Sandbox018B
{
    public static void main(String[] args)
    {
        new Sandbox018C_Centralized().run();
    }
    
    public Sandbox018C_Centralized()
    {
        super();
        
        try
        {
            out = new PrintStream(new TeeOutputStream(new FileOutputStream("/tmp/Sandbox018C.csv"), System.out));
        }
        catch (FileNotFoundException e)
        {
            throw new RuntimeException(e);
        }
        
        //// Override configuration ////
        
                    PROJECT_X = true;
               FORMAT_NUMBERS = false;
            
            VALIDATE_BACKTRACKING = false;
       VALIDATE_LAGRANGE_DECREASE = false;
            
           START_WITH_TRUE_VOTLAGES  = false;
           START_WITH_OPTIMAL_POWERS = true;
            USE_SIGMOID_ALPHA_UPDATE = false;
                          USE_TRUE_G = false;
            
                ROUND_TO_ZERO = false;
                    MIN_VALUE = 1e-12;
                  MIN_G_VALUE = 1e-6;
                  MIN_H_VALUE = 1e-6;
            
                           KEEP_G_AND_H_CLOSE = true;
            DONT_UPDATE_ALPHA_WHILE_IMPROVING = true;
                   MIN_CONSTRAINT_IMPROVEMENT = 0.99;
                            ALPHA_UPDATE_RATE = 1;
            
            APPROXIMATE_CONSTRAINTS = true;
                    CONSTRAINT_BASE = 1.0;
                  CONSTRAINT_TARGET = 0.001;
            
                         DISABLE_COST = false;
             DISABLE_VOTLAGES_UPDATES = false;
            DISABLE_POWER_FLOW_UPDATE = false;
                DISABLE_LAMBDA_UPDATE = false;
            
            DISABLE_CONSENSUS_UPDATE = false;
                   DISABLE_MU_UPDATE = false;
                   DISABLE_LINE_LOSS = false;
            
            INITIAL_G_AUG_SCALE = 0.001;
               G_AUG_SCALE_STEP = 1.001;
                G_MAX_AUG_SCALE = 0.01;//1e6;
            
            INITIAL_H_AUG_SCALE = 8;
               H_AUG_SCALE_STEP = 1.001;
                H_MAX_AUG_SCALE = 8;//1e6;
            
                        ETA_G = 3;
                        ETA_H = 0.005;
                       
                           XI = 0.2;
                   
                                      K = 20000;
                             DEBUG_RATE = K / 1000;
            AGENT_SELECTION_PROBABILITY = 1.0;
            
            INITIAL_STEP_SIZE = 1;
                MIN_STEP_SIZE = 1e-50;
                  MAX_X_STEPS = 1;
                  
              EPSILON_BASE_PQ = 100;//0.1;
              EPSILON_BASE_EF = 1000;//100;
               EPSILON_TARGET = 1;//10e-2;
    }
    
    @Override
    public void loop()
    {
        double min = minimise();
        out.println("min = "+min);
    }
    
    public double minimise()
    {
        RealVector x = p
               .append(q)
               .append(e)
               .append(f);
        
        RealVector lambda = lp
                    .append(lq);
        
        RealVector mu = vector(2, 0.0);
        
        double a_g = INITIAL_G_AUG_SCALE;
        double a_h = INITIAL_H_AUG_SCALE;
        
        for(int k = 0; k < K; ++k)
        {
if(k > 10)
    Util.nullop();
            // 1. Projected gradient descent:
            double lagBefore = L(x, lambda, mu, a_g, a_h);
//            for(int i = 0; i < 100; ++i)
            {
                RealVector grad = gradL();
                int dimension = grad.getDimension();
                for(int j = 0; j < dimension; ++j)
                {
                    RealVector _grad = new ArrayRealVector(dimension);
                    _grad.setEntry(j, grad.getEntry(j));
                    double gamma = backtrack(_grad, x, lambda, mu, a_g, a_h);
                    x = project(x.subtract(_grad.mapMultiply(gamma)));
                }
            }
            double lagAfter = L(x, lambda, mu, a_g, a_h);
            if(lagAfter > lagBefore)
                System.err.println("L() increased from "+lagBefore+" to "+lagAfter);
            
            // 2. Projected lambda step:
            lambda = lambda.add(g(x).mapMultiply(a_g));
            
            // 3. Projected mu step:
            mu = mu.add(h(x).mapMultiply(a_h));
            
            // 4. Increment alpha:
            if(a_g < G_MAX_AUG_SCALE)
                a_g *= G_AUG_SCALE_STEP;
            if(a_h < H_MAX_AUG_SCALE)
                a_h *= H_AUG_SCALE_STEP;
            
            // Update values:
            {
                int dimension = p.getDimension();
                RealVector p = x.getSubVector(0, dimension);
                RealVector q = x.getSubVector(dimension, dimension);
                RealVector e = x.getSubVector(2*dimension, dimension);
                RealVector f = x.getSubVector(3*dimension, dimension);
                RealVector lp = lambda.getSubVector(0, dimension);
                RealVector lq = lambda.getSubVector(dimension, dimension);
                double mup = mu.getEntry(0);
                double muq = mu.getEntry(1);
                for(int i = 0; i < dimension; ++i)
                {
                    this.p.setEntry(i, p.getEntry(i));
                    this.q.setEntry(i, q.getEntry(i));
                    this.e.setEntry(i, e.getEntry(i));
                    this.f.setEntry(i, f.getEntry(i));
                    this.lp.setEntry(i, lp.getEntry(i));
                    this.lq.setEntry(i, lq.getEntry(i));
                    this.mup.setEntry(i, mup);
                    this.muq.setEntry(i, muq);
                    this.hp.setEntry(i, hp(p, q, e, f));
                    this.hq.setEntry(i, hq(p, q, e, f));
                    this.alpha_g.setEntry(i, a_g);
                    this.alpha_h.setEntry(i, a_h);
                }
            }

            // Log:
            if(k%DEBUG_RATE == 0)
                debug(k);
        }
        
        return L(x, lambda, mu, a_g, a_h);
    }

    protected RealVector gradL()
    {
        int dimension = p.getDimension();
        return          vector(dimension, this::gradL_p, slackIndex)
                .append(vector(dimension, this::gradL_q, slackIndex))
                .append(vector(dimension, this::gradL_e, slackIndex))
                .append(vector(dimension, this::gradL_f, slackIndex));
    }

    private double backtrack(RealVector grad, RealVector x, RealVector lambda, RealVector mu, double a_g, double a_h)
    {
        double t = 1.0;
        double alpha = 0.5; // Step update size.
        double beta = 0.5;
        
        double slope = alpha*grad.dotProduct(grad);
        double L0 = L(x, lambda, mu, a_g, a_h);
        while( L(x.subtract(grad.mapMultiply(t)), lambda, mu, a_g, a_h) > L0 - t*slope )
            t *= beta;
        
        return t;
    }

    private double L(RealVector x, RealVector lambda, RealVector mu, double a_g, double a_h)
    {
        int dimension = x.getDimension()/4;
        double h_norm = h(x).getNorm();
        double g_norm = g(x).getNorm();
        RealVector p = x.getSubVector(0, dimension);
        RealVector q = x.getSubVector(dimension, 2*dimension);
        
        return   cost(p, q)
                +lambda.dotProduct(g(x))
                +mu.dotProduct(h(x))
                +(a_g/2)*g_norm*g_norm
                +(a_h/2)*h_norm*h_norm;
    }

    private RealVector project(RealVector vin)
    {
        int dimension = p.getDimension();
        RealVector vout = new ArrayRealVector(vin.getDimension());
        
        for(int i = 0; i < dimension; ++i)
        {
            // Power:
            {
                double p = vin.getEntry(i);
                double q = vin.getEntry(i+dimension);
                
                if(generatorMask.getEntry(i) == 1)
                {
                    if(p < 0)
                        p = 0;
                    
                    Complex s = new Complex(p, q);
                    double abs = s.abs();
                    double p_max_i = p_max.getEntry(i);
                    if(abs > p_max_i)
                    {
                        s = s.multiply(p_max_i/abs);
                        p = s.getReal();
                        q = s.getImaginary();
                    }
                }
                
                vout.setEntry(i, p);
                vout.setEntry(i+dimension, q);
            }
            
            // Voltage:
            {
                double e = vin.getEntry(i+2*dimension);
                double f = vin.getEntry(i+3*dimension);
                if(e < 0)
                    e = 0;
                
                // Clamp |v|:
                Complex v = new Complex(e, f);
                double abs = v.abs();
                while(abs > V_MAX || abs < V_MIN)
                {
                    if(abs > V_MAX)
                    {
                        v = v.multiply(V_MAX/abs);
                    }
                    else if(abs < V_MIN)
                    {
                        v = v.multiply(V_MIN/abs);
                    }
        
                    abs = v.abs();
                }
                
                // Clamp angle:
                double arg = v.getArgument();
                if(arg < V_ARG_MIN)
                {
                    v = ComplexUtils.polar2Complex(abs, V_ARG_MIN);
                }
                else if(arg > V_ARG_MAX)
                {
                    v = ComplexUtils.polar2Complex(abs, V_ARG_MAX);
                }
                
                e = v.getReal();
                f = v.getImaginary();
                
                vout.setEntry(i+2*dimension, e);
                vout.setEntry(i+3*dimension, f);
            }
        }
        
        return vout;
    }

    private RealVector g(RealVector x)
    {
        int dimension = x.getDimension()/4;
        RealVector p = x.getSubVector(0, dimension);
        RealVector q = x.getSubVector(dimension, dimension);
        RealVector e = x.getSubVector(2*dimension, dimension);
        RealVector f = x.getSubVector(3*dimension, dimension);
        return gp(p, q, e, f).append(gq(p, q, e, f));
    }

    private RealVector h(RealVector x)
    {
        int dimension = x.getDimension()/4;
        RealVector p = x.getSubVector(0, dimension);
        RealVector q = x.getSubVector(dimension, dimension);
        RealVector e = x.getSubVector(2*dimension, dimension);
        RealVector f = x.getSubVector(3*dimension, dimension);
        return vector(hp(p, q, e, f), hq(p, q, e, f));
    }
    
    
    //// Debug ////
    
    @SuppressWarnings("unused")
    private Object zz__showBacktrack = new Object()
    {
        public String toString() 
        {
            StringBuffer sb = new StringBuffer();
            sb.append("t,delta L(), target,lower,delta part(), true delta part()\n");
            
            RealVector x = p.append(q).append(e).append(f);
            RealVector xOld = x;
            RealVector lambda = lp.append(lq);
            RealVector mu = vector(mup.getEntry(0), muq.getEntry(0));
            double a_g = alpha_g.getEntry(0);
            double a_h = alpha_h.getEntry(0);
            
            RealVector grad = gradL();
            double d = -0.5*grad.dotProduct(grad);
            double before = L(x, lambda, mu, a_g, a_h);
            double partBefore = 0;//zz__LagrangeLocal_inequalitySquared(); // WARNING: Don't forget to change all three
            
            double tmin = -2e-7;
            double tmax = 2e-6;
            double tinc = (tmax-tmin)/50.0; // Note that more than this truncates results
            for(double t = tmin; t < tmax; t+=tinc)
            {
                RealVector step = grad.mapMultiply(-t);
                
                x = xOld.add(step);
                
                double after = L(x, lambda, mu, a_g, a_h);
                double lagrangeDelta = after - before;
                double truePartDelta = 0;//zz__LagrangeLocal_inequalitySquared() - partBefore;
            
                // Reset:
                x = xOld;
                
                double partDelta = 0;//inequalitySquaredDelta(step);
                
                sb.append(t+","+lagrangeDelta+","+(t*d)+","+(t*d/0.5)+","+partDelta+","+truePartDelta);
                sb.append('\n');
            }
            return sb.toString();
        };
    };
    
    private Object zz__showLagrangeRegion = new Object()
    {
        public String toString() 
        {
            StringBuffer sb = new StringBuffer();
            sb.append("i,L(+delta),L(x),'L(-delta),grad\n");
            
            RealVector x = p.append(q).append(e).append(f);
            RealVector lambda = lp.append(lq);
            RealVector mu = vector(mup.getEntry(0), muq.getEntry(0));
            double a_g = alpha_g.getEntry(0);
            double a_h = alpha_h.getEntry(0);
            RealVector grad = gradL();
            double l = L(x, lambda, mu, a_g, a_h);
            
            int dimension = x.getDimension();
            for (int i = 0; i < dimension; i++)
            {
                sb.append(i);
                sb.append(",");
                RealVector step = new ArrayRealVector(dimension);
                step.setEntry(i, 1e-3);
                sb.append(L(x.add(step), lambda, mu, a_g, a_h));
                sb.append(",");
                sb.append(l);
                sb.append(",");
                sb.append(L(x.subtract(step), lambda, mu, a_g, a_h));
                sb.append(",");
                sb.append(grad.getEntry(i));
                sb.append("\n");
            }
            
            return sb.toString();
        };
    };
}
