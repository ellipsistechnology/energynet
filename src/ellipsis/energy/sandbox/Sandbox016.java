package ellipsis.energy.sandbox;

import static ellipsis.energy.sandbox.GridConstants.SLACK_SOURCE;
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Timer;
import com.mls.util.Util;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.GridGenerator;
import ellipsis.energy.grid.Load;
import ellipsis.energy.grid.SlackSource;
import ellipsis.energy.grid.Source;
import ellipsis.util.MatrixHelper;
import ellipsis.util.TeeOutputStream;
import static ellipsis.util.Sum.*;
import static ellipsis.util.VectorHelper.*;

/**
 * Asynchronous version of {@link Sandbox014} - copied from {@link Sandbox015} - with larger network.
 * @author bmillar
 *
 */
public class Sandbox016
{
    public static void main(String[] args)
	{
		new Sandbox016().run();
	}

    
    public static interface NeighbourFunction
    {
        double value(Agent neighbour);
    }
    
    public static double sumA(Map<Agent, Complex> keys, NeighbourFunction f)
    {
        double sum = 0;
        for (Agent agent : keys.keySet())
        {
            double v = f.value(agent);
            sum += v;
        }
        
        return sum;
    }
	
    public static PrintStream out;
    static
    {
        try
        {
            out = new PrintStream(new TeeOutputStream(new FileOutputStream("/tmp/Sandbox016.csv"), System.out));
        } 
        catch (FileNotFoundException e)
        {
            throw new RuntimeException(e);
        }
    }
    
    
    //// Global variables and constants ////
    
    private static final double BASE_VOLTAGE = 120;
    private static final double SLACK_VOLTAGE = 1.00;
    private static final double BASE_POWER = 10e3;
    
    private static final double V_MAX = 1.05;
    private static final double V_MIN = 0.95;
    private static final double V_ARG_MIN = -Math.PI/4;
    private static final double V_ARG_MAX = Math.PI/4;

    private static final double MIN_CONSTRAINT_IMPROVEMENT = 0.99; // Minimum improvement percentage of g(x) constraint
    private static final double INITIAL_AUG_SCALE = 1e-3; // 1e-3; FIXME DO NOT COMMIT
    private static final double AUG_SCALE_STEP = 1.0025;//1.0001; //1.0025; FIXME DO NOT COMMIT
    private static final double MAX_AUG_SCALE = 1e6; //1e6; FIXME DO NOT COMMIT
    
    private static final double COST_MULTIPLIER = 1;
    
    private static final int X_STEPS_PER_ITERATION = 1;
    
    private static final int K = 3000;//50000; // FIXME DO NOT COMMIT
    private static final int DEBUG_RATE = K/1000;
    
    private static final double INITIAL_STEP_SIZE = 1;
    private static final double MIN_STEP_SIZE = 1e-24;
    
    private static final double MAX_LAMBDA = 1e14;
    
    // Switches (not final to avoid warnings when set to false):
    protected static boolean PROJECT_X = true;
    protected static boolean USE_NEWTON = false;
    protected static boolean USE_NESTEROV = false;
    protected static boolean FORMAT_NUMBERS = false;
    protected static boolean START_WITH_TRUE_VOTLAGES = false;
    protected static boolean DISABLE_VOTLAGES_UPDATES = false;
    
    protected RealVector p, q, e, f, lp, lq;
    protected RealMatrix G, B;
    protected RealVector p_max;
    protected RealVector p_costMin;
    protected int slackIndex;
    protected RealVector generatorMask;
    protected Set<Agent> agents;
    protected Grid grid;
    protected RealVector c;
    protected RealVector t;
    protected double cpuTime = 0;
    
    protected void run()
    {
        AnalysisResults results = init();
        loop();
        testResults(results);
        finish();
    }
    
    protected void finish()
    {
        out.println();
        out.println(SimpleDateFormat.getDateTimeInstance().format((Calendar.getInstance().getTime())));
    }
    
    protected void testResults(AnalysisResults results)
    {
        // Set values to results:
        Map<String, Integer> busNumbers = results.getBusNumbers();
        for (String busName : busNumbers.keySet())
        {
            Bus bus = grid.getBus(busName);
            Source dg = null;
            for (Source source : bus.getSources())
            {
                dg = source;
                break;
            }
            if(dg != null)
            {
                int i = busNumbers.get(busName);
                dg.setPowerOutput(BASE_POWER*p(i), BASE_POWER*q(i));
            }
        }
        
        // Log solution:
        out.println("\nFinal Values:");
        out.println("index,bus,e,f,p,q,,|v|,|s|");
        for (String busName : busNumbers.keySet())
        {
            int i = busNumbers.get(busName);
            double f_i = f(i);
            double e_i = e(i);
            String v_abs = format(Math.sqrt(e_i*e_i+f_i*f_i));
            double p_i = p(i);
            double q_i = q(i);
            String s_abs = format(Math.sqrt(p_i*p_i+q_i*q_i));
            out.println(i+","+busName+','+format(e_i)+','+format(f_i)+','+format(p_i)+','+format(q_i)+",,"+v_abs+','+s_abs);
        }
        
        // Analyse grid to get parameters:
        LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        lfa.setIterations(100);
        lfa.setTargetError(1e-6);
        AnalysisResults results2 = lfa.analyse();//complexVoltages()); // WARNING Passing in voltages is changing the load powers somehow in the results
        
        out.println("\nAnalysis Results:");
//        GridDisplay.showInFrame(grid, grid.getBasePower(), grid.getBaseVoltage()).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        logBusses(results2);
    }
    
    
    //// Agents ////
    
    public class Agent
    {   
        public double p, q, e, f; // primal variables
        double p_new, q_new, e_new, f_new;
        public double p_costMin;
        public double p_max;
        public double lp, lq; // dual variables
        double lp_new, lq_new;
        public double gp, gq; // derived values
        protected double c;// = 1e-6;
        public boolean isSlack;
        public boolean isGenerator;
        public int index;
        protected Map<Agent, Complex> neighbours = new HashMap<>();
        public double G, B; // Self admittances.
        public double lastT = 1;
        
        public Agent(int i)
        {
            this.index = i;
        }
        
        public void addNeighbour(Agent neighbour, Complex admittance)
        {
            neighbours.put(neighbour, admittance);
        }

        protected boolean debug_backtrack = false;
        
        public void step(int k)
        {
            if(isSlack)
                return;
    
if(k >= 800*3 && index == 17)
{
    Util.nullop();
//  debug_backtrack = true;
}
            
            // Remember current constraint values:
            double gp_old = gp;//(p, q, e, f);
            double gq_old = gq;//(p, q, e, f);
            
            for(int i = 0; i < X_STEPS_PER_ITERATION; ++i)
            {
//                if(isGenerator)
//                {
                    stepX(PQ);
                    stepX(EF);
//                }
//                else
//                {
//                    stepX(EF);
//                }
            }
            
            stepLambda();
            stepAugScale(k, gp_old, gq_old);
        }

        Random rand = new Random(0);
        int stepCount = 0;
        public static final int PQ = 1;
        public static final int EF = 2;
        public static final int PQEF = 3;
        protected void stepX(int pqef)
        {
            // Find x step direction:
            double dL_p = isGenerator && (pqef&PQ)!=0 ? gradL_p() : 0;
            double dL_q = isGenerator && (pqef&PQ)!=0 ? gradL_q() : 0;
            double dL_e = !isSlack    && !DISABLE_VOTLAGES_UPDATES && (pqef&EF)!=0 ? gradL_e() : 0;
            double dL_f = !isSlack    && !DISABLE_VOTLAGES_UPDATES && (pqef&EF)!=0 ? gradL_f() : 0;

            RealMatrix Hx = USE_NEWTON ? hessian() : MatrixUtils.createRealIdentityMatrix(4);
            RealVector x = new ArrayRealVector(new double[]{p, q, e, f});
            RealVector dL_x = new ArrayRealVector(new double[]{dL_p, dL_q, dL_e, dL_f});
            RealVector dx = MatrixHelper.solve(Hx, dL_x);
            if(dx == null)
                dx = dL_x; // Revert to gradient decent.
            dx = dx.mapMultiply(-1); // delta x = -H^-1 * grad_x L(x,lambda)
            
            // Debug:
            if(debug_backtrack)
            {
                double cOld = Sandbox016.this.c.getEntry(index);
                Sandbox016.this.c.setEntry(index, c);
                
                compareGlobalLocal(dL_p, dL_q, dL_e, dL_f);
                System.out.println();
                checkGradient(dL_p, dL_q, dL_e, dL_f);
                System.out.println();
                plotImprovements();
                System.out.println();
                plotComparison();
                
                Sandbox016.this.c.setEntry(index, cOld);
            }
            
            // Backtrack to find best step size:
            double t = backTrack(dL_x, dx, p, q, e, f, INITIAL_STEP_SIZE);
            
            // Update state:
            RealVector x_new = x.add(dx.mapMultiply(t));
            p_new = x_new.getEntry(0);
            q_new = x_new.getEntry(1);
            e_new = x_new.getEntry(2);
            f_new = x_new.getEntry(3);
            
            // Nesterov�s Accelerated Method - momentum step:
            if(USE_NESTEROV) // FIXME causes null pointer, but only when not debugging
            {
                // Ref. Nesterov 1998.
                // beta = ||f'(x_k)||^2 / ||f'(x_{k+1})||^2:
                double dL_p2 = isGenerator && (pqef&PQ)!=0 ? gradL_p(p_new, q_new, e_new, f_new) : 0;
                double dL_q2 = isGenerator && (pqef&PQ)!=0 ? gradL_q(p_new, q_new, e_new, f_new) : 0;
                double dL_e2 = !isSlack    && (pqef&EF)!=0 ? gradL_e(p_new, q_new, e_new, f_new) : 0;
                double dL_f2 = !isSlack    && (pqef&EF)!=0 ? gradL_f(p_new, q_new, e_new, f_new) : 0;
                RealVector dx2 = vector(dL_p2, dL_q2, dL_e2, dL_f2);
                double nes_beta = normSquared(dx)/normSquared(dx2);
                
                // p = f'(x_{k+1}) - beta*f'(x_k):
                RealVector p = dx2.subtract(dx.mapMultiply(nes_beta));
                
                // x_{k+1} = x_new + gamma*p:
                double t2 = backTrack(dx2, dx2, p_new, q_new, e_new, f_new, INITIAL_STEP_SIZE);
                x_new = x_new.add(p.mapMultiply(t2));
                p_new = x_new.getEntry(0);
                q_new = x_new.getEntry(1);
                e_new = x_new.getEntry(2);
                f_new = x_new.getEntry(3);
            }
            
            // Check for NaN:
            if(x_new.isNaN())
            {
                System.err.println("x_new is NaN");
//                System.out.println();
                Util.nullop();
            }
            
            // Project onto constrained set X:
            if(PROJECT_X)
            {
                Complex ef = projectVoltage(e_new, f_new);
                e_new = ef.getReal();
                f_new = ef.getImaginary();
                if(isGenerator)
                {
                    Complex pq = projectPowers(p_new, q_new);
                    p_new = pq.getReal();
                    q_new = pq.getImaginary();
                }
            }
            
            // Update x:
            p = isGenerator ? p_new : p;
            q = isGenerator ? q_new : q;
            if(!DISABLE_VOTLAGES_UPDATES)
            {
                e = e_new;
                f = f_new;
            }
            
            // Update cost function values:
            gp = gp();
            gq = gq();
            
            // FIXME DO NOT COMMIT
//            if(gp < 1e-6 && gq < 1e-6)
//            {
////                lp = 0;
////                lq = 0;
////                c = INITIAL_AUG_SCALE;
//                lp *= 0.99;
//                lq *= 0.99;
//            }
            
            lastT = t;
            
            ++stepCount;
        }

        protected double backTrack(RealVector dL_x, RealVector dx, double p, double q, double e, double f, double t)
        {
            double alpha = 0.5; // Step update size.
            double beta = 0.5;
            double L_improvement;
            double step = alpha*dL_x.dotProduct(dx);
            double minImprovement;

            double dp = dx.getEntry(0);
            double dq = dx.getEntry(1);
            double de = dx.getEntry(2);
            double df = dx.getEntry(3);
            
            // Check for NaN:
            if(dx.isNaN())
                Util.nullop();
            
            do
            {
                double dp2 = dp*t;
                double dq2 = dq*t;
                double de2 = de*t;
                double df2 = df*t;
                
                // Project onto constrained set X:
                if(PROJECT_X)
                {
                    Complex ef = projectVoltage(e+de2, f+df2);
                    de2 = ef.getReal() - e;
                    df2 = ef.getImaginary() - f;
                    if(isGenerator)
                    {
                        Complex pq = projectPowers(p+dp2, q+dq2);
                        dp2 = pq.getReal() - p;
                        dq2 = pq.getImaginary() - q;
                    }
                }
                
                L_improvement = lagrangianImprovement(dp2, dq2, de2, df2);
                minImprovement = t*step;

                if(debug_backtrack)
                {
                    // Check improvement by L(x-t*dx)-L(x)
//                    compareImprovement(dp2, dq2, de2, df2, 1, L_improvement);
                    double lowerBound = t*step/alpha;
                    System.out.println(t+","+L_improvement+","+minImprovement+","+lowerBound);
                }
                
                t *= beta;
            } while(L_improvement > minImprovement && t > MIN_STEP_SIZE);// || debug_backtrack);
            
            // Undo last step update:
            t /= beta;
            
            // Limit step size:
            if(t < MIN_STEP_SIZE)
                t = MIN_STEP_SIZE;

            // Debug:
            if(debug_backtrack)
            {
                System.out.println("0,0,0,0");
                double _t = -1e-2;
                L_improvement = lagrangianImprovement(
                    /*_t,*/
                    dp*_t, 
                    dq*_t, 
                    de*_t,
                    df*_t);
                double lowerBound = _t*step/alpha;
                System.out.println(_t+","+L_improvement+","+(_t*step)+","+lowerBound);
            }

            return t;
        }
        
        public void showPVsE()
        {
            RealVector _p = new ArrayRealVector(Sandbox016.this.p);
            RealVector _q = Sandbox016.this.q;
            RealVector _e = new ArrayRealVector(Sandbox016.this.e);
            RealVector _f = Sandbox016.this.f;
            RealVector _lp = Sandbox016.this.lp;
            RealVector _lq = Sandbox016.this.lq;
            for(double _de = -10; _de <= 10; _de += 0.2)
            {
                _e.setEntry(index, Sandbox016.this.e.getEntry(index)+_de);
                for(double _dp = -10; _dp <= 10; _dp += 0.1)
                {
                    _p.setEntry(index, Sandbox016.this.p.getEntry(index)+_dp);
                    System.out.print(lagrangian(_p, _q, _e, _f, _lp, _lq));
                    System.out.print(',');
                }
                System.out.println(',');
            }
        }

        protected static final int P = 0;
        protected static final int Q = 1;
        protected static final int E = 2;
        protected static final int F = 3;
        protected RealMatrix hessian()
        {
            RealMatrix H = new Array2DRowRealMatrix(4, 4);
            H.setEntry(P, P, grad2L_pp());
            H.setEntry(P, Q, grad2L_pq());
            H.setEntry(P, E, grad2L_pe());
            H.setEntry(P, F, grad2L_pf());

            H.setEntry(Q, P, grad2L_qp());
            H.setEntry(Q, Q, grad2L_qq());
            H.setEntry(Q, E, grad2L_qe());
            H.setEntry(Q, F, grad2L_qf());

            H.setEntry(E, P, grad2L_ep());
            H.setEntry(E, Q, grad2L_eq());
            H.setEntry(E, E, grad2L_ee());
            H.setEntry(E, F, grad2L_ef());

            H.setEntry(F, P, grad2L_fp());
            H.setEntry(F, Q, grad2L_fq());
            H.setEntry(F, E, grad2L_fe());
            H.setEntry(F, F, grad2L_ff());
            
            return H;
        }

        protected double grad2L_pp()
        {
            return 1 + c;
        }

        protected double grad2L_pq()
        {
            return 0;
        }

        protected double grad2L_pe()
        {
            return -c*gradGp_e();
        }

        protected double grad2L_pf()
        {
            return -c*gradGp_f();
        }

        protected double grad2L_qp()
        {
            return 0;
        }

        protected double grad2L_qq()
        {
            return 1 + c;
        }

        protected double grad2L_qe()
        {
            return -c*gradGq_e();
        }

        protected double grad2L_qf()
        {
            return -c*gradGq_f();
        }

        protected double grad2L_ep()
        {
            return -c*sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) + 2*G*e;
        }

        protected double grad2L_eq()
        {
            return -c*sumA(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent)) - 2*B*e;
        }

        protected double grad2L_ee()
        {
            return
                    2*lp*G - 2*lq*B
                    +c*sumA(neighbours, agent -> 
                              gradGp_e(agent)*(agent.e*G(agent) + agent.f*B(agent))
                            + gradGq_e(agent)*(agent.f*G(agent) - agent.e*B(agent))
                        )
                    +c*gradGp_e()*sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
                    +c*(gradGp_e()*2*G*e + 2*G*gp())
                    +c*gradGq_e()*sumA(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent))
                    -c*(gradGq_e()*2*B*e + 2*B*gq())
                    ;
        }

        protected double gradGp_e(Agent agent)
        {
            return agent.e*G(agent) + agent.f*B(agent);
        }

        protected double gradGp_f(Agent agent)
        {
            return agent.f*G(agent) - agent.e*B(agent);
        }

        protected double gradGq_e(Agent agent)
        {
            return agent.f*G(agent) - agent.e*B(agent);
        }

        protected double gradGq_f(Agent agent)
        {
            return -agent.e*G(agent) - agent.f*B(agent);
        }

        protected double grad2L_ef()
        {
            return
                    c*sumA(neighbours, agent -> gradGp_f(agent)*(agent.e*G(agent) + agent.f*B(agent)) + gradGq_f(agent)*(agent.f*G(agent) - agent.e*B(agent)))
                    +c*gradGp_f()*(sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) + 2*G*e)
                    +c*gradGq_f()*(sumA(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent)) - 2*B*e);
        }

        protected double grad2L_fp()
        {
            return -c*sumA(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent)) + 2*G*f;
        }

        protected double grad2L_fq()
        {
            return -c*sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) - 2*B*f;
        }

        protected double grad2L_fe()
        {
            return
                    c*sumA(neighbours, agent -> gradGp_e(agent)*(agent.f*G(agent) - agent.e*B(agent)) + gradGq_e(agent)*(-agent.e*G(agent) - agent.f*B(agent)))
                    +c*gradGp_e()*(sumA(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent)) + 2*G*f)
                    +c*gradGq_e()*(sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) - 2*B*f);
        }

        protected double grad2L_ff()
        {
            return
                    2*lp*G - 2*lq*B
                    +c*sumA(neighbours, agent -> 
                              gradGp_f(agent)*(agent.f*G(agent) - agent.e*B(agent))
                            + gradGq_f(agent)*(-agent.e*G(agent) - agent.f*B(agent))
                        )
                    +c*gradGp_f()*sumA(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent))
                    +c*(gradGp_f()*2*G*f + 2*G*gp())
                    +c*gradGq_f()*sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
                    -c*(gradGq_f()*2*B*f + 2*B*gq())
                    ;
        }
        
        protected double gradGp_e()
        {
            return sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) + 2*G*e;
        }
        
        protected double gradGp_f()
        {
            return sumA(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent)) + 2*G*f;
        }
        
        protected double gradGq_e()
        {
            return sumA(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent)) - 2*B*e;
        }
        
        protected double gradGq_f()
        {
            return sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) - 2*B*f;
        }

        protected void stepLambda()
        {
            // Maximise (l := l - c*g(x), ref http://en.wikipedia.org/wiki/Augmented_Lagrangian_method
            //           and Bertsekas' book):
            double dlp = isSlack ? 0 : gp;
            lp = lp + dlp*c;
            double dlq = isSlack ? 0 : gq;
            lq = lq + dlq*c;
            
            lp = project(lp, -MAX_LAMBDA, MAX_LAMBDA);
            lq = project(lq, -MAX_LAMBDA, MAX_LAMBDA);
        }

        protected void stepAugScale(int k, double gp_old, double gq_old)
        {
            // Update augmentation scale:
            if(Math.abs(gp) > MIN_CONSTRAINT_IMPROVEMENT*Math.abs(gp_old) || Math.abs(gq) > MIN_CONSTRAINT_IMPROVEMENT*Math.abs(gq_old))
            {
                if(c < MAX_AUG_SCALE)
                    c *= AUG_SCALE_STEP;
            }
            
            // Experimental:
//            double gpImp = abs(gp) - MIN_CONSTRAINT_IMPROVEMENT*abs(gp_old);
//            double gqImp = abs(gq) - MIN_CONSTRAINT_IMPROVEMENT*abs(gq_old);
//            
//            averageGImprovement = (1-G_IMPROVEMENT_AVERAGE_RATE)*averageGImprovement + G_IMPROVEMENT_AVERAGE_RATE*(gpImp+gqImp)/2;
//            if(averageGImprovement >= 0)
//            {
//                if(c < MAX_AUG_SCALE)
//                    c *= AUG_SCALE_STEP;
//            }
        }

        protected Complex projectVoltage(double e, double f)
        {
            if(e < 0)
                e = 0;
            
            // Clamp |v|:
            Complex v = new Complex(e, f);
            double abs = v.abs();
            if(abs > V_MAX)
            {
                v = v.multiply(V_MAX/abs);
            }
            else if(abs < V_MIN)
            {
                v = v.multiply(V_MIN/abs);
            }
            
            // Clamp angle:
            abs = v.abs();
            double arg = v.getArgument();
            if(arg < V_ARG_MIN)
            {
                v = ComplexUtils.polar2Complex(abs, V_ARG_MIN);
            }
            else if(arg > V_ARG_MAX)
            {
                v = ComplexUtils.polar2Complex(abs, V_ARG_MAX);
            }
            
            return v;
        }

        protected Complex projectPowers(double p, double q)
        {
            if(p < 0)
                p = 0;
            
            Complex s = new Complex(p, q);
            double abs = s.abs();
            double p_max_i = p_max;
            if(abs > p_max_i)
            {
                s = s.multiply(p_max_i/abs);
            }
            
            return s;
        }

        /**
         * This is equivalent to L(X+tdx) - L(x) but is calculable within the agent's neighbourhood.
         * @param t
         * @param dp
         * @param dq
         * @param de
         * @param df
         * @return
         */
        private double lagrangianImprovement(/*double t, */double dp, double dq, double de, double df)
        {
            double dc = cost(p+/*t**/dp, q+/*t**/dq) - cost(p, q);
            
            double constraintTerms = lagrangianImprovement_constraintTerms(/*t, */dp, dq, de, df);
            double augTermP = lagrangianImprovement_augTermP(/*t, */dp, dq, de, df);
            double augTermQ = lagrangianImprovement_augTermQ(/*t, */dp, dq, de, df);
            double augTerm = augTermP + augTermQ;
            return 
                    dc
                    +constraintTerms
                    +augTerm;
        }

        public double lagrangianImprovement_augTermP(/*double t, */double dp, double dq, double de, double df)
        {
            return 0.5*c*sumA(neighbours, agent -> agent.isSlack ? 0 : dgp(agent, /*t, */dp, dq, de, df)*(dgp(agent, /*t, */dp, dq, de, df) + 2*agent.gp()))
                                +0.5*c*( dgpSelf(/*t, */dp, dq, de, df)*(dgpSelf(/*t, */dp, dq, de, df) + 2*gp()) );
        }

        public double lagrangianImprovement_augTermQ(/*double t, */double dp, double dq, double de, double df)
        {
            return 
                     0.5*c*sumA(neighbours, agent -> agent.isSlack ? 0 : dgq(agent, /*t, */dp, dq, de, df)*(dgq(agent, /*t, */dp, dq, de, df) + 2*agent.gq()))
                    +0.5*c*(dgqSelf(/*t, */dp, dq, de, df)*(dgqSelf(/*t, */dp, dq, de, df) + 2*gq()));
        }

        public double lagrangianImprovement_constraintTerms(/*double t,*/
                double dp, double dq, double de, double df)
        {
            return 
                    sumA(neighbours, agent -> agent.isSlack ? 0 : agent.lp*dgp(agent, /*t, */dp, dq, de, df))
                    +lp*dgpSelf(/*t, */dp, dq, de, df)
                    +sumA(neighbours, agent -> agent.isSlack ? 0 : agent.lq*dgq(agent, /*t, */dp, dq, de, df))
                    +lq*dgqSelf(/*t, */dp, dq, de, df);
        }

        private double dgpSelf(/*double t, */double dp, double dq, double de, double df)
        {
            return
                    /*t*t**/(de*de + df*df)*G - /*t**/dp
                    +/*t**/(e*de*G + f*df*G + f*de*B - e*df*B)
                    +/*t**/sumA(neighbours, agent -> de*agent.e*G(agent) + df*agent.f*G(agent) + df*agent.e*B(agent) - de*agent.f*B(agent))
                    +/*t**/(de*e*G + df*f*G + df*e*B - de*f*B); // This last term is needed since neighbours does not include this agent
        }

        private double dgqSelf(/*double t, */double dp, double dq, double de, double df)
        {
            return
                    -/*t*t**/(de*de + df*df)*B - /*t**/dq
                    +/*t**/(f*de*G - e*df*G - e*de*B - f*df*B)
                    +/*t**/sumA(neighbours, agent -> df*agent.e*G(agent) - de*agent.f*G(agent) - de*agent.e*B(agent) - df*agent.f*B(agent))
                    +/*t**/(df*e*G - de*f*G - de*e*B - df*f*B); // This last term is needed since neighbours does not include this agent
        }

        private double dgp(Agent agent, /*double t, */double dp, double dq, double de, double df)
        {
            double G_in = G(agent);
            double B_in = B(agent);
            return /*t**/(agent.e*de*G_in + agent.f*df*G_in + agent.f*de*B_in - agent.e*df*B_in);
        }

        private double dgq(Agent agent, /*double t, */double dp, double dq, double de, double df)
        {
            double G_in = G(agent);
            double B_in = B(agent);
            return /*t**/(agent.f*de*G_in - agent.e*df*G_in - agent.e*de*B_in - agent.f*df*B_in);
        }

        protected double project(double v, double low, double high)
        {
            return min(max(low, v), high);
        }

        protected double gp()
        {
            return gp(p, q, e, f);
        }
        
        /**
         * WARNING: This still uses neighbour values, only this agent's values are passed in.
         * @param p
         * @param q
         * @param e
         * @param f
         * @return
         */
        protected double gp(double p, double q, double e, double f)
        {
            double sum = 0;
            for (Agent agent : neighbours.keySet())
            {
                sum += e*agent.e*G(agent) + f*agent.f*G(agent) + f*agent.e*B(agent) - e*agent.f*B(agent);
            }
            sum += e*e*G + f*f*G; // self
            return  sum - p;
        }
        
        protected double gq()
        {
            return gq(p, q, e, f);
        }
        
        protected double gq(double p, double q, double e, double f)
        {
            double sum1 = sumA(neighbours, agent -> f*agent.e*G(agent) - e*agent.f*G(agent) - e*agent.e*B(agent) - f*agent.f*B(agent));
            double sum = sum1
                    - e*e*B - f*f*B;
            return 
                sum
                - q;
        }
        
        /**
         * c(x) = 0.5*(p - p_max)^2 + 0.5*q^2
         * @param p Real power output.
         * @param q Reactive power output.
         * @return
         */
        protected double cost(double p, double q)
        {
            double d = p-p_costMin;
            return COST_MULTIPLIER*(0.5*d*d + 0.5*q*q);
        }

        protected double gradL_p()
        {
            return gradC_p() - lp - c*gp();
        }
        
        protected double gradL_p(double p, double q, double e, double f)
        {
            return gradC_p(p, q, e, f) - lp - c*gp(p, q, e, f);
        }
        
        /**
         * Assuming C_i(p_i) = 0.5*(p_i - p_max)^2
         * => dC_i/dp_i = p_i - p_max
         * @return
         */
        protected double gradC_p()
        {
            return COST_MULTIPLIER*(p - p_costMin);
        }
        
        private double gradC_p(double p, double q, double e, double f)
        {
            return COST_MULTIPLIER*(p - p_costMin);
        }

        protected double gradL_q()
        {
            return gradC_q() - lq - c*gq();
        }

        protected double gradL_q(double p, double q, double e, double f)
        {
            return gradC_q(p, q, e, f) - lq - c*gq(p, q, e, f);
        }
        
        /**
         * Assuming C_i(q_i) = 0.5*q_i^2
         * => dC_i/dp_i = q_i
         * @param i
         * @return
         */
        protected double gradC_q()
        {
            return COST_MULTIPLIER*q;
        }
        
        private double gradC_q(double p, double q, double e, double f)
        {
            return COST_MULTIPLIER*q;
        }

        protected double gradL_e()
        {
            return gradL_e(p, q, e, f);
        }
        
        protected double gradL_e(double p, double q, double e, double f)
        {
            double constraintTerms = 
                         sumA(neighbours, agent -> agent.lp*(agent.e*G(agent) + agent.f*B(agent)))
                        +sumA(neighbours, agent -> agent.lq*(agent.f*G(agent) - agent.e*B(agent)))
                        +lp*(sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
                             +2*G*e)
                        +lq*(sumA(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent))
                             -2*B*e);

            double pt1 = c*sumA(neighbours, agent -> agent.isSlack ? 0 : 
                                 gp(agent)*(agent.e*G(agent) + agent.f*B(agent))
                                +gq(agent)*(agent.f*G(agent) - agent.e*B(agent)));
            double pt2 = c*gp()*(sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) 
                               +2*e*G);
            double pt3 = c*gq()*(sumA(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent)) 
                               -2*e*B);
            double penaltyTerms = pt1
                                +
                                (isSlack ? 0 : 
                                (
                                    +pt2
                                    +pt3
                                ));
            
            return
                constraintTerms
                +penaltyTerms
                ;
        }

        protected double gq(Agent agent)
        {
            return agent.gq;
        }

        protected double gp(Agent agent)
        {
            return agent.gp;
        }

        protected double gradL_f()
        {
            return gradL_f(p, q, e, f);
        }
        
        protected double gradL_f(double p, double q, double e, double f)
        {
            double constraintTerms = 
                         sumA(neighbours, agent -> agent.lp*(agent.f*G(agent) - agent.e*B(agent)))
                        +sumA(neighbours, agent -> agent.lq*(-agent.e*G(agent) - agent.f*B(agent)))
                        +lp*(sumA(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent))
                             +2*G*f)
                        +lq*(sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
                             -2*B*f);

            double pt1 = c*sumA(neighbours, agent -> agent.isSlack ? 0 : 
                               gp(agent)*(agent.f*G(agent) - agent.e*B(agent))
                               +gq(agent)*(-agent.e*G(agent) - agent.f*B(agent)));
            double pt2 = c*gp()*(sumA(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent)) 
                               +2*f*G);
            double pt3 = c*gq()*(sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) 
                               -2*f*B);
            double penaltyTerms = pt1
                                +
                                (isSlack ? 0 : 
                                (
                                    +pt2
                                    +pt3
                                ));

            return
                constraintTerms
                +penaltyTerms
                ;
        }

        protected double G(Agent agent)
        {
            double G_ij = neighbours.get(agent).getReal();
            return G_ij;
        }

        protected double B(Agent agent)
        {
            return neighbours.get(agent).getImaginary();
        }

        public void init()
        {
            gp = gp();
            gq = gq();
            
            // Find unconstrained optimum of x as the starting point:
            if(isGenerator)
            {
                p = p_costMin;
                q = 0;
            }
        }
        
        public String toString()
        {
            StringBuffer buf = new StringBuffer();
            buf.append("Agent ");
            buf.append(index);
            buf.append(": \n");
            buf.append("\tp=");
            buf.append(p);
            buf.append("\n\tq=");
            buf.append(q);
            buf.append("\n\te=");
            buf.append(e);
            buf.append("\n\tf=");
            buf.append(f);
            buf.append("\n\tlp=");
            buf.append(lp);
            buf.append("\n\tlq=");
            buf.append(lq);
            buf.append("\n\tc=");
            buf.append(c);
            buf.append("\n\tt=");
            buf.append(lastT);
            buf.append("\n");
                
            return buf.toString();
        }
        
        
        //// Debug ////
        
        public Object showNeighbours = new Object()
        {
            public String toString()
            {
                StringBuffer buf = new StringBuffer();

                for (Agent n : neighbours.keySet())
                {
                    buf.append(n);
                    buf.append("\tG_ij=");
                    buf.append(neighbours.get(n).getReal());
                    buf.append("\n\tB_ij=");
                    buf.append(neighbours.get(n).getImaginary());
                    buf.append("\n\n");
                }
                
                return buf.toString();
            }
        };

        protected void compareGlobalLocal(double dp, double dq, double de,
                double df)
        {
            System.out.println("Compare Global/Local:");
            
            double globalDp = Sandbox016.this.gradL_p(index);
            System.out.println(globalDp + ",?=," + dp);

            double globalDq = Sandbox016.this.gradL_q(index);
            System.out.println(globalDq + ",?=," + dq);
            
            double globalDe = Sandbox016.this.gradL_e(index);
            System.out.println(globalDe + ",?=," + de);

            double globalDf = Sandbox016.this.gradL_f(index);
            System.out.println(globalDf + ",?=," + df);
        }

        protected void compareImprovement(double dp, double dq, double de,
                double df, double t, double L_improvement)
        {
            RealVector _p = new ArrayRealVector(Sandbox016.this.p);
            RealVector _q = new ArrayRealVector(Sandbox016.this.q);
            RealVector _e = new ArrayRealVector(Sandbox016.this.e);
            RealVector _f = new ArrayRealVector(Sandbox016.this.f);
            double Lx = lagrangian(_p, _q, _e, _f, Sandbox016.this.lp, Sandbox016.this.lq);
            _p.setEntry(index, p-t*dp);
            _q.setEntry(index, q-t*dq);
            _e.setEntry(index, e-t*de);
            _f.setEntry(index, f-t*df);
            double Ldx = lagrangian(_p, _q, _e, _f, Sandbox016.this.lp, Sandbox016.this.lq);
            System.out.println("L(x-t*dx)-L(x) =, "+(Ldx-Lx)+", L_imp =, "+L_improvement);
        }

        protected void plotComparison()
        {
            System.out.println(
                    "Comparison:\ndx"
//                    + ",imp_c(x)_e,imp_c(x)_f,c(x+de)-c(x),c(x_df)-c(x)"
                    + ",imp_g(x)_e,imp_g(x)_f,g(x+de)-g(x),g(x_df)-g(x)"
                    + ",imp_augP(x)_e,imp_augP(x)_f,augP(x+de)-augP(x),augP(x_df)-augP(x)"
                    + ",imp_augQ(x)_e,imp_augQ(x)_f,augQ(x+de)-augQ(x),augQ(x_df)-augQ(x)"
                    + ",imp_sum(x)_e,imp_sum(x)_f,L_sum(x+de)-L_sum(x),L_sum(x_df)-L_sum(x)"
                    + ",imp(x)_e,imp(x)_f,L(x+de)-L(x),L(x_df)-L(x)");
            for(double d = -1e-1; d <= 1.0000001e-1; d += 1e-2)
            {
                // Local calcs:
                double _de_cons = lagrangianImprovement_constraintTerms(/*1, */0, 0, -d, 0);
                double _df_cons = lagrangianImprovement_constraintTerms(/*1, */0, 0, 0, -d);
                
                double _de_augP = lagrangianImprovement_augTermP(/*1, */0, 0, -d, 0);
                double _df_augP = lagrangianImprovement_augTermP(/*1, */0, 0, 0, -d);
                
                double _de_augQ = lagrangianImprovement_augTermQ(/*1, */0, 0, -d, 0);
                double _df_augQ = lagrangianImprovement_augTermQ(/*1, */0, 0, 0, -d);

                double _de = lagrangianImprovement(/*1, */0, 0, -d, 0);
                double _df = lagrangianImprovement(/*1, */0, 0, 0, -d);
                
                // Global calcs:
                RealVector _p = new ArrayRealVector(Sandbox016.this.p);
                RealVector _q = new ArrayRealVector(Sandbox016.this.q);
                RealVector _e = new ArrayRealVector(Sandbox016.this.e);
                RealVector _f = new ArrayRealVector(Sandbox016.this.f);
                RealVector gp_global = Sandbox016.this.gp(_p, _q, _e, _f);
                RealVector gq_global = Sandbox016.this.gq(_p, _q, _e, _f);
                gp_global.setEntry(slackIndex, 0);
                gq_global.setEntry(slackIndex, 0);
                double Lx_augP = 0.5*c*gp_global.dotProduct(gp_global);
                double Lx_augQ = 0.5*c*gq_global.dotProduct(gq_global);
                double Lx_lg = Sandbox016.this.lp.dotProduct(gp_global) + Sandbox016.this.lq.dotProduct(gq_global);
                double Lx = lagrangian(_p, _q, _e, _f, Sandbox016.this.lp, Sandbox016.this.lq);

                _e.setEntry(index, e-d);
                gp_global = Sandbox016.this.gp(_p, _q, _e, _f);
                gq_global = Sandbox016.this.gq(_p, _q, _e, _f);
                gp_global.setEntry(slackIndex, 0);
                gq_global.setEntry(slackIndex, 0);
                double Lx_augP_e = 0.5*c*gp_global.dotProduct(gp_global);
                double Lx_augQ_e = 0.5*c*gq_global.dotProduct(gq_global);
                double Lx_lg_e = Sandbox016.this.lp.dotProduct(gp_global) + Sandbox016.this.lq.dotProduct(gq_global);
                double Lx_e = lagrangian(_p, _q, _e, _f, Sandbox016.this.lp, Sandbox016.this.lq);
                
                _e.setEntry(index, e);
                _f.setEntry(index, f-d);
                gp_global = Sandbox016.this.gp(_p, _q, _e, _f);
                gq_global = Sandbox016.this.gq(_p, _q, _e, _f);
                gp_global.setEntry(slackIndex, 0);
                gq_global.setEntry(slackIndex, 0);
                double Lx_augP_f = 0.5*c*gp_global.dotProduct(gp_global);
                double Lx_augQ_f = 0.5*c*gq_global.dotProduct(gq_global);
                double Lx_lg_f = Sandbox016.this.lp.dotProduct(gp_global) + Sandbox016.this.lq.dotProduct(gq_global);
                double Lx_f = lagrangian(_p, _q, _e, _f, Sandbox016.this.lp, Sandbox016.this.lq);
                
                System.out.println(d+","+_de_cons+","+_df_cons+","+(Lx_lg_e-Lx_lg)+","+(Lx_lg_f-Lx_lg)
                        +","+_de_augP+","+_df_augP+","+(Lx_augP_e-Lx_augP)+","+(Lx_augP_f-Lx_augP)
                        +","+_de_augQ+","+_df_augQ+","+(Lx_augQ_e-Lx_augQ)+","+(Lx_augQ_f-Lx_augQ)
                        +","+(_de_cons+_de_augP+_de_augQ)+","+(_df_cons+_df_augP+_df_augQ)+","
                            +((Lx_lg_e+Lx_augP_e+Lx_augQ_e)-(Lx_lg+Lx_augP+Lx_augQ))+","+((Lx_lg_f+Lx_augP_f+Lx_augQ_f)-(Lx_lg+Lx_augP+Lx_augQ))
                        +","+_de+","+_df+","+(Lx_e-Lx)+","+(Lx_f-Lx));
            }
        }

        protected void plotImprovements()
        {
            System.out.println("Improvements:\ndx,imp_e,imp_f,L(x+de)-L(x),L(x_df)-L(x)");
            for(double d = -1e-12; d <= 1.0000001e-12; d += 1e-13)
            {
                double _de = lagrangianImprovement(/*1, */0, 0, -d, 0);
                double _df = lagrangianImprovement(/*1, */0, 0, 0, -d);
                
                RealVector _p = new ArrayRealVector(Sandbox016.this.p);
                RealVector _q = new ArrayRealVector(Sandbox016.this.q);
                RealVector _e = new ArrayRealVector(Sandbox016.this.e);
                RealVector _f = new ArrayRealVector(Sandbox016.this.f);
                double Lx = lagrangian(_p, _q, _e, _f, Sandbox016.this.lp, Sandbox016.this.lq);
                _e.setEntry(index, e-d);
                double Ldx_e = lagrangian(_p, _q, _e, _f, Sandbox016.this.lp, Sandbox016.this.lq);
                
                _e.setEntry(index, e);
                _f.setEntry(index, f-d);
                double Ldx_f = lagrangian(_p, _q, _e, _f, Sandbox016.this.lp, Sandbox016.this.lq);
                
                System.out.println(d+","+_de+","+_df+","+(Ldx_e-Lx)+","+(Ldx_f-Lx));
            }
            
//          for(double d = -1e-12; d < 1e-12; d += 1e-13)
//          {
//              double _imp = lagrangianImprovement(1, 0, 0, -d, -d);
//              RealVector _p = new ArrayRealVector(Sandbox016.this.p);
//              RealVector _q = new ArrayRealVector(Sandbox016.this.q);
//              RealVector _e = new ArrayRealVector(Sandbox016.this.e);
//              RealVector _f = new ArrayRealVector(Sandbox016.this.f);
//              double Lx = lagrangian(_p, _q, _e, _f, Sandbox016.this.lp, Sandbox016.this.lq);
//              _e.setEntry(index, e-d);
//              _f.setEntry(index, f-d);
//              double Ldx = lagrangian(_p, _q, _e, _f, Sandbox016.this.lp, Sandbox016.this.lq);
//              out.println(d+","+_imp+","+(Ldx-Lx));
//          }
            
            Util.nullop();
        }

        protected void checkGradient(double dp, double dq, double de, double df)
        {
            double delta = 1e-24;
            
            System.out.println("Check Gradient:\nindex,dp,,dq,,de,,df");
            
            // Estimate p slope:
            double left = lagrangianImprovement(/*1, */-delta, 0, 0, 0);
            double right = lagrangianImprovement(/*1, */delta, 0, 0, 0);
            double pSlope = (right-left)/(2*delta);

            // Estimate q slope:
            left = lagrangianImprovement(/*1, */0, -delta, 0, 0);
            right = lagrangianImprovement(/*1, */0, delta, 0, 0);
            double qSlope = (right-left)/(2*delta);
            
            // Estimate e slope:
            left = lagrangianImprovement(/*1, */0, 0, -delta, 0);
            right = lagrangianImprovement(/*1, */0, 0, delta, 0);
            double eSlope = (right-left)/(2*delta);

            // Estimate f slope:
            left = lagrangianImprovement(/*1, */0, 0, 0, -delta);
            right = lagrangianImprovement(/*1, */0, 0, 0, delta);
            double fSlope = (right-left)/(2*delta);
            
            // Log:
            System.out.println(index
                    +","+dp+","+pSlope
                    +","+dq+","+qSlope
                    +","+de+","+eSlope
                    +","+df+","+fSlope);
        }
    } // End Agent
    
    
    //// Setup ////

    protected AnalysisResults init()
    {
        initGrid();
        
        // Analyse grid to get parameters:
        LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        lfa.setIterations(40);
//        lfa.setTargetError();
        AnalysisResults results = lfa.analyse();
        if(!results.getDidConverge())
            throw new RuntimeException("Initial analysis did not converge.");
        
        logBusses(results);
        
        // Initialise variables:
        initVariables(results);
        
        // Create agents:
        initAgents();
        
        // Debug:
        debugHeader(results);
        
        return results;
    }
    
    protected void initGrid()
    {
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox016-lineConfig");
        InputStream lineData = getClass().getResourceAsStream("Sandbox016-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox016-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream("Sandbox016-loadData.csv");
        grid = GridGenerator.loadGrid("IEEE123", lineConfig, lineData, switchData, loadData);
        
        // Add slack bus:
        SlackSource slack = new SlackSource();
        slack.setVoltage(new Complex(BASE_VOLTAGE*SLACK_VOLTAGE));
        slack.setName(SLACK_SOURCE);
        grid.getBus(GridGenerator.BUS_PREFIX+149).addChild(slack);
        
        // Transformer:
//        SlackSource slack = new SlackSource();
//        slack.setVoltage(new Complex(BASE_VOLTAGE));
//        slack.setName(SLACK_SOURCE);
//        
//        Bus slackBus = new Bus(grid);
//        slackBus.setName(SLACK_BUS);
//        slackBus.addChild(slack);
//        
//        Transformer slackXFMR = new Transformer();
//        slackXFMR.setName(SLACK_BUS+GridGenerator.LINE_JOIN+149);
//        slackXFMR.setLength(1);
//        Complex impedancePerMetre = grid.getLine(149+GridGenerator.LINE_JOIN+1).impedancePerMetre();
//        slackXFMR.setImpedencePerMetre(impedancePerMetre);
//        slackXFMR.setImpedance(slackXFMR.impedance()); // Slight mismatch between lines and transformers here
//        slackXFMR.setMinimumT(0.9);
//        slackXFMR.setMaximumT(1.1);
//        slackXFMR.setT(1.01);
//        slackXFMR.setBaseVoltageRatio(1);
//        slackXFMR.setFromBus(slackBus);
//        slackXFMR.setToBus(grid.getBus(GridGenerator.BUS_PREFIX+149));
        
        grid.setBaseVoltage(BASE_VOLTAGE);
        grid.setBasePower(BASE_POWER);
        
        // Add DG:
//        addDG(1,  500e3, 0, 0);
//        addDG(3,  500e3, 0, 0);
//        addDG(8,  500e3, 0, 0);
//        addDG(14, 500e3, 0, 0);
//        addDG(18, 500e3, 0, 0);
//        addDG(26, 500e3, 0, 0);
//        addDG(29, 500e3, 0, 0);
        addDG(1, 500e3, 100e3, 50e3);
        addDG(3, 500e3, 100e3, 50e3);
        addDG(8, 500e3, 100e3, 100e3);
        addDG(14, 500e3, 100e3, 50e3);
        addDG(18, 500e3, 100e3, 100e3);
        addDG(26, 500e3, 100e3, 50e3);
        addDG(29, 500e3, 100e3, 50e3);
        
//      GridDisplay.showInFrame(grid).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    /**
     * Creates and adds a DG, and removes any loads.
     * @param i
     * @param pMax
     * @param p
     * @param q
     */
    protected void addDG(int i, double pMax, double p, double q)
    {
        Bus bus = grid.getBus(GridGenerator.BUS_PREFIX+i);
        for (Load load : bus.getLoads())
        {
            bus.removeChild(load);
        }
        bus.addChild(makeDG("DG"+i, pMax, p, q));
    }

    public static DistributedSource makeDG(String name, double pMax, double p, double q)
    {
        DistributedSource dg = new DistributedSource();
        dg.setName(name);
        dg.setPmax(pMax);
        dg.setPowerOutput(p, q);
        return dg;
    }

    protected void initVariables(AnalysisResults results)
    {
        slackIndex = results.getSlackIndex();
        
        // Get admittances:
        FieldMatrix<Complex> Y = results.getAdmittanceMatrix().Y;
        int dimension = Y.getRowDimension();
        B = new Array2DRowRealMatrix(dimension, dimension);
        G = new Array2DRowRealMatrix(dimension, dimension);
        for(int i = 0; i < dimension; ++i)
        {
            for(int j = 0; j < dimension; ++j)
            {
                Complex Y_ij = Y.getEntry(i, j);
                G.setEntry(i, j, Y_ij.getReal());
                B.setEntry(i, j, Y_ij.getImaginary());
            }
        }
        
        // Prepare generator mask (1 if generator, 0 otherwise):
        generatorMask = new ArrayRealVector(dimension);
        Map<String, Integer> busNumbers = results.getBusNumbers();
        for (String name : busNumbers.keySet())
        {
            int index = busNumbers.get(name);
            Bus bus = grid.getBus(name);
            Collection<Source> sources = bus.getSources();
            if(!sources.isEmpty() && sources.iterator().next() instanceof DistributedSource)
            //if(!bus.getGeneratedPower().equals(Complex.ZERO))
                generatorMask.setEntry(index, 1);
        }
        
        // Init power vector:
        p = p0(results);
        q = q0(results);
        p.setEntry(slackIndex, 0);
        q.setEntry(slackIndex, 0);
        
        p_max = maxPower(results, grid);
        p_costMin = new ArrayRealVector(dimension, 0.0);
        
        // Init voltage vectors:
        e = new ArrayRealVector(dimension, SLACK_VOLTAGE);
        f = new ArrayRealVector(dimension, 0);
        
        if(START_WITH_TRUE_VOTLAGES)
        {
            for (String name : busNumbers.keySet())
            {
                Complex v = results.getBusVoltage(name);
                int i = busNumbers.get(name);
                e.setEntry(i, v.getReal());
                f.setEntry(i, v.getImaginary());
            }
        }
        
        // Init Lagrange multipliers:
        lp = new ArrayRealVector(dimension, 0);
        lq = new ArrayRealVector(dimension, 0);
        
        // Init debug variables:
        c = new ArrayRealVector(dimension, INITIAL_AUG_SCALE);
        t = new ArrayRealVector(dimension, 0);

        // FIXME DO NOT COMMIT
//        p = p.add(generatorMask.ebeMultiply(new ArrayRealVector(p.getDimension(), 10*Math.random())));
//        q = q.add(generatorMask.ebeMultiply(new ArrayRealVector(p.getDimension(), 10*Math.random())));
//        e = new ArrayRealVector(p.getDimension(), Math.random());
//        f = new ArrayRealVector(p.getDimension(), Math.random());
    }

    public static RealVector p0(AnalysisResults results)
    {
        Map<String, Integer> numbers = results.getBusNumbers();
        double[] p0 = new double[numbers.size()];
        
        for (String name : numbers.keySet())
        {
            int i = numbers.get(name);
            p0[i] = results.getBusPower(name).getReal();
        }
        
        return new ArrayRealVector(p0);
    }

    public static RealVector q0(AnalysisResults results)
    {
        Map<String, Integer> numbers = results.getBusNumbers();
        double[] q0 = new double[numbers.size()];
        
        for (String name : numbers.keySet())
        {
            int i = numbers.get(name);
            q0[i] = results.getBusPower(name).getImaginary();
        }
        
        return new ArrayRealVector(q0);
    }

    public static RealVector maxPower(AnalysisResults results, Grid grid)
    {
        Map<String, Integer> numbers = results.getBusNumbers();
        RealVector powers = new ArrayRealVector(numbers.size());
        for (String name : numbers.keySet())
        {
            int index = numbers.get(name);
            Bus bus = grid.getBus(name);
            Collection<Source> sources = bus.getSources();
            if(sources.size() > 0)
            {
                Source s = sources.iterator().next();
                powers.setEntry(index, s.getPmax()/grid.getBasePower());
            }
        }
        return powers;
    }
    
    protected void initAgents()
    {
        // Make agents:
        agents = new LinkedHashSet<>();
        int dimension = p.getDimension();
        for(int i = 0; i < dimension; ++i)
        {
            Agent agent = new Agent(i);
            agent.c = this.c.getEntry(i);
            agent.p = p(i);
            agent.q = q(i);
            agent.e = e(i);
            agent.f = f(i);
            agent.p_max = p_max.getEntry(i);
            agent.p_costMin = p_costMin.getEntry(i);
            agent.lp = lp(i);
            agent.lq = lq(i);
            agent.isGenerator = generatorMask.getEntry(i) == 1;
            agent.isSlack = i == slackIndex;
            agents.add(agent);
        }
        
        // Setup Neighbours:
        for (Agent agent_i : agents)
        {
            int i = agent_i.index;
            for (Agent agent_j : agents)
            {
                int j = agent_j.index;
                Complex Y_ij = new Complex(G(i,j), B(i,j));
                if(!Y_ij.equals(Complex.ZERO))
                {
                    if(i == j)
                    {
                        agent_i.G = Y_ij.getReal();
                        agent_i.B = Y_ij.getImaginary();
                    }
                    else
                    {
                        agent_i.addNeighbour(agent_j, Y_ij);
                    }
                }
            }
        }
        
        for (Agent agent : agents)
            agent.init();
    }
    
    
    //// Iterations ////

    public void loop()
    {
        Timer timer = Timer.startNewTimer();
        for(int k = 0; k < K; ++k)
        {
            for (Agent agent : agents)
            {
                agent.step(k);
                update();
            }
            
            cpuTime = timer.lapNano()/1e6/agents.size();

            // Log:
            if(k%DEBUG_RATE == 0)
                debug(k);
        }
    }
    
    /**
     * Update global variables from agent values.
     */
    protected void update()
    {
        for (Agent agent : agents)
        {
            int i = agent.index;
            p.setEntry(i, agent.p);
            q.setEntry(i, agent.q);
            e.setEntry(i, agent.e);
            f.setEntry(i, agent.f);
            lp.setEntry(i, agent.lp);
            lq.setEntry(i, agent.lq);

            c.setEntry(i, agent.c);
            t.setEntry(i, agent.lastT);
        }
        
        if(p.isNaN() || q.isNaN() || e.isNaN() || f.isNaN())
            Util.nullop();
    }
    
    
    //// Convenience methods ////

    public static String format(double d)
    {
        if(!FORMAT_NUMBERS)
            return ""+d;
        
        double abs = Math.abs(d);
        if(abs < 1e-24)
            return "0.0";
        else if(abs > 1e-3 && abs < 1e6)
            return String.format("%.4f", d);
        else
            return String.format("%.4e", d);
    }

    public static double norm(double... ds)
    {
        return Math.sqrt(sum(ds, d -> d*d));
    }

    protected double p(int i)
    {
        return p.getEntry(i);
    }

    protected double q(int i)
    {
        return q.getEntry(i);
    }

    protected double lp(int i)
    {
        return lp.getEntry(i);
    }

    protected double lq(int i)
    {
        return lq.getEntry(i);
    }

    protected double B(int i, int j)
    {
        return B.getEntry(i, j);
    }

    protected double G(int i, int j)
    {
        return G.getEntry(i, j);
    }

    protected double e(int i)
    {
        return e.getEntry(i);
    }

    protected double f(int i)
    {
        return f.getEntry(i);
    }
    
    
    //// Debug ////

    protected void debugHeader(AnalysisResults results)
    {
        Map<String, Integer> numbers = results.getBusNumbers();
        out.print("k,");
        
        int dimension = p.getDimension();
        debugLogNames(numbers, dimension, "p");
        debugLogNames(numbers, dimension, "q");
        debugLogNames(numbers, dimension, "e");
        debugLogNames(numbers, dimension, "f");
        debugLogNames(numbers, dimension, "l_p");
        debugLogNames(numbers, dimension, "l_q");
        debugLogNames(numbers, dimension, "gp");
        debugLogNames(numbers, dimension, "gq");
        debugLogNames(numbers, dimension, "gradL_p");
        debugLogNames(numbers, dimension, "gradL_q");
        debugLogNames(numbers, dimension, "gradL_e");
        debugLogNames(numbers, dimension, "gradL_f");
        out.print("C(pq),L(pql),");
        debugLogNames(numbers, dimension, "t");
        debugLogNames(numbers, dimension, "c");
        debugLogNames(numbers, dimension, "inf");
        out.print("||grad_p||,||grad_q||,||grad_e||,||grad_f||,CPU Time,beta_p_max,betakg");
        out.println();
    }

    public static void debugLogNames(Map<String, Integer> numbers, int dimension, String type)
    {
        for(int i = 0; i < dimension; ++i)
        {
            for (String name : numbers.keySet())
            {
                int index = numbers.get(name);
                if(index == i)
                {
                    out.print(name);
                    out.print('(');
                    out.print(type);
                    out.print("),");
                    break;
                }
            }
        }
    }
    
    /**
     * 
     * @param k Iteration counter.
     * @param t 
     */
    protected void debug(int k)
    {
        out.print(k);
        out.print(',');
        int dimension = p.getDimension();
        
        // Power:
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(p(i)));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(q(i)));
            out.print(',');
        }
        
        // Voltage:
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(e(i)));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(f(i)));
            out.print(',');
        }
        
        // Lambda:
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(lp(i)));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(lq.getEntry(i)));
            out.print(',');
        }
        
        // Calculated Power Error:
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gp(i)));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gq(i)));
            out.print(',');
        }
        
        // Power derivatives:
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gradL_p(i)));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gradL_q(i)));
            out.print(',');
        }
        
        // Voltage derivatives:
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gradL_e(i)));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gradL_f(i)));
            out.print(',');
        }
        
        // Slack power:
//      out.print(format(pSlack()));
//      out.print(',');
//      out.print(format(qSlack()));
//      out.print(',');
        
        // Costs:
        out.print(format(cost()));
        out.print(',');
        out.print(format(L()));
        out.print(',');

        // Step size:
        for(int i = 0; i < dimension; ++i)
        {
            out.print(t.getEntry(i));
            out.print(',');
        }

        // Augmentation scale:
        for(int i = 0; i < dimension; ++i)
        {
            out.print(c.getEntry(i));
            out.print(',');
        }

        // Controllability:
        for(int i = 0; i < dimension; ++i)
        {
            Agent agent = findAgent(i);
            out.print(format(controllability(agent)));
            out.print(',');
        }
        
        out.print(format(gradLNorm(this::gradL_p, dimension)));
        out.print(',');
        out.print(format(gradLNorm(this::gradL_q, dimension)));
        out.print(',');
        out.print(format(gradLNorm(this::gradL_e, dimension)));
        out.print(',');
        out.print(format(gradLNorm(this::gradL_f, dimension)));
        out.print(',');
        out.print(format(cpuTime));
        out.print(',');
        out.print(betaMax(k, dimension));
        out.print(',');
        out.print(betakg(k, dimension));
        out.print(',');
        
        out.println();
    }

    private double betakg(int k, int dimension)
    {
        double gnorm = norm(this::gp, dimension);
        return Math.pow(AUG_SCALE_STEP, k)*gnorm;
    }

    final double a = 1;
    final double r = 0.5;
    private double betaMax(int k, int dimension)
    {
        double gnorm = norm(this::gp, dimension);
        return Math.pow(a/gnorm, 1.0/k)*r;
    }

    private double norm(IndexedFunction f, int dimension)
    {
        RealVector v = vector(dimension, f);
        
        for(int i = 0; i < dimension; ++i)
        {
            if(i == slackIndex)
                v.setEntry(i, 0);
        }
        return v.getNorm();
    }

    private double gradLNorm(IndexedFunction gradFn, int dimension)
    {
         RealVector v = vector(dimension, gradFn);
        
        for(int i = 0; i < dimension; ++i)
        {
            if(i == 17 || i == slackIndex)
                v.setEntry(i, 0);
        }
        return v.getNorm();
    }

    protected Agent findAgent(int i)
    {
        for (Agent agent : agents)
        {
            if(agent.index == i)
                return agent;
        }
        return null;
    }

    protected double controllability(Agent agent)
    {
        double gradCNorm = norm(agent.gradC_p(), agent.gradC_q());
        double gradGNorm = norm(agent.gradGp_e(), agent.gradGp_f(), agent.gradGq_e(), agent.gradGq_f());
        double gNorm = norm(agent.gp(), agent.gq());
        double lNorm = norm(agent.lp, agent.lq);
        return gradCNorm/(gradGNorm*(lNorm+agent.c*gNorm));
    }

    protected double gradL_p(int i)
    {
        if(generatorMask.getEntry(i) == 1)
            return  gradC_p(i) - lp(i) - c.getEntry(i)*gp(i);
        else
            return 0;
    }
    
    /**
     * Assuming F_i(p_i) = 0.5*(p_i - p_max)^2
     * => dF_i/dp_i = p_i - p_max
     * @param i
     * @return
     */
    protected double gradC_p(int i)
    {
        return gradC_p(p, q, i);
    }
    
    protected double gradC_p(RealVector p, RealVector q, int i)
    {
        return COST_MULTIPLIER*(p.getEntry(i) - p_costMin.getEntry(i));
    }

    protected double gradL_q(int i)
    {
        if(generatorMask.getEntry(i) == 1)
            return gradC_q(i) - lq(i) - c.getEntry(i)*gq(i);
        else
            return 0;
    }
    
    /**
     * Assuming F_i(q_i) = 0.5*q_i^2
     * => dF_i/dp_i = q_i
     * @param i
     * @return
     */
    protected double gradC_q(int i)
    {
        return gradC_q(p, q, i);
    }
    
    protected double gradC_q(RealVector p, RealVector q, int i)
    {
        return COST_MULTIPLIER*q.getEntry(i);
    }
    
    protected double gp(int i)
    {
        int dimension = p.getDimension();
        double d = sum(n -> e(i)*e(n)*G(i,n) + f(i)*f(n)*G(i,n) + f(i)*e(n)*B(i,n) - e(i)*f(n)*B(i,n), dimension) - p(i);
        return d;
    }
    
    protected double gq(int i)
    {
        int dimension = q.getDimension();
        double d = sum(n -> f(i)*e(n)*G(i,n) - e(i)*f(n)*G(i,n) - e(i)*e(n)*B(i,n) - f(i)*f(n)*B(i,n), dimension) - q(i);
        return d;
    }

    protected double gradL_e(int j)
    {
        int dimension = p.getDimension();
        double constraintTerms = 
                 sum(i -> lp(i)*(e(i)*G(i,j) + f(i)*B(i,j)), dimension, j) // i != j
                +sum(i -> lq(i)*(f(i)*G(i,j) - e(i)*B(i,j)), dimension, j) // i != j
                +lp(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
                        +2*G(j,j)*e(j))
                +lq(j)*(sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) // n != j
                        -2*B(j,j)*e(j));
        
        double pt1 = c.getEntry(j)*sum(i -> i == slackIndex ? 0 : 
                             gp(i)*(e(i)*G(i,j) + f(i)*B(i,j))
                            +gq(i)*(f(i)*G(i,j) - e(i)*B(i,j)),
                            dimension, j);
        double pt2 = c.getEntry(j)*gp(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) 
                  +2*e(j)*G(j,j));
        double pt3 = c.getEntry(j)*gq(j)*(sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) 
                  -2*e(j)*B(j,j));
        double penaltyTerms = pt1
                            +
                            (j == slackIndex ? 0 : 
                            (
                                +pt2
                                +pt3
                            ));
        return
            constraintTerms
            +penaltyTerms
            ;
    }

    protected double gradL_f(int j)
    {
        int dimension = p.getDimension();
        double constraintTerms = sum(i -> lp(i)*(f(i)*G(i,j) - e(i)*B(i,j)), dimension, j) // i != j
                    +sum(i -> lq(i)*(-e(i)*G(i,j) - f(i)*B(i,j)), dimension, j) // i != j
                    +lp(j)*(sum(n -> f(n)*G(j,n) + e(n)*B(j,n), dimension, j) // n != j
                            +2*G(j,j)*f(j))
                    +lq(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
                            -2*B(j,j)*f(j));
        
        double pt1 = c.getEntry(j)*sum(i -> i == slackIndex ? 0 : 
                        gp(i)*( f(i)*G(i,j) - e(i)*B(i,j))
                       +gq(i)*(-e(i)*G(i,j) - f(i)*B(i,j)),
                       dimension, j);
        double pt2 = c.getEntry(j)*gp(j)*(sum(n -> f(n)*G(j,n) + e(n)*B(j,n), dimension, j) 
                        +2*f(j)*G(j,j));
        double pt3 = c.getEntry(j)*gq(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) 
                        -2*f(j)*B(j,j));
        double penaltyTerms = pt1
                            + 
                            (j == slackIndex ? 0 : 
                            (
                                +pt2
                                +pt3
                            ));
        
//out.print(c*gq(j)+",*(,");
//for(int i = 0; i < dimension; ++i)
//{
//  if(i == 1)
//      out.print(e(i)*G(j,i) - f(i)*B(j,i));
//}
//out.println(",) + "+(-2*f(j)*B(j,j)));
        
        return 
            constraintTerms
            +penaltyTerms
            ;
    }
    
    /**
     * Lagrangian.
     * @return
     */
    protected double L()
    {
        int dimension = p.getDimension();
        
        return 
                cost()
                +sum(i -> lp(i)*gp(i), dimension)
                +sum(i -> lq(i)*gq(i), dimension)
                +0.5*sum(i -> c.getEntry(i)*gp(i)*gp(i), dimension, slackIndex)
                +0.5*sum(i -> c.getEntry(i)*gq(i)*gq(i), dimension, slackIndex);
    }
    
    protected double lagrangian(RealVector p, RealVector q, RealVector e, RealVector f, RealVector lp, RealVector lq)
    {
        RealVector gp = gp(p, q, e, f);
        RealVector gq = gq(p, q, e, f);
        gp.setEntry(slackIndex, 0);
        gq.setEntry(slackIndex, 0);
        
        @SuppressWarnings("deprecation")
        double aug = c.ebeMultiply(gp).dotProduct(gp) + c.ebeMultiply(gq).dotProduct(gq);
        double lgp = lp.dotProduct(gp);
        double lgq = lq.dotProduct(gq);
        
        return cost(p, q) + lgp + lgq + 0.5*aug;
    }
    
    public double cost(RealVector p, RealVector q)
    {
        int dimension = p.getDimension();
        return 
               sum(i -> generatorMask.getEntry(i) == 1 ? cost_i(p.getEntry(i), q.getEntry(i), p_costMin.getEntry(i)) : 0, dimension);
    //              sum(i -> (1-e(i))*(1-e(i)) + f(i)*f(i), dimension);
    }
    
    public RealVector gp(RealVector p, RealVector q, RealVector e, RealVector f)
    {
        int dimension = p.getDimension();
        return vector(dimension, i -> i == slackIndex ? 0 : gp(p, q, e, f, i));
    }
    
    protected double gp(RealVector p, RealVector q, RealVector e, RealVector f, int i)
    {
        int dimension = p.getDimension();
        double p_i = p.getEntry(i);
        double e_i = e.getEntry(i);
        double f_i = f.getEntry(i);
        double d = sum(n -> e_i*e.getEntry(n)*G(i,n) + f_i*f.getEntry(n)*G(i,n) + f_i*e.getEntry(n)*B(i,n) - e_i*f.getEntry(n)*B(i,n), dimension);
        return d - p_i;
    }
    
    protected RealVector gq(RealVector p, RealVector q, RealVector e, RealVector f)
    {
        int dimension = q.getDimension();
        return vector(dimension, i -> i == slackIndex ? 0 : gq(p, q, e, f, i));
    }
    
    protected double gq(RealVector p, RealVector q, RealVector e, RealVector f, int i)
    {
        int dimension = p.getDimension();
        double q_i = q.getEntry(i);
        double e_i = e.getEntry(i);
        double f_i = f.getEntry(i);
        double d = sum(j -> f_i*e.getEntry(j)*G(i,j) - e_i*f.getEntry(j)*G(i,j) - e_i*e.getEntry(j)*B(i,j) - f_i*f.getEntry(j)*B(i,j), dimension);
        return d - q_i;
    }

    protected double cost()
    {
    	return cost(p, q);
//        int dimension = p.getDimension();
//        return 
//                sum(i -> generatorMask.getEntry(i) == 1 ? cost_i(p(i), q(i), p_costMin.getEntry(i)) : 0, dimension);
////              sum(i -> (1-e(i))*(1-e(i)) + f(i)*f(i), dimension);
    }

//    protected double cost(int i)
//    {
//        return cost_i(p(i), q(i), p_costMin.getEntry(i));
//    }

    protected double cost_i(double p_i, double q_i, double p_max_i)
    {
        double d = p_i - p_max_i;
        return COST_MULTIPLIER*(0.5*d*d + 0.5*q_i*q_i);
    }
    
    public static void logBusses(AnalysisResults results)
    {
        out.println("index,bus,e,f,p,q");
        Map<String, Integer> numbers = results.getBusNumbers();
        for (String name : numbers.keySet())
        {
            int index = numbers.get(name);
            out.print(index);
            out.print(',');
            out.print(name);
            out.print(',');
            Complex v = results.getBusVoltage(name);
            out.print(format(v.getReal()));
            out.print(',');
            out.print(format(v.getImaginary()));
            out.print(',');
            Complex s = results.getBusPower(name);
            out.print(format(s.getReal()));
            out.print(',');
            out.print(format(s.getImaginary()));
            out.println();
        }
    }
    
    public void logSubGradients()
    {
        double eMin = 0.999;
        double eRange = 0.002;
        double fMin = -0.001;
        double fRange = 0.002;
        
        RealVector _p = new ArrayRealVector(p);
        RealVector _q = new ArrayRealVector(q);
        RealVector _e = new ArrayRealVector(e);
        RealVector _f = new ArrayRealVector(f);
        
        System.out.print("e,f,");
        for (Agent agent : agents)
            if(!agent.isSlack)
                System.out.print("L(e_"+agent.index+"),"+"L(f_"+agent.index+"),");
        System.out.println();
        for(double x = 0; x <= 1; x += 0.01)
        {
            System.out.print((eMin+eRange*x)+","+(fMin+fRange*x)+",");
            for (Agent agent : agents)
            {
                if(agent.isSlack)
                    continue;
                _e.setEntry(agent.index, eMin+x*eRange);
                double L_e = lagrangian(_p, _q, _e, _f, lp, lq);
                _e.setEntry(agent.index, e.getEntry(agent.index));

                _f.setEntry(agent.index, fMin+x*fRange);
                double L_f = lagrangian(_p, _q, _e, _f, lp, lq);
                _f.setEntry(agent.index, f.getEntry(agent.index));
                
                System.out.print(L_e+","+L_f+",");
            }
            System.out.println();
        }
    }
    
    public void loge0e1()
    {
        RealVector _p = new ArrayRealVector(p);
        RealVector _q = new ArrayRealVector(q);
        RealVector _e = new ArrayRealVector(e);
        RealVector _f = new ArrayRealVector(f);
        
        double e0Min = 0.97;
        double e0Max = 1.01;
        double e0Inc = 0.001;
        double e1Min = 0.97;
        double e1Max = 1.01;
        double e1Inc = 0.001;
        for(double e1 = e1Min; e1 < e1Max; e1 += e1Inc)
            System.out.print(e1+",");
        System.out.println();
        for(double e0 = e0Min; e0 < e0Max; e0 += e0Inc)
        {
            System.out.print(e0+",");
            _e.setEntry(0, e0);
            for(double e1 = e1Min; e1 < e1Max; e1 += e1Inc)
            {
                _e.setEntry(1, e1);
                double L_e = lagrangian(_p, _q, _e, _f, lp, lq);
                System.out.print(L_e+",");
            }
            System.out.println();
        }
    }
}