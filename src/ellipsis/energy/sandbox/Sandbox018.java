package ellipsis.energy.sandbox;

import static ellipsis.energy.sandbox.GridConstants.SLACK_SOURCE;
import static ellipsis.util.MatrixHelper.invert;
import static ellipsis.util.Sum.sum;
import static ellipsis.util.VectorHelper.vector;
import static java.lang.Math.abs;
import static java.lang.Math.max;

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
import ellipsis.util.Sum.IndexedFunction;
import ellipsis.util.TeeOutputStream;

/**
 * Extension to {@link Sandbox016} to use asynchronous consensus.
 * Uses Sandbox016 data.
 * @author bmillar
 * FIXME The gradients etc. assume h is the sum of powers whereas it's actually the average power.
 */
public class Sandbox018
{
	protected static final int DISPLAY_DIMENSION = 35;
	
    public static void main(String[] args)
	{
		new Sandbox018().run();
	}

    public static interface NeighbourFunction
    {
        double value(Agent neighbour);
    }
    
    public static double sumA(Map<Agent, Complex> keys, NeighbourFunction f)
    {
        return sumA(keys.keySet(), f);
    }
    
    public static double sumA(Set<Agent> keys, NeighbourFunction f)
    {
        double sum = 0;
        for (Agent agent : keys)
        {
            double v = f.value(agent);
            sum += v;
        }
        
        return sum;
    }
    
    
    //// Global variables ////
    
    protected final double BASE_VOLTAGE = 120;
    protected final double SLACK_VOLTAGE = 1.00;
    protected final double BASE_POWER = 10e3;
    
    protected final double V_MAX = 1.05;
    protected final double V_MIN = 0.95;
    protected final double V_ARG_MIN = -Math.PI/4;
    protected final double V_ARG_MAX = Math.PI/4;

    protected double MIN_CONSTRAINT_IMPROVEMENT = 0.99; // Minimum improvement percentage of g(x) constraint
    protected double INITIAL_G_AUG_SCALE = 1;// 1e-3; 
    protected double G_AUG_SCALE_STEP = 1; //1.0025;
    protected double G_MAX_AUG_SCALE = 1e6;

    protected double INITIAL_H_AUG_SCALE = 10; // 1e-0; // FIXME 1e-3
    protected double H_AUG_SCALE_STEP = 1;//1.0025; // FIXME 1.0025
    protected double H_MAX_AUG_SCALE = 10; // FIXME 1e6;
    
    private final double COST_MULTIPLIER = 1;
    
    protected int X_STEPS_PER_ITERATION = 1;
    protected int MAX_X_STEPS = 1;//100;
    
    protected int K = 6000;
    protected int DEBUG_RATE = K/1000;
    protected double AGENT_SELECTION_PROBABILITY = 1.0; // FIXME 0.5;
    
    protected double INITIAL_STEP_SIZE = 1;
    protected double MIN_STEP_SIZE = 1e-50;

    protected double XI = 0.05; // Consensus step size. FIXME changing this may help with non-zero grad convergence
    protected double ETA_G = 1; // 1.5; // 2.5; // Multiplier update step size.
    protected double ETA_H = 0.005;
    
    protected double MOMENTUM_RATE = 0.0;
    
    protected double EPSILON_BASE_PQ = 0.1;
    protected double EPSILON_BASE_EF = 100;
    protected double EPSILON_TARGET = 10e-2;
    
    protected double CONSTRAINT_BASE = 1;
    protected double CONSTRAINT_TARGET = 1e-4;
    
    // Switches 
    // (not final to avoid warnings when set to false and to 
    //  allow override in subclass or during debug):
    protected boolean PROJECT_X = true; // FIXME
    
    protected boolean FORMAT_NUMBERS = false;
    protected boolean VALIDATE_BACKTRACKING = false;
    protected boolean VALIDATE_LAGRANGE_DECREASE = true; // Use with USE_TRUE_G = true. FIXME
    protected boolean VALIDATE_CONSENSUS = false; // FIXME
    
    protected boolean START_WITH_TRUE_VOTLAGES = true;
    protected boolean START_WITH_OPTIMAL_POWERS = true; // WARNING This has no effect for small initial epsilon since the first iteration will set powers to their unconstrained optimum anyway.
    protected boolean USE_TRUE_G = true;
    protected boolean USE_TRUE_H = false; // FIXME
    
    protected boolean DISABLE_VOTLAGES_UPDATES = false;
    protected boolean DISABLE_LAMBDA_UPDATE = false;
    protected boolean DISABLE_POWER_FLOW_UPDATE = false;
    protected boolean DISABLE_COST = false;
    
    protected boolean DISABLE_CONSENSUS_UPDATE = false;
    protected boolean DISABLE_MU_UPDATE = false; // FIXME even though ETA=0 this is still having an effect
    protected boolean DISABLE_LINE_LOSS = false;
    
    protected boolean USE_SIGMOID_ALPHA_UPDATE = false;
    
    protected boolean ROUND_TO_ZERO = false; // FIXME
    protected double MIN_VALUE = 1e-12;
    protected double MIN_G_VALUE = 1e-6;
    protected double MIN_H_VALUE = 1e-6;
    
    protected boolean DONT_UPDATE_ALPHA_WHILE_IMPROVING = false;
    protected int ALPHA_UPDATE_RATE = 1;
    protected boolean KEEP_G_AND_H_CLOSE = true; // FIXME
    protected boolean APPROXIMATE_CONSTRAINTS = true; // FIXME
    
    protected boolean DEBUG = true;
    
    public PrintStream out;
    public Sandbox018()
    {
        try
        {
            out = new PrintStream(new TeeOutputStream(new FileOutputStream("/tmp/Sandbox018-"+K+".csv"), System.out));
        }
        catch (FileNotFoundException e)
        {
            throw new RuntimeException(e);
        }
    }
    
    protected double round(double d)
    {
        return round(d, MIN_VALUE);
    }
    
    protected double round(double d, double min) 
    { 
        if(ROUND_TO_ZERO)
        {
            if(d < 0)
                return -round(-d, min);
            else
                return d < min ? 0 : d;
        }
        else
        {
            return d;
        }
    }

    public final int P = 1;
    public final int Q = 2;
    public final int E = 4;
    public final int F = 8;
    
    protected RealVector p, q, e, f, lp, lq;
    protected RealMatrix G, B;
    protected RealVector p_max;
    protected RealVector p_costMin;
    protected int slackIndex;
    protected RealVector generatorMask;
    protected Set<Agent> agents;
    protected Grid grid;
    protected RealVector alpha_g, alpha_h;
    protected RealVector t;
    protected double cpuTime = 0;
    protected RealVector hp, hq, wp, wq, mup, muq;
    protected double epsilon_pq, epsilon_ef;
    
//    protected MeanVarEstimator variance_pNorm =  new MeanVarEstimator();
//    protected MeanVarEstimator variance_qNorm =  new MeanVarEstimator();
//    protected MeanVarEstimator variance_eNorm =  new MeanVarEstimator();
//    protected MeanVarEstimator variance_fNorm =  new MeanVarEstimator();
//    protected MeanVarEstimator variance_gpNorm = new MeanVarEstimator();
//    protected MeanVarEstimator variance_gqNorm = new MeanVarEstimator();
//    protected MeanVarEstimator variance_hpNorm = new MeanVarEstimator();
//    protected MeanVarEstimator variance_hqNorm = new MeanVarEstimator();
    
    // Debug:
    double diff = 0;
    
    protected void run()
    {
        start();
        AnalysisResults results = init();
        loop();
        testResults(results);
        finish();
    }

    protected Timer timer;
    protected void start()
    {
        timer = Timer.startNewTimer();
    }

    protected void finish()
    {
        double t = timer.stop()/1000.0;
        int m = (int)(t/60.0);
        double s = t - m*60;
        if(DEBUG)
        {
            out.printf("\nTotal time %d:%.3f\n", m, s);
            out.println(SimpleDateFormat.getDateTimeInstance().format((Calendar.getInstance().getTime())));
        }
        
        out.flush();
    }
    
    protected void testResults(AnalysisResults results)
    {
        if(!DEBUG)
            return;
        
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
        out.println("index,bus,e,f,p,q,loss_p,loss_q,,|v|,|s|");
        for (String busName : busNumbers.keySet())
        {
            int i = busNumbers.get(busName);
            double f_i = f(i);
            double e_i = e(i);
            String v_abs = format(Math.sqrt(e_i*e_i+f_i*f_i));
            double p_i = p(i);
            double q_i = q(i);
            String s_abs = format(Math.sqrt(p_i*p_i+q_i*q_i));
            double loss_p_i = 0.5*lineLoss_p(i);
            double loss_q_i = 0.5*lineLoss_q(i);
            out.println(i+","+busName+','+format(e_i)+','+format(f_i)+','+format(p_i)+','+format(q_i)+','+format(loss_p_i)+','+format(loss_q_i)+",,"+v_abs+','+s_abs);
        }
        
        // Analyse grid to get parameters:
        LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        lfa.setIterations(100);
        lfa.setTargetError(1e-6);
        AnalysisResults results2 = lfa.analyse();//complexVoltages()); // WARNING Passing in voltages is changing the load powers somehow in the results
        
        out.println("\nAnalysis Results:");
        
        logBusses(results2);
        
//        GridDisplay.showInFrame(grid, grid.getBasePower(), grid.getBaseVoltage()).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
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
        protected double alpha_g, alpha_h;
        public boolean isSlack;
        public boolean isGenerator;
        public int index;
        protected Map<Agent, Complex> neighbours = new HashMap<>();
        public double G, B; // Self admittances.
        public double lastT = 1;
        public double epsilon_pq, epsilon_ef;
        
        // Momentum parameters:
        private double previousDp, previousDq, previousDe, previousDf;
        
        // Consensus parameters:
        private double hp, hq; // Average network power estimate.
        private double mup, muq; // Average network power constraint multiplier.
        private double wp, wq; // Balancing value from neighbours.
        
        public Agent(int i)
        {
            this.index = i;
        }
        
        public void addNeighbour(Agent neighbour, Complex admittance)
        {
            neighbours.put(neighbour, admittance);
        }

        protected boolean debug_backtrack = false;
        
        @SuppressWarnings("unused")
        public void step(int k)
        {
boolean debug = false;
if(k >= 2000 && index == 0)
    Util.nullop();

            if(!isSlack)
            {
                nextEpsilon(k);
                
                // Update constraint function values:
                gp = gp();
                gq = gq();
            }
                
            // Remember current constraint values:
            double gp_old = gp;
            double gq_old = gq;
            double hp_old = hp;
            double hq_old = hq;
            double p_old = p;
            double q_old = q;
            double e_old = e;
            double f_old = f;
            
            if(!isSlack)
            {
                //for(int i = 0; i < X_STEPS_PER_ITERATION; ++i)
                boolean steppedPQ, steppedEF;
if(debug)
    System.out.println("initial,"+lagrangian()+","+normGradL_p()+","+lastT);
int i = 1;
                do
                {
                    double agentLagBefore = lagrangian();
if(i >= 99)
    Util.nullop();
                    if(isGenerator)
                        steppedPQ = stepX(P|Q);
                    else
                        steppedPQ = false;
                    double t = lastT;
if(debug)
    System.out.println("stepX(P|Q),"+lagrangian()+","+normGradL_p()+","+lastT);
                    
                    if(!DISABLE_VOTLAGES_UPDATES)
                    {
                        steppedEF = stepX(E|F);
                        t = max(t, lastT);
                    }
                    else
                    {
                        steppedEF = false;
                    }
if(debug)
    System.out.println("stepX(E|F),"+lagrangian()+","+normGradL_p()+","+lastT);
                    
                    lastT = max(t, lastT);
++i;
                    // Before and after check of L(x,l,mu):
                    if(VALIDATE_LAGRANGE_DECREASE)
                    {
                        double agentLagAfter = lagrangian();
                        if(agentLagAfter > agentLagBefore && Math.abs((agentLagAfter - agentLagBefore)/max(agentLagAfter, agentLagBefore)) > 1e-6) 
                            throw new RuntimeException("L() increased after x step(s): "+agentLagAfter +" > "+ agentLagBefore);
                    }
                } while(   i < MAX_X_STEPS 
                           && 
                           (
                               (normGradL_p() > epsilon_pq || normGradL_q() > epsilon_pq)&steppedPQ
                               ||
                               (normGradL_e() > epsilon_ef || normGradL_f() > epsilon_ef)&steppedEF
                           )
                       );
                
                if(!DISABLE_LAMBDA_UPDATE)
                    stepLambda();
            }
            
            /* reset(p_old, q_old, e_old, f_old, gp_old, gq_old, hp_old, hq_old) */

            stepAugScale(k, gp_old, gq_old, hp_old, hq_old);
            
            // Consensus step:
            stepConsensus(k);//, p_old, q_old, e_old, f_old);
        }
        
        protected void reset(
                double p_old, double q_old, double e_old, double f_old, 
                double gp_old, double gq_old,
                double hp_old, double hq_old)
        {
            p = p_old;
            q = q_old;
            e = e_old;
            f = f_old;
            gp = gp_old;
            gq = gq_old;
            hp = hp_old;
            hq = hq_old;
        }

        /**
         * If x is on the surface of X then return 0.
         * Otherwise return the norm of the gradient.
         * @return
         */
        private double normGradL_p()
        {
            if(!isGenerator)
                return 0;
            
            if(0 < p && abs(p-p_max) > 1e-12) // 0 <= p <= p_max 
                return abs(gradL_p());
            else
                return 0;
        }

        /**
         * If x is on the surface of X then return 0.
         * Otherwise return the norm of the gradient.
         * @return
         */
        private double normGradL_q()
        {
            if(!isGenerator)
                return 0;
            
            return abs(gradL_q());
        }

        /**
         * If x is on the surface of X then return 0.
         * Otherwise return the norm of the gradient.
         * @return
         */
        private double normGradL_e()
        {
            if(!onXSurface_V())
            {
                return abs(gradL_e());
            }
            else
                return 0;
        }

        /**
         * If x is on the surface of X then return 0.
         * Otherwise return the norm of the gradient.
         * @return
         */
        private double normGradL_f()
        {
            if(!onXSurface_V())
                return abs(gradL_f());
            else
                return 0;
        }

        protected boolean onXSurface_V()
        {
            Complex v = new Complex(e, f);
            double abs = v.abs();
            double arg = v.getArgument();
            return abs(abs-V_MAX) < 1e-12     || abs(abs-V_MIN) < 1e-12 ||
                   abs(arg-V_ARG_MAX) < 1e-12 || abs(arg-V_ARG_MIN) < 1e-12;
        }

        protected void nextEpsilon(int k)
        {
            double epsilonGeometric = Math.pow(EPSILON_TARGET, k/(double)K);
            epsilon_pq = EPSILON_BASE_PQ*epsilonGeometric;
            epsilon_ef = EPSILON_BASE_EF*epsilonGeometric;
        }

        /**
         * Uses agent's p, q, e and f values.
         * @return
         */
        protected double lagrangian()
        {
            return lagrangian(p, q, e, f);
        }
        
        /**
         * Assumes hp and hq values correspond to current p and q values.
         * @param p Test value.
         * @param q
         * @param e
         * @param f
         * @return
         */
        protected double lagrangian(double p, double q, double e, double f)//, double p_old, double q_old)
        {
            RealVector _p = new ArrayRealVector(Sandbox018.this.p);
            _p.setEntry(index, p);
            RealVector _q = new ArrayRealVector(Sandbox018.this.q);
            _q.setEntry(index, q);
            
            RealVector _e = new ArrayRealVector(Sandbox018.this.e);
            _e.setEntry(index, e);
            RealVector _f = new ArrayRealVector(Sandbox018.this.f);
            _f.setEntry(index, f);

            RealVector _hp = new ArrayRealVector(Sandbox018.this.hp);
            double dhp = p-this.p +  sumA(neighbours, agent -> lineLoss_p(agent, e, f) - lineLoss_p(agent, this.e, this.f));
            _hp.setEntry(index, hp+dhp);
            RealVector _hq = new ArrayRealVector(Sandbox018.this.hq);
            double dhq = q-this.q + sumA(neighbours, agent -> lineLoss_q(agent, e, f) - lineLoss_q(agent, this.e, this.f));
            _hq.setEntry(index, hq+dhq);
            
            double L = Sandbox018.this.lagrangian(
                    _p, _q, _e, _f, 
                    _hp, 
                    _hq, 
                    Sandbox018.this.mup, 
                    Sandbox018.this.muq, 
                    Sandbox018.this.lp, 
                    Sandbox018.this.lq);
            return L;
        }

        private void stepConsensus(int k)//, double p_old, double q_old, double e_old, double f_old)
        {
            // Bias: Adjust for change in power: Moved to stepX
//            if(!DISABLE_CONSENSUS_UPDATE)
//            {
//                hp += p - p_old;
//                hq += q - q_old;
//                if(!DISABLE_LINE_LOSS)
//                {
//                    hp += sumA(neighbours, agent -> lineLoss_p(agent, e, f) - lineLoss_p(agent, e_old, f_old));
//                    hq += sumA(neighbours, agent -> lineLoss_q(agent, e, f) - lineLoss_q(agent, e_old, f_old));
//                }
//            }
            
            // Consensus: Power mismatch estimation:
            double pStep = 0;
            double qStep = 0;
            for (Agent agent : neighbours.keySet())
            {
                if(agent.isSlack)
                    continue;
                
                double pStep_j = XI*(agent.hp - hp);
                double qStep_j = XI*(agent.hq - hq);
                
                pStep += pStep_j;
                qStep += qStep_j;
                
                agent.addPOffset(pStep_j, qStep_j);
            }
            
            if(!DISABLE_CONSENSUS_UPDATE)
            {
                hp += pStep - wp;
                hq += qStep - wq;
                wp = 0;
                wq = 0;
            }
            
            // FIXME Experimental:
            // Round h to zero through a bias:
//            if(ROUND_TO_ZERO)
//            {
//                if(abs(hp) < MIN_H_VALUE)
//                {
//                    double offset_p;
//                    offset_p = -hp/neighbours.size();
//                    hp = 0;
//                    
//                    for (Agent agent : neighbours.keySet())
//                    {
//                        if(agent.isSlack)
//                        {
//                            this.addPOffset(offset_p, 0); // Must keep global average.
//                        }
//                        else
//                        {
//                            agent.addPOffset(offset_p, 0);
//                        }
//                    }
//                }
//                
//                if(abs(hq) < MIN_H_VALUE)
//                {
//                    double offset_q;
//                    offset_q = -hq/neighbours.size();
//                    hq = 0;
//                    
//                    for (Agent agent : neighbours.keySet())
//                    {
//                        if(agent.isSlack)
//                        {
//                            this.addPOffset(0, offset_q); // Must keep global average.
//                        }
//                        else
//                        {
//                            agent.addPOffset(0, offset_q);
//                        }
//                    }
//                }
//            }
            
            // Multiplier estimation:
            if(!DISABLE_MU_UPDATE)
            {
                double muStepSize = muStepSize();
                mup += XI*sumA(neighbours, agent -> agent.mup - mup) + muStepSize*hp;
                muq += XI*sumA(neighbours, agent -> agent.muq - muq) + muStepSize*hq;
            }
        }
        
        private double muStepSize()
        {
            if(alpha_h == 0)
                return ETA_H;
            else
                return ETA_H*alpha_h;
        }

        private double lineLoss_p(Agent agent, double e, double f)
        {
            double de = e - agent.e;
            double df = f - agent.f;
            return round(G(agent)*(de*de + df*df));
        }

        private double lineLoss_q(Agent agent, double e, double f)
        {
            double de = e - agent.e;
            double df = f - agent.f;
            return round(-B(agent)*(de*de + df*df));
        }

        private void addPOffset(double deltaP, double deltaQ)
        {
            wp += deltaP;
            wq += deltaQ;
        }
        
        public Complex getAverageNetworkPowerEstimate()
        {
            return new Complex(hp(), hq());
        }
        
        public double hp()
        {
            if(USE_TRUE_H)
            {
                Sandbox018.this.p.setEntry(index, p);
                Sandbox018.this.q.setEntry(index, q);
                Sandbox018.this.e.setEntry(index, e);
                Sandbox018.this.f.setEntry(index, f);
                return averagePowerMismatch_p();
            }
            else
            {
                return hp;
            }
        }
        
        public double hq()
        {
            if(USE_TRUE_H)
            {
                Sandbox018.this.p.setEntry(index, p);
                Sandbox018.this.q.setEntry(index, q);
                Sandbox018.this.e.setEntry(index, e);
                Sandbox018.this.f.setEntry(index, f);
                return averagePowerMismatch_q();
            }
            else
            {
                return hq;
            }
        }

//        Random rand = new Random(0);
        int stepCount = 0;
        /**
         * 
         * @param pqef
         * @return true iff a step could be taken
         */
        private boolean stepX(int pqef)
        {
            double p_old = p;
            double q_old = q;
            double e_old = e;
            double f_old = f;
            
            // Find x step direction:
            double dL_p = isGenerator && (pqef&P)!=0 ? gradL_p() : 0;
            double dL_q = isGenerator && (pqef&Q)!=0 ? gradL_q() : 0;
            double dL_e = !isSlack    && !DISABLE_VOTLAGES_UPDATES && (pqef&E)!=0 ? gradL_e() : 0;
            double dL_f = !isSlack    && !DISABLE_VOTLAGES_UPDATES && (pqef&F)!=0 ? gradL_f() : 0;
            
//            if(MOMENTUM_RATE > 0)
//            {
//                dL_p = MOMENTUM_RATE*previousDp + (1-MOMENTUM_RATE)*dL_p;
//                dL_q = MOMENTUM_RATE*previousDq + (1-MOMENTUM_RATE)*dL_q;
//                dL_e = MOMENTUM_RATE*previousDe + (1-MOMENTUM_RATE)*dL_e;
//                dL_f = MOMENTUM_RATE*previousDf + (1-MOMENTUM_RATE)*dL_f;
//                
//                previousDp = dL_p;
//                previousDq = dL_q;
//                previousDe = dL_e;
//                previousDf = dL_f;
//            }

            RealVector x = new ArrayRealVector(new double[]{p, q, e, f});
            RealVector dL_x = new ArrayRealVector(new double[]{dL_p, dL_q, dL_e, dL_f});
            RealVector dx = dL_x.mapMultiply(-1); // delta x = -grad_x L(x,lambda)
            
            // Debug:
            if(debug_backtrack)
            {
                update(); // Global values are used in some testing calculations, so make sure they are up to date here.
                
                double cOld_g = Sandbox018.this.alpha_g.getEntry(index);
                double cOld_h = Sandbox018.this.alpha_h.getEntry(index);
                Sandbox018.this.alpha_g.setEntry(index, alpha_g);
                Sandbox018.this.alpha_h.setEntry(index, alpha_h);
                
                compareGlobalLocal(dL_p, dL_q, dL_e, dL_f);
                System.out.println();
                checkGradient(dL_p, dL_q, dL_e, dL_f);
                System.out.println();
                plotImprovements();
                System.out.println();
                plotComparison();
                
                Sandbox018.this.alpha_g.setEntry(index, cOld_g);
                Sandbox018.this.alpha_h.setEntry(index, cOld_h);
            }
            
            // Backtrack to find best step size:
            double t = backTrack(dL_x, dx, p, q, e, f, INITIAL_STEP_SIZE);
            if(t == 0.0)
            {
                lastT = 0.0;
                return false;
            }
            
            // Update state:
            RealVector step = dx.mapMultiply(t);
            
            if(MOMENTUM_RATE > 0)
            {
                step.setEntry(0, MOMENTUM_RATE*previousDp + (1-MOMENTUM_RATE)*step.getEntry(0));
                step.setEntry(1, MOMENTUM_RATE*previousDq + (1-MOMENTUM_RATE)*step.getEntry(1));
                step.setEntry(2, MOMENTUM_RATE*previousDe + (1-MOMENTUM_RATE)*step.getEntry(2));
                step.setEntry(3, MOMENTUM_RATE*previousDf + (1-MOMENTUM_RATE)*step.getEntry(3));
                
                previousDp = step.getEntry(0);
                previousDq = step.getEntry(1);
                previousDe = step.getEntry(2);
                previousDf = step.getEntry(3);
            }
            
            RealVector x_new = x.add(step);
            p_new = x_new.getEntry(0);
            q_new = x_new.getEntry(1);
            e_new = x_new.getEntry(2);
            f_new = x_new.getEntry(3);
            
            // Check for NaN:
            if(x_new.isNaN())
            {
                System.err.println("x_new is NaN");
                Util.breakpoint();
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

            // Debug step:
            boolean debug = false;
            if(debug)
            {
                System.out.println(lagrangianChange(
                        step.getEntry(0),
                        step.getEntry(1),
                        step.getEntry(2),
                        step.getEntry(3)
                        ));
            }
            
            // Update x:
            p = isGenerator ? p_new : p;
            q = isGenerator ? q_new : q;
            if(!DISABLE_VOTLAGES_UPDATES)
            {
                e = e_new;
                f = f_new;
            }
            
            // Consensus Bias: Adjust for change in power:
            if(!DISABLE_CONSENSUS_UPDATE)
            {
                double step_p = p - p_old;
                double step_q = q - q_old;
                
                hp += step_p;
                hq += step_q;
                
                if(!DISABLE_LINE_LOSS)
                {
                    hp += sumA(neighbours, agent -> lineLoss_p(agent, e, f) - lineLoss_p(agent, e_old, f_old));
                    hq += sumA(neighbours, agent -> lineLoss_q(agent, e, f) - lineLoss_q(agent, e_old, f_old));
                }
            }
            
            // Update constraint function values:
            gp = gp();
            gq = gq();
            
            lastT = t;
            
            ++stepCount;
            
            return true;
        }

        protected String doubleToLongBits(RealVector v)
        {
            StringBuffer sb = new StringBuffer();
            for(int i = 0; i < v.getDimension(); ++i)
            {
                sb.append(Double.doubleToLongBits(v.getEntry(i)));
                sb.append("\n");
            }
            return sb.toString();
        }

        private double backTrack(RealVector dL_x, RealVector dx, double p, double q, double e, double f, double t)
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
            
            boolean debug1 = false;
            if(debug1)
            {
                FORMAT_NUMBERS = true;
                
                // Plot improvement for e, f:
                double st = 2e-3;
                double dp_max = 10*st;
                double dq_max = 10*st;
                for(double _dp = -10*st; _dp < dp_max; _dp += st)
                {
                    System.out.print(_dp+",");
                    for(double _dq = -10*st; _dq < dq_max; _dq += st)
                    {
                        double val = lagrangianChange(_dp, _dq, 0, 0);
                        System.out.print(format(val)+",");
                    }
                    System.out.println(",");
                }
            }

            double dp2, dq2, de2, df2;
            do
            {
                dp2 = dp*t;
                dq2 = dq*t;
                de2 = de*t;
                df2 = df*t;
                
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
                
                L_improvement = lagrangianChange(dp2, dq2, de2, df2);
                minImprovement = t*step;
                
                if(debug_backtrack)
                {
                    double lowerBound = t*step/alpha;
                    double Lx = lagrangian(p, q, e, f);
                    double Lx2 = lagrangian(p+dp2, q+dq2, e+de2, f+df2);
                    double dLx = Lx2 - Lx;
                    System.out.println(t+","+L_improvement+","+minImprovement+","+lowerBound+","+dLx+","+Lx+","+Lx2);
                }

                t *= beta;
            } while((Double.isNaN(L_improvement) || L_improvement > minImprovement) && t > MIN_STEP_SIZE); // || debug_backtrack);
            
            if(t <= MIN_STEP_SIZE)
                return 0.0;
            
            if(L_improvement > 0)
                throw new RuntimeException("Improvement was positive.");
            
            // Undo last step update:
            t /= beta;

            // Debug:
            if(debug_backtrack)
            {
                debugBacktrackingChange(p, q, e, f, L_improvement, dp2, dq2, de2, df2);
                
                // Output two more points for analysing change curve:
                System.out.println("0,0,0,0,0");
                double _t = -t;
                L_improvement = lagrangianChange(
                    dp*_t, 
                    dq*_t, 
                    de*_t,
                    df*_t);
                double lowerBound = _t*step/alpha;
                double Lx = lagrangian(p, q, e, f);
                double Lx2 = lagrangian(p+dp*_t, q+dq*_t, e+de*_t, f+df*_t);
                double dLx = Lx2 - Lx;
                System.out.println(_t+","+L_improvement+","+(_t*step)+","+lowerBound+","+dLx);
            }
            
            // Limit step size:
            if(t < MIN_STEP_SIZE)
                t = MIN_STEP_SIZE;

            return t;
        }

        protected void debugBacktrackingChange(double p, double q, double e,
                double f, double L_improvement, double dp2, double dq2,
                double de2, double df2)
        {
            RealVector _p = new ArrayRealVector(Sandbox018.this.p);
            RealVector _q = new ArrayRealVector(Sandbox018.this.q);
            RealVector _e = new ArrayRealVector(Sandbox018.this.e);
            RealVector _f = new ArrayRealVector(Sandbox018.this.f);
            RealVector _hp = new ArrayRealVector(Sandbox018.this.hp);
            RealVector _hq = new ArrayRealVector(Sandbox018.this.hq);
            
            double Lx = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                    Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);

            final double _de2 = de2;
            final double _df2 = df2;
            double dhp = sumA(neighbours, agent -> lineLoss_p(agent, e+_de2, f+_df2) - lineLoss_p(agent, e, f)) + dp2;
            double dhq = sumA(neighbours, agent -> lineLoss_q(agent, e+_de2, f+_df2) - lineLoss_q(agent, e, f)) + dq2;
            
            _p.setEntry(index, p+dp2);
            _q.setEntry(index, q+dq2);
            _e.setEntry(index, e+de2);
            _f.setEntry(index, f+df2);
            _hp.setEntry(index, hp+dhp);
            _hq.setEntry(index, hq+dhq);
            
            double Lx2 = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                    Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
            
            double change = Lx2 - Lx;
            if(abs(change - L_improvement) > 1e-9)
                System.err.println("Lagrange change disagreement: Local change = "+ L_improvement +", global change = "+change);
        }
        
        public void showPVsE()
        {
            RealVector _p = new ArrayRealVector(Sandbox018.this.p);
            RealVector _q = Sandbox018.this.q;
            RealVector _e = new ArrayRealVector(Sandbox018.this.e);
            RealVector _f = Sandbox018.this.f;
            RealVector _lp = Sandbox018.this.lp;
            RealVector _lq = Sandbox018.this.lq;
            RealVector _mup = Sandbox018.this.mup;
            RealVector _muq = Sandbox018.this.muq;
            RealVector _hp = Sandbox018.this.hp;
            RealVector _hq = Sandbox018.this.hq;
            
            for(double _de = -10; _de <= 10; _de += 0.2)
            {
                _e.setEntry(index, Sandbox018.this.e.getEntry(index)+_de);
                for(double _dp = -10; _dp <= 10; _dp += 0.1)
                {
                    _p.setEntry(index, Sandbox018.this.p.getEntry(index)+_dp);
                    System.out.print(Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq, _mup, _muq, _lp, _lq));
                    System.out.print(',');
                }
                System.out.println(',');
            }
        }

        protected void stepLambda()
        {
            // Maximise (l := l - c*g(x), ref http://en.wikipedia.org/wiki/Augmented_Lagrangian_method
            //           and Bertsekas' book):
            double dlp = isSlack ? 0 : gp;
            lp = lp + ETA_G*dlp*alpha_g;
            double dlq = isSlack ? 0 : gq;
            lq = lq + ETA_G*dlq*alpha_g;
        }

        private double a_g = (2.0/K)*Math.log(-1+G_MAX_AUG_SCALE/INITIAL_G_AUG_SCALE);
        private double a_h = (2.0/K)*Math.log(-1+H_MAX_AUG_SCALE/INITIAL_H_AUG_SCALE);
        private void stepAugScale(int k, double gp_old, double gq_old, double hp_old, double hq_old)
        {
            if(k%ALPHA_UPDATE_RATE != 0)
                return;
            
            if(APPROXIMATE_CONSTRAINTS)
            {
                double geometric = Math.pow(CONSTRAINT_TARGET, k/(double)K);
                double target = CONSTRAINT_BASE*geometric;
                
                if(abs(gp) < target && abs(gq) < target
                   &&
                   abs(hp) < target && abs(hq) < target)
                {
                    return;
                }
            }
            
            if(DONT_UPDATE_ALPHA_WHILE_IMPROVING)
            {
                if(Math.abs(gp) < MIN_CONSTRAINT_IMPROVEMENT*Math.abs(gp_old) && 
                   Math.abs(gq) < MIN_CONSTRAINT_IMPROVEMENT*Math.abs(gq_old) &&
                   Math.abs(hp) < MIN_CONSTRAINT_IMPROVEMENT*Math.abs(hp_old) && 
                   Math.abs(hq) < MIN_CONSTRAINT_IMPROVEMENT*Math.abs(hq_old))
                    return;
            }
            
            if(KEEP_G_AND_H_CLOSE)
            {
                // Only update alpha_g OR alpha_h:
                if(max(abs(gp), abs(gq)) > max(abs(hp), abs(hq)))
                {
                    if(alpha_g < G_MAX_AUG_SCALE)
                        alpha_g *= G_AUG_SCALE_STEP;
                }
                else
                {
                    if(alpha_h < H_MAX_AUG_SCALE)
                        alpha_h *= H_AUG_SCALE_STEP;
                }
            }
            else if(USE_SIGMOID_ALPHA_UPDATE)
            {
                alpha_g = G_MAX_AUG_SCALE/(1+Math.exp(-a_g*(k-K/2.0)));
                alpha_h = H_MAX_AUG_SCALE/(1+Math.exp(-a_h*(k-K/2.0)));
            }
            else
            {
                if(alpha_g < G_MAX_AUG_SCALE)
                    alpha_g *= G_AUG_SCALE_STEP;
                if(alpha_h < H_MAX_AUG_SCALE)
                    alpha_h *= H_AUG_SCALE_STEP;
            }
        }

        private Complex projectVoltage(double e, double f)
        {
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
            
            // Clamp angle:)
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

        private Complex projectPowers(double p, double q)
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
         * This is equivalent to L(x+dx) - L(x) but is calculable within the agent's neighbourhood.
         * @param t
         * @param dp
         * @param dq
         * @param de
         * @param df
         * @return
         */
        private double lagrangianChange(double dp, double dq, double de, double df)
        {
            double dc = isGenerator ? cost(p+dp, q+dq) - cost(p, q) : 0;
            
            double constraintTerms = lagrangianChange_constraintTerms(dp, dq, de, df);
            double augTermP = lagrangianChange_augTermP(dp, dq, de, df);
            double augTermQ = lagrangianChange_augTermQ(dp, dq, de, df);
            double augTerm = DISABLE_POWER_FLOW_UPDATE ? 0 : augTermP + augTermQ;
            
            double consensusConstraintTerms = lagrangeChange_consensusConstraintTerms(dp, dq, de, df);
            double consensusAugTerms = lagrangeChange_consensusAugTerms(dp, dq, de, df);

            double L_change = 
            		dc
					+constraintTerms
					+augTerm
					+consensusConstraintTerms
					+consensusAugTerms;

            // Compare improvement against L(x-t*dx)-L(x)
            if(VALIDATE_BACKTRACKING)
            	compareChange(
            		dp, dq, de, df, 
            		L_change,
            		dc, constraintTerms, consensusConstraintTerms, augTermP, augTermQ, consensusAugTerms);
            
			return round(L_change);
        }

        public double lagrangeChange_consensusConstraintTerms(double dp, double dq, double de, double df)
        {
            return DISABLE_MU_UPDATE ? 0 : 
                mup*dhp(dp, de, df) + muq*dhq(dq, de, df);
        }

        protected double dhp(double dp, double de, double df)
        {
            if(DISABLE_LINE_LOSS)
                return dp;
            
            return dp + sumA(neighbours, agent -> 
                G(agent)*(de*de + df*df +2*de*(e - agent.e) + 2*df*(f - agent.f))
            );
        }

        protected double dhq(double dq, double de, double df)
        {
            if(DISABLE_LINE_LOSS)
                return dq;
            
            return dq + sumA(neighbours, agent -> 
                -B(agent)*(de*de + df*df +2*de*(e - agent.e) + 2*df*(f - agent.f))
            );
        }

        public double lagrangeChange_consensusAugTerms(double dp, double dq, double de, double df)
        {
            double dhp = dhp(dp, de, df);
            double dhq = dhq(dq, de, df);
            return DISABLE_CONSENSUS_UPDATE ? 0 : 
                0.5*alpha_h*(dhp*(dhp + 2*hp()) + dhq*(dhq + 2*hq()));
        }

        public double lagrangianChange_augTermP(double dp, double dq, double de, double df)
        {
            double sumA = sumA(neighbours, agent -> agent.isSlack ? 0 : agent.alpha_g*dgp(agent, dp, dq, de, df)*(dgp(agent, /*t, */dp, dq, de, df) + 2*gp(agent)));
            return 
            		0.5*sumA
                   +0.5*alpha_g*( dgpSelf(dp, dq, de, df)*(dgpSelf(dp, dq, de, df) + 2*gp()) );
        }

        public double lagrangianChange_augTermQ(double dp, double dq, double de, double df)
        {
            return 
                     0.5*sumA(neighbours, agent -> agent.isSlack ? 0 : agent.alpha_g*dgq(agent, dp, dq, de, df)*(dgq(agent, dp, dq, de, df) + 2*gq(agent)))
                    +0.5*alpha_g*(dgqSelf(dp, dq, de, df)*(dgqSelf(dp, dq, de, df) + 2*gq()));
        }

        public double lagrangianChange_constraintTerms(double dp, double dq, double de, double df)
        {
            if(DISABLE_LAMBDA_UPDATE)
                return 0;
            
            double lgp = sumA(neighbours, agent -> agent.isSlack ? 0 : agent.lp*dgp(agent, /*t, */dp, dq, de, df))
                                +lp*dgpSelf(dp, dq, de, df);
            double lgq = sumA(neighbours, agent -> agent.isSlack ? 0 : agent.lq*dgq(agent, /*t, */dp, dq, de, df))
                                +lq*dgqSelf(dp, dq, de, df);
            return lgp + lgq;
        }

        private double dgpSelf(double dp, double dq, double de, double df)
        {
            return
                     (de*de + df*df)*G - dp
                    +(e*de*G + f*df*G + f*de*B - e*df*B)
                    +sumA(neighbours, agent -> de*agent.e*G(agent) + df*agent.f*G(agent) + df*agent.e*B(agent) - de*agent.f*B(agent))
                    +(de*e*G + df*f*G + df*e*B - de*f*B); // This last term is needed since neighbours does not include this agent
        }

        private double dgqSelf(double dp, double dq, double de, double df)
        {
            return
                    -(de*de + df*df)*B - dq
                    +(f*de*G - e*df*G - e*de*B - f*df*B)
                    +sumA(neighbours, agent -> df*agent.e*G(agent) - de*agent.f*G(agent) - de*agent.e*B(agent) - df*agent.f*B(agent))
                    +(df*e*G - de*f*G - de*e*B - df*f*B); // This last term is needed since neighbours does not include this agent
        }

        private double dgp(Agent agent, double dp, double dq, double de, double df)
        {
            double G_in = G(agent);
            double B_in = B(agent);
            return (agent.e*de*G_in + agent.f*df*G_in + agent.f*de*B_in - agent.e*df*B_in);
        }

        private double dgq(Agent agent, double dp, double dq, double de, double df)
        {
            double G_in = G(agent);
            double B_in = B(agent);
            return (agent.f*de*G_in - agent.e*df*G_in - agent.e*de*B_in - agent.f*df*B_in);
        }

        private double gp()
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
        private double gp(double p, double q, double e, double f)
        {
            if(DISABLE_POWER_FLOW_UPDATE)
                return 0;
            
            double sum = 0;
            for (Agent agent : neighbours.keySet())
            {
                sum += e*agent.e*G(agent) + f*agent.f*G(agent) + f*agent.e*B(agent) - e*agent.f*B(agent);
            }
            sum += e*e*G + f*f*G; // self
            
            return  round(sum - p, MIN_G_VALUE);
        }
        
        private double gq()
        {
            return gq(p, q, e, f);
        }
        
        private double gq(double p, double q, double e, double f)
        {
            if(DISABLE_POWER_FLOW_UPDATE)
                return 0;
            double sum1 = sumA(neighbours, agent -> f*agent.e*G(agent) - e*agent.f*G(agent) - e*agent.e*B(agent) - f*agent.f*B(agent));
            double sum = sum1 - e*e*B - f*f*B;
            
            return round(sum - q, MIN_G_VALUE);
        }
        
        /**
         * c(x) = 0.5*(p - p_max)^2 + 0.5*q^2
         * @param p Real power output.
         * @param q Reactive power output.
         * @return
         */
        private double cost(double p, double q)
        {
            if(DISABLE_COST)
                return 0;
            double d = p-p_costMin;
            
            return round(COST_MULTIPLIER*(0.5*d*d + 0.5*q*q));
        }

        private double gradL_p()
        {
            return gradC_p() - lp - alpha_g*gp() + mup + alpha_h*hp();
        }
        
//        private double gradL_p(double p, double q, double e, double f)
//        {
//            return gradC_p(p, q, e, f) - lp - alpha*gp(p, q, e, f);
//        }
        
        /**
         * Assuming C_i(p_i) = 0.5*(p_i - p_max)^2
         * => dC_i/dp_i = p_i - p_max
         * @return
         */
        private double gradC_p()
        {
            return DISABLE_COST ? 0 : COST_MULTIPLIER*(p - p_costMin);
        }
        
//        private double gradC_p(double p, double q, double e, double f)
//        {
//            return COST_MULTIPLIER*(p - p_costMin);
//        }

        private double gradL_q()
        {
            return gradC_q() - lq - alpha_g*gq() + muq + alpha_h*hq();
        }

//        private double gradL_q(double p, double q, double e, double f)
//        {
//            return gradC_q(p, q, e, f) - lq - alpha*gq(p, q, e, f);
//        }
        
        /**
         * Assuming C_i(q_i) = 0.5*q_i^2
         * => dC_i/dp_i = q_i
         * @param i
         * @return
         */
        private double gradC_q()
        {
            return DISABLE_COST ? 0 : COST_MULTIPLIER*q;
        }
        
//        private double gradC_q(double p, double q, double e, double f)
//        {
//            return COST_MULTIPLIER*q;
//        }

        private double gradL_e()
        {
            return gradL_e(p, q, e, f);
        }
        
        private double gradL_e(double p, double q, double e, double f)
        {
            double constraintTerms = gradL_e_constraintTerms(e);

            double penaltyTerms = gradL_e_augTerms(e);
            
            // grad_e h(x)(mu + alpha*h(x))
            // grad_e hp(x) = 2*sum{G_ij(e_i - e_j)}
            // grad_e hq(x) = 2*sum{-B_ij(e_i - e_j)}
            double consensusTerms = DISABLE_LINE_LOSS ? 0 :
                     2*sumA(neighbours, agent -> G(agent)*(e - agent.e))*(mup + alpha_h*hp())
                    +2*sumA(neighbours, agent -> -B(agent)*(e - agent.e))*(muq + alpha_h*hq());
            
            return
                constraintTerms
                +penaltyTerms
                +consensusTerms
                ;
        }

        protected double gradL_e_augTerms(double e)
        {
            double pt1a = sumA(neighbours, agent -> agent.isSlack ? 0 :
                               agent.alpha_g*gp(agent)*(agent.e*G(agent) + agent.f*B(agent)));
            
            
            double pt1b = sumA(neighbours, agent -> agent.isSlack ? 0 : 
                               agent.alpha_g*gq(agent)*(agent.f*G(agent) - agent.e*B(agent)));
            double sum2 = sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent));
            double pt2 = alpha_g*gp()*(sum2+2*e*G);
            double sum3 = sumA(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent));
            double pt3 = alpha_g*gq()*(sum3-2*e*B);
            double penaltyTerms = pt1a+pt1b
                                +
                                (isSlack ? 0 : 
                                (
                                    +pt2
                                    +pt3
                                ));
            return penaltyTerms;
        }

        protected double gradL_e_constraintTerms(double e)
        {
            return sumA(neighbours, agent -> agent.lp*(agent.e*G(agent) + agent.f*B(agent)))
            +sumA(neighbours, agent -> agent.lq*(agent.f*G(agent) - agent.e*B(agent)))
            +lp*(sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
                 +2*G*e)
            +lq*(sumA(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent))
                 -2*B*e);
        }

        private double gp(Agent agent)
        {
            if(USE_TRUE_G)
                return agent.gp();
            else
                return agent.gp;
        }

        private double gq(Agent agent)
        {
            if(USE_TRUE_G)
                return agent.gq();
            else
                return agent.gq;
        }

        private double gradL_f()
        {
            return gradL_f(p, q, e, f);
        }
        
        private double gradL_f(double p, double q, double e, double f)
        {
            double constraintTerms = 
                         sumA(neighbours, agent -> agent.lp*(agent.f*G(agent) - agent.e*B(agent)))
                        +sumA(neighbours, agent -> agent.lq*(-agent.e*G(agent) - agent.f*B(agent)))
                        +lp*(sumA(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent))
                             +2*G*f)
                        +lq*(sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
                             -2*B*f);

            double pt1 = sumA(neighbours, agent -> agent.isSlack ? 0 :
                               agent.alpha_g*gp(agent)*(agent.f*G(agent) - agent.e*B(agent))
                              +agent.alpha_g*gq(agent)*(-agent.e*G(agent) - agent.f*B(agent)));
            double pt2 = alpha_g*gp()*(sumA(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent)) 
                               +2*f*G);
            double pt3 = alpha_g*gq()*(sumA(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) 
                               -2*f*B);
            double penaltyTerms = pt1
                                +
                                (isSlack ? 0 : 
                                (
                                    +pt2
                                    +pt3
                                ));
            
            // grad_f h(x)(mu + alpha*h(x))
            // grad_f hp(x) = 2*sum{G_ij(f_i - f_j)}
            // grad_f hq(x) = 2*sum{-B_ij(f_i - f_j)}
            double consensusTerms = DISABLE_LINE_LOSS ? 0 :
                     2*sumA(neighbours, agent -> G(agent)*(f - agent.f))*(mup + alpha_h*hp())
                    +2*sumA(neighbours, agent -> -B(agent)*(f - agent.f))*(muq + alpha_h*hq());

            return
                constraintTerms
                +penaltyTerms
                +consensusTerms
                ;
        }

        private double G(Agent agent)
        {
            double G_ij = neighbours.get(agent).getReal();
            return G_ij;
        }

        private double B(Agent agent)
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
                if(START_WITH_OPTIMAL_POWERS)
                {
                    p = p_costMin;
                    q = 0;
                }
            }
            
            // Consensus values:
            hp = DISABLE_CONSENSUS_UPDATE ? 0 : p 
                 +
                 (DISABLE_LINE_LOSS ? 0 : 0.5*sumA(neighbours, agent -> lineLoss_p(agent, e, f)));
            hq = DISABLE_CONSENSUS_UPDATE ? 0 : q 
                 +
                 (DISABLE_LINE_LOSS ? 0 : 0.5*sumA(neighbours, agent -> lineLoss_q(agent, e, f)));
            mup = 0;
            muq = 0;
            wp = 0;
            wq = 0;
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
            buf.append("\n\talpha_g=");
            buf.append(alpha_g);
            buf.append("\n\talpha_h=");
            buf.append(alpha_h);
            buf.append("\n\tt=");
            buf.append(lastT);
            buf.append("\n");
            buf.append("\th=("+hp+","+hq+")\n");
            buf.append("\tw=("+wp+","+wq+")\n");
            buf.append("\tmu=("+mup+","+muq+")\n");
                
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

        private void compareGlobalLocal(double dp, double dq, double de,
                double df)
        {
            System.out.println("Compare Global/Local:");
            
            double globalDp = Sandbox018.this.gradL_p(index);
            System.out.println(globalDp + ",?=," + dp);

            double globalDq = Sandbox018.this.gradL_q(index);
            System.out.println(globalDq + ",?=," + dq);
            
            double globalDe = Sandbox018.this.gradL_e(index);
            System.out.println(globalDe + ",?=," + de);

            double globalDf = Sandbox018.this.gradL_f(index);
            System.out.println(globalDf + ",?=," + df);
        }

        private void compareChange(
        		double dp, double dq, double de, double df, // Change in x
        		double L_improvement, // Change in L(x,l)
        		double costImp, double lgImp, double muhImp, double gAugpImp, double gAugqImp, double hAugImp // Changes in parts of L()
                )
        {
            final double SMALL = 1e-3;
            
            RealVector _p = new ArrayRealVector(Sandbox018.this.p);
            RealVector _q = new ArrayRealVector(Sandbox018.this.q);
            RealVector _e = new ArrayRealVector(Sandbox018.this.e);
            RealVector _f = new ArrayRealVector(Sandbox018.this.f);
            RealVector _hp = new ArrayRealVector(Sandbox018.this.hp);
            RealVector _hq = new ArrayRealVector(Sandbox018.this.hq);

            _p.setEntry(index, p);
            _q.setEntry(index, q);
            _e.setEntry(index, e);
            _f.setEntry(index, f);
            _hp.setEntry(index, hp);
            _hq.setEntry(index, hq);
            
//            double Lx = lagrangian(_p, _q, _e, _f, 
//                    _hp, _hq,
//                    Sandbox018.this.mup, Sandbox018.this.muq, 
//                    Sandbox018.this.lp, Sandbox018.this.lq);
            double gAugp_x = lagrangian_gAug_p(_p, _q, _e, _f);
            double gAugq_x = lagrangian_gAug_q(_p, _q, _e, _f);
            double lg_x = lagrangian_lg(_p, _q, _e, _f, Sandbox018.this.lp, Sandbox018.this.lq);
            double muh_x = lagrangian_muh(_hp, _hq, Sandbox018.this.mup, Sandbox018.this.muq);
            double hAug_x = lagrangian_hAug(_hp, _hq);
            double cost_x = Sandbox018.this.cost(_p, _q);
			double Lx = cost_x + lg_x + (gAugp_x + gAugq_x)
                    + muh_x
                    + hAug_x;
            
            if(Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq, Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq) != Lx)
                throw new RuntimeException("Sum of Lagrange parts != lagrangian()");
            
            _p.setEntry(index, p+dp);
            _q.setEntry(index, q+dq);
            _e.setEntry(index, e+de);
            _f.setEntry(index, f+df);
            _hp.setEntry(index, hp+dhp(dp, de, df));
            _hq.setEntry(index, hq+dhq(dq, de, df));
//            double Ldx = lagrangian(_p, _q, _e, _f, 
//                    _hp, _hq,
//                    Sandbox018.this.mup, Sandbox018.this.muq, 
//                    Sandbox018.this.lp, Sandbox018.this.lq);
            
            double gAugp_dx = lagrangian_gAug_p(_p, _q, _e, _f);
            double gAugq_dx = lagrangian_gAug_q(_p, _q, _e, _f);
            double lg_dx = lagrangian_lg(_p, _q, _e, _f, Sandbox018.this.lp, Sandbox018.this.lq);
            double muh_dx = lagrangian_muh(_hp, _hq, Sandbox018.this.mup, Sandbox018.this.muq);
            double hAug_dx = lagrangian_hAug(_hp, _hq);
            double cost_dx = Sandbox018.this.cost(_p, _q);
			double Ldx = cost_dx + lg_dx + (gAugp_dx + gAugq_dx)  
                    + muh_dx
                    + hAug_dx;
            
            if(Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq, Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq) != Ldx)
                throw new RuntimeException("Sum of Lagrange parts != lagrangian()");
            
			double diff_gAugp = gAugp_dx - gAugp_x;
			double diff_gAugq = gAugq_dx - gAugq_x;
			double diff_lg = lg_dx - lg_x;
			double diff_muh = muh_dx - muh_x;
			double diff_hAug = hAug_dx - hAug_x;
			double diff_cost = cost_dx - cost_x;
            double diff = Ldx-Lx;
            
            boolean failed = false;
            if(abs(diff_gAugp-gAugpImp)/max(abs(gAugp_x), abs(gAugp_dx)) > SMALL)
            {
            	System.err.println("diff_gAugp != 0, "+diff_gAugp+","+gAugpImp);
            	failed = true;
            }
            
            if(abs(diff_gAugq-gAugqImp)/max(abs(gAugq_x), abs(gAugq_dx)) > SMALL)
            {
            	System.err.println("diff_gAugq != 0, "+diff_gAugq+","+gAugqImp);
            	failed = true;
            }

            if(abs(diff_lg-lgImp)/max(abs(lg_x), abs(lg_dx)) > SMALL)
            {
            	System.err.println("diff_lg != 0, "+diff_lg+","+lgImp);
            	failed = true;
            }

            if(abs(diff_muh-muhImp)/max(abs(muh_x), abs(muh_dx)) > SMALL)
            {
            	System.err.println("diff_muh != 0, "+diff_muh+","+muhImp);
            	failed = true;
            }

            if(abs(diff_hAug-hAugImp)/max(abs(hAug_x), abs(hAug_dx)) > SMALL)
            {
            	System.err.println("diff_hAug != 0, "+diff_hAug+","+hAugImp);
            	failed = true;
            }

            if(abs(diff_cost-costImp)/max(abs(cost_x), abs(cost_dx)) > SMALL)
            {
            	System.err.println("diff_cost != 0, "+diff_cost+","+costImp);
            	failed = true;
            }

			if(abs(diff-L_improvement)/max(abs(Lx), abs(Ldx)) > SMALL)
			{
			    System.err.println("diff_L != 0, diff,L_improvement , Normalised, abs(diff)/abs(Lx), abs(L_improvement)/abs(Lx)" );
				System.err.println("diff_L != 0, "+diff+","+L_improvement
						+", Normalised, "+(abs(diff)/abs(Lx))+", "+(abs(L_improvement)/abs(Lx)) );
            	Util.breakpoint();
            	failed = true;
			}
			
			if(failed)
				throw new RuntimeException("Backtracking improvement check failed: L(x-t*dx)-L(x) =, "+diff+", L_imp =, "+L_improvement);
        }

        @SuppressWarnings("deprecation")
        private void plotComparison()
        {
            System.out.println(
                    "Comparison (ef):\ndx"
//                    + ",imp_c(x)_e,imp_c(x)_f,c(x+de)-c(x),c(x_df)-c(x)"
                    + ",imp_g(x)_e,imp_g(x)_f,g(x+de)-g(x),g(x_df)-g(x)"
                    + ",imp_augP(x)_e,imp_augP(x)_f,augP(x+de)-augP(x),augP(x_df)-augP(x)"
                    + ",imp_augQ(x)_e,imp_augQ(x)_f,augQ(x+de)-augQ(x),augQ(x_df)-augQ(x)"
                    + ",imp_sum(x)_e,imp_sum(x)_f,L_sum(x+de)-L_sum(x),L_sum(x_df)-L_sum(x)"
                    + ",imp(x)_e,imp(x)_f,L(x+de)-L(x),L(x_df)-L(x)");
            for(double d = -2e-4; d <= 2.0000001e-4; d += 2e-5)
            {
                final double _d = d;
                
                // Local calcs:
                double _de_cons = lagrangianChange_constraintTerms(/*1, */0, 0, d, 0);
                double _df_cons = lagrangianChange_constraintTerms(/*1, */0, 0, 0, d);
                
                double _de_augP = lagrangianChange_augTermP(/*1, */0, 0, d, 0);
                double _df_augP = lagrangianChange_augTermP(/*1, */0, 0, 0, d);
                
                double _de_augQ = lagrangianChange_augTermQ(/*1, */0, 0, d, 0);
                double _df_augQ = lagrangianChange_augTermQ(/*1, */0, 0, 0, d);
                
                double _de_consAug =  lagrangeChange_consensusAugTerms(0, 0, d, 0);
                double _df_consAug =  lagrangeChange_consensusAugTerms(0, 0, 0, d);
                
                double _de_consCons =  lagrangeChange_consensusConstraintTerms(0, 0, d, 0);
                double _df_consCons =  lagrangeChange_consensusConstraintTerms(0, 0, 0, d);

                double _de = lagrangianChange(/*1, */0, 0, d, 0);
                double _df = lagrangianChange(/*1, */0, 0, 0, d);
                
                // Global calcs:
                RealVector _p = new ArrayRealVector(Sandbox018.this.p);
                RealVector _q = new ArrayRealVector(Sandbox018.this.q);
                RealVector _e = new ArrayRealVector(Sandbox018.this.e);
                RealVector _f = new ArrayRealVector(Sandbox018.this.f);
                RealVector _hp = new ArrayRealVector(Sandbox018.this.hp);
                RealVector _hq = new ArrayRealVector(Sandbox018.this.hq);
                RealVector gp_global = Sandbox018.this.gp(_p, _q, _e, _f);
                RealVector gq_global = Sandbox018.this.gq(_p, _q, _e, _f);
                gp_global.setEntry(slackIndex, 0);
                gq_global.setEntry(slackIndex, 0);
                double Lx_augP = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gp_global).dotProduct(gp_global);
                double Lx_augQ = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gq_global).dotProduct(gq_global);
                double Lx_augCons =  0.5*Sandbox018.this.alpha_h.ebeMultiply(_hp).dotProduct(_hp)
                                    +0.5*Sandbox018.this.alpha_h.ebeMultiply(_hq).dotProduct(_hq);
                double Lx_muh = Sandbox018.this.mup.dotProduct(_hp) + Sandbox018.this.muq.dotProduct(_hq);
                double Lx_lg = Sandbox018.this.lp.dotProduct(gp_global) + Sandbox018.this.lq.dotProduct(gq_global);
                double Lx = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);

                _e.setEntry(index, e+d);
                double dhp = sumA(neighbours, agent -> lineLoss_p(agent, e+_d, f) - lineLoss_p(agent, e, f));
                double dhq = sumA(neighbours, agent -> lineLoss_q(agent, e+_d, f) - lineLoss_q(agent, e, f));
                _hp.setEntry(index, hp+dhp);
                _hq.setEntry(index, hq+dhq);
                gp_global = Sandbox018.this.gp(_p, _q, _e, _f);
                gq_global = Sandbox018.this.gq(_p, _q, _e, _f);
                gp_global.setEntry(slackIndex, 0);
                gq_global.setEntry(slackIndex, 0);
                double Lx_augP_e = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gp_global).dotProduct(gp_global);
                double Lx_augQ_e = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gq_global).dotProduct(gq_global);
                double Lx_augCons_e =  0.5*Sandbox018.this.alpha_h.ebeMultiply(_hp).dotProduct(_hp)
                                    +0.5*Sandbox018.this.alpha_h.ebeMultiply(_hq).dotProduct(_hq);
                double Lx_muh_e = Sandbox018.this.mup.dotProduct(_hp) + Sandbox018.this.muq.dotProduct(_hq);
                double Lx_lg_e = Sandbox018.this.lp.dotProduct(gp_global) + Sandbox018.this.lq.dotProduct(gq_global);
                double Lx_e = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
                
                _e.setEntry(index, e);
                _f.setEntry(index, f+d);
                dhp = sumA(neighbours, agent -> lineLoss_p(agent, e, f+_d) - lineLoss_p(agent, e, f));
                dhq = sumA(neighbours, agent -> lineLoss_q(agent, e, f+_d) - lineLoss_q(agent, e, f));
                _hp.setEntry(index, hp+dhp);
                _hq.setEntry(index, hq+dhq);
                gp_global = Sandbox018.this.gp(_p, _q, _e, _f);
                gq_global = Sandbox018.this.gq(_p, _q, _e, _f);
                gp_global.setEntry(slackIndex, 0);
                gq_global.setEntry(slackIndex, 0);
                double Lx_augP_f = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gp_global).dotProduct(gp_global);
                double Lx_augQ_f = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gq_global).dotProduct(gq_global);
                double Lx_augCons_f =  0.5*Sandbox018.this.alpha_h.ebeMultiply(_hp).dotProduct(_hp)
                                      +0.5*Sandbox018.this.alpha_h.ebeMultiply(_hq).dotProduct(_hq);
                double Lx_muh_f = Sandbox018.this.mup.dotProduct(_hp) + Sandbox018.this.muq.dotProduct(_hq);
                double Lx_lg_f = Sandbox018.this.lp.dotProduct(gp_global) + Sandbox018.this.lq.dotProduct(gq_global);
                double Lx_f = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
                
                System.out.println(d+","+_de_cons+","+_df_cons+","+(Lx_lg_e-Lx_lg)+","+(Lx_lg_f-Lx_lg)
                        +","+_de_augP+","+_df_augP+","+(Lx_augP_e-Lx_augP)+","+(Lx_augP_f-Lx_augP)
                        +","+_de_augQ+","+_df_augQ+","+(Lx_augQ_e-Lx_augQ)+","+(Lx_augQ_f-Lx_augQ)
                        +","+(_de_cons+_de_augP+_de_augQ+_de_consAug+_de_consCons)+","+(_df_cons+_df_augP+_df_augQ+_df_consAug+_df_consCons)+","
                            +((Lx_lg_e+Lx_augP_e+Lx_augQ_e+Lx_augCons_e+Lx_muh_e)-(Lx_lg+Lx_augP+Lx_augQ+Lx_augCons+Lx_muh))+","+((Lx_lg_f+Lx_augP_f+Lx_augQ_f+Lx_augCons_f+Lx_muh_f)-(Lx_lg+Lx_augP+Lx_augQ+Lx_augCons+Lx_muh))
                        +","+_de+","+_df+","+(Lx_e-Lx)+","+(Lx_f-Lx));
            }
            
            System.out.println(
                    "\n\n\n\n\n\n\n\n\n\nComparison (pq):\ndx"
//                    + ",imp_c(x)_e,imp_c(x)_f,c(x+de)-c(x),c(x_df)-c(x)"
                    + ",imp_g(x)_p,imp_g(x)_q,g(x+dp)-g(x),g(x_dq)-g(x)"
                    + ",imp_augP(x)_p,imp_augP(x)_q,augP(x+dp)-augP(x),augP(x_dq)-augP(x)"
                    + ",imp_augQ(x)_p,imp_augQ(x)_q,augQ(x+dp)-augQ(x),augQ(x_dq)-augQ(x)"
                    + ",imp_cost_p,imp_cost_q,cost(x+dp)-cost(x),cost(x+dq)-cost(x)"
                    + ",imp_sum(x)_p,imp_sum(x)_q,L_sum(x+dp)-L_sum(x),L_sum(x_dq)-L_sum(x)"
                    + ",imp(x)_p,imp(x)_q,L(x+dp)-L(x),L(x_dq)-L(x)");
            for(double d = -1e-7; d <= 1.0000001e-7; d += 1e-8)
            {
                // Local calcs:
                double _dp_cost, _dq_cost;
                if(isGenerator)
                {
                    _dp_cost = cost(p+d, q) - cost(p, q);
                    _dq_cost = cost(p, q+d) - cost(p, q);
                }
                else
                {
                    _dp_cost = 0;
                    _dq_cost = 0;
                }
                
                double _dp_cons =  lagrangianChange_constraintTerms(d, 0, 0, 0);
                double _dq_cons =  lagrangianChange_constraintTerms(0, d, 0, 0);
                
                double _dp_augP = lagrangianChange_augTermP(d, 0, 0, 0);
                double _dq_augP = lagrangianChange_augTermP(0, d, 0, 0);
                
                double _dp_augQ = lagrangianChange_augTermQ(d, 0, 0, 0);
                double _dq_augQ = lagrangianChange_augTermQ(0, d, 0, 0);
                
                double _dp_consAug =  lagrangeChange_consensusAugTerms(d, 0, 0, 0);
                double _dq_consAug =  lagrangeChange_consensusAugTerms(0, d, 0, 0);
                
                double _dp_consCons =  lagrangeChange_consensusConstraintTerms(d, 0, 0, 0);
                double _dq_consCons =  lagrangeChange_consensusConstraintTerms(0, d, 0, 0);

                double _dp = lagrangianChange(d, 0, 0, 0);
                double _dq = lagrangianChange(0, d, 0, 0);
                
                // Global calcs:
                RealVector _p = new ArrayRealVector(Sandbox018.this.p);
                RealVector _q = new ArrayRealVector(Sandbox018.this.q);
                RealVector _e = new ArrayRealVector(Sandbox018.this.e);
                RealVector _f = new ArrayRealVector(Sandbox018.this.f);
                RealVector _hp = new ArrayRealVector(Sandbox018.this.hp);
                RealVector _hq = new ArrayRealVector(Sandbox018.this.hq);
                RealVector gp_global = Sandbox018.this.gp(_p, _q, _e, _f);
                RealVector gq_global = Sandbox018.this.gq(_p, _q, _e, _f);
                gp_global.setEntry(slackIndex, 0);
                gq_global.setEntry(slackIndex, 0);
                double Lx_augP = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gp_global).dotProduct(gp_global);
                double Lx_augQ = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gq_global).dotProduct(gq_global);
                double Lx_augCons =  0.5*Sandbox018.this.alpha_h.ebeMultiply(_hp).dotProduct(_hp)
                                      +0.5*Sandbox018.this.alpha_h.ebeMultiply(_hq).dotProduct(_hq);
                double Lx_muh = Sandbox018.this.mup.dotProduct(_hp) + Sandbox018.this.muq.dotProduct(_hq);
                double Lx_lg = Sandbox018.this.lp.dotProduct(gp_global) + Sandbox018.this.lq.dotProduct(gq_global);
                double Lx_cost = Sandbox018.this.cost(_p, _q);
                double Lx = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);

                _p.setEntry(index, p+d);
                _hp.setEntry(index, hp+d);
                gp_global = Sandbox018.this.gp(_p, _q, _e, _f);
                gq_global = Sandbox018.this.gq(_p, _q, _e, _f);
                gp_global.setEntry(slackIndex, 0);
                gq_global.setEntry(slackIndex, 0);
                double Lx_augP_p = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gp_global).dotProduct(gp_global);
                double Lx_augQ_p = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gq_global).dotProduct(gq_global);
                double Lx_augCons_p =  0.5*Sandbox018.this.alpha_h.ebeMultiply(_hp).dotProduct(_hp)
                                      +0.5*Sandbox018.this.alpha_h.ebeMultiply(_hq).dotProduct(_hq);
                double Lx_muh_p = Sandbox018.this.mup.dotProduct(_hp) + Sandbox018.this.muq.dotProduct(_hq);
                double Lx_lg_p = Sandbox018.this.lp.dotProduct(gp_global) + Sandbox018.this.lq.dotProduct(gq_global);
                double Lx_cost_p = Sandbox018.this.cost(_p, _q);
                double Lx_p = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
                
                _p.setEntry(index, p);
                _hp.setEntry(index, hp);
                _q.setEntry(index, q+d);
                _hq.setEntry(index, hq+d);
                gp_global = Sandbox018.this.gp(_p, _q, _e, _f);
                gq_global = Sandbox018.this.gq(_p, _q, _e, _f);
                gp_global.setEntry(slackIndex, 0);
                gq_global.setEntry(slackIndex, 0);
                double Lx_augP_q = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gp_global).dotProduct(gp_global);
                double Lx_augQ_q = 0.5*Sandbox018.this.alpha_g.ebeMultiply(gq_global).dotProduct(gq_global);
                double Lx_augCons_q =  0.5*Sandbox018.this.alpha_h.ebeMultiply(_hp).dotProduct(_hp)
                                      +0.5*Sandbox018.this.alpha_h.ebeMultiply(_hq).dotProduct(_hq);
                double Lx_muh_q = Sandbox018.this.mup.dotProduct(_hp) + Sandbox018.this.muq.dotProduct(_hq);
                double Lx_lg_q = Sandbox018.this.lp.dotProduct(gp_global) + Sandbox018.this.lq.dotProduct(gq_global);
                double Lx_cost_q = Sandbox018.this.cost(_p, _q);
                double Lx_q = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
                
                System.out.println(d+","+_dp_cons+","+_dq_cons+","+(Lx_lg_p-Lx_lg)+","+(Lx_lg_q-Lx_lg)
                        +","+_dp_augP+","+_dq_augP+","+(Lx_augP_p-Lx_augP)+","+(Lx_augP_q-Lx_augP)
                        +","+_dp_augQ+","+_dq_augQ+","+(Lx_augQ_p-Lx_augQ)+","+(Lx_augQ_q-Lx_augQ)
                        +","+_dp_cost+","+_dq_cost+","+(Lx_cost_p-Lx_cost)+","+(Lx_cost_q-Lx_cost)
                        +","+(_dp_cons+_dp_augP+_dp_augQ+_dp_cost+_dp_consAug+_dp_consCons)+","+(_dq_cons+_dq_augP+_dq_augQ+_dq_cost+_dq_consAug+_dq_consCons)+","
                            +((Lx_lg_p+Lx_augP_p+Lx_augQ_p+Lx_cost_p+Lx_augCons_p+Lx_muh_p)-(Lx_lg+Lx_augP+Lx_augQ+Lx_cost+Lx_augCons+Lx_muh))+","+((Lx_lg_q+Lx_augP_q+Lx_augQ_q+Lx_cost_q+Lx_augCons_q+Lx_muh_q)-(Lx_lg+Lx_augP+Lx_augQ+Lx_cost+Lx_augCons+Lx_muh))
                        +","+_dp+","+_dq+","+(Lx_p-Lx)+","+(Lx_q-Lx));
            }
            
            Util.nullop();
        }

        private void plotImprovements()
        {
            System.out.println("Improvements:\ndx,imp_p,imp_q,imp_e,imp_f,L(x+dp)-L(x),L(x+dq)-L(x),L(x+de)-L(x),L(x_df)-L(x)");
            double low = -2e-4;
            double high = 2.0000001e-4;
            double step = 2e-5;
            for(double d = low; d <= high; d += step)
            {
                final double _d = d;
                
                double _dp = lagrangianChange(d, 0, 0, 0);
                double _dq = lagrangianChange(0, d, 0, 0);
                double _de = lagrangianChange(0, 0, d, 0);
                double _df = lagrangianChange(0, 0, 0, d);
                
                RealVector _p = new ArrayRealVector(Sandbox018.this.p);
                RealVector _q = new ArrayRealVector(Sandbox018.this.q);
                RealVector _e = new ArrayRealVector(Sandbox018.this.e);
                RealVector _f = new ArrayRealVector(Sandbox018.this.f);
                RealVector _hp = new ArrayRealVector(Sandbox018.this.hp);
                RealVector _hq = new ArrayRealVector(Sandbox018.this.hq);
                
                double Lx = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
                _p.setEntry(index, p+d);
                _hp.setEntry(index, hp+d);
                double Ldx_p = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
                
                _p.setEntry(index, p);
                _hp.setEntry(index, hp);
                _q.setEntry(index, q+d);
                _hq.setEntry(index, hq+d);
                double Ldx_q = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);

                _q.setEntry(index, q);
                double dhp = sumA(neighbours, agent -> lineLoss_p(agent, e+_d, f) - lineLoss_p(agent, e, f));
                double dhq = sumA(neighbours, agent -> lineLoss_q(agent, e+_d, f) - lineLoss_q(agent, e, f));
                _hp.setEntry(index, hp+dhp);
                _hq.setEntry(index, hq+dhq);
                _e.setEntry(index, e+d);
                double Ldx_e = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
                
                _e.setEntry(index, e);
                dhp = sumA(neighbours, agent -> lineLoss_p(agent, e, f+_d) - lineLoss_p(agent, e, f));
                dhq = sumA(neighbours, agent -> lineLoss_q(agent, e, f+_d) - lineLoss_q(agent, e, f));
                _hp.setEntry(index, hp+dhp);
                _hq.setEntry(index, hq+dhq);
                _f.setEntry(index, f+d);
                double Ldx_f = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
                        Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
                
                System.out.println(d+","+_dp+","+_dq+","+_de+","+_df+","
                                   +(Ldx_p-Lx)+","+(Ldx_q-Lx)+","+(Ldx_e-Lx)+","+(Ldx_f-Lx));
            }
            
            double maxDiff = 0;
            for(double dp = low; dp <= high; dp += step)
            {
                for(double dq = low; dq <= high; dq += step)
                {
                    double change = lagrangianChange(/*1, */dp, dq, 0, 0);
                    
//                    RealVector _p = new ArrayRealVector(Sandbox018.this.p);
//                    RealVector _q = new ArrayRealVector(Sandbox018.this.q);
//                    RealVector _e = new ArrayRealVector(Sandbox018.this.e);
//                    RealVector _f = new ArrayRealVector(Sandbox018.this.f);
//                    RealVector _hp = new ArrayRealVector(Sandbox018.this.hp);
//                    RealVector _hq = new ArrayRealVector(Sandbox018.this.hq);
                    
//                    double Lx = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
//                            Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
                    double Lx = lagrangian(p, q, e, f);
//                    _p.setEntry(index, p+dp);
//                    _hp.setEntry(index, hp+dp);
//                    _q.setEntry(index, q+dq);
//                    _hq.setEntry(index, hq+dq);
//                    double Ldx_p = Sandbox018.this.lagrangian(_p, _q, _e, _f, _hp, _hq,
//                            Sandbox018.this.mup, Sandbox018.this.muq, Sandbox018.this.lp, Sandbox018.this.lq);
                    double Ldx_p = lagrangian(p+dp, q+dq, e, f);
                    
                    maxDiff = max(maxDiff, abs(change-(Ldx_p-Lx)));
                }
            }
            
            System.out.println("maxDiff,"+maxDiff);
            
            Util.nullop();
        }

        private void checkGradient(double dp, double dq, double de, double df)
        {
            double delta = 1e-7;
            
            System.out.println("Check Gradient:\nindex,dp,,dq,,de,,df,,de_cons,,de_aug");
            
            // Estimate p slope:
            double left = lagrangianChange(/*1, */-delta, 0, 0, 0);
            double right = lagrangianChange(/*1, */delta, 0, 0, 0);
            double pSlope = (right-left)/(2*delta);

            // Estimate q slope:
            left = lagrangianChange(/*1, */0, -delta, 0, 0);
            right = lagrangianChange(/*1, */0, delta, 0, 0);
            double qSlope = (right-left)/(2*delta);
            
            // Estimate e slope:
            left = lagrangianChange(/*1, */0, 0, -delta, 0);
            right = lagrangianChange(/*1, */0, 0, delta, 0);
            double eSlope = (right-left)/(2*delta);

            // Estimate f slope:
            left = lagrangianChange(/*1, */0, 0, 0, -delta);
            right = lagrangianChange(/*1, */0, 0, 0, delta);
            double fSlope = (right-left)/(2*delta);
            
            // Log:
            System.out.print(index
                    +","+dp+","+pSlope
                    +","+dq+","+qSlope
                    +","+de+","+eSlope
                    +","+df+","+fSlope);
            
            // Check components of gradient:
            left = lagrangianChange_constraintTerms(0, 0, -delta, 0);
            right = lagrangianChange_constraintTerms(0, 0, delta, 0);
            double lg_eSlope = (right-left)/(2*delta);
            double de_constraintTerms = gradL_e_constraintTerms(e);
            
            left = lagrangianChange_augTermP(0, 0, -delta, 0) + lagrangianChange_augTermQ(0, 0, -delta, 0);
            right = lagrangianChange_augTermP(0, 0, delta, 0) + lagrangianChange_augTermQ(0, 0, delta, 0);
            double aug_eSlope = (right-left)/(2*delta);
            double de_augTerms = gradL_e_augTerms(e);
            
            System.out.println(","+de_constraintTerms+","+lg_eSlope+","+de_augTerms+","+aug_eSlope);
        }

        public Complex getAverageNetworkPowerEstimateOffset()
        {
            return new Complex(wp, wq);
        }

        public Complex getAverageNetworkPowerEstimateMultiplier()
        {
            return new Complex(mup, muq);
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
        
        logConfig();
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
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox016-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox016-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox016-loadData.csv");
        grid = GridGenerator.loadGrid("IEEE123-35bus", lineConfig, lineData, switchData, loadData);
        
        // Add slack bus:
        SlackSource slack = new SlackSource();
        slack.setVoltage(new Complex(BASE_VOLTAGE*SLACK_VOLTAGE));
        slack.setName(SLACK_SOURCE);
        grid.getBus(GridGenerator.BUS_PREFIX+149).addChild(slack);
        
        // Transformer:
        grid.setBaseVoltage(BASE_VOLTAGE);
        grid.setBasePower(BASE_POWER);
        
        // Add DG:
        int[] dgIndices = new int[]{1,3,8,14,18,26,29};
        for (int i = 0; i < dgIndices.length; i++)
        {
            /*DistributedSource dg = */addDG(dgIndices[i], 500e3, 100e3, 50e3);
            
//            // Init with optimal powers:
//            if(START_WITH_OPTIMAL_POWERS)
//            {
//                dg.setPowerOutput(0, 0); // WARNING These are hard coded since optimal isn't calculated until initVariables.
//            }
        }
//        addDG(1,  500e3, 0, 0);
//        addDG(3,  500e3, 0, 0);
//        addDG(8,  500e3, 0, 0);
//        addDG(14, 500e3, 0, 0);
//        addDG(18, 500e3, 0, 0);
//        addDG(26, 500e3, 0, 0);
//        addDG(29, 500e3, 0, 0);
        
//        addDG(1, 500e3, 100e3, 50e3);
//        addDG(3, 500e3, 100e3, 50e3);
//        addDG(8, 500e3, 100e3, 100e3);
//        addDG(14, 500e3, 100e3, 50e3);
//        addDG(18, 500e3, 100e3, 100e3);
//        addDG(26, 500e3, 100e3, 50e3);
//        addDG(29, 500e3, 100e3, 50e3);
        
//      GridDisplay.showInFrame(grid).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    /**
     * Creates and adds a DG, and removes any loads.
     * @param i
     * @param pMax
     * @param p
     * @param q
     * @return 
     */
    protected DistributedSource addDG(int i, double pMax, double p, double q)
    {
        Bus bus = grid.getBus(GridGenerator.BUS_PREFIX+i);
        for (Load load : bus.getLoads())
        {
            bus.removeChild(load);
        }
        DistributedSource dg = makeDG("DG"+i, pMax, p, q);
        bus.addChild(dg);
        return dg;
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
        alpha_g = new ArrayRealVector(dimension, INITIAL_G_AUG_SCALE);
        alpha_h = new ArrayRealVector(dimension, INITIAL_H_AUG_SCALE);
        t = new ArrayRealVector(dimension, 0);
        
        // Init consensus vectors:
        hp = new ArrayRealVector(dimension, 0);
        hq = new ArrayRealVector(dimension, 0);
        wp = new ArrayRealVector(dimension, 0);
        wq = new ArrayRealVector(dimension, 0);
        mup = new ArrayRealVector(dimension, 0);
        muq = new ArrayRealVector(dimension, 0);
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
            agent.alpha_g = this.alpha_g.getEntry(i);
            agent.alpha_h = this.alpha_h.getEntry(i);
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
        
        // Update global values according to initialised agent values:
        update();
        
        // Validate initialisation:
        validateConsensus();
    }
    
    
    //// Iterations ////

    public void loop()
    {
        Random agentSelector = new Random(0);
        Timer timer = Timer.startNewTimer();
        for(int k = 0; k < K; ++k)
        {
            int count = 0;
if(k == 2000)
    Util.nullop();
            for (Agent agent : agents)
            {
                if(agentSelector.nextDouble() > AGENT_SELECTION_PROBABILITY)
                    continue;

                agent.step(k);
                update();
                ++count;
                
                debugStepEnd(k, agent);
            }

            // Validate consensus average:
            validateConsensus();
            
            cpuTime = timer.lapNano()/1e6/count;
            
            // Estimate convergence:
//            addSamples();

            // Log:
            if(k%DEBUG_RATE == 0)
                debug(k);
        }
    }

//    double pscale = 1;
//    double qscale = 1;
//    double escale = 1;
//    double fscale = 1;
//    double gscale = 1;
//    double hscale = 1;
//    private void addSamples()
//    {
//        int dimension = p.getDimension();
//        
//        variance_pNorm.addSample(p.getNorm()/pscale);
//        variance_qNorm.addSample(q.getNorm()/qscale);
//        variance_eNorm.addSample(e.getNorm()/escale);
//        variance_fNorm.addSample(f.getNorm()/fscale);
//        
//        variance_gpNorm.addSample(vector(dimension, this::gp).getNorm()/gscale);
//        variance_gqNorm.addSample(vector(dimension, this::gq).getNorm()/gscale);
//        
//        variance_hpNorm.addSample(hp.getNorm()/hscale);
//        variance_hqNorm.addSample(hq.getNorm()/hscale);
//    }

    protected void validateConsensus()
    {
        if(!VALIDATE_CONSENSUS)
            return;
        
        double mismatch_p = powerMismatch_p();
        double hSum_p = hsum_p();
        if(abs(mismatch_p - hSum_p) > 1e-9)
            throw new RuntimeException("Consensus power mismatch (p) incorrect:, mismatch=,"+mismatch_p+", sum=,"+hSum_p);
        
        double mismatch_q = powerMismatch_q();
        double hSum_q = hsum_q();
        if(abs(mismatch_q - hSum_q) > 1e-9)
            throw new RuntimeException("Consensus power mismatch (q) incorrect:, mismatch=,"+mismatch_q+", sum=,"+hSum_q);
    }

    protected double hsum_q()
    {
        return sumA(agents, agent -> agent.hq-agent.wq);
    }

    protected double hsum_p()
    {
        return sumA(agents, agent -> agent.hp-agent.wp);
    }

    protected void debugStepEnd(int k, Agent agent)
    {
        if(agent.index == 2 && k >= 2800)
        {
            Util.nullop();
//            double global = gradL_e(agent.index);
//            double local = agent.normGradL_e();
//            System.out.println(global+","+local);
        }
    }
    
//private void checkConsensusAverage()
//{
//    for (Agent agent : agents)
//    {
//        Complex h = agent.getAverageNetworkPowerEstimate();
//        System.out.print(h.getReal());
//        System.out.print(",");
//        System.out.print(h.getImaginary());
//        System.out.print(",");
//    }
//    System.out.print(",");
//    for (Agent agent : agents)
//    {
//        Complex w = agent.getAverageNetworkPowerEstimateOffset();
//        System.out.print(-w.getReal());
//        System.out.print(",");
//        System.out.print(-w.getImaginary());
//        System.out.print(",");
//    }
//    System.out.println();
//}

    /**
     * Update global variables from agent values.
     */
    protected void update()
    {
        epsilon_pq = 0;
        epsilon_ef = 0;
        
        for (Agent agent : agents)
        {
            int i = agent.index;
            p.setEntry(i, agent.p);
            q.setEntry(i, agent.q);
            e.setEntry(i, agent.e);
            f.setEntry(i, agent.f);
            lp.setEntry(i, agent.lp);
            lq.setEntry(i, agent.lq);

            alpha_g.setEntry(i, agent.alpha_g);
            alpha_h.setEntry(i, agent.alpha_h);
            t.setEntry(i, agent.lastT);
            
            Complex h = agent.getAverageNetworkPowerEstimate();
            hp.setEntry(i, h.getReal());
            hq.setEntry(i, h.getImaginary());
            Complex w = agent.getAverageNetworkPowerEstimateOffset();
            wp.setEntry(i, w.getReal());
            wq.setEntry(i, w.getImaginary());
            Complex mu = agent.getAverageNetworkPowerEstimateMultiplier();
            mup.setEntry(i, mu.getReal());
            muq.setEntry(i, mu.getImaginary());
            
            epsilon_pq = max(epsilon_pq, agent.epsilon_pq);
            epsilon_ef = max(epsilon_ef, agent.epsilon_ef);
        }
        
        if(p.isNaN() || q.isNaN() || e.isNaN() || f.isNaN())
            Util.nullop();
    }
    
    
    //// Convenience methods ////

    public String format(double d)
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
        if(!DEBUG)
            return;
        
        Map<String, Integer> numbers = results.getBusNumbers();
        debugHeaderRow(numbers);
        out.println();
    }

	protected void debugHeaderRow(Map<String, Integer> numbers)
	{
		out.print("k,");
        
        int dimension = debugDimension();
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
        debugLogNames(numbers, dimension, "gradP_p");
        debugLogNames(numbers, dimension, "gradP_q");
        debugLogNames(numbers, dimension, "gradP_e");
        debugLogNames(numbers, dimension, "gradP_f");
        out.print("C(pq),L_lg,L_gAug,L_muh,L_hAug,L(pql),");
        debugLogNames(numbers, dimension, "t");
        debugLogNames(numbers, dimension, "alpha_g");
        debugLogNames(numbers, dimension, "alpha_h");
        out.print("||grad_p||,||grad_q||,||grad_e||,||grad_f||,CPUTime,");
        debugLogNames(numbers, dimension, "hp");
        out.print("hpSum,h_p(x),");
        debugLogNames(numbers, dimension, "wp");
        out.print("wpSum,");
        debugLogNames(numbers, dimension, "hq");
        out.print("hqSum,h_q(x),");
        debugLogNames(numbers, dimension, "wq");
        out.print("wqSum,");
        debugLogNames(numbers, dimension, "mup");
        debugLogNames(numbers, dimension, "muq");
        out.print("epsilon_pq,epsilon_ef,");
        out.print("gNorm,hNorm,gradNorm,");
        out.print("lambda*,");
	}

    public void debugLogNames(Map<String, Integer> numbers, int dimension, String type)
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
        if(!DEBUG)
            return;
        
        debugRow(k);
        
        out.println();
    }

	protected void debugRow(int k)
	{
		out.print(k);
        out.print(',');
        int dimension = debugDimension();
        
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
        
        // Power Aug Lag derivatives:
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
        
        // Voltage Aug Lag derivatives:
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
        
        // Power Non-aug derivatives:
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gradP_p(i)));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gradP_q(i)));
            out.print(',');
        }
        
        // Voltage Non-aug derivatives:
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gradP_e(i)));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gradP_f(i)));
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
        out.print(format(lagrangian_lg(p, q, e, f, lp, lq)));
        out.print(',');
        out.print(format(lagrangian_gAug(p, q, e, f)));
        out.print(',');
        out.print(format(lagrangian_muh(hp, hq, mup, muq)));
        out.print(',');
        out.print(format(lagrangian_hAug(hp, hq)));
        out.print(',');
        out.print(format(lagrangian()));
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
            out.print(alpha_g.getEntry(i));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(alpha_h.getEntry(i));
            out.print(',');
        }
        
        // Augmented Lagrange Gradients:
        out.print(format(gradLNorm(this::gradL_p, dimension)));
        out.print(',');
        out.print(format(gradLNorm(this::gradL_q, dimension)));
        out.print(',');
        out.print(format(gradLNorm(i -> onXSurface_V(i) ? 0 : gradL_e(i), dimension)));
        out.print(',');
        out.print(format(gradLNorm(i -> onXSurface_V(i) ? 0 : gradL_f(i), dimension)));
        out.print(',');
        
        // Processing time:
        out.print(format(cpuTime));
        out.print(',');

        // Consensus variables:
        double hpSum = 0;
        for(int i = 0; i < dimension; ++i)
        {
            double hp_i = hp.getEntry(i);
            hpSum += hp_i;
            out.print(format(hp_i));
            out.print(',');
        }
        out.print(format(hpSum));
        out.print(',');
        out.print(format(averagePowerMismatch_p()));
        out.print(',');

        double wpSum = 0;
        for(int i = 0; i < dimension; ++i)
        {
            double wp_i = wp.getEntry(i);
            wpSum += wp_i;
            out.print(format(wp_i));
            out.print(',');
        }
        out.print(format(wpSum));
        out.print(',');
        
        double hqSum = 0;
        for(int i = 0; i < dimension; ++i)
        {
            double hq_i = hq.getEntry(i);
            hqSum += hq_i;
            out.print(format(hq_i));
            out.print(',');
        }
        out.print(format(hqSum));
        out.print(',');
        out.print(format(averagePowerMismatch_q()));
        out.print(',');

        double wqSum = 0;
        for(int i = 0; i < dimension; ++i)
        {
            double wq_i = wq.getEntry(i);
            wqSum += wq_i;
            out.print(format(wq_i));
            out.print(',');
        }
        out.print(format(wqSum));
        out.print(',');

        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(mup.getEntry(i)));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(muq.getEntry(i)));
            out.print(',');
        }
        
        // Epsilon:
        out.print(format(epsilon_pq));
        out.print(',');
        out.print(format(epsilon_ef));
        out.print(',');
        
        // Convergence measures:
        out.print(format(gNorm()));
        out.print(',');
        out.print(format(hNorm()));
        out.print(',');
        out.print(format(gradNorm()));
        out.print(',');
        
        // lambda*
        out.print(format(lambdaStar().getMaxValue()));
        out.print(',');
	}

	protected int debugDimension()
	{
		int dimension = p.getDimension();
		return Math.min(dimension, DISPLAY_DIMENSION);
	}
    
    protected RealVector lambdaStar()
    {
        // lambda* = -(gradG^{T}gradG)^{-1}gradG^{T}gradC
        RealMatrix grad_xG = grad_xG();
        RealMatrix grad_xGtrans = grad_xG.transpose();
        RealMatrix squareG = grad_xGtrans.multiply(grad_xG);
        RealVector costGradient = costGradient();
        RealVector lambdaStar = invert(squareG).multiply(grad_xGtrans).operate(costGradient).mapMultiply(-1);
        
//        // Crosscheck gradC + gradG.lambda:
//        RealVector grad = costGradient.add(grad_xG.operate(lambdaStar));
//        double gradNorm = grad.getNorm();
//        int rank = rank(grad_xG);
        
        return lambdaStar;
    }
    
    /**
     * [ dc(x)/dp ]
     * [ dc(x)/dq ]
     * [ dc(x)/de ]
     * [ dc(x)/df ]
     * 
     * @return Gradient of c(x) w.r.t. x.
     */
    public RealVector costGradient()
    {
        int dimension = p.getDimension()-1;
        return 
                        vector(dimension, i -> i < slackIndex ? gradC_p(i) : gradC_p(i+1))  // dc(x)/dp
                .append(vector(dimension, i -> i < slackIndex ? gradC_q(i) : gradC_q(i+1))) // dc(x)/dp
                .append(new ArrayRealVector(2*dimension)) // dc(x)/e = dc(x)/df = 0
                ;
    }
    
    /**
     * [ dg_p(x)/dp   dg_q(x)/dp ]
     * [ dg_p(x)/dq   dg_q(x)/dq ]
     * [ dg_p(x)/de   dg_q(x)/de ]
     * [ dg_p(x)/df   dg_q(x)/df ]
     * 
     * @return Gradient of g(x) w.r.t. x.
     */
    protected RealMatrix grad_xG()
    {
        int dimension = p.getDimension()-1; // ignore slack
        RealMatrix gradG = new Array2DRowRealMatrix(4*dimension, 2*dimension);
        
        // dg_p(x)/dp:
        for(int _i = 0; _i < dimension; ++_i)
        {
            int i = _i < slackIndex ? _i : _i+1;

            for(int _j = 0; _j < dimension; ++_j)
            {
                int j = _j < slackIndex ? _j : _j+1;

                if(i == j)
                    gradG.setEntry(i, j, -1.0);
            }
        }
        
        // dg_p(x)/dq = 0.
        
        // dg_p(x)/de:
        for(int _i = 0; _i < dimension; ++_i)
        {
            int i = _i < slackIndex ? _i : _i+1;

            for(int _j = 0; _j < dimension; ++_j)
            {
                final int j = _j < slackIndex ? _j : _j+1;

                double dg_pjdei = 
                        i != j ? e(j)*G(j,i) + f(j)*B(j,i)
                               : sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
                                 +2*G(j,j)*e(j);
                gradG.setEntry(i+2*dimension, j, dg_pjdei);
            }
        }
        
        // dg_p(x)/df:
        for(int _i = 0; _i < dimension; ++_i)
        {
            int i = _i < slackIndex ? _i : _i+1;

            for(int _j = 0; _j < dimension; ++_j)
            {
                final int j = _j < slackIndex ? _j : _j+1;

                double dg_pjdfi = 
                        i != j ? f(j)*G(j,i) - e(j)*B(j,i)
                               : sum(n -> f(n)*G(j,n) + e(n)*B(j,n), dimension, j) // n != j
                                 +2*G(j,j)*f(j);
                gradG.setEntry(i+3*dimension, j, dg_pjdfi);
            }
        }
        
        // dg_q(x)/dp = 0.
        
        // dg_q(x)/dq:
        for(int _i = 0; _i < dimension; ++_i)
        {
            int i = _i < slackIndex ? _i : _i+1;

            for(int _j = 0; _j < dimension; ++_j)
            {
                final int j = _j < slackIndex ? _j : _j+1;

                if(i == j)
                    gradG.setEntry(i+dimension, j+dimension, -1.0);
            }
        }
        
        // dg_q(x)/de:
        for(int _i = 0; _i < dimension; ++_i)
        {
            int i = _i < slackIndex ? _i : _i+1;

            for(int _j = 0; _j < dimension; ++_j)
            {
                final int j = _j < slackIndex ? _j : _j+1;

                double dg_pjdei = 
                        i != j ? f(j)*G(j,i) - e(j)*B(j,i)
                               : sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) // n != j
                                 -2*B(j,j)*e(j);
                gradG.setEntry(i+2*dimension, j+dimension, dg_pjdei);
            }
        }
        
        // dg_q(x)/df:
        for(int _i = 0; _i < dimension; ++_i)
        {
            int i = _i < slackIndex ? _i : _i+1;

            for(int _j = 0; _j < dimension; ++_j)
            {
                final int j = _j < slackIndex ? _j : _j+1;

                double dg_pjdfi = 
                        i != j ? -e(j)*G(j,i) - f(j)*B(j,i)
                               : sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
                                 -2*B(j,j)*f(j);
                gradG.setEntry(i+3*dimension, j+dimension, dg_pjdfi);
            }
        }
        
        return gradG;
    }

    protected double gNorm()
    {
        return gp(p, q, e, f).append(gq(p, q, e, f)).getNorm();
    }

    protected double hNorm()
    {
        return hp.append(hq).getNorm();
    }

    protected double gradNorm()
    {
        int dimension = p.getDimension();
        return      vector(dimension, this::gradL_p, slackIndex)
            .append(vector(dimension, this::gradL_q, slackIndex))
            .append(vector(dimension, i -> onXSurface_V(i) ? 0 : gradL_e(i), slackIndex))
            .append(vector(dimension, i -> onXSurface_V(i) ? 0 : gradL_f(i), slackIndex))
            .getNorm();
    }

    protected boolean onXSurface_V(int i)
    {
        Complex v = new Complex(e(i), f(i));
        return onSurface(v);
    }

    protected boolean onSurface(Complex v)
    {
        double abs = v.abs();
        double arg = v.getArgument();
        return abs(abs-V_MAX) < 1e-12     || abs(abs-V_MIN) < 1e-12 ||
               abs(arg-V_ARG_MAX) < 1e-12 || abs(arg-V_ARG_MIN) < 1e-12;
    }

    private double powerMismatch_p()
    {
        int dimension = p.getDimension();
        return sum(i -> (i == slackIndex ? 0 : p(i)) + 0.5*lineLoss_p(i), dimension);
    }
    
    private double averagePowerMismatch_p()
    {
        return powerMismatch_p()/p.getDimension();
    }

    protected double lineLoss_p(int i)
    {
        int dimension = p.getDimension();
        return sum(j -> lineLoss_p(i, j), dimension);
    }

    private double lineLoss_p(int i, int j)
    {
        if(DISABLE_LINE_LOSS)
            return 0;
        
        double e2 = e(i)-e(j);
        double f2 = f(i)-f(j);
        return G(i,j)*(e2*e2 + f2*f2);
    }

    private double powerMismatch_q()
    {
        int dimension = q.getDimension();
        return sum(i -> (i == slackIndex ? 0 : q(i)) + 0.5*lineLoss_q(i), dimension);
    }
    
    private double averagePowerMismatch_q()
    {
        return powerMismatch_q()/q.getDimension();
    }

    protected double lineLoss_q(int i)
    {
        return sum(j -> lineLoss_q(i, j), q.getDimension());
    }

    private double lineLoss_q(int i, int j)
    {
        if(DISABLE_LINE_LOSS)
            return 0;
        
        double e2 = e(i)-e(j);
        double f2 = f(i)-f(j);
        return -B(i,j)*(e2*e2 + f2*f2);
    }

    protected double gradLNorm(IndexedFunction gradFn, int dimension)
    {
         RealVector v = vector(dimension, gradFn);
        
        for(int i = 0; i < dimension; ++i)
        {
            if(/*i == 17 || */i == slackIndex)
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

    protected double gradL_p(int i)
    {
        if(generatorMask.getEntry(i) == 1)
        {
            double hp_i = hp.getEntry(i);
            double alpha_g_i = alpha_g.getEntry(i);
            double alpha_h_i = alpha_h.getEntry(i);
            double p = gradP_p(i);
            return  p - alpha_g_i*gp(i) + alpha_h_i*hp_i;
        } 
        else
            return 0;
    }

    protected double gradP_p(int i)
    {
        if(generatorMask.getEntry(i) == 1)
        {
            double mup_i = mup.getEntry(i);
            double p = gradC_p(i) - lp(i) + mup_i;
            return p;
        }
        else
        {
            return 0;
        }
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
        if(generatorMask.getEntry(i) == 0)
            return 0;
        return DISABLE_COST ? 0 : COST_MULTIPLIER*(p.getEntry(i) - p_costMin.getEntry(i));
    }

    protected double gradL_q(int i)
    {
        if(generatorMask.getEntry(i) == 1)
        {
            double alpha_g_i = alpha_g.getEntry(i);
            double alpha_h_i = alpha_h.getEntry(i);
            double hq_i = hq.getEntry(i);
            double p = gradP_q(i);
            return p - alpha_g_i*gq(i) + alpha_h_i*hq_i;
        } 
        else
            return 0;
    }

    protected double gradP_q(int i)
    {
        if(generatorMask.getEntry(i) == 1)
        {
            double muq_i = muq.getEntry(i);
            double p = gradC_q(i) - lq(i) + muq_i;
            return p;
        }
        else
        {
            return 0;
        }
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
        if(generatorMask.getEntry(i) == 0)
            return 0;
        return DISABLE_COST ? 0 : COST_MULTIPLIER*q.getEntry(i);
    }
    
    protected double gp(int i)
    {
        if(DISABLE_POWER_FLOW_UPDATE)
            return 0;
        int dimension = p.getDimension();
        double d = sum(n -> e(i)*e(n)*G(i,n) + f(i)*f(n)*G(i,n) + f(i)*e(n)*B(i,n) - e(i)*f(n)*B(i,n), dimension) - p(i);
        return round(d, MIN_G_VALUE);
    }
    
    protected double gq(int i)
    {
        if(DISABLE_POWER_FLOW_UPDATE)
            return 0;
        int dimension = q.getDimension();
        double sum = sum(n -> f(i)*e(n)*G(i,n) - e(i)*f(n)*G(i,n) - e(i)*e(n)*B(i,n) - f(i)*f(n)*B(i,n), dimension);
        double d = sum - q(i);
        
        return round(d, MIN_G_VALUE);
    }

    protected double gradL_e(int j)
    {
        int dimension = p.getDimension();
        double constraintTerms = gradL_e_constraints(j);
        
        double alpha_g_j = alpha_g.getEntry(j);
        double pt1a = sum(i -> i == slackIndex ? 0 : alpha_g.getEntry(i)*gp(i)*(e(i)*G(i,j) + f(i)*B(i,j)), dimension, j);
        double pt1b = sum(i -> i == slackIndex ? 0 : alpha_g.getEntry(i)*gq(i)*(f(i)*G(i,j) - e(i)*B(i,j)), dimension, j);
        double sum2 = sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j);
        double pt2 = alpha_g_j*gp(j)*(sum2 + 2*e(j)*G(j,j));
        double sum3 = sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j);
        double pt3 = alpha_g_j*gq(j)*(sum3 - 2*e(j)*B(j,j));
        double penaltyTerms = pt1a+pt1b
                            +
                            (j == slackIndex ? 0 : 
                            (
                                +pt2
                                +pt3
                            ));
        
        // grad_e h(x)(mu + alpha*h(x))
        // grad_e hp(x) = 2*sum{G_ij(e_i - e_j)}
        // grad_e hq(x) = 2*sum{-B_ij(e_i - e_j)}
        double consensusTerms;
        if(!DISABLE_LINE_LOSS)
        {
            double alpha_h_j = alpha_h.getEntry(j);
            double hp_j = hp.getEntry(j);
            double hq_j = hq.getEntry(j);
            double mup_j = mup.getEntry(j);
            double muq_j = muq.getEntry(j);
            double sum_p = sum(n -> G(n,j)*(e(j) - e(n)), dimension, j);
            double sum_q = sum(n -> -B(n,j)*(e(j) - e(n)), dimension, j);
            consensusTerms = 
                     2*sum_p*mup_j
                    +2*sum_q*muq_j
                    +2*sum_p*alpha_h_j*hp_j
                    +2*sum_q*alpha_h_j*hq_j;
        }
        else
        {
            consensusTerms = 0;
        }
        
        return
            constraintTerms
            +penaltyTerms
            +consensusTerms
            ;
    }
    
    protected double gradP_e(int j)
    {
        double constraintsTerms = gradL_e_constraints(j);
        
        double consensusTerms;
        int dimension = p.getDimension();
        if(!DISABLE_LINE_LOSS)
        {
            double mup_j = mup.getEntry(j);
            double muq_j = muq.getEntry(j);
            double sum_p = sum(n -> G(n,j)*(e(j) - e(n)), dimension, j);
            double sum_q = sum(n -> -B(n,j)*(e(j) - e(n)), dimension, j);
            consensusTerms = 
                     2*sum_p*mup_j
                    +2*sum_q*muq_j;
        }
        else
        {
            consensusTerms = 0;
        }
        
        return constraintsTerms + consensusTerms;
    }

    protected double gradL_e_constraints(int j)
    {
        int dimension = p.getDimension();
        
        return   sum(i -> lp(i)*(e(i)*G(i,j) + f(i)*B(i,j)), dimension, j) // i != j
                +sum(i -> lq(i)*(f(i)*G(i,j) - e(i)*B(i,j)), dimension, j) // i != j
                +lp(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
                        +2*G(j,j)*e(j))
                +lq(j)*(sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) // n != j
                        -2*B(j,j)*e(j));
    }

    protected double gradL_f(int j)
    {
        int dimension = p.getDimension();
        double constraintTerms = gradL_f_constraints(j);
        
        double alpha_g_j = alpha_g.getEntry(j);
        double pt1 = sum(i -> i == slackIndex ? 0 : 
                        alpha_g.getEntry(i)*gp(i)*( f(i)*G(i,j) - e(i)*B(i,j))
                       +alpha_g.getEntry(i)*gq(i)*(-e(i)*G(i,j) - f(i)*B(i,j)),
                       dimension, j);
        double pt2 = alpha_g_j*gp(j)*(sum(n -> f(n)*G(j,n) + e(n)*B(j,n), dimension, j) 
                        +2*f(j)*G(j,j));
        double pt3 = alpha_g_j*gq(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) 
                        -2*f(j)*B(j,j));
        double penaltyTerms = pt1
                            + 
                            (j == slackIndex ? 0 : 
                            (
                                +pt2
                                +pt3
                            ));
        
        // grad_e h(x)(mu + alpha*h(x))
        // grad_e hp(x) = 2*sum{G_ij(e_i - e_j)}
        // grad_e hq(x) = 2*sum{-B_ij(e_i - e_j)}
        double consensusTerms;
        if(!DISABLE_LINE_LOSS)
        {
            double hp_j = hp.getEntry(j);
            double hq_j = hq.getEntry(j);
            double mup_j = mup.getEntry(j);
            double muq_j = muq.getEntry(j);
            double alpha_h_j = alpha_h.getEntry(j);
            double sum_p = sum(n -> G(n,j)*(f(j) - f(n)), dimension, j);
            double sum_q = sum(n -> -B(n,j)*(f(j) - f(n)), dimension, j);
            consensusTerms =
                     2*sum_p*mup_j
                    +2*sum_p*alpha_h_j*hp_j
                    +2*sum_q*muq_j
                    +2*sum_q*alpha_h_j*hq_j;
        }
        else
        {
            consensusTerms = 0;
        }
        
        return 
            constraintTerms
            +penaltyTerms
            +consensusTerms
            ;
    }
    
    protected double gradP_f(int j)
    {
        int dimension = p.getDimension();
        double constraintTerms = gradL_f_constraints(j);
        
        double consensusTerms;
        if(!DISABLE_LINE_LOSS)
        {
            double mup_j = mup.getEntry(j);
            double muq_j = muq.getEntry(j);
            double sum_p = sum(n -> G(n,j)*(f(j) - f(n)), dimension, j);
            double sum_q = sum(n -> -B(n,j)*(f(j) - f(n)), dimension, j);
            consensusTerms =
                     2*sum_p*mup_j
                    +2*sum_q*muq_j;
        }
        else
        {
            consensusTerms = 0;
        }
        
        return constraintTerms + consensusTerms;
    }

    protected double gradL_f_constraints(int j)
    {
        int dimension = p.getDimension();
        return sum(i -> lp(i)*(f(i)*G(i,j) - e(i)*B(i,j)), dimension, j) // i != j
                    +sum(i -> lq(i)*(-e(i)*G(i,j) - f(i)*B(i,j)), dimension, j) // i != j
                    +lp(j)*(sum(n -> f(n)*G(j,n) + e(n)*B(j,n), dimension, j) // n != j
                            +2*G(j,j)*f(j))
                    +lq(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
                            -2*B(j,j)*f(j));
    }
    
    protected double lagrangian()
    {
        return lagrangian(p, q, e, f, hp, hq, mup, muq, lp, lq);
    }
    
    protected double lagrangian(RealVector p, RealVector q, RealVector e, RealVector f, RealVector hp, RealVector hq, RealVector mup, RealVector muq, RealVector lp, RealVector lq)
    {
		double gAug = lagrangian_gAug(p, q, e, f);
        double lg = lagrangian_lg(p, q, e, f, lp, lq);
        double muh = lagrangian_muh(hp, hq, mup, muq);
        double hAug = lagrangian_hAug(hp, hq);
        return cost(p, q) + lg + gAug 
                + muh
                + hAug;
    }

	public double lagrangian_hAug(RealVector hp, RealVector hq)
	{
    	if(DISABLE_CONSENSUS_UPDATE)
    		return 0;
    	double _hp = hp(p, q, e, f);
    	double _hq = hq(p, q, e, f);
		return  0.5*alpha_h.getEntry(0)*_hp*_hp 
		       +0.5*alpha_h.getEntry(0)*_hq*_hq;
	}

	protected double hp(RealVector p, RealVector q, RealVector e, RealVector f)
    {
        int dimension = p.getDimension();
        return sum(i -> (i == slackIndex ? 0 : p.getEntry(i)) + 0.5*lineLoss_p(p, q, e, f, i), dimension);
    }

	protected double hq(RealVector p, RealVector q, RealVector e, RealVector f)
    {
        int dimension = p.getDimension();
        return sum(i -> (i == slackIndex ? 0 : q.getEntry(i)) + 0.5*lineLoss_q(p, q, e, f, i), dimension);
    }

    protected double lineLoss_p(RealVector p, RealVector q, RealVector e, RealVector f, int i)
    {
        int dimension = p.getDimension();
        return sum(j -> lineLoss_p(e.getEntry(i), f.getEntry(i), e.getEntry(j), f.getEntry(j), i, j), dimension);
    }

    protected double lineLoss_q(RealVector p, RealVector q, RealVector e, RealVector f, int i)
    {
        int dimension = p.getDimension();
        return sum(j -> lineLoss_q(e.getEntry(i), f.getEntry(i), e.getEntry(j), f.getEntry(j), i, j), dimension);
    }

    private double lineLoss_p(double ei, double fi, double ej, double fj, int i, int j)
    {
        double de = ei-ej;
        double df = fi-fj;
        return G(i,j)*(de*de + df*df);
    }

    private double lineLoss_q(double ei, double fi, double ej, double fj, int i, int j)
    {
        double de = ei - ej;
        double df = fi - fj;
        return round(-B(i, j)*(de*de + df*df));
    }

	public double lagrangian_muh(RealVector hp, RealVector hq, RealVector mup, RealVector muq)
	{
    	if(DISABLE_CONSENSUS_UPDATE)
    		return 0;
		return mup.dotProduct(hp) + muq.dotProduct(hq);
	}

	public double lagrangian_lg(RealVector p, RealVector q, RealVector e, RealVector f, RealVector lp, RealVector lq)
	{
		RealVector gp = gp(p, q, e, f);
        RealVector gq = gq(p, q, e, f);
        gp.setEntry(slackIndex, 0);
        gq.setEntry(slackIndex, 0);
		return lp.dotProduct(gp) + lq.dotProduct(gq);
	}

	public double lagrangian_gAug(RealVector p, RealVector q, RealVector e, RealVector f)
	{
		return lagrangian_gAug_p(p, q, e, f) + lagrangian_gAug_q(p, q, e, f);
	}

    @SuppressWarnings("deprecation")
	public double lagrangian_gAug_q(RealVector p, RealVector q, RealVector e, RealVector f)
	{
        RealVector gq = gq(p, q, e, f);
        gq.setEntry(slackIndex, 0);
		return 0.5*alpha_g.ebeMultiply(gq).dotProduct(gq);
	}

    @SuppressWarnings("deprecation")
	public double lagrangian_gAug_p(RealVector p, RealVector q, RealVector e, RealVector f)
	{
		RealVector gp = gp(p, q, e, f);
        gp.setEntry(slackIndex, 0);
		return 0.5*alpha_g.ebeMultiply(gp).dotProduct(gp);
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
        if(DISABLE_POWER_FLOW_UPDATE)
            return 0;
        int dimension = p.getDimension();
        double p_i = p.getEntry(i);
        double e_i = e.getEntry(i);
        double f_i = f.getEntry(i);
        double d = sum(n -> e_i*e.getEntry(n)*G(i,n) + f_i*f.getEntry(n)*G(i,n) + f_i*e.getEntry(n)*B(i,n) - e_i*f.getEntry(n)*B(i,n), dimension);
        return round(d - p_i, MIN_G_VALUE);
    }
    
    protected RealVector gq(RealVector p, RealVector q, RealVector e, RealVector f)
    {
        int dimension = q.getDimension();
        return vector(dimension, i -> i == slackIndex ? 0 : gq(p, q, e, f, i));
    }
    
    protected double gq(RealVector p, RealVector q, RealVector e, RealVector f, int i)
    {
        if(DISABLE_POWER_FLOW_UPDATE)
            return 0;
        int dimension = p.getDimension();
        double q_i = q.getEntry(i);
        double e_i = e.getEntry(i);
        double f_i = f.getEntry(i);
        double d = sum(j -> f_i*e.getEntry(j)*G(i,j) - e_i*f.getEntry(j)*G(i,j) - e_i*e.getEntry(j)*B(i,j) - f_i*f.getEntry(j)*B(i,j), dimension);
        return round(d - q_i, MIN_G_VALUE);
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
        if(DISABLE_COST)
            return 0;
        double d = p_i - p_max_i;
        return round(COST_MULTIPLIER*(0.5*d*d + 0.5*q_i*q_i));
    }
    
    public void logConfig()
    {
        if(!DEBUG)
            return;
        
        // Switches:
        out.print("Switches:,");
        if(PROJECT_X                 ) out.print("PROJECT_X                 ,");
        if(FORMAT_NUMBERS            ) out.print("FORMAT_NUMBERS            ,");
        if(VALIDATE_BACKTRACKING     ) out.print("VALIDATE_BACKTRACKING     ,");
        if(VALIDATE_LAGRANGE_DECREASE) out.print("VALIDATE_LAGRANGE_DECREASE,");
        if(VALIDATE_CONSENSUS        ) out.print("VALIDATE_CONSENSUS        ,");
        if(START_WITH_TRUE_VOTLAGES  ) out.print("START_WITH_TRUE_VOTLAGES  ,");
        if(START_WITH_OPTIMAL_POWERS ) out.print("START_WITH_OPTIMAL_POWERS ,");
        if(USE_TRUE_G                ) out.print("USE_TRUE_G                ,");
        if(USE_TRUE_H                ) out.print("USE_TRUE_H                ,");
        if(DISABLE_VOTLAGES_UPDATES  ) out.print("DISABLE_VOTLAGES_UPDATES  ,");
        if(DISABLE_LAMBDA_UPDATE     ) out.print("DISABLE_LAMBDA_UPDATE     ,");
        if(DISABLE_POWER_FLOW_UPDATE ) out.print("DISABLE_POWER_FLOW_UPDATE ,");
        if(DISABLE_COST              ) out.print("DISABLE_COST              ,");
        if(DISABLE_CONSENSUS_UPDATE  ) out.print("DISABLE_CONSENSUS_UPDATE  ,");
        if(DISABLE_MU_UPDATE         ) out.print("DISABLE_MU_UPDATE         ,");
        if(DISABLE_LINE_LOSS         ) out.print("DISABLE_LINE_LOSS         ,");
        if(USE_SIGMOID_ALPHA_UPDATE  ) out.print("USE_SIGMOID_ALPHA_UPDATE  ,");
        if(ROUND_TO_ZERO             ) out.print("ROUND_TO_ZERO             ,");
        if(KEEP_G_AND_H_CLOSE        ) out.print("KEEP_G_AND_H_CLOSE        ,");
        if(APPROXIMATE_CONSTRAINTS    ) out.print("APPROXIMATE_CONSTRAINTS  ,");
        out.println();
        
        // Constant values:
        out.print("BASE_VOLTAGE               ,");
        out.print("SLACK_VOLTAGE              ,");
        out.print("BASE_POWER                 ,");
        out.print("V_MAX                      ,");
        out.print("V_MIN                      ,");
        out.print("V_ARG_MIN                  ,");
        out.print("V_ARG_MAX                  ,");
        out.print("INITIAL_G_AUG_SCALE        ,");
        out.print("G_AUG_SCALE_STEP           ,");
        out.print("G_MAX_AUG_SCALE            ,");
        out.print("INITIAL_H_AUG_SCALE        ,");
        out.print("H_AUG_SCALE_STEP           ,");
        out.print("H_MAX_AUG_SCALE            ,");
        out.print("COST_MULTIPLIER            ,");
        out.print("X_STEPS_PER_ITERATION      ,");
        out.print("MAX_X_STEPS                ,");
        out.print("K                          ,");
        out.print("DEBUG_RATE                 ,");
        out.print("AGENT_SELECTION_PROBABILITY,");
        out.print("INITIAL_STEP_SIZE          ,");
        out.print("MIN_STEP_SIZE              ,");
        out.print("XI                         ,");
        out.print("ETA_G                      ,");
        out.print("ETA_H                      ,");
        out.print("MOMENTUM_RATE              ,");
        out.print("EPSILON_BASE_PQ            ,");
        out.print("EPSILON_BASE_EF            ,");
        out.print("EPSILON_TARGET             ,");
        out.print("CONSTRAINT_BASE            ,");
        out.print("CONSTRAINT_TARGET          ,");
        out.println();
        
        out.print(BASE_VOLTAGE               +",");
        out.print(SLACK_VOLTAGE              +",");
        out.print(BASE_POWER                 +",");
        out.print(V_MAX                      +",");
        out.print(V_MIN                      +",");
        out.print(V_ARG_MIN                  +",");
        out.print(V_ARG_MAX                  +",");
        out.print(INITIAL_G_AUG_SCALE        +",");
        out.print(G_AUG_SCALE_STEP           +",");
        out.print(G_MAX_AUG_SCALE            +",");
        out.print(INITIAL_H_AUG_SCALE        +",");
        out.print(H_AUG_SCALE_STEP           +",");
        out.print(H_MAX_AUG_SCALE            +",");
        out.print(COST_MULTIPLIER            +",");
        out.print(X_STEPS_PER_ITERATION      +",");
        out.print(MAX_X_STEPS                +",");
        out.print(K                          +",");
        out.print(DEBUG_RATE                 +",");
        out.print(AGENT_SELECTION_PROBABILITY+",");
        out.print(INITIAL_STEP_SIZE          +",");
        out.print(MIN_STEP_SIZE              +",");
        out.print(XI                         +",");
        out.print(ETA_G                      +",");
        out.print(ETA_H                      +",");
        out.print(MOMENTUM_RATE              +",");
        out.print(EPSILON_BASE_PQ            +",");
        out.print(EPSILON_BASE_EF            +",");
        out.print(EPSILON_TARGET             +",");
        out.print(CONSTRAINT_BASE            +",");
        out.print(CONSTRAINT_TARGET          +",");
        out.println();
    }
    
    public void logBusses(AnalysisResults results)
    {
        if(!DEBUG)
            return;
        
        out.println("index,bus,e,f,p,q");
        Map<String, Integer> numbers = results.getBusNumbers();
        int i = 0;
        for (String name : numbers.keySet())
        {
        	if(i >= DISPLAY_DIMENSION)
        		break;
        	
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
            
            ++i;
        }
    }
    
    public void logApproximateGradients()
    {
        double eMin = 0.999;
        double eRange = 0.002;
        double fMin = -0.001;
        double fRange = 0.002;
        
        RealVector _p = new ArrayRealVector(p);
        RealVector _q = new ArrayRealVector(q);
        RealVector _e = new ArrayRealVector(e);
        RealVector _f = new ArrayRealVector(f);
        RealVector _hp = new ArrayRealVector(hp);
        RealVector _hq = new ArrayRealVector(hq);
        
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
                double L_e = lagrangian(_p, _q, _e, _f, _hp, _hq, mup, muq, lp, lq);
                _e.setEntry(agent.index, e.getEntry(agent.index));

                _f.setEntry(agent.index, fMin+x*fRange);
                double L_f = lagrangian(_p, _q, _e, _f, _hp, _hq, mup, muq, lp, lq);
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
        RealVector _hp = new ArrayRealVector(hp);
        RealVector _hq = new ArrayRealVector(hq);
        
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
                double L_e = lagrangian(_p, _q, _e, _f, _hp, _hq, mup, muq, lp, lq);
                System.out.print(L_e+",");
            }
            System.out.println();
        }
    }
}