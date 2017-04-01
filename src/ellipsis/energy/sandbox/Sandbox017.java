package ellipsis.energy.sandbox;

import static ellipsis.energy.sandbox.GridConstants.BASE_IMPEDANCE;
import static ellipsis.energy.sandbox.GridConstants.BASE_POWER;
import static ellipsis.energy.sandbox.GridConstants.BASE_VOLTAGE;
import static ellipsis.energy.sandbox.GridConstants.BUS_2;
import static ellipsis.energy.sandbox.GridConstants.BUS_3;
import static ellipsis.energy.sandbox.GridConstants.BUS_4;
import static ellipsis.energy.sandbox.GridConstants.BUS_5;
import static ellipsis.energy.sandbox.GridConstants.DG_3;
import static ellipsis.energy.sandbox.GridConstants.DG_4;
import static ellipsis.energy.sandbox.GridConstants.LINE_1_2;
import static ellipsis.energy.sandbox.GridConstants.LINE_1_3;
import static ellipsis.energy.sandbox.GridConstants.LINE_2_3;
import static ellipsis.energy.sandbox.GridConstants.LINE_2_4;
import static ellipsis.energy.sandbox.GridConstants.LINE_3_5;
import static ellipsis.energy.sandbox.GridConstants.LOAD_2;
import static ellipsis.energy.sandbox.GridConstants.LOAD_5;
import static ellipsis.energy.sandbox.GridConstants.SLACK_BUS;
import static ellipsis.energy.sandbox.GridConstants.SLACK_SOURCE;
import static ellipsis.util.Sum.sum;
import static ellipsis.util.VectorHelper.vector;
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Collection;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

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
import ellipsis.util.TeeOutputStream;

/**
 * Test Implicit Recurrent Neural Network solver for the same problem as
 * {@link Sandbox016}.
 * @author bmillar
 *
 */
public class Sandbox017
{
    public static void main(String[] args)
    {
        new Sandbox017().run(); 
    }
    
    public static PrintStream out;
    static
    {
        try
        {
            out = new PrintStream(new TeeOutputStream(new FileOutputStream("/tmp/Sandbox017.csv"), System.out));
        } 
        catch (FileNotFoundException e)
        {
            throw new RuntimeException(e);
        }
    }
    
    protected static int K = 1000;
    protected static int M = 1000*K;
    protected static boolean FORMAT_NUMBERS = true;
    protected static boolean START_WITH_TRUE_VOTLAGES = false;
    protected static boolean _35BUS = true;
    protected static final int T = 20*M;
    protected static final int DEBUG_RATE = T/1200;
    protected static final double COST_MULTIPLIER = 10;
    protected static final double MAX_LAMBDA = 1e14;
    
    private static final double BASE_VOLTAGE_35 = 120;
    private static final double BASE_POWER_35 = 10e3;
    private static final double SLACK_VOLTAGE = 1.0;
    
    private static final double V_MAX = 1.05;
    private static final double V_MIN = 0.95;
    private static final double V_ARG_MIN = -Math.PI/4;
    private static final double V_ARG_MAX = Math.PI/4;

    protected Grid grid;
    protected RealVector p, q, e, f, lp, lq;
    protected double alpha = 1e-1;
    protected double beta = 1;//4; // aug step multiplier per second
    protected double timeStep = 100/T;
    protected double timeStepRate = 1.0 - 1e-6;
    protected double minTimeStep = 10.0/T;
    protected double maxTime = 24.0;
    protected int slackIndex;
    protected RealMatrix G, B;
    protected RealVector generatorMask;
    protected RealVector p_max;
    protected RealVector p_costMin;

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
    
    
    //// Simulation/Solution ////
    
    @SuppressWarnings("deprecation")
    private void loop()
    {
        double time = 0;
        //for(int t = 0; t < T; ++t)
        int t = 0;
        while(time < maxTime)
        {
            // Log:
            if(t%DEBUG_RATE == 0)
                debug(time);

            // Calculate ODE:
            RealVector gp = zeroIfSmall(gp());
            RealVector gq = zeroIfSmall(gq());
            
            RealVector lpDot = zeroIfSmall(gp.mapMultiply(alpha));
            RealVector lqDot = zeroIfSmall(gq.mapMultiply(alpha));
            RealVector pDot =  zeroIfSmall(grad_pC().add(lp.add(lpDot).mapMultiply(-1)).mapMultiply(-1).ebeMultiply(generatorMask));
            RealVector qDot =  zeroIfSmall(grad_qC().add(lq.add(lqDot).mapMultiply(-1)).mapMultiply(-1).ebeMultiply(generatorMask));
            RealVector l = lp.append(lq);
            RealVector g = gp.append(gq);
            RealVector lPlusAlphaG = l.add(g.mapMultiply(alpha));
            RealVector eDot = zeroIfSmall(grad_eG().operate(lPlusAlphaG).mapMultiply(-1));
            RealVector fDot = zeroIfSmall(grad_fG().operate(lPlusAlphaG).mapMultiply(-1));
            double alphaDot = (beta-1)*alpha;
            
            // Integrate:
            p = p.add(pDot.mapMultiply(timeStep));
            q = q.add(qDot.mapMultiply(timeStep));
            e = e.add(eDot.mapMultiply(timeStep));
            f = f.add(fDot.mapMultiply(timeStep));
            lp = lp.add(lpDot.mapMultiply(timeStep));
            lq = lq.add(lqDot.mapMultiply(timeStep));
            alpha += alphaDot*timeStep;
            
            if(e.isNaN())
                Util.breakpoint();
            
            // Constrain:
            projectPowers();
            projectVoltages();
            project(lp, -MAX_LAMBDA, MAX_LAMBDA);
            project(lq, -MAX_LAMBDA, MAX_LAMBDA);
            
            if(e.isNaN())
                Util.breakpoint();
            
            time += timeStep;
            timeStep = Math.max(timeStepRate*timeStep, minTimeStep);
            ++t;
        }
        
        debug(time);
    }

    protected RealVector zeroIfSmall(RealVector v)
    {
        if(v.getNorm() < 1e-6)
            v = vector(v.getDimension(), 0.0);
        return v;
    }

    private void projectPowers()
    {
        int dimension = p.getDimension();
        for(int i = 0; i < dimension; ++i)
        {
            if(generatorMask.getEntry(i) != 1)
                continue;
            
            double p_i = p.getEntry(i);
            double q_i = q.getEntry(i);
            Complex s_i = projectPower(p_i, q_i, p_max.getEntry(i));
            p.setEntry(i, s_i.getReal());
            q.setEntry(i, s_i.getImaginary());
        }
    }

    private void projectVoltages()
    {
        int dimension = e.getDimension();
        for(int i = 0; i < dimension; ++i)
        {
            if(i == slackIndex)
                continue;
            
            double e_i = e.getEntry(i);
            double f_i = f.getEntry(i);
            Complex v_i = projectVoltage(e_i, f_i);
            e.setEntry(i, v_i.getReal());
            f.setEntry(i, v_i.getImaginary());
        }
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

    protected Complex projectPower(double p, double q, double p_max)
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
    
    protected void project(RealVector v, double low, double high)
    {
        int dimension = v.getDimension();
        for (int i = 0; i < dimension; i++)
        {
            double v_i = v.getEntry(i);
            v_i = project(v_i, low, high);
            v.setEntry(i, v_i);
        }
    }

    protected double project(double d, double low, double high)
    {
        return min(max(low, d), high);
    }
    
    public RealMatrix grad_eG()
    {
        RealMatrix grad_eGp = gradG(e, f, this::gradGpi_dej);
        RealMatrix grad_eGq = gradG(e, f, this::gradGqi_dej);
        
        int dimension = e.getDimension();
        RealMatrix grad = new Array2DRowRealMatrix(dimension, 2*dimension);
        
        // [ grad_eGp  grad_eGq ]
        for(int i = 0; i < dimension; ++i)
        {
            for(int j = 0; j < dimension; ++j)
            {
                grad.setEntry(i, j, grad_eGp.getEntry(i, j));
                grad.setEntry(i, j+dimension, grad_eGq.getEntry(i, j));
            }
        }
        
        return grad;
    }
    
    public RealMatrix grad_fG()
    {
        RealMatrix grad_fGp = gradG(e, f, this::gradGpi_dfj);
        RealMatrix grad_fGq = gradG(e, f, this::gradGqi_dfj);
        
        int dimension = e.getDimension();
        RealMatrix grad = new Array2DRowRealMatrix(dimension, 2*dimension);
        
        // [ grad_fGp  grad_fGq ]
        for(int i = 0; i < dimension; ++i)
        {
            for(int j = 0; j < dimension; ++j)
            {
                grad.setEntry(i, j, grad_fGp.getEntry(i, j));
                grad.setEntry(i, j+dimension, grad_fGq.getEntry(i, j));
            }
        }
        
        return grad;
    }
    
    protected static interface Grad
    {
        /**
         * Gradient of constraint i W.R.T. index j.
         */
        double value(RealVector e, RealVector f, int i, int j);
    }
    protected RealMatrix gradG(RealVector e, RealVector f, Grad grad)
    {
        int dimension = e.getDimension();
        RealMatrix jacobian = new Array2DRowRealMatrix(dimension, dimension);
        for(int i = 0; i < jacobian.getRowDimension(); ++i)
        {
            if(i == slackIndex)
                continue;
            
            for(int j = 0; j < jacobian.getColumnDimension(); ++j)
            {
                if(j == slackIndex)
                    continue;
                
                jacobian.setEntry(i, j, grad.value(e, f, j, i));
            }
        }
        return jacobian;
    }
    
    private double gradGpi_dej(RealVector e, RealVector f, int i, int j)
    {
        double e_i = e.getEntry(i);
        double f_i = f.getEntry(i);
        
        if(i != j)
            return e_i*G(i, j) + f_i*B(i, j);
        
        int dimension = p.getDimension();
        return sum(n -> e.getEntry(n)*G(i, n) - f.getEntry(n)*B(i, n), dimension, i) + 2*G(i, i)*e_i;
    }
    
    private double gradGpi_dfj(RealVector e, RealVector f, int i, int j)
    {
        double e_i = e.getEntry(i);
        double f_i = f.getEntry(i);
        
        if(i != j)
            return f_i*G(i, j) - e_i*B(i, j);
        
        int dimension = p.getDimension();
        return sum(n -> f.getEntry(n)*G(i, n) + e.getEntry(n)*B(i, n), dimension, i) + 2*G(i, i)*f_i;
    }
    
    private double gradGqi_dej(RealVector e, RealVector f, int i, int j)
    {
        double e_i = e.getEntry(i);
        double f_i = f.getEntry(i);
        
        if(i != j)
            return f_i*G(i, j) - e_i*B(i, j);
        
        int dimension = p.getDimension();
        return sum(n -> -f.getEntry(n)*G(i, n) - e.getEntry(n)*B(i, n), dimension, i) - 2*B(i, i)*e_i;
    }
    
    private double gradGqi_dfj(RealVector e, RealVector f, int i, int j)
    {
        double e_i = e.getEntry(i);
        double f_i = f.getEntry(i);
        
        if(i != j)
            return -e_i*G(i, j) - f_i*B(i, j);
        
        int dimension = p.getDimension();
        return sum(n -> e.getEntry(n)*G(i, n) - f.getEntry(n)*B(i, n), dimension, i) - 2*B(i, i)*f_i;
    }
    
    public RealVector grad_pC()
    {
        int length = p.getDimension();
        RealVector grad = new ArrayRealVector(length, 0.0);
        for(int i = 0; i < length; ++i)
        {
            if(generatorMask.getEntry(i) == 1)
            {
                double gradC_p_i = gradC_p(p, q, i);
                grad.setEntry(i, gradC_p_i);
            }
        }
        return grad;
    }
    
    protected double gradC_p(RealVector p, RealVector q, int i)
    {
        return COST_MULTIPLIER*(p.getEntry(i) - p_costMin.getEntry(i));
    }
    
    public RealVector grad_qC()
    {
        int length = q.getDimension();
        RealVector grad = new ArrayRealVector(length, 0.0);
        for(int i = 0; i < length; ++i)
        {
            if(generatorMask.getEntry(i) == 1)
            {
                double gradC_q_i = gradC_q(p, q, i);
                grad.setEntry(i, gradC_q_i);
            }
        }
        return grad;
    }
    
    protected double gradC_q(RealVector p, RealVector q, int i)
    {
        return COST_MULTIPLIER*q.getEntry(i);
    }
    
    public RealVector gp()
    {
        return gp(p, q, e, f);
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
    
    public RealVector gq()
    {
        return gq(p, q, e, f);
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
    
    
    //// Accessors ////

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
    
    
    //// Setup ////

    protected AnalysisResults init()
    {
        initGrid();
        
        // Analyse grid to get parameters:
        LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        lfa.setIterations(40);
        AnalysisResults results = lfa.analyse();
        if(!results.getDidConverge())
            throw new RuntimeException("Initial analysis did not converge.");
        
        logBusses(results);
        
        initVariables(results);
        
        debugHeader(results);
        
        return results;
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
    
    protected void initGrid()
    {
        if(_35BUS)
            initGrid_35Bus();
        else
            initGrid_5Bus();
    }

    protected void initGrid_5Bus()
    {
        // Copied from Sandbox008:
        grid = 
              Grid.grid().
              
                  Bus(SLACK_BUS).
                  
                      SlackSource(SLACK_SOURCE, 1.00*BASE_VOLTAGE, 0, 0, 0).
                      
                      Line(LINE_1_2, 1, 0.02*BASE_IMPEDANCE, 0.04*BASE_IMPEDANCE).
                          Bus(BUS_2).
                              Load(LOAD_2).
                              Line(LINE_2_3, BUS_3, 1, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
                              Line(LINE_2_4, 1.5, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
                                Bus(BUS_4).
                                    DistributedSource(DG_4).
                                terminate().
                            terminate().
                          terminate().
                      terminate().
                      
                      Line(LINE_1_3, 1, 0.01*BASE_IMPEDANCE, 0.03*BASE_IMPEDANCE).
                          Bus(BUS_3).
                            DistributedSource(DG_3).
                            Line(LINE_3_5, 1.5, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
                                Bus(BUS_5).
                                    Load(LOAD_5).
                                terminate().
                            terminate().
                          terminate().
                      terminate().
                      
                  terminate().
                  
              grid();
        
        grid.setBaseVoltage(BASE_VOLTAGE);
        grid.setBasePower(BASE_POWER);
        
        grid.getLoad(LOAD_2).setLoad(new Complex(150e6, 50e6));
        grid.getLoad(LOAD_5).setLoad(new Complex(80e6, 25e6));
      
        grid.getSource(DG_3).setPmax(650e6);
        grid.getSource(DG_3).setPowerOutput(100e6, 0);
        grid.getSource(DG_4).setPmax(650e6);
        grid.getSource(DG_4).setPowerOutput(100e6, 0);
  }
    
    protected void initGrid_35Bus()
    {
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox017-lineConfig");
        InputStream lineData = getClass().getResourceAsStream("Sandbox017-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox017-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream("Sandbox017-loadData.csv");
        grid = GridGenerator.loadGrid("IEEE123", lineConfig, lineData, switchData, loadData);
        
        grid.setBaseVoltage(BASE_VOLTAGE_35);
        grid.setBasePower(BASE_POWER_35);
        
        // Add slack bus:
        SlackSource slack = new SlackSource();
        slack.setVoltage(new Complex(grid.getBaseVoltage()*SLACK_VOLTAGE));
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
    
    
    //// Debug and Logging ////

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

    protected void debugHeader(AnalysisResults results)
    {
        Map<String, Integer> numbers = results.getBusNumbers();
        out.print("t,");
        
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
//        out.print("C(pq),L(pql),");
//        debugLogNames(numbers, dimension, "t");
        out.print("alpha,");
//        debugLogNames(numbers, dimension, "inf");
//        out.print("||grad_p||,||grad_q||,||grad_e||,||grad_f||,CPU Time,");
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
    
    protected void debug(double t)
    {
        out.print(t);
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
            out.print(format(gp(p, q, e, f, i)));
            out.print(',');
        }
        for(int i = 0; i < dimension; ++i)
        {
            out.print(format(gq(p, q, e, f, i)));
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
        
        // Costs:
//        out.print(format(cost()));
//        out.print(',');
//        out.print(format(L()));
//        out.print(',');

        // Step size:
//        for(int i = 0; i < dimension; ++i)
//        {
//            out.print(t.getEntry(i));
//            out.print(',');
//        }

        // Augmentation scale:
        out.print(alpha);
        out.print(',');
        
//        out.print(format(gradLNorm(this::gradL_p, dimension)));
//        out.print(',');
//        out.print(format(gradLNorm(this::gradL_q, dimension)));
//        out.print(',');
//        out.print(format(gradLNorm(this::gradL_e, dimension)));
//        out.print(',');
//        out.print(format(gradLNorm(this::gradL_f, dimension)));
//        out.print(',');
//        out.print(format(cpuTime));
//        out.print(',');
        
        out.println();
    }

    protected double gradL_p(int i)
    {
        if(generatorMask.getEntry(i) == 1)
            return  gradC_p(i) - lp(i) - alpha*gp(i);
        else
            return 0;
    }
    
    public double gp(int i)
    {
        return gp(p, q, e, f, i);
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

    protected double gradL_q(int i)
    {
        if(generatorMask.getEntry(i) == 1)
            return gradC_q(i) - lq(i) - alpha*gq(i);
        else
            return 0;
    }
    
    public double gq(int i)
    {
        return gq(p, q, e, f, i);
    }

    protected double gradC_q(int i)
    {
        return gradC_q(p, q, i);
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
        
        double pt1 = alpha*sum(i -> i == slackIndex ? 0 : 
                             gp(i)*(e(i)*G(i,j) + f(i)*B(i,j))
                            +gq(i)*(f(i)*G(i,j) - e(i)*B(i,j)),
                            dimension, j);
        double pt2 = alpha*gp(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) 
                  +2*e(j)*G(j,j));
        double pt3 = alpha*gq(j)*(sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) 
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
        
        double pt1 = alpha*sum(i -> i == slackIndex ? 0 : 
                        gp(i)*( f(i)*G(i,j) - e(i)*B(i,j))
                       +gq(i)*(-e(i)*G(i,j) - f(i)*B(i,j)),
                       dimension, j);
        double pt2 = alpha*gp(j)*(sum(n -> f(n)*G(j,n) + e(n)*B(j,n), dimension, j) 
                        +2*f(j)*G(j,j));
        double pt3 = alpha*gq(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) 
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
                dg.setPowerOutput(grid.getBasePower()*p(i), grid.getBasePower()*q(i));
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
        lfa.setTargetError(1e-9);
        AnalysisResults results2 = lfa.analyse();
        if(!results2.getDidConverge())
            throw new RuntimeException("Results did not converge.");
        
        out.println("\nAnalysis Results:");
//        GridDisplay.showInFrame(grid, grid.getBasePower(), grid.getBaseVoltage()).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        logBusses(results2);
    }

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
}