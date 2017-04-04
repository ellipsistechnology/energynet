package ellipsis.energy.sandbox;

import static ellipsis.energy.sandbox.GridConstants.SLACK_SOURCE;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Random;

import javax.swing.JFrame;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;

import com.mls.util.ThreadPool;

import ellipsis.energy.grid.GridDisplay;
import ellipsis.energy.grid.GridGenerator;
import ellipsis.energy.grid.SlackSource;
import ellipsis.util.EmptyOutputStream;
import ellipsis.util.TeeOutputStream;

/**
 * Alternate configureation to {@link Sandbox018}.
 * @author bmillar
 *
 */
public class Sandbox018B extends Sandbox018
{
    static boolean tune = false;
    static boolean analyseLambdaStar = false;
    
    public Sandbox018B()
    {
        super();
        
        if(tune)
        {
            // No output from test cases:
            out = new PrintStream(new EmptyOutputStream());
        }
        else
        {
            try
            {
                out = new PrintStream(new TeeOutputStream(new FileOutputStream("/tmp/Sandbox018B.csv"), System.out));
            }
            catch (FileNotFoundException e)
            {
                throw new RuntimeException(e);
            }
        }
              
        
        //// Override configuration ////
        
                                PROJECT_X = true;
                           FORMAT_NUMBERS = false;
              
                    VALIDATE_BACKTRACKING = false;
               VALIDATE_LAGRANGE_DECREASE = false;
              
                 START_WITH_TRUE_VOTLAGES = false;
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
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 8;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 3;
                                    ETA_H = 0.005;
                                   
                                       XI = 0.2;
                               
                                        K = 3000;
                               DEBUG_RATE = K / 1000;
              AGENT_SELECTION_PROBABILITY = 1.0;
             
                        INITIAL_STEP_SIZE = 1;
                            MIN_STEP_SIZE = 1e-50;
                              MAX_X_STEPS = 1;
                              
                          EPSILON_BASE_PQ = 100;//0.1;
                          EPSILON_BASE_EF = 1000;//100;
                           EPSILON_TARGET = 1;//10e-2;
    }
    
    public static void main(String[] args)
    {
        if(tune)
        {
            // Tune parameters and then run again to get results:
            System.out.println("Tuning...");
            tune();
            System.out.println("Tuning complete.\nRunning...");
            new Sandbox018B().run();
            System.out.println("Running complete.");
        }
        else if(analyseLambdaStar)
        {
            Sandbox018B sandbox018b = new Sandbox018B();
            sandbox018b.start();
            sandbox018b.init();
            sandbox018b.analyseLambdaStar();
            sandbox018b.finish();
        }
        else
        {
            new Sandbox018B().run();
        }
    }
    
    /** 
     * WARNING: incomplete implementation.
     */
    protected void analyseLambdaStar()
    {
        // Randomise powers and voltages:
        Random rand = new Random(0);
        int dimension = p.getDimension();
        for (int i = 0; i < dimension; i++)
        {
            // P & Q:
            Complex s_i = new Complex(p.getEntry(i), q.getEntry(i)); // start with p&q to preserve angle
            s_i = s_i.multiply(rand.nextDouble()*p_max.getEntry(i)/s_i.abs());
            p.setEntry(i, s_i.getReal());
            q.setEntry(i, s_i.getImaginary());
            
            // E & F:
            Complex v_i = ComplexUtils.polar2Complex(
                rand.nextDouble()*(V_MAX-V_MIN)+V_MIN,
                rand.nextDouble()*(V_ARG_MAX-V_ARG_MIN)+V_ARG_MIN
            );
            
            e.setEntry(i, v_i.getReal());
            f.setEntry(i, v_i.getImaginary());
        }
        
        // Analyse lambda*:
        try
        {
            RealVector lambdaStar = lambdaStar();
            for (int i = 0; i < lambdaStar.getDimension(); i++)
            {
                out.print(lambdaStar.getEntry(i));
                out.print(",");
            }
            out.println();
        }
        catch(SingularMatrixException e){}
    }

    static class TuneResult
    {
        double gNorm, hNorm, gradNorm;
        Sandbox018B sb18B;
        
        public TuneResult(Sandbox018B sb)
        {
            sb18B = sb;
        }
        
        public TuneResult run()
        {
            try
            {
                sb18B.run();
            }
            catch(MathIllegalArgumentException e)
            {
                gNorm = Double.MAX_VALUE;
                hNorm = Double.MAX_VALUE;
                gradNorm = Double.MAX_VALUE;
                
                return this;
            }
            
            gNorm = sb18B.gNorm();
            hNorm = sb18B.hNorm();
            gradNorm = sb18B.gradNorm();
            
            return this;
        }

        protected double cost()
        {
            return Math.max(1.0, gNorm/targetG) + Math.max(1.0, hNorm/targetH) + Math.max(1.0, gradNorm/targetGrad);
        }
    }
    
    static double targetG = 0.001;
    static double targetH = 0.001;
    static double targetGrad = 0.1;
    
    static double augGStep = 0.001;
    static double augHStep = 0.1;
    static double etaGStep = 0.1;
    static double etaHStep = 0.0001;
    
    static double randomRange = 10;
    
    static int MAX_ITERATIONS = 100;
    static double INITIAL_STEP = 10;
    
    static boolean changed;
    static TuneResult result_augG = new TuneResult(new Sandbox018B());
    static TuneResult result_augH = new TuneResult(new Sandbox018B());
    static TuneResult result_etaG = new TuneResult(new Sandbox018B());
    static TuneResult result_etaH = new TuneResult(new Sandbox018B());
    
    static int SET_AUG_G = 1;
    static int SET_AUG_H = 2;
    static int SET_ETA_G = 3;
    static int SET_ETA_H = 4;
    
    private static synchronized void setChanged(boolean c)
    {
        changed = c; 
    }

    private static void updateAugG(double iNITIAL_G_AUG_SCALE)
    {
        result_augG.sb18B.INITIAL_G_AUG_SCALE = iNITIAL_G_AUG_SCALE;
        result_augH.sb18B.INITIAL_G_AUG_SCALE = iNITIAL_G_AUG_SCALE;
        result_etaG.sb18B.INITIAL_G_AUG_SCALE = iNITIAL_G_AUG_SCALE;
        result_etaH.sb18B.INITIAL_G_AUG_SCALE = iNITIAL_G_AUG_SCALE;
    }

    private static void updateAugH(double iNITIAL_H_AUG_SCALE)
    {
        result_augG.sb18B.INITIAL_H_AUG_SCALE = iNITIAL_H_AUG_SCALE;
        result_augH.sb18B.INITIAL_H_AUG_SCALE = iNITIAL_H_AUG_SCALE;
        result_etaG.sb18B.INITIAL_H_AUG_SCALE = iNITIAL_H_AUG_SCALE;
        result_etaH.sb18B.INITIAL_H_AUG_SCALE = iNITIAL_H_AUG_SCALE;
    }

    private static void updateEtaG(double eTA_G)
    {
        result_augG.sb18B.ETA_G = eTA_G;
        result_augH.sb18B.ETA_G = eTA_G;
        result_etaG.sb18B.ETA_G = eTA_G;
        result_etaH.sb18B.ETA_G = eTA_G;
    }

    private static void updateEtaH(double eTA_H)
    {
        result_augG.sb18B.ETA_H = eTA_H;
        result_augH.sb18B.ETA_H = eTA_H;
        result_etaG.sb18B.ETA_H = eTA_H;
        result_etaH.sb18B.ETA_H = eTA_H;
    }
    
    protected static void tune()
    {
        System.out.println("Iter.,INITIAL_G_AUG_SCALE,INITIAL_H_AUG_SCALE,ETA_G,ETA_H,gNorm,hNorm,gradNorm,cost");

        int i = 0;
        ThreadPool pool = new ThreadPool(4);
//        pool.setDisallowed(true);
        do
        {
            setChanged(false);
            
//            augGStep = 0.1/(i+1);
//            augHStep = 10/(i+1);
//            etaGStep = 10/(i+1);
//            etaHStep = 0.01/(i+1);
            
            // Check for better g aug scales:
            pool.queueTask(new Runnable()
            {
                @Override
                public void run()
                {
                    result_augG.run();
                    double baseCost = result_augG.cost();
                    
                    result_augG.sb18B.INITIAL_G_AUG_SCALE -= augGStep;
                    result_augG.run();
                    double leftCost = result_augG.cost();
                    
                    result_augG.sb18B.INITIAL_G_AUG_SCALE += 2*augGStep;
                    result_augG.run();
                    double rightCost = result_augG.cost();
                    
                    if(leftCost < baseCost)
                    {
                        result_augG.sb18B.INITIAL_G_AUG_SCALE -= 2*augGStep;
                        setChanged(true);
                    }
                    else if(baseCost < rightCost)
                    {
                        result_augG.sb18B.INITIAL_G_AUG_SCALE -= augGStep;
                    }
                    else
                    {
                        setChanged(true);
                    }
                }
            });
            
            // Check for better h aug scales:
            pool.queueTask(new Runnable()
            {
                @Override
                public void run()
                {
                    result_augH.run();
                    double baseCost = result_augH.cost();
                    
                    result_augH.sb18B.INITIAL_H_AUG_SCALE -= augHStep;
                    result_augH.run();
                    double leftCost = result_augH.cost();
                    
                    result_augH.sb18B.INITIAL_H_AUG_SCALE += 2*augHStep;
                    result_augH.run();
                    double rightCost = result_augH.cost();
                    
                    if(leftCost < baseCost)
                    {
                        result_augH.sb18B.INITIAL_H_AUG_SCALE -= 2*augHStep;
                        setChanged(true);
                    }
                    else if(baseCost < rightCost)
                    {
                        result_augH.sb18B.INITIAL_H_AUG_SCALE -= augHStep;
                    }
                    else
                    {
                        setChanged(true);
                    }
                }
            });
            
            // Check for better eta_g:
            pool.queueTask(new Runnable()
            {
                @Override
                public void run()
                {
                    result_etaG.run();
                    double baseCost = result_etaG.cost();
                    
                    result_etaG.sb18B.ETA_G -= etaGStep;
                    result_etaG.run();
                    double leftCost = result_etaG.cost();
                    
                    result_etaG.sb18B.ETA_G += 2*etaGStep;
                    result_etaG.run();
                    double rightCost = result_etaG.cost();
                    
                    if(leftCost < baseCost)
                    {
                        result_etaG.sb18B.ETA_G -= 2*etaGStep;
                        setChanged(true);
                    }
                    else if(baseCost < rightCost)
                    {
                        result_etaG.sb18B.ETA_G -= etaGStep;
                    }
                    else
                    {
                        setChanged(true);
                    }
                }
            });
            
            // Check for better eta_h:
            pool.queueTask(new Runnable()
            {
                @Override
                public void run()
                {
                    result_etaH.run();
                    double baseCost = result_etaH.cost();
                    
                    result_etaH.sb18B.ETA_H -= etaHStep;
                    result_etaH.run();
                    double leftCost = result_etaH.cost();
                    
                    result_etaH.sb18B.ETA_H += 2*etaHStep;
                    result_etaH.run();
                    double rightCost = result_etaH.cost();
                    
                    
                    if(leftCost < baseCost)
                    {
                        result_etaH.sb18B.ETA_H -= 2*etaHStep;
                        setChanged(true);
                    }
                    else if(baseCost < rightCost)
                    {
                        result_etaH.sb18B.ETA_H -= etaHStep;
                    }
                    else
                    {
                        setChanged(true);
                    }
                }
            });
            
            pool.waitForAll();
            
            System.out.println(
                    i+","
                    +result_augG.sb18B.INITIAL_G_AUG_SCALE+","+result_augG.sb18B.INITIAL_H_AUG_SCALE+","+result_augG.sb18B.ETA_G+","+result_augG.sb18B.ETA_H+","
                    +result_augG.gNorm+","+result_augG.hNorm+","+result_augG.gradNorm+","+result_augG.cost());
            System.out.println(
                    i+","
                    +result_augH.sb18B.INITIAL_G_AUG_SCALE+","+result_augH.sb18B.INITIAL_H_AUG_SCALE+","+result_augH.sb18B.ETA_G+","+result_augH.sb18B.ETA_H+","
                    +result_augH.gNorm+","+result_augH.hNorm+","+result_augH.gradNorm+","+result_augH.cost());
            System.out.println(
                    i+","
                    +result_etaG.sb18B.INITIAL_G_AUG_SCALE+","+result_etaG.sb18B.INITIAL_H_AUG_SCALE+","+result_etaG.sb18B.ETA_G+","+result_etaG.sb18B.ETA_H+","
                    +result_etaG.gNorm+","+result_etaG.hNorm+","+result_etaG.gradNorm+","+result_etaG.cost());
            System.out.println(
                    i+","
                    +result_etaH.sb18B.INITIAL_G_AUG_SCALE+","+result_etaH.sb18B.INITIAL_H_AUG_SCALE+","+result_etaH.sb18B.ETA_G+","+result_etaH.sb18B.ETA_H+","
                    +result_etaH.gNorm+","+result_etaH.hNorm+","+result_etaH.gradNorm+","+result_etaH.cost());
            
            updateAugG(result_augG.sb18B.INITIAL_G_AUG_SCALE);
            updateAugH(result_augH.sb18B.INITIAL_H_AUG_SCALE);
            updateEtaG(result_etaG.sb18B.ETA_G);
            updateEtaH(result_etaH.sb18B.ETA_H);
            
            // Random step:
//            double range = randomRange/(i+1);
//            INITIAL_G_AUG_SCALE += (1-2*Math.random())*augGStep*range;
//            INITIAL_H_AUG_SCALE += (1-2*Math.random())*augHStep*range;
//            ETA_G               += (1-2*Math.random())*etaGStep*range;
//            ETA_H               += (1-2*Math.random())*etaHStep*range;
            
            ++i;
        } while(/*changed &&*/ 
                   (result_augG.gNorm > targetG || result_augG.hNorm > targetH || result_augG.gradNorm > targetGrad)
                && (result_augH.gNorm > targetG || result_augH.hNorm > targetH || result_augH.gradNorm > targetGrad)
                && (result_etaG.gNorm > targetG || result_etaG.hNorm > targetH || result_etaG.gradNorm > targetGrad)
                && (result_etaH.gNorm > targetG || result_etaH.hNorm > targetH || result_etaH.gradNorm > targetGrad)
                && i < MAX_ITERATIONS);
    }

    @Override
    protected void initGrid()
    {
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018B-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018B-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018B-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018B-loadData.csv");
        grid = GridGenerator.loadGrid("IEEE123-nbus", lineConfig, lineData, switchData, loadData);
        
        // Add slack bus:
        SlackSource slack = new SlackSource();
        slack.setVoltage(new Complex(BASE_VOLTAGE*SLACK_VOLTAGE));
        slack.setName(SLACK_SOURCE);
        grid.getBus(GridGenerator.BUS_PREFIX+149).addChild(slack);
        
        // Transformer:
        grid.setBaseVoltage(BASE_VOLTAGE);
        grid.setBasePower(BASE_POWER);
        
        // Add DG:
        addDG(1, 500e3, 100e3, 50e3);
        addDG(3, 500e3, 100e3, 50e3);
        addDG(8, 500e3, 100e3, 100e3);
        addDG(14, 500e3, 100e3, 50e3);
        addDG(18, 500e3, 100e3, 100e3);
        addDG(26, 500e3, 100e3, 50e3);
        addDG(29, 500e3, 100e3, 50e3);
        
        GridDisplay.showInFrame(grid).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }
}
