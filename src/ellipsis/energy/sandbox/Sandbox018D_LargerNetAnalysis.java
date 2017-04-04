package ellipsis.energy.sandbox;

import static ellipsis.common.math.VectorHelper.maxAbs;
import static ellipsis.energy.sandbox.GridConstants.SLACK_SOURCE;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.PrintStream;

import org.apache.commons.math3.complex.Complex;

import com.mls.util.Util;

import ellipsis.energy.grid.GridGenerator;
import ellipsis.energy.grid.SlackSource;
import ellipsis.util.EmptyOutputStream;
import ellipsis.util.TeeOutputStream;

public class Sandbox018D_LargerNetAnalysis extends Sandbox018B
{
	private int gConvergedIteration = -1;
	private double gConvergeThreshold = 5e-1;
	private int gradConvergedIteration = -1;
	private double gradConvergeThreshold = 5e-1;
	private int hConvergedIteration = -1;
	private double hConvergeThreshold = 1e-1;
	private int consensusConvergedIteration = -1;
	private double consensusConvergeThreshold = 1e-1;
	
	public Sandbox018D_LargerNetAnalysis()
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
                out = new PrintStream(new TeeOutputStream(new FileOutputStream("/tmp/Sandbox018D.csv"), System.out));
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
                 
                      INITIAL_G_AUG_SCALE = 0.0037;//0.00075;//0.001;
                         G_AUG_SCALE_STEP = 1.001;
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 8;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 3.2;//3;
                                    ETA_H = 0.0051;//0.005;
                                   
                                       XI = 0.2;
                               
                                        K = 5000;//3000;
                               DEBUG_RATE = K / 1000;
              AGENT_SELECTION_PROBABILITY = 1.0;
             
                        INITIAL_STEP_SIZE = 1;
                            MIN_STEP_SIZE = 1e-50;
                              MAX_X_STEPS = 1;
                              
                          EPSILON_BASE_PQ = 100;
                          EPSILON_BASE_EF = 1000;
                           EPSILON_TARGET = 1;
    }
	
	@Override
	protected void debug(int k)
	{
		super.debug(k);
		
if(k > 2900)
	Util.nullop();
		
		// Check for convergence:
		int dimension = p.getDimension();
		double maxG = Math.max(maxAbs(dimension, this::gp), maxAbs(dimension, this::gq));
		double hp = hp(p, q, e, f);
		double hq = hq(p, q, e, f);
		double maxH = Math.max(Math.abs(hp), Math.abs(hq));
		double gradNorm = gradNorm();
		Agent[] as = agents.toArray(new Agent[agents.size()]);
		double consensusError = Math.max(maxAbs(dimension, i -> as[i].hp() - hp), maxAbs(dimension, i -> as[i].hq() - hq));
		
		if(maxG < gConvergeThreshold && gConvergedIteration == -1)
			gConvergedIteration = k;
		if(maxH < hConvergeThreshold && hConvergedIteration == -1)
			hConvergedIteration = k;
		if(gradNorm < gradConvergeThreshold && gradConvergedIteration == -1)
			gradConvergedIteration = k;
		if(maxG < gConvergeThreshold && gConvergedIteration == -1)
			gConvergedIteration = k;
		if(consensusError < consensusConvergeThreshold && consensusConvergedIteration == -1)
			consensusConvergedIteration = k;
	}

	public static void main(String[] args)
    {
		tune = false;
		
		if(tune)
        {
			result_augG = new TuneResult(new Sandbox018D_LargerNetAnalysis());
		    result_augH = new TuneResult(new Sandbox018D_LargerNetAnalysis());
		    result_etaG = new TuneResult(new Sandbox018D_LargerNetAnalysis());
		    result_etaH = new TuneResult(new Sandbox018D_LargerNetAnalysis());
		    
            // Tune parameters and then run again to get results:
            System.out.println("Tuning...");
            tune();
            System.out.println("Tuning complete.\nRunning...");
            new Sandbox018B().run();
            System.out.println("Running complete.");
        }
        else if(analyseLambdaStar)
        {
            Sandbox018B sandbox018b = new Sandbox018D_LargerNetAnalysis();
            sandbox018b.start();
            sandbox018b.init();
            sandbox018b.analyseLambdaStar();
            sandbox018b.finish();
        }
        else
        {
            new Sandbox018D_LargerNetAnalysis().run();
        }
    }
	
	@Override
	protected void finish()
	{
		out.println(agents.size()+" bus network:");
		out.println("g() convergence achieved to threshold of "+gConvergeThreshold+" at iteration "+gConvergedIteration);
		out.println("h() convergence achieved to threshold of "+hConvergeThreshold+" at iteration "+hConvergedIteration);
		out.println("||grad|| convergence achieved to threshold of "+gradConvergeThreshold+" at iteration "+gradConvergedIteration);
		out.println("Consensus convergence achieved to threshold of "+consensusConvergeThreshold+" at iteration "+consensusConvergedIteration);
		
		super.finish();
	}
	
	@Override
    protected void initGrid()
    {
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018D-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018D-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018D-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018D-loadData.csv");
        grid = GridGenerator.loadGrid("IEEE123-subset", lineConfig, lineData, switchData, loadData);
        
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
        
        addDG(36, 500e3, 100e3, 50e3);
        
//        GridDisplay.showInFrame(grid).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }
}
