package ellipsis.energy.sandbox;

import static ellipsis.common.math.VectorHelper.maxAbs;
import static ellipsis.energy.sandbox.GridConstants.SLACK_SOURCE;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Timer;
import com.mls.util.Util;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.grid.GridGenerator;
import ellipsis.energy.grid.SlackSource;
import ellipsis.util.EmptyOutputStream;
import ellipsis.util.TeeOutputStream;

public class Sandbox018D_LargerNetAnalysis extends Sandbox018B
{
	static enum Network
	{
		_35BUS,
		_44BUS,
		_52BUS,
		_60BUS,
		_67BUS,
		_77BUS,
		_88BUS,
		_97BUS,
		_107BUS,
		_117BUS
	}
	
	public Network testNetwork = Network._107BUS;
	public boolean convergenceAnalysis = false;
	
	private int gConvergedIteration = -1;
	private double gConvergeThreshold = 5e-1;
	private int gradConvergedIteration = -1;
	private double gradConvergeThreshold = 5e-1;
	private int hConvergedIteration = -1;
	private double hConvergeThreshold = 1e-1;
	private int consensusConvergedIteration = -1;
	private double consensusConvergeThreshold = 1e-1;
	
	public Sandbox018D_LargerNetAnalysis(Network testNetwork, boolean convergenceAnalysis, PrintStream out)
    {
        this.testNetwork = testNetwork;
        this.convergenceAnalysis = convergenceAnalysis;
        
        initD(out);
    }
	
	public Sandbox018D_LargerNetAnalysis() 
	{
		initD(null);
	}
	
	
	private void initD(PrintStream out) 
	{
		this.out = out;
		
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
		
        if(this.out == null)
        {
			if(tune)
	        {
	            // No output from test cases:
				this.out = new PrintStream(new EmptyOutputStream());
	        }
	        else
	        {
	            try
	            {
	            	this.out = new PrintStream(new TeeOutputStream(new FileOutputStream("/tmp/Sandbox018D.csv"), System.out));
	            }
	            catch (FileNotFoundException e)
	            {
	                throw new RuntimeException(e);
	            }
	        }
        }
        
        //// Override configuration ////
        switch (testNetwork)
		{
		case _35BUS:
			// Config already performed above
			break;
		case _44BUS:
			config44();
			break;
		case _52BUS:
			config52();
			break;
		case _60BUS:
			config60();
			break;
		case _67BUS:
			config67();
			break;
		case _77BUS:
			config77();
			break;
		case _88BUS:
			config88();
			break;
		case _97BUS:
			config97();
			break;
		case _107BUS:
			config107();
			break;
		case _117BUS:
			config117();
			break;
		default:
			break;
		}
	}
	
	@Override
	protected void testResults(AnalysisResults results) 
	{
		if(!convergenceAnalysis)
			super.testResults(results);
	}
	
	@Override
	public void loop()
    {
		if(convergenceAnalysis)
		{
			// Set up all test cases to run concurrently:
			List<Sandbox018D_LargerNetAnalysis> testCases = new ArrayList<Sandbox018D_LargerNetAnalysis>();
			for (Network net : Network.values()) 
			{
				Sandbox018D_LargerNetAnalysis testCase = new Sandbox018D_LargerNetAnalysis(net, true, out)
				{
					@Override
					protected void debugHeader(AnalysisResults results) {} // Avoid logging headers.
				};
				testCase.init();
				testCases.add(testCase);
			}
			
			// Iterations:
	        Random agentSelector = new Random(0);
	        Timer timer = Timer.startNewTimer();
	        for(int k = 0; k < 10000; ++k)
	        {
	        	for (Sandbox018D_LargerNetAnalysis testCase : testCases) 
	        	{
	        		testCase.iteration(agentSelector, timer, k);
				}
	        	storeConvergenceRow(k, testCases);
	        }
	        
	        // Prepare and log convergence data:
	        for (Sandbox018D_LargerNetAnalysis testCase : testCases)
	        {
	        	testCase.finalState = testCase.states.get(testCase.states.size()-1);
	        }
	        
	        for(int k = 0; k < K; ++k)
	        {
	        	logConvergenceRow(k, testCases);
	        }
		}
		else
		{
			super.loop();
		}
    }

    private void storeConvergenceRow(int k, List<Sandbox018D_LargerNetAnalysis> testCases) 
    {
    	for (Sandbox018D_LargerNetAnalysis testCase : testCases)
    	{
    		testCase.storeConvergenceRow();
    	}
	}

    private ArrayList<Double> consensusErrors = new ArrayList<>();
    private ArrayList<RealVector> states = new ArrayList<>();
    private RealVector finalState;
	private void storeConvergenceRow() 
	{
		consensusErrors.add(consensusError()); 
		states.add(state());
	}

	protected RealVector state() 
	{
		return p
		  .append(q)
		  .append(e)
		  .append(f);
	}

	protected void logConvergenceRow(int k, List<Sandbox018D_LargerNetAnalysis> testCases) 
    {
    	out.print(k);
    	out.print(",");
		for (Sandbox018D_LargerNetAnalysis testCase : testCases) 
		{
			// Consensus error:
			out.print(testCase.consensusErrors.get(k));
			out.print(",");
			
			// State norm:
			RealVector state = testCase.states.get(k);
			out.print(state.subtract(testCase.finalState).getNorm());
			out.print(",");
		}
		out.println();
	}

	private void config()
    {
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
    }
    
    private void config44()
    {
    	config();
                 
                      INITIAL_G_AUG_SCALE = 0.0077;//0.0037;//0.00075;//0.001;
                         G_AUG_SCALE_STEP = 1.001;
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 7.7;//8;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 3.6;//3.2;//3;
                                    ETA_H = 0.005;//0.0051;//0.005;
                                   
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
    
    private void config52()
    {
    	config();
                 
                      INITIAL_G_AUG_SCALE = 0.0097;
                         G_AUG_SCALE_STEP = 1.001;
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 7.8;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 3.7;
                                    ETA_H = 0.0051;
                                   
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
    
    private void config60()
    {
    	config();
                 
                      INITIAL_G_AUG_SCALE = 0.0107;
                         G_AUG_SCALE_STEP = 1.001;
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 7.9;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 3.8;
                                    ETA_H = 0.0052;
                                   
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
    
    private void config67()
    {
    	config();
                 
                      INITIAL_G_AUG_SCALE = 0.0127;
                         G_AUG_SCALE_STEP = 1.001;
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 8.2;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 4.1;
                                    ETA_H = 0.0056;
                                   
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
    
    private void config77()
    {
    	config();
                 
                      INITIAL_G_AUG_SCALE = 0.0127;
                         G_AUG_SCALE_STEP = 1.001;
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 8.8;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 4.1;
                                    ETA_H = 0.0063;
                                   
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
    
    private void config88()
    {
    	config();
                 
                      INITIAL_G_AUG_SCALE = 0.0127;
                         G_AUG_SCALE_STEP = 1.001;
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 7;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 4.2;
                                    ETA_H = 0.0085;
                                   
                                       XI = 0.2;
                               
                                        K = 7000;//3000;
                               DEBUG_RATE = K / 1000;
              AGENT_SELECTION_PROBABILITY = 1.0;
             
                        INITIAL_STEP_SIZE = 1;
                            MIN_STEP_SIZE = 1e-50;
                              MAX_X_STEPS = 1;
                              
                          EPSILON_BASE_PQ = 100;
                          EPSILON_BASE_EF = 1000;
                           EPSILON_TARGET = 1;
    }
    
    private void config97()
    {
    	config();
                 
                      INITIAL_G_AUG_SCALE = 0.0127;
                         G_AUG_SCALE_STEP = 1.001;
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 7;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 4.2;
                                    ETA_H = 0.0085;
                                   
                                       XI = 0.2;
                               
                                        K = 7000;//3000;
                               DEBUG_RATE = K / 1000;
              AGENT_SELECTION_PROBABILITY = 1.0;
             
                        INITIAL_STEP_SIZE = 1;
                            MIN_STEP_SIZE = 1e-50;
                              MAX_X_STEPS = 1;
                              
                          EPSILON_BASE_PQ = 100;
                          EPSILON_BASE_EF = 1000;
                           EPSILON_TARGET = 1;
    }
    
    private void config107()
    {
    	config();
                 
                      INITIAL_G_AUG_SCALE = 0.0127;
                         G_AUG_SCALE_STEP = 1.001;
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 7;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 4.2;
                                    ETA_H = 0.0085;
                                   
                                       XI = 0.2;
                               
                                        K = 8000;//3000;
                               DEBUG_RATE = K / 1000;
              AGENT_SELECTION_PROBABILITY = 1.0;
             
                        INITIAL_STEP_SIZE = 1;
                            MIN_STEP_SIZE = 1e-50;
                              MAX_X_STEPS = 1;
                              
                          EPSILON_BASE_PQ = 100;
                          EPSILON_BASE_EF = 1000;
                           EPSILON_TARGET = 1;
    }
    
    private void config117()
    {
    	config();
                 
                      INITIAL_G_AUG_SCALE = 0.0127;
                         G_AUG_SCALE_STEP = 1.001;
                          G_MAX_AUG_SCALE = 1e6;
                      
                      INITIAL_H_AUG_SCALE = 7;
                         H_AUG_SCALE_STEP = 1.001;
                          H_MAX_AUG_SCALE = 1e6;
               
                                    ETA_G = 4.2;
                                    ETA_H = 0.0085;
                                   
                                       XI = 0.2;
                               
                                        K = 9000;//3000;
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
		if(convergenceAnalysis)
			return;
		
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
            new Sandbox018D_LargerNetAnalysis().run();
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
		switch (testNetwork)
		{
		case _35BUS:
			super.initGrid();
			break;
		case _44BUS:
			initGrid44();
			break;
		case _52BUS:
			initGrid52();
			break;
		case _60BUS:
			initGrid60();
			break;
		case _67BUS:
			initGrid67();
			break;
		case _77BUS:
			initGrid77();
			break;
		case _88BUS:
			initGrid88();
			break;
		case _97BUS:
			initGrid97();
			break;
		case _107BUS:
			initGrid107();
			break;
		case _117BUS:
			initGrid117();
			break;
		default:
			break;
		}
    }
	
	protected void initGrid44()
	{
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018D44-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018D44-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018D44-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018D44-loadData.csv");
        initGridWith(lineConfig, lineData, switchData, loadData);
    }
	
	protected void initGrid52()
	{
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018D-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018D52-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018D-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018D52-loadData.csv");
        initGridWith(lineConfig, lineData, switchData, loadData);
        
//        GridDisplay.showInFrame(grid).setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }
	
	protected void initGrid60()
	{
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018D-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018D60-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018D-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018D60-loadData.csv");
        initGridWith(lineConfig, lineData, switchData, loadData);
    }
	
	protected void initGrid67()
	{
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018D-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018D67-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018D-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018D67-loadData.csv");
        initGridWith(lineConfig, lineData, switchData, loadData);
    }
	
	protected void initGrid77()
	{
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018D-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018D77-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018D-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018D77-loadData.csv");
        initGridWith(lineConfig, lineData, switchData, loadData);
    }
	
	protected void initGrid88()
	{
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018D-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018D88-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018D-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018D88-loadData.csv");
        initGridWith(lineConfig, lineData, switchData, loadData);
    }
	
	protected void initGrid97()
	{
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018D-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018D97-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018D-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018D97-loadData.csv");
        initGridWith(lineConfig, lineData, switchData, loadData);
    }
	
	protected void initGrid107()
	{
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018D-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018D107-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018D-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018D107-loadData.csv");
        initGridWith(lineConfig, lineData, switchData, loadData);
    }
	
	protected void initGrid117()
	{
        InputStream lineConfig = getClass().getResourceAsStream("Sandbox018D-lineConfig");
        InputStream lineData = getClass().getResourceAsStream(  "Sandbox018D117-lineData.csv");
        InputStream switchData = getClass().getResourceAsStream("Sandbox018D-switchData.csv");
        InputStream loadData = getClass().getResourceAsStream(  "Sandbox018D117-loadData.csv");
        initGridWith(lineConfig, lineData, switchData, loadData);
    }

	protected void initGridWith(InputStream lineConfig, InputStream lineData, InputStream switchData, InputStream loadData)
	{
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
        
        switch (testNetwork)
		{
        case _117BUS:
        	addDG(110, 500e3, 100e3, 50e3);
        case _107BUS:
        	addDG(102, 500e3, 100e3, 50e3);
        	addDG(99, 500e3, 100e3, 50e3);
        case _97BUS:
        	addDG(74, 500e3, 100e3, 50e3);
        	addDG(68, 500e3, 100e3, 50e3);
        case _88BUS:
        	addDG(76, 500e3, 100e3, 50e3);
        	addDG(81, 500e3, 100e3, 50e3);
        case _77BUS:
        	addDG(93, 500e3, 100e3, 50e3);
        	addDG(87, 500e3, 100e3, 50e3);
        case _67BUS:
        	addDG(60, 500e3, 100e3, 50e3);
        	addDG(64, 500e3, 100e3, 50e3);
        case _60BUS:
            addDG(54, 500e3, 100e3, 50e3);
            addDG(57, 500e3, 100e3, 50e3);
        case _52BUS:
        	addDG(42, 500e3, 100e3, 50e3);
        	addDG(47, 500e3, 100e3, 50e3);
		case _44BUS:
			addDG(36, 500e3, 100e3, 50e3);
			break;

		default:
			break;
		}
	}
	
	@Override
	public void logConfig() 
	{
		if(!convergenceAnalysis)
			super.logConfig();
	}
	
	@Override
	public void logBusses(AnalysisResults results) 
	{
		if(!convergenceAnalysis)
			super.logBusses(results);
	}
	
	@Override
	protected void debugHeaderRow(Map<String, Integer> numbers) 
	{
		if(convergenceAnalysis)
		{
			out.println("Consensus convergence: Iteration at which the error permanently drops below the threshold percentage of total system DG power.");
			out.println("OPF convergence: Iteration at which the norm of the mean difference between final state value and iteration state value drops below the threshold percentage.");
			out.println();
			out.print("Test case:,");
			for (Network net : Network.values()) 
			{
				out.print(net);
				out.print(",,");
			}
			out.println();
			
			out.print("Iteration,");
			for (@SuppressWarnings("unused") Network net : Network.values()) 
			{
				out.print("Con. Error,State Norm,");
			}
		}
		else
		{
			super.debugHeaderRow(numbers);
		}
	}
}
