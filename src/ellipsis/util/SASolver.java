package ellipsis.util;

import java.util.Random;

/**
 * Simulated Annealing Solver.
 * @author bmillar
 *
 */
public class SASolver
{
    private static final boolean TEST = false;
    
	private double minT = 0.0001;
	private double startT = 1;
	private double TMultiplier = 0.5;
	private int maxAttempts = 10000;
	
	private static final int SOLUTIONS = 4;
	
	@SuppressWarnings("unchecked")
    public <X> X solve(SAProblem<X> problem)
	{
		Random rand = new Random(1);
		Object xs[] = new Object[SOLUTIONS];
		double costs[] = new double[SOLUTIONS];
		double minCost = Double.MAX_VALUE;
		X minSolution = null;
		
		if(TEST) System.out.println("A,B,C,D");

		for (int i = 0; i < SOLUTIONS; i++)
        {
		    xs[i] = problem.getInitialSolution();
		    costs[i] = problem.cost((X) xs[i]);
		    
            if(TEST) System.out.print(costs[i]);
            if(TEST) System.out.print(',');
        }
		
		if(TEST) System.out.println();
		
		for(double T = startT; T > minT; T = TMultiplier*T)
		{
            problem. prepareNextSolution(T);
            
		    for (int s = 0; s < SOLUTIONS; s++)
		    {
    			for(int i = 0; i < maxAttempts; ++i) // Loop will end if we find an acceptable next solution.
    			{
    				// Choose a random solution:
    				X nextX = problem.nextSolution((X) xs[s], i);
    					
    				// Calculate solution cost:
    				double nextCost = problem.cost(nextX);
    					
    				// If cost is lower than previous cost accept it:
    				if(nextCost < costs[s])
    				{
    					xs[s] = nextX;
    					costs[s] = nextCost;
    					break;
    				}
    				// Otherwise randomly decide if we keep it:
    				else
    				{
    					double probability = Math.exp((costs[s]-nextCost)/T); // e^(delta C/T) => [0-1]; delta C is negative
//System.out.println(costs[s]+","+(costs[s]-nextCost));
    					double r = rand.nextDouble();
    					if(r < probability)
    					{
    						xs[s] = nextX;
    						costs[s] = nextCost;
    						break;
    					}
    					// else try again with next solution (start at beginning of while loop without reducing T)
    				}
    			}
    			
    			if(costs[s] < minCost)
    			{
    			    minCost = costs[s];
    			    minSolution = (X) xs[s];
    			}

                if(TEST) System.out.print(costs[s]);
                if(TEST) System.out.print(',');
		    }
	        if(TEST) System.out.println();
		}
        if(TEST) System.out.print("Final cost, "+minCost+",");
		
		return problem.finish(minSolution);
	}
	
	public static void main(String[] args)
	{
		SASolver solver = new SASolver();
		double[] prefixes = solver.solve(new SAProblem<double[]>()
		{
			double neighbourhoodRadius = 50;
			Random rand = new Random(1);
			double[] sampleX = new double[]{0, 6,    19,     35,      47};
			double[] sampleY = new double[]{5, 1.16, -68.91, -108.75, 32.73};
			
			@Override
			public void prepareNextSolution(double T)
			{
			    neighbourhoodRadius *= 0.999;
			}
			
			@Override
			public double[] nextSolution(double[] prefixes, int attempt)
			{
				double[] newPrefixes = new double[4];
				newPrefixes[0] = nextPrefix(prefixes[0]);
				newPrefixes[1] = nextPrefix(prefixes[1]);
				newPrefixes[2] = nextPrefix(prefixes[2]);
				newPrefixes[3] = nextPrefix(prefixes[3]);

				return newPrefixes;
			}

			private double nextPrefix(double d)
			{
				return d + (2*rand.nextDouble() - 1)*neighbourhoodRadius;
			}
			
			@Override
			public double[] getInitialSolution()
			{
				return new double[]{0, 0, 0, 0};
			}
			
			@Override
			public double cost(double[] prefixes)
			{
				double cost = 0;
				
				for (int i = 0; i < sampleX.length; i++)
				{
					double x = sampleX[i];
					double y = prefixes[0]*x*x*x + prefixes[1]*x*x + prefixes[2]*x + prefixes[3];
					double diff = y - sampleY[i];
					cost += diff*diff;
				}
				
				return cost/5e6;
			}

            @Override
            public double[] finish(double[] x)
            {
                // Nothing to do
                return null;
            }
		});
		
		System.out.println();
		System.out.println("a,b,c,d");
		System.out.println(prefixes[0]+"," + prefixes[1]+"," + prefixes[2]+"," + prefixes[3]);
	}

    public double getMinT()
    {
        return minT;
    }

    public void setMinT(double minT)
    {
        this.minT = minT;
    }

    public double getStartT()
    {
        return startT;
    }

    public void setStartT(double startT)
    {
        this.startT = startT;
    }

    public double getTMultiplier()
    {
        return TMultiplier;
    }

    public void setTMultiplier(double tMultiplier)
    {
        TMultiplier = tMultiplier;
    }

    public int getMaxAttempts()
    {
        return maxAttempts;
    }

    public void setMaxAttempts(int maxAttempts)
    {
        this.maxAttempts = maxAttempts;
    }
}
