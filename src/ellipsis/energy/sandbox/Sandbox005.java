package ellipsis.energy.sandbox;

/**
 * Max. P
 * S.T. 0 < P < P_max
 * 
 * Max. P - l1(P - P_max) - l2(-P)
 */
public class Sandbox005
{
	public static void main(String[] args)
	{
		test01();
	}
	
	private static final double epsilon = 0.01;

	private static void test01()
	{
		double p = 0;
		double l1 = 1;
		double l2 = 1;
		
		for(int k = 0; k < 1000; ++k)
		{
			double nextP = nextP(p, l1, l2);
			double nextL2 = nextL2(l1, l2, p);
			double nextL1 = //1 + nextL2; // From KKT conditions
							nextL1(l1, l2, p);
			
			p = 
//				nextP < 0 ? 0 :
//				nextP > p_max ? p_max :
						nextP;
			l1 = Math.max(0, nextL1);
			l2 = Math.max(0, nextL2);
			
			log(k, p, l1, l2);
		}
	}

	private static void log(int k, double p, double l1, double l2)
	{
		double f = p + l2*p - l1*(p - p_max);
		System.out.println(k+","+p+","+l1+","+l2+","+f);
	}

	private static final double epsilon_l1 = epsilon;
	private static double nextL1(double l1, double l2, double p)
	{
		return l1 - epsilon_l1*dfdl1(p, l1, l2);
	}

	private static final double p_max = 1;
	private static double dfdl1(double p, double l1, double l2)
	{
		return p_max - p;
	}

	private static final double epsilon_l2 = epsilon;
	private static double nextL2(double l1, double l2, double p)
	{
		return l2 - epsilon_l2*dfdl2(p, l1, l2);
	}

	private static double dfdl2(double p, double l1, double l2)
	{
		return p;
	}

	private static final double epsilon_p = epsilon;
	private static double nextP(double p, double l1, double l2)
	{
		return p + epsilon_p*dfdp(p, l1, l2);
	}

	private static double dfdp(double p, double l1, double l2)
	{
		return 1 - l1 + l2;
	}
}