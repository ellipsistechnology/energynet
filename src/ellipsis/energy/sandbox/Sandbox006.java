package ellipsis.energy.sandbox;

/**
 * Max. P + Q
 * S.T. P^2 + Q^2 = 1
 * 
 * Max. P + Q + l(P^2 + Q^2 - 1)
 */
public class Sandbox006
{
	public static void main(String[] args)
	{
		test01();
	}
	
	private static final double epsilon = 0.1;

	private static void test01()
	{
		double p = 1;
		double q = -1;
		double l = 1;
		
		for(int k = 0; k < 100; ++k)
		{
			p = nextP(p, q, l);
			q = nextQ(p, q, l);
			l = nextL(p, q, l);
			
			log(k, p, q, l);
		}
	}

	private static void log(int k, double p, double q, double l)
	{
		double f = p + q + l*(p*p + q*q -1);
		System.out.println(k+","+p+","+q+","+l+","+f);
	}

	private static final double epsilon_l = epsilon;
	private static double nextL(double p, double q, double l)
	{
		return l - epsilon_l*dfdl(p, q, l);
	}

	private static double dfdl(double p, double q, double l)
	{
		return p*p + q*q -1;
	}

	private static final double epsilon_p = epsilon;
	private static double nextP(double p, double q, double l)
	{
		return p + epsilon_p*dfdp(p, q, l);
	}

	private static double dfdp(double p, double q, double l)
	{
		return 1 + 2*l*p;
	}

	private static final double epsilon_q = epsilon;
	private static double nextQ(double p, double q, double l)
	{
		return q + epsilon_q*dfdq(p, q, l);
	}

	private static double dfdq(double p, double q, double l)
	{
		return 1 + 2*l*q;
	}
}