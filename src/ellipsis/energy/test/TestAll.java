package ellipsis.energy.test;

public class TestAll
{
	public static void main(String[] args) throws Exception
	{
		System.out.println("GridlessADPTests");
		GridlessADPTests.main(args);
		System.out.println("\nGridlessADPFactoryTest");
		GridlessADPFactoryTest.main(args);
		System.out.println("\nGridlessADPTests2");
		GridlessADPTests2.main(args);
		System.out.println("\nGridlessADPIEEE13BusGridTest");
		GridlessADPIEEE13BusGridTest.main(args);
		System.out.println("\nGridlessCentralControllerTest");
		GridlessCentralControllerTest.main(args);
		System.out.println("\nGridlessTransformerRegulatorTest");
		GridlessTrasformerRegulatorTest.main(args);
	}
}