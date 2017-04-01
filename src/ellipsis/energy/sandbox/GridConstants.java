package ellipsis.energy.sandbox;

public class GridConstants
{
	public static final String SLACK_SOURCE = "Slack Source";
    public static final String SLACK_BUS = "Slack Bus";
    public static final double BASE_POWER = 100e6;
    public static final double BASE_VOLTAGE = 1e3;
    public static final double BASE_IMPEDANCE = BASE_VOLTAGE*BASE_VOLTAGE/BASE_POWER;
    
    public static final String LINE_2_3 = "2-3";
    public static final String LINE_2_4 = "2-4";
    public static final String LINE_3_5 = "3-5";
    public static final String LINE_1_2 = "1-2";
    public static final String LINE_1_3 = "1-3";
    
    public static final String LOAD_2 = "Load 2";
    public static final String LOAD_5 = "Load 5";
    
    public static final String DG_3 = "DG 3";
    public static final String DG_4 = "DG 4";
    
    public static final String BUS_2 = "Bus 2";
    public static final String BUS_3 = "Bus 3";
    public static final String BUS_4 = "Bus 4";
    public static final String BUS_5 = "Bus 5";
}