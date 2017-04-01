package ellipsis.energy.test;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.GridDisplay;
import ellipsis.energy.grid.Line;
import ellipsis.energy.grid.Load;
import ellipsis.energy.grid.SlackSource;
import ellipsis.energy.grid.Switch;
import ellipsis.energy.grid.Transformer;

public class IEEE13BusGrid extends Grid
{
	public static final double BASE_VOLTAGE = 4.16e3;
	public static final double BASE_POWER = 1e6;
    public static final double BASE_IMPEDANCE = BASE_VOLTAGE*BASE_VOLTAGE/BASE_POWER;
	
	// Busses:
	public Bus bus611 = new Bus(this);
	public Bus bus632 = new Bus(this);
	public Bus bus633 = new Bus(this);
	public Bus bus634 = new Bus(this);
	public Bus bus645 = new Bus(this);
	public Bus bus646 = new Bus(this);
	public Bus bus650 = new Bus(this);
	public Bus bus650xfm = new Bus(this);
	public Bus bus652 = new Bus(this);
	public Bus bus671 = new Bus(this);
	public Bus bus675 = new Bus(this);
	public Bus bus680 = new Bus(this);
	public Bus bus684 = new Bus(this);
	public Bus bus692 = new Bus(this);
	
	// Lines:
	public Line line650xfm_632 = new Line();
	public Line line632_633 = new Line();
	public Line line632_645 = new Line();
	public Line line632_671 = new Line();
	public Line line645_646 = new Line();
	public Line line671_684 = new Line();
	public Line line671_680 = new Line();
	public Line line692_675 = new Line();
	public Line line684_611 = new Line();
	public Line line684_652 = new Line();
	
	// Transformers:
	public Transformer xfm650 = new Transformer();
	public Transformer xfm633_634 = new Transformer();
	
	// Switch:
	public Switch sw671_692 = new Switch();
	
	// Loads:
	public Load load634 = new Load();
	public Load load645 = new Load();
	public Load load646 = new Load();
	public Load load652 = new Load();
	public Load load671 = new Load();
	public Load load675 = new Load();
	public Load load692 = new Load();
	public Load load611 = new Load();

	public IEEE13BusGrid()
	{
		setBasePower(BASE_POWER);
		setBaseVoltage(BASE_VOLTAGE);
		
		// Set names:
		bus611.setName("611");
		bus632.setName("632");
		bus633.setName("633");
		bus634.setName("634");
		bus645.setName("645");
		bus646.setName("646");
		bus650.setName("650");
		bus650xfm.setName("650xfm");
		bus652.setName("652");
		bus671.setName("671");
		bus675.setName("675");
		bus680.setName("680");
		bus684.setName("684");
		bus692.setName("692");
		
		line632_633.setName("632-633");
		line632_645.setName("632-645");
		line632_671.setName("632-672");
		line645_646.setName("645-646");
		line650xfm_632.setName("650-632");
		line671_680.setName("671-680");
		line671_684.setName("671-684");
		line684_611.setName("684-611");
		line684_652.setName("684-652");
		line692_675.setName("692-675");
		
		sw671_692.setName("671-692");
		
		xfm633_634.setName("633-634");
		xfm650.setName("650-");
		
		load611.setName("L611");
		load634.setName("L634");
		load645.setName("L645");
		load646.setName("L646");
		load652.setName("L652");
		load671.setName("L671");
		load675.setName("L675");
		load692.setName("L692");
		
		// Set slack:
		Complex slackVoltage = new Complex(BASE_VOLTAGE);
		SlackSource slack = new SlackSource();
		slack.setName("slack");
		slack.setVoltage(slackVoltage);
		bus650.addChild(slack);
		
		// Set switch:
		sw671_692.setResistance(1e-3);
		
		// Set lines:
		Map<Integer, Complex> impedences = new HashMap<>();
		
		//             Config.          Res. (ohms/km)
		impedences.put(601, new Complex(0.11625, 0.01));
		impedences.put(602, new Complex(0.37, 0.04));
		impedences.put(603, new Complex(0.7, 0.07));
		impedences.put(604, new Complex(0.7, 0.07));
		impedences.put(605, new Complex(0.7, 0.07));
		impedences.put(606, new Complex(0.25625, 0.02));
		impedences.put(607, new Complex(0.379375, 0.03));
		
		line632_645.setImpedencePerMetre(impedences.get(603));
		line632_633.setImpedencePerMetre(impedences.get(602));
		line645_646.setImpedencePerMetre(impedences.get(603));
		line650xfm_632.setImpedencePerMetre(impedences.get(601));
		line684_652.setImpedencePerMetre(impedences.get(607));
		line632_671.setImpedencePerMetre(impedences.get(601));
		line671_684.setImpedencePerMetre(impedences.get(604));
		line671_680.setImpedencePerMetre(impedences.get(601));
		line684_611.setImpedencePerMetre(impedences.get(605));
		line692_675.setImpedencePerMetre(impedences.get(606));

		line632_645.setLength(152.4390243902e-3);
		line632_633.setLength(152.4390243902e-3);
		line645_646.setLength(91.4634146341e-3);
		line650xfm_632.setLength(609.756097561e-3);
		line684_652.setLength(243.9024390244e-3);
		line632_671.setLength(609.756097561e-3);
		line671_684.setLength(91.4634146341e-3);
		line671_680.setLength(304.8780487805e-3);
		line684_611.setLength(91.4634146341e-3);
		line692_675.setLength(152.4390243902e-3);
		
		// Set loads:
		load634.setLoad(400e3,290e3);
		load645.setLoad(170e3,125e3);
		load646.setLoad(230e3,132e3);
		load652.setLoad(128e3,86e3);
		load671.setLoad(1155e3,660e3);
		load675.setLoad(843e3,462e3);
		load692.setLoad(170e3,151e3);
		load611.setLoad(170e3,80e3);
		
		bus634.addChild(load634);
		bus645.addChild(load645);
		bus646.addChild(load646);
		bus652.addChild(load652);
		bus671.addChild(load671);
		bus675.addChild(load675);
		bus692.addChild(load692);
		bus611.addChild(load611);
		
		// Setup transformers:
//		xfm650.setBaseVoltageRatio(1);
//		xfm650.setMinimumT(0.9);
//		xfm650.setMaximumT(1.1);
//		xfm650.setImpedance(new Complex(0.01*BASE_IMPEDANCE, 0.08*BASE_IMPEDANCE));
		xfm650.setBaseVoltageRatio(1);
		xfm650.setImpedance(new Complex(0.001, 0.008));//new Complex(0.03, 0.27));
		xfm650.setMinimumT(0.9);
		xfm650.setMaximumT(1.1);
		xfm650.setT(1);
		xfm650.setTStepSize(0.01);
		xfm633_634.setBaseVoltageRatio(1);
		xfm633_634.setImpedance(new Complex(0.035, 0.06));
		
		// Construct network:
		xfm650.setFromBus(bus650);
		xfm650.setToBus(bus650xfm);
		
		line650xfm_632.setFromBus(bus650xfm);
		line650xfm_632.setToBus(bus632);

		line632_633.setFromBus(bus632);
		line632_633.setToBus(bus633);
		line632_645.setFromBus(bus632);
		line632_645.setToBus(bus645);
		line632_671.setFromBus(bus632);
		line632_671.setToBus(bus671);
		
		line645_646.setFromBus(bus645);
		line645_646.setToBus(bus646);
		
		xfm633_634.setFromBus(bus633);
		xfm633_634.setToBus(bus634);
		
		line671_680.setFromBus(bus671);
		line671_680.setToBus(bus680);
		line671_684.setFromBus(bus671);
		line671_684.setToBus(bus684);
		
		sw671_692.setFromBus(bus671);
		sw671_692.setToBus(bus692);
		line692_675.setFromBus(bus692);
		line692_675.setToBus(bus675);
		
		line684_611.setFromBus(bus684);
		line684_611.setToBus(bus611);
		line684_652.setFromBus(bus684);
		line684_652.setToBus(bus652);
	}
	
	public static void main(String[] args)
	{
        double basePower = 1e6;
        
		Grid grid = new IEEE13BusGrid();
		
		GridDisplay.showInFrame(grid, basePower, BASE_VOLTAGE);
	}
}
