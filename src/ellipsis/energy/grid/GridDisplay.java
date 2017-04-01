package ellipsis.energy.grid;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.LayoutManager;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Collection;
import java.util.HashSet;
import java.util.Locale;
import java.util.Set;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexFormat;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.smartgrid.ControllableDemand;

public class GridDisplay extends JPanel implements ActionListener
{
	private static final String ANALYSE_ACTION_NAME = "Analyse";
    private static final long serialVersionUID = 8717271841183976878L;

	private static class BusUI extends JPanel
	{
		private class BusMouseAdapter extends MouseAdapter
		{
			private Point clickPos;

			@Override
			public void mousePressed(MouseEvent e)
			{
				clickPos = e.getPoint();
			}

			@Override
			public void mouseDragged(MouseEvent e)
			{
				Point drag = e.getPoint();
				Point delta = new Point(drag.x-clickPos.x, drag.y-clickPos.y);
				Point location = getLocation();
				setLocation(new Point(location.x+delta.x, location.y+delta.y));
				
				getParent().repaint();
			}

			@Override
			public void mouseEntered(MouseEvent e)
			{
				setBackground(new Color(0.0f, 0.4f, 0.4f));
			}

			@Override
			public void mouseExited(MouseEvent e)
			{
				setBackground(Color.BLACK);
			}

			@Override
			public void mouseClicked(MouseEvent e)
			{
				for (Unit child : BusUI.this.bus.getChildren())
				{
					if(child instanceof Transformer)
					{
						showChangeTapPositionDialog((Transformer)child);
						break;
					}
				}
			}
		}
		
		private static final long serialVersionUID = -5011064104121122220L;
		
		Bus bus;
		Set<BusUI> children = new HashSet<BusUI>();
		public BusUI(Bus bus)
		{
			this.bus = bus;
			setSize(20, 4);
			setOpaque(true);
			setBackground(Color.BLACK);
			
			// Draggable component:
			MouseAdapter mouse = new BusMouseAdapter();
			addMouseListener(mouse);
			addMouseMotionListener(mouse);
		}

		private void showChangeTapPositionDialog(Transformer child)
		{
			String newTap = JOptionPane.showInputDialog(this, "Please enter new tap position.");
			if(newTap == null) // Cancelled.
				return;
			double pos = Double.parseDouble(newTap);
			child.setT(pos);
			repaint();
		}
	}
	
	private static class GridLayout implements LayoutManager
	{
		static int MARGIN = 0;
		
		@Override
		public void addLayoutComponent(String name, Component comp) {}

		@Override
		public void removeLayoutComponent(Component comp) {}

		@Override
		public Dimension preferredLayoutSize(Container parent)
		{
			return parent.getParent().getSize();
		}

		@Override
		public Dimension minimumLayoutSize(Container parent)
		{
			return parent.getParent().getSize();
		}

		@Override
		public void layoutContainer(Container parent)
		{
			BusUI slack = ((GridDisplay)parent).slackBusUI;
			int width = parent.getWidth();
			int slackWidth = slack.getWidth();
			slack.setLocation(width/2-slackWidth/2, 20);
			layoutBus((GridDisplay)parent, slack, width-MARGIN*2);
		}

		private void layoutBus(GridDisplay parent, BusUI parentBus, int width)
		{
			int count = parentBus.children.size();
			if(count == 0)
				return;
			
			int gap = count == 1 ? 0 : width/count;

			int left = count == 1 ? 
					parentBus.getLocation().x+parentBus.getWidth()/2 : 
					parentBus.getLocation().x+parentBus.getWidth()/2-width/2;

			int x = left+gap/2;
			int y = parentBus.getLocation().y+100;

			for (BusUI bus : parentBus.children)
			{
				bus.setLocation(x-bus.getWidth()/2, y);
				x += gap;
				layoutBus(parent, bus, gap == 0 ? width : gap-MARGIN*2);
			}
		}
	}
	
	private BusUI slackBusUI;
	private Grid grid;
	public AnalysisResults results;
	private double basePower;
	private double baseVoltage;
	
	public GridDisplay(Grid grid, double basePower, double baseVoltage)
	{
		this.grid = grid;
		this.basePower = basePower;
		this.baseVoltage = baseVoltage;
		
		setLayout(new GridLayout());
		
		// Add components for each bus - start from the slack bus:
		Bus slack = grid.getSlackBus();
		slackBusUI = addBus(slack, null, new HashSet<Bus>());
		
		new Thread() {
			public void run() 
			{
				while(true)
				{
					repaint();
					try
					{
						Thread.sleep(1000);
					} catch (InterruptedException e)
					{
						e.printStackTrace();
					}
				}
			};
		}.start();
	}
	
	private BusUI addBus(Bus bus, BusUI parent, Set<Bus> added)
	{
		// Only add once:
		if(added.contains(bus))
			return null;
		added.add(bus);
		
		// Add to parent container:
		BusUI b = new BusUI(bus);
		if(parent != null)
			parent.children.add(b);
		add(b);
		
		// Add children:
		for (Line line : bus.getLines())
		{
			Bus childBus = line.getFromBus();
			if(childBus == b.bus)
				childBus = line.getToBus();
			addBus(childBus, b, added);
		}
		
		return b;
	}

	@Override
	protected void paintComponent(Graphics g)
	{
		super.paintComponent(g);
		
		g.setColor(Color.BLACK);
		
		BusUI bus = slackBusUI;
		paintBus(g, bus, new HashSet<BusUI>());
	}

	protected void paintBus(Graphics g, BusUI bus, Set<BusUI> hashSet)
	{
		// Only draw for each bus once:
		if(hashSet.contains(bus))
			return;
		hashSet.add(bus);
		
		// Draw labels:
		int x = bus.getX()+bus.getWidth()+1;
		int y = bus.getY()+bus.getHeight();
		StringBuffer label = new StringBuffer(bus.bus.getName());
		Collection<Unit> children = bus.bus.getChildren();
		if(!children.isEmpty())
		{
			label.append("<");
			for (Unit child : children)
			{
				if(
						child instanceof ControllableDemand || 
						(child instanceof DistributedSource && !(child instanceof ControllableDemand.CDSource)) || 
						(child instanceof Load && !(child instanceof ControllableDemand.CDLoad))
				  )
				{
					label.append(child.getName());
					label.append(',');
				}
				else if(child instanceof Transformer)
				{
					label.append("XFMR[m=");
					label.append(((Transformer)child).getT());
					label.append("],");
				}
			}
//			label.deleteCharAt(label.length()-1); // remove last comma
			label.append(">");
		}
		g.drawString(label.toString(), x, y);
		
		if(results != null && bus.bus.getSlackVoltage().abs() > 0)
			g.drawString(results.getBusPower(bus.bus.getName())+"W", x, y+14);
		else
			g.drawString(ComplexFormat.getInstance("j", Locale.getDefault()).format(bus.bus.netPower())+"W", x, y+14);
		
		if(results != null)
		{
			Complex v = results.getBusVoltage(bus.bus.getName());
			g.drawString(ComplexFormat.getInstance("j", Locale.getDefault()).format(v)+"V", x, y+28);
		}
		
		// Draw lines:
		Point loc1 = bus.getLocation();
		loc1.x += bus.getWidth()/2;
		for (BusUI b : bus.children)
		{
			Point loc2 = b.getLocation();
			loc2.x += b.getWidth()/2;
			g.drawLine(loc1.x, loc1.y, loc2.x, loc2.y);
			paintBus(g, b, hashSet);
		}
	}
	
	public static JFrame showInFrame(Grid grid)
	{
	    return showInFrame(grid, grid.getBasePower(), grid.getBaseVoltage());
	}

	public static JFrame showInFrame(Grid grid, double basePower, double baseVoltage)
	{
		GridDisplay gridDisplay = new GridDisplay(grid, basePower, baseVoltage);
		return gridDisplay.showInFrame();
	}

	public JFrame showInFrame()
	{
		JFrame frame = new JFrame(grid.getName());
		JScrollPane sp = new JScrollPane(this);
		frame.getContentPane().add(sp);
		frame.getContentPane().add(analysisButton(this), BorderLayout.SOUTH);
		frame.setExtendedState(JFrame.MAXIMIZED_BOTH);
		frame.setSize(400, 400);
		frame.setVisible(true);
		
		return frame;
	}

	private static Component analysisButton(GridDisplay gridDisplay)
	{
		JButton b = new JButton(ANALYSE_ACTION_NAME);
		b.setActionCommand(ANALYSE_ACTION_NAME);
		b.addActionListener(gridDisplay);
		return b;
	}

	@Override
	public void actionPerformed(ActionEvent action)
	{
		switch (action.getActionCommand())
		{
		case ANALYSE_ACTION_NAME:
			analyse();
			break;

		default:
			break;
		}
	}

	public void analyse()
	{
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
		lfa.setBasePower(basePower);
		lfa.setBaseVoltage(baseVoltage);
		lfa.setIterations(100);
		lfa.setTargetError(1e-6);
		
		results = lfa.analyse();
		if(!results.getDidConverge())
			System.err.println("Did not converge!");
	}
}