using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace BipedRobot
{
    public partial class graph : Form
    {
        private Biped _biped;
        public graph(Biped biped)
        {
            _biped = biped;
            InitializeComponent();
            ShowDialog();
        }

        private void graph_Load(object sender, EventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e)
        {

            for (int i = 0; i < _biped.simulationData.RES.Length; i++)
            {
                if (_biped.simulationData.RES[i] == null)
                {
                    break;
                }
                double q1 = _biped.simulationData.RES[i].Item1[0];
                double q2 = _biped.simulationData.RES[i].Item1[1];
                double q3 = _biped.simulationData.RES[i].Item1[2];
                double time = _biped.simulationData.RES[i].Item2;
                dynamicsChart.Series["q1"].Points.AddXY(time, q1);
                dynamicsChart.Series["q2"].Points.AddXY(time, q2);
                dynamicsChart.Series["q3"].Points.AddXY(time, q3);
            }
            dynamicsChart.Series["q1"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            dynamicsChart.Series["q1"].Color = Color.Red;

            dynamicsChart.Series["q2"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            dynamicsChart.Series["q2"].Color = Color.Black;

            dynamicsChart.Series["q3"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            dynamicsChart.Series["q3"].Color = Color.Blue;
        }

        private void dynamicsChart_Click(object sender, EventArgs e)
        {

        }

        private void button2_Click(object sender, EventArgs e)
        {
            for (int i = 0; i < _biped.reducedSimulationData.RES.Count; i++)
            {
                double theta = _biped.reducedSimulationData.RES[i].Item1[0];
                double dtheta = _biped.reducedSimulationData.RES[i].Item1[1];
                double ddtheta = _biped.reducedSimulationData.RES[i].Item1[2];
                double time = _biped.reducedSimulationData.RES[i].Item2;
                zeroDynamics.Series["phaseportrait"].Points.AddXY(theta, dtheta);
                //zeroDynamics.Series["dtheta"].Points.AddXY(time, dtheta);
                //zeroDynamics.Series["ddtheta"].Points.AddXY(time, ddtheta);
            }
            zeroDynamics.Series["phaseportrait"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            zeroDynamics.Series["phaseportrait"].Color = Color.Red;

            //zeroDynamics.Series["dtheta"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            //zeroDynamics.Series["dtheta"].Color = Color.Black;

            //zeroDynamics.Series["ddtheta"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            //zeroDynamics.Series["ddtheta"].Color = Color.Blue;
        }
    }
}
