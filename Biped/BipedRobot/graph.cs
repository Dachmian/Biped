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

        private void button1_Click_1(object sender, EventArgs e)
        {
            for (int i = 0; i < _biped.reducedSimulationData.RES.Count; i++)
            {
                double theta = _biped.reducedSimulationData.RES[i].Item1[0];
                double dtheta = _biped.reducedSimulationData.RES[i].Item1[1];
                double ddtheta = _biped.reducedSimulationData.RES[i].Item1[2];
                double time = _biped.reducedSimulationData.RES[i].Item2;
                zeroDynamics.Series["Phaseportrait"].Points.AddXY(theta, dtheta);
                //zeroDynamics.Series["dtheta"].Points.AddXY(time, dtheta);
                //zeroDynamics.Series["ddtheta"].Points.AddXY(time, ddtheta);
            }
            zeroDynamics.Series["Phaseportrait"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            zeroDynamics.Series["Phaseportrait"].Color = Color.Red;

            //zeroDynamics.Series["dtheta"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            //zeroDynamics.Series["dtheta"].Color = Color.Black;

            //zeroDynamics.Series["ddtheta"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            //zeroDynamics.Series["ddtheta"].Color = Color.Blue;
            zeroDynamics.SaveImage(@"../../../../pictures/phaseabg.png", System.Drawing.Imaging.ImageFormat.Png);
        }
    }
}
