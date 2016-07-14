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
    public partial class TrajectoryGraph : Form
    {
        private BRgait _gait;
        private double[,] _THETA;

        public TrajectoryGraph(BRgait gait, double[,] THETA)
        {
            _gait = gait;
            _THETA = THETA;
            InitializeComponent();
            ShowDialog();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            for (int i = 0; i < _THETA.Length / _THETA.Rank; i++)
            {
                double q1 = _gait.vhc.evalPhi1(_THETA[0, i]);
                double q2 = _gait.vhc.evalPhi2(_THETA[0, i]);
                double q3 = _gait.vhc.evalPhi3(_THETA[0, i]);
                double theta = _THETA[0, i];
                trajectoryChart.Series["q1"].Points.AddXY(theta, q1);
                trajectoryChart.Series["q2"].Points.AddXY(theta, q2);
                trajectoryChart.Series["q3"].Points.AddXY(theta, q3);
            }
            trajectoryChart.Series["q1"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            trajectoryChart.Series["q1"].Color = Color.Red;
            trajectoryChart.Series["q2"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            trajectoryChart.Series["q2"].Color = Color.Blue;
            trajectoryChart.Series["q3"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            trajectoryChart.Series["q3"].Color = Color.Green;
            trajectoryChart.SaveImage(@"../../../../pictures/trajectory.png", System.Drawing.Imaging.ImageFormat.Png);
        }
    }
}
