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
    public partial class Phaseportrait : Form
    {
        private double[,] _THETA;
        public Phaseportrait(double[,] THETA)
        {
            _THETA = THETA;
            InitializeComponent();
            ShowDialog();
        }

        private void Phaseportrait_Load(object sender, EventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e)
        {
            for (int i = 0; i < _THETA.Length/_THETA.Rank; i++)
            {
                double theta = _THETA[i, 0];
                double dtheta = _THETA[i, 1];
                dynamicsChart.Series["dtheta"].Points.AddXY(theta, Math.Sqrt(dtheta));
            }
            dynamicsChart.Series["dtheta"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            dynamicsChart.Series["dtheta"].Color = Color.Red;
        }
    }
}
