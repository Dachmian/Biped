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

    public partial class ErrorGraph : Form
    {
        private BRgait _gait;
        private double[,] _THETA;
        private BRReducedSimulationData _data;

        public ErrorGraph(BRgait gait, double[,] THETA, BRReducedSimulationData data)
        {
            _gait = gait;
            _THETA = THETA;
            _data = data;
            InitializeComponent();
            ShowDialog();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            double error = 0;
            for (int i = 0; i < _THETA.Length / _THETA.Rank; i++)
            {
                double theta = _THETA[0, i];
                for (int j = 0; j < _data.RES.Count; j++)
                {
                    if (_data.RES[j].Item1[0] > theta)
                    {
                        error = _data.RES[j-1].Item1[1] - _THETA[1, i];
                        break;
                    }
                }
                errorChart.Series["Error"].Points.AddXY(theta, error);
            }
            errorChart.Series["Error"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            errorChart.Series["Error"].Color = Color.Blue;
            errorChart.SaveImage(@"../../../../pictures/phaseError.png", System.Drawing.Imaging.ImageFormat.Png);
        }
    }
}
