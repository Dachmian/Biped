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
    public partial class TorquesGraph : Form
    {
        private BRTorques _torques;
        public TorquesGraph(BRTorques torques)
        {
            _torques = torques;
            InitializeComponent();
            ShowDialog();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            for (int i = 0; i < _torques.torque1.Length; i++)
            {
                double torque1 = _torques.torque1[i];
                double torque2 = _torques.torque2[i];
                double theta = _torques.theta[i];
                torquesChart.Series["Torque1"].Points.AddXY(theta, torque1);
                torquesChart.Series["Torque2"].Points.AddXY(theta, torque2);
            }
            torquesChart.Series["Torque1"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            torquesChart.Series["Torque1"].Color = Color.Red;
            torquesChart.Series["Torque2"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            torquesChart.Series["Torque2"].Color = Color.Blue;
            torquesChart.SaveImage(@"../../../../pictures/torques.png", System.Drawing.Imaging.ImageFormat.Png);
        }
    }
}
