﻿using System;
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
    public partial class TrajectoryPortrait : Form
    {
        private BRgait _gait;
        private double[,] _THETA;

        public TrajectoryPortrait(BRgait gait, double[,] THETA)
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
                double dq1 = _gait.vhc.evalDphi1(_THETA[0, i]) * _THETA[1,i];
                double dq2 = _gait.vhc.evalDphi2(_THETA[0, i]) * _THETA[1, i];
                double dq3 = _gait.vhc.evalDphi3(_THETA[0, i]) * _THETA[1, i];
                double theta = _THETA[0, i];
                trajectoryPortraitChart.Series["q1"].Points.AddXY(q1, dq1);
                trajectoryPortraitChart.Series["q2"].Points.AddXY(q2, dq2);
                trajectoryPortraitChart.Series["q3"].Points.AddXY(q3, dq3);
            }
            trajectoryPortraitChart.Series["q1"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            trajectoryPortraitChart.Series["q1"].Color = Color.Red;
            trajectoryPortraitChart.Series["q2"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            trajectoryPortraitChart.Series["q2"].Color = Color.Blue;
            trajectoryPortraitChart.Series["q3"].ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.FastLine;
            trajectoryPortraitChart.Series["q3"].Color = Color.Green;
            trajectoryPortraitChart.SaveImage(@"../../../../pictures/trajectoryPortrait.png", System.Drawing.Imaging.ImageFormat.Png);
        }
    }
}
