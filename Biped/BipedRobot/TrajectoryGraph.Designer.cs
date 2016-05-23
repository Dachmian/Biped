namespace BipedRobot
{
    partial class TrajectoryGraph
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea1 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend1 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series1 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series2 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series3 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.trajectoryChart = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.button1 = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.trajectoryChart)).BeginInit();
            this.SuspendLayout();
            // 
            // trajectoryChart
            // 
            chartArea1.Name = "ChartArea1";
            this.trajectoryChart.ChartAreas.Add(chartArea1);
            legend1.Name = "Legend1";
            this.trajectoryChart.Legends.Add(legend1);
            this.trajectoryChart.Location = new System.Drawing.Point(12, 2);
            this.trajectoryChart.Name = "trajectoryChart";
            series1.ChartArea = "ChartArea1";
            series1.Legend = "Legend1";
            series1.Name = "q1";
            series2.ChartArea = "ChartArea1";
            series2.Legend = "Legend1";
            series2.Name = "q2";
            series3.ChartArea = "ChartArea1";
            series3.Legend = "Legend1";
            series3.Name = "q3";
            this.trajectoryChart.Series.Add(series1);
            this.trajectoryChart.Series.Add(series2);
            this.trajectoryChart.Series.Add(series3);
            this.trajectoryChart.Size = new System.Drawing.Size(704, 488);
            this.trajectoryChart.TabIndex = 3;
            this.trajectoryChart.Text = "chart1";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(609, 519);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(107, 36);
            this.button1.TabIndex = 4;
            this.button1.Text = "Plot";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // Trajectory
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(784, 567);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.trajectoryChart);
            this.Name = "Trajectory";
            this.Text = "Trajectory";
            ((System.ComponentModel.ISupportInitialize)(this.trajectoryChart)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart trajectoryChart;
        private System.Windows.Forms.Button button1;
    }
}