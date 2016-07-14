namespace BipedRobot
{
    partial class TrajectoryPortrait
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
            this.trajectoryPortraitChart = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.button1 = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.trajectoryPortraitChart)).BeginInit();
            this.SuspendLayout();
            // 
            // trajectoryPortraitChart
            // 
            chartArea1.Name = "ChartArea1";
            this.trajectoryPortraitChart.ChartAreas.Add(chartArea1);
            legend1.Name = "Legend1";
            this.trajectoryPortraitChart.Legends.Add(legend1);
            this.trajectoryPortraitChart.Location = new System.Drawing.Point(67, 3);
            this.trajectoryPortraitChart.Margin = new System.Windows.Forms.Padding(4);
            this.trajectoryPortraitChart.Name = "trajectoryPortraitChart";
            series1.ChartArea = "ChartArea1";
            series1.Legend = "Legend1";
            series1.Name = "q1";
            series2.ChartArea = "ChartArea1";
            series2.Legend = "Legend1";
            series2.Name = "q2";
            series3.ChartArea = "ChartArea1";
            series3.Legend = "Legend1";
            series3.Name = "q3";
            this.trajectoryPortraitChart.Series.Add(series1);
            this.trajectoryPortraitChart.Series.Add(series2);
            this.trajectoryPortraitChart.Series.Add(series3);
            this.trajectoryPortraitChart.Size = new System.Drawing.Size(939, 601);
            this.trajectoryPortraitChart.TabIndex = 4;
            this.trajectoryPortraitChart.Text = "chart1";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(900, 612);
            this.button1.Margin = new System.Windows.Forms.Padding(4);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(143, 44);
            this.button1.TabIndex = 5;
            this.button1.Text = "Plot";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // TrajectoryPortrait
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1080, 680);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.trajectoryPortraitChart);
            this.Name = "TrajectoryPortrait";
            this.Text = "TrajectoryPortrait";
            ((System.ComponentModel.ISupportInitialize)(this.trajectoryPortraitChart)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart trajectoryPortraitChart;
        private System.Windows.Forms.Button button1;
    }
}