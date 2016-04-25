namespace BipedRobot
{
    partial class graph
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
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea2 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend2 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series4 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.dynamicsChart = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.button1 = new System.Windows.Forms.Button();
            this.button2 = new System.Windows.Forms.Button();
            this.zeroDynamics = new System.Windows.Forms.DataVisualization.Charting.Chart();
            ((System.ComponentModel.ISupportInitialize)(this.dynamicsChart)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.zeroDynamics)).BeginInit();
            this.SuspendLayout();
            // 
            // dynamicsChart
            // 
            chartArea1.Name = "ChartArea1";
            this.dynamicsChart.ChartAreas.Add(chartArea1);
            legend1.Name = "Legend1";
            this.dynamicsChart.Legends.Add(legend1);
            this.dynamicsChart.Location = new System.Drawing.Point(12, 12);
            this.dynamicsChart.Name = "dynamicsChart";
            series1.ChartArea = "ChartArea1";
            series1.Legend = "Legend1";
            series1.Name = "q1";
            series2.ChartArea = "ChartArea1";
            series2.Legend = "Legend1";
            series2.Name = "q2";
            series3.ChartArea = "ChartArea1";
            series3.Legend = "Legend1";
            series3.Name = "q3";
            this.dynamicsChart.Series.Add(series1);
            this.dynamicsChart.Series.Add(series2);
            this.dynamicsChart.Series.Add(series3);
            this.dynamicsChart.Size = new System.Drawing.Size(704, 488);
            this.dynamicsChart.TabIndex = 0;
            this.dynamicsChart.Text = "chart1";
            this.dynamicsChart.Click += new System.EventHandler(this.dynamicsChart_Click);
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(609, 506);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(107, 36);
            this.button1.TabIndex = 1;
            this.button1.Text = "Plot";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // button2
            // 
            this.button2.Location = new System.Drawing.Point(1165, 506);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(107, 36);
            this.button2.TabIndex = 2;
            this.button2.Text = "Plot";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.button2_Click);
            // 
            // zeroDynamics
            // 
            chartArea2.Name = "ChartArea1";
            this.zeroDynamics.ChartAreas.Add(chartArea2);
            legend2.Name = "Legend1";
            this.zeroDynamics.Legends.Add(legend2);
            this.zeroDynamics.Location = new System.Drawing.Point(722, 12);
            this.zeroDynamics.Name = "zeroDynamics";
            series4.ChartArea = "ChartArea1";
            series4.Legend = "Legend1";
            series4.Name = "phaseportrait";
            this.zeroDynamics.Series.Add(series4);
            this.zeroDynamics.Size = new System.Drawing.Size(563, 488);
            this.zeroDynamics.TabIndex = 3;
            this.zeroDynamics.Text = "chart1";
            // 
            // graph
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1297, 588);
            this.Controls.Add(this.zeroDynamics);
            this.Controls.Add(this.button2);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.dynamicsChart);
            this.Name = "graph";
            this.Text = "graph";
            this.Load += new System.EventHandler(this.graph_Load);
            ((System.ComponentModel.ISupportInitialize)(this.dynamicsChart)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.zeroDynamics)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart dynamicsChart;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.Button button2;
        private System.Windows.Forms.DataVisualization.Charting.Chart zeroDynamics;
    }
}