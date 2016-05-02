namespace BipedRobot
{
    partial class TorquesGraph
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
            this.torquesChart = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.button1 = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.torquesChart)).BeginInit();
            this.SuspendLayout();
            // 
            // torquesChart
            // 
            chartArea1.Name = "ChartArea1";
            this.torquesChart.ChartAreas.Add(chartArea1);
            legend1.Name = "Legend1";
            this.torquesChart.Legends.Add(legend1);
            this.torquesChart.Location = new System.Drawing.Point(13, 13);
            this.torquesChart.Margin = new System.Windows.Forms.Padding(4);
            this.torquesChart.Name = "torquesChart";
            series1.ChartArea = "ChartArea1";
            series1.Legend = "Legend1";
            series1.Name = "Torque1";
            series2.ChartArea = "ChartArea1";
            series2.Legend = "Legend1";
            series2.Name = "Torque2";
            this.torquesChart.Series.Add(series1);
            this.torquesChart.Series.Add(series2);
            this.torquesChart.Size = new System.Drawing.Size(939, 601);
            this.torquesChart.TabIndex = 2;
            this.torquesChart.Text = "chart1";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(809, 622);
            this.button1.Margin = new System.Windows.Forms.Padding(4);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(143, 44);
            this.button1.TabIndex = 3;
            this.button1.Text = "Plot";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // TorquesGraph
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1079, 716);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.torquesChart);
            this.Name = "TorquesGraph";
            this.Text = "TorquesGraph";
            ((System.ComponentModel.ISupportInitialize)(this.torquesChart)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart torquesChart;
        private System.Windows.Forms.Button button1;
    }
}