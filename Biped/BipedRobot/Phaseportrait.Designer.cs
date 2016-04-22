namespace BipedRobot
{
    partial class Phaseportrait
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
            this.dynamicsChart = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.button1 = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.dynamicsChart)).BeginInit();
            this.SuspendLayout();
            // 
            // dynamicsChart
            // 
            chartArea1.Name = "ChartArea1";
            this.dynamicsChart.ChartAreas.Add(chartArea1);
            legend1.Name = "Legend1";
            this.dynamicsChart.Legends.Add(legend1);
            this.dynamicsChart.Location = new System.Drawing.Point(64, 23);
            this.dynamicsChart.Name = "dynamicsChart";
            series1.ChartArea = "ChartArea1";
            series1.Legend = "Legend1";
            series1.Name = "dtheta";
            this.dynamicsChart.Series.Add(series1);
            this.dynamicsChart.Size = new System.Drawing.Size(704, 488);
            this.dynamicsChart.TabIndex = 1;
            this.dynamicsChart.Text = "chart1";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(661, 526);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(107, 36);
            this.button1.TabIndex = 2;
            this.button1.Text = "Plot";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // Phaseportrait
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(874, 624);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.dynamicsChart);
            this.Name = "Phaseportrait";
            this.Text = "Phaseportrait";
            this.Load += new System.EventHandler(this.Phaseportrait_Load);
            ((System.ComponentModel.ISupportInitialize)(this.dynamicsChart)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart dynamicsChart;
        private System.Windows.Forms.Button button1;
    }
}