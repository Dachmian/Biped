namespace BipedRobot
{
    partial class ErrorGraph
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
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea2 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend2 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series2 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.errorChart = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.button1 = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.errorChart)).BeginInit();
            this.SuspendLayout();
            // 
            // errorChart
            // 
            chartArea2.Name = "ChartArea1";
            this.errorChart.ChartAreas.Add(chartArea2);
            legend2.Name = "Legend1";
            this.errorChart.Legends.Add(legend2);
            this.errorChart.Location = new System.Drawing.Point(13, 13);
            this.errorChart.Margin = new System.Windows.Forms.Padding(4);
            this.errorChart.Name = "errorChart";
            series2.ChartArea = "ChartArea1";
            series2.Legend = "Legend1";
            series2.Name = "Error";
            this.errorChart.Series.Add(series2);
            this.errorChart.Size = new System.Drawing.Size(939, 601);
            this.errorChart.TabIndex = 4;
            this.errorChart.Text = "chart1";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(809, 640);
            this.button1.Margin = new System.Windows.Forms.Padding(4);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(143, 44);
            this.button1.TabIndex = 5;
            this.button1.Text = "Plot";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // ErrorGraph
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1047, 697);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.errorChart);
            this.Name = "ErrorGraph";
            this.Text = "ErrorGraph";
            ((System.ComponentModel.ISupportInitialize)(this.errorChart)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart errorChart;
        private System.Windows.Forms.Button button1;
    }
}