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
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea2 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend2 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series4 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series5 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series6 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.dynamicsChart = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.button1 = new System.Windows.Forms.Button();
            ((System.ComponentModel.ISupportInitialize)(this.dynamicsChart)).BeginInit();
            this.SuspendLayout();
            // 
            // dynamicsChart
            // 
            chartArea2.Name = "ChartArea1";
            this.dynamicsChart.ChartAreas.Add(chartArea2);
            legend2.Name = "Legend1";
            this.dynamicsChart.Legends.Add(legend2);
            this.dynamicsChart.Location = new System.Drawing.Point(12, 12);
            this.dynamicsChart.Name = "dynamicsChart";
            series4.ChartArea = "ChartArea1";
            series4.Legend = "Legend1";
            series4.Name = "q1";
            series5.ChartArea = "ChartArea1";
            series5.Legend = "Legend1";
            series5.Name = "q2";
            series6.ChartArea = "ChartArea1";
            series6.Legend = "Legend1";
            series6.Name = "q3";
            this.dynamicsChart.Series.Add(series4);
            this.dynamicsChart.Series.Add(series5);
            this.dynamicsChart.Series.Add(series6);
            this.dynamicsChart.Size = new System.Drawing.Size(704, 488);
            this.dynamicsChart.TabIndex = 0;
            this.dynamicsChart.Text = "chart1";
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
            // graph
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(728, 588);
            this.Controls.Add(this.button1);
            this.Controls.Add(this.dynamicsChart);
            this.Name = "graph";
            this.Text = "graph";
            this.Load += new System.EventHandler(this.graph_Load);
            ((System.ComponentModel.ISupportInitialize)(this.dynamicsChart)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.DataVisualization.Charting.Chart dynamicsChart;
        private System.Windows.Forms.Button button1;
    }
}