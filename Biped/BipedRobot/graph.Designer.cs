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
            this.button1 = new System.Windows.Forms.Button();
            this.zeroDynamics = new System.Windows.Forms.DataVisualization.Charting.Chart();
            ((System.ComponentModel.ISupportInitialize)(this.zeroDynamics)).BeginInit();
            this.SuspendLayout();
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(806, 667);
            this.button1.Margin = new System.Windows.Forms.Padding(4);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(143, 44);
            this.button1.TabIndex = 3;
            this.button1.Text = "Plot";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click_1);
            // 
            // zeroDynamics
            // 
            chartArea1.Name = "ChartArea1";
            this.zeroDynamics.ChartAreas.Add(chartArea1);
            legend1.Name = "Legend1";
            this.zeroDynamics.Legends.Add(legend1);
            this.zeroDynamics.Location = new System.Drawing.Point(22, 62);
            this.zeroDynamics.Margin = new System.Windows.Forms.Padding(4);
            this.zeroDynamics.Name = "zeroDynamics";
            series1.ChartArea = "ChartArea1";
            series1.Legend = "Legend1";
            series1.Name = "Phaseportrait";
            this.zeroDynamics.Series.Add(series1);
            this.zeroDynamics.Size = new System.Drawing.Size(939, 601);
            this.zeroDynamics.TabIndex = 4;
            this.zeroDynamics.Text = "chart1";
            // 
            // graph
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(8F, 16F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(983, 724);
            this.Controls.Add(this.zeroDynamics);
            this.Controls.Add(this.button1);
            this.Margin = new System.Windows.Forms.Padding(4, 4, 4, 4);
            this.Name = "graph";
            this.Text = "graph";
            this.Load += new System.EventHandler(this.graph_Load);
            ((System.ComponentModel.ISupportInitialize)(this.zeroDynamics)).EndInit();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.DataVisualization.Charting.Chart zeroDynamics;


    }
}