namespace BipedRobot
{
    partial class MainForm
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
            this.txtParameters = new System.Windows.Forms.TextBox();
            this.btnParameters = new System.Windows.Forms.Button();
            this.dlgParameters = new System.Windows.Forms.OpenFileDialog();
            this.grpParameters = new System.Windows.Forms.GroupBox();
            this.btnInitialize = new System.Windows.Forms.Button();
            this.grpParameters.SuspendLayout();
            this.SuspendLayout();
            // 
            // txtParameters
            // 
            this.txtParameters.Location = new System.Drawing.Point(6, 19);
            this.txtParameters.Name = "txtParameters";
            this.txtParameters.Size = new System.Drawing.Size(184, 20);
            this.txtParameters.TabIndex = 0;
            // 
            // btnParameters
            // 
            this.btnParameters.Location = new System.Drawing.Point(196, 17);
            this.btnParameters.Name = "btnParameters";
            this.btnParameters.Size = new System.Drawing.Size(75, 23);
            this.btnParameters.TabIndex = 1;
            this.btnParameters.Text = "...";
            this.btnParameters.UseVisualStyleBackColor = true;
            this.btnParameters.Click += new System.EventHandler(this.btnParameters_Click);
            // 
            // dlgParameters
            // 
            this.dlgParameters.FileName = "openFileDialog1";
            // 
            // grpParameters
            // 
            this.grpParameters.Controls.Add(this.txtParameters);
            this.grpParameters.Controls.Add(this.btnParameters);
            this.grpParameters.Location = new System.Drawing.Point(12, 12);
            this.grpParameters.Name = "grpParameters";
            this.grpParameters.Size = new System.Drawing.Size(277, 64);
            this.grpParameters.TabIndex = 2;
            this.grpParameters.TabStop = false;
            this.grpParameters.Text = "Choose parameter set";
            // 
            // btnInitialize
            // 
            this.btnInitialize.Location = new System.Drawing.Point(12, 104);
            this.btnInitialize.Name = "btnInitialize";
            this.btnInitialize.Size = new System.Drawing.Size(277, 23);
            this.btnInitialize.TabIndex = 3;
            this.btnInitialize.Text = "Initialize Biped";
            this.btnInitialize.UseVisualStyleBackColor = true;
            this.btnInitialize.Click += new System.EventHandler(this.btnInitialize_Click);
            // 
            // MainForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1094, 537);
            this.Controls.Add(this.btnInitialize);
            this.Controls.Add(this.grpParameters);
            this.Name = "MainForm";
            this.Text = "MainForm";
            this.grpParameters.ResumeLayout(false);
            this.grpParameters.PerformLayout();
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.TextBox txtParameters;
        private System.Windows.Forms.Button btnParameters;
        private System.Windows.Forms.OpenFileDialog dlgParameters;
        private System.Windows.Forms.GroupBox grpParameters;
        private System.Windows.Forms.Button btnInitialize;
    }
}