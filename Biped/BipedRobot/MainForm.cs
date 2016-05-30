using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;

namespace BipedRobot
{
    public partial class MainForm : Form
    {
        private Biped _biped;
        public MainForm()
        {
            InitializeComponent();
            ShowDialog();
            _biped = null;
        }

        private void btnParameters_Click(object sender, EventArgs e)
        {
            Stream filestream = null;
            OpenFileDialog dlgParameters = new OpenFileDialog();
            dlgParameters.InitialDirectory = @"C:\\Users\damira\Documents\Visual Studio 2015\Projects\BipedRobot\Biped";
            dlgParameters.Filter = "XML files (*.xml)|*.xml";
            dlgParameters.RestoreDirectory = true;
            if (dlgParameters.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    if ((filestream = dlgParameters.OpenFile()) != null)
                    {
                        using (filestream)
                        {
                            txtParameters.Text = dlgParameters.FileName;
                        }
                    }
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Error: Could not read file from disk. Original error: " + ex.Message);
                }
            }
        }

        private void btnInitialize_Click(object sender, EventArgs e)
        {
            _biped = new Biped(@txtParameters.Text);
        }
    }
}
