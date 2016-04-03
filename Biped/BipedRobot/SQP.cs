using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Symbolics;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using System.IO;

namespace BipedRobot
{

    //the performanceindex is to be integrated from theta0 to thetaT 
    public class SQP
    {
        private Expression _performanceIndex;
        private Expression _ineqConstraint1;
        private Expression _ineqConstraint2;
        private Expression _eqConstraint1;
        private Expression _eqConstraint2;

        private double _performanceIndexVal;

        private Dictionary<string, FloatingPoint> _parameters;

        public SQP(BRVHC vhc)
        {

            Expression ddtheta = Expression.Symbol("ddtheta");
            Expression dthetaSquared = Expression.Symbol("dtheta^2");
            Expression theta = Expression.Symbol("theta");

            _parameters = new Dictionary<string, FloatingPoint>();
            _parameters.Add("theta", 0);
            _parameters.Add("dtheta", 0);
            _parameters.Add("ddtheta", 0);

            foreach (KeyValuePair<string, FloatingPoint> entry in vhc.phi1Parameters)
            {
                _parameters.Add(entry.Key, entry.Value);
            }
            foreach (KeyValuePair<string, FloatingPoint> entry in vhc.phi2Parameters)
            {
                _parameters.Add(entry.Key, entry.Value);
            }
            foreach (KeyValuePair<string, FloatingPoint> entry in vhc.phi3Parameters)
            {
                _parameters.Add(entry.Key, entry.Value);
            }

            StreamReader fs = null;
            fs = new StreamReader(@"../../../alpha1.txt");
            string temp = fs.ReadLine();
            Expression alpha1 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            alpha1 = Structure.Substitute("dphi1", vhc.dphi1, alpha1);
            alpha1 = Structure.Substitute("dphi2", vhc.dphi2, alpha1);
            alpha1 = Structure.Substitute("dphi3", vhc.dphi3, alpha1);
            fs.Close();

            fs = new StreamReader(@"../../../beta1.txt");
            temp = fs.ReadLine();
            Expression beta1 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            beta1 = Structure.Substitute("dphi2", vhc.dphi2, beta1);
            beta1 = Structure.Substitute("dphi3", vhc.dphi3, beta1);
            beta1 = Structure.Substitute("ddphi1", vhc.ddphi1, beta1);
            beta1 = Structure.Substitute("ddphi2", vhc.ddphi2, beta1);
            fs.Close();

            fs = new StreamReader(@"../../../gamma1.txt");
            temp = fs.ReadLine();
            Expression gamma1 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            _ineqConstraint1 = alpha1 * ddtheta + beta1 * dthetaSquared + gamma1 - 150;


            fs = new StreamReader(@"../../../alpha3.txt");
            temp = fs.ReadLine();
            Expression alpha3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            alpha3 = Structure.Substitute("dphi1", vhc.dphi1, alpha3);
            alpha3 = Structure.Substitute("dphi3", vhc.dphi3, alpha3);
            fs.Close();

            fs = new StreamReader(@"../../../beta3.txt");
            temp = fs.ReadLine();
            Expression beta3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            beta3 = Structure.Substitute("dphi1", vhc.dphi1, beta3);
            beta3 = Structure.Substitute("dphi3", vhc.dphi3, beta3);
            beta3 = Structure.Substitute("ddphi1", vhc.ddphi1, beta3);
            fs.Close();

            fs = new StreamReader(@"../../../gamma3.txt");
            temp = fs.ReadLine();
            Expression gamma3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            _ineqConstraint2 = alpha3 * ddtheta + beta3 * dthetaSquared + gamma3 - 150;

            _performanceIndex = Expression.Abs((alpha1 * ddtheta + beta1 * dthetaSquared + gamma1) * vhc.dphi1 + (alpha3 * ddtheta + beta3 * dthetaSquared + gamma3) * vhc.dphi3);

            _eqConstraint1 = vhc.impactNegFirstLine / vhc.impactPosFirstLine - vhc.impactNegSecondLine / vhc.impactPosSecondLine;
            _eqConstraint2 = vhc.impactNegFirstLine / vhc.impactPosFirstLine - vhc.impactNegThirdLine / vhc.impactPosThirdLine;
        }


        public double EvalPerformanceIndex(ref BRgait gait)
        {
            gait.data = integrationReducedDynamics.run(gait.vhc,
                Vector<double>.Build.Dense(new[] { gait.gaitParam.theta0, gait.gaitParam.dtheta0 }),
                Vector<double>.Build.Dense(new[] { gait.gaitParam.thetaT, gait.gaitParam.dthetaT }));
            double sum = 0;
            for (int i = 0; i < gait.data.RES.Count; i++)
            {

            }
            return 0;
        }

        public void run(BRgait gait)
        {

            _performanceIndexVal = EvalPerformanceIndex(ref gait);


        }


    }
}
