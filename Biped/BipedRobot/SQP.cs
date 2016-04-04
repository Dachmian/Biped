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
        private Expression _ineqConstraint3;
        private Expression _ineqConstraint4;
        private Expression _lagrangian;
        private Expression[] _gradientArray;
        private Expression[,] _hessianMatrix;
        private Vector<double> _gradient;
        private Matrix<double> _hessian;
        private double _performanceIndexVal;
        private Dictionary<string, FloatingPoint> _parameters;

        public delegate double func(Expression exp);

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

            _ineqConstraint3 = vhc.impactNegFirstLine / vhc.impactPosFirstLine - vhc.impactNegSecondLine / vhc.impactPosSecondLine;
            _ineqConstraint4 = vhc.impactNegFirstLine / vhc.impactPosFirstLine - vhc.impactNegThirdLine / vhc.impactPosThirdLine;


            _lagrangian = _performanceIndex - "lambda1" * _ineqConstraint1 - "lambda2" * _ineqConstraint2 - "lambda3" * _ineqConstraint3 - "lambda4" * _ineqConstraint4;
            _gradientArray = new Expression[_parameters.Count - 3];
            _hessianMatrix = new Expression[_parameters.Count - 3, _parameters.Count - 3];
            for (int i = 0; i < _gradientArray.Length; i++)
            {
                _gradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _performanceIndex);
                for(int j = 0; j <_gradientArray.Length; j++)
                {
                    _hessianMatrix[i, j] = Calculus.Differentiate("P" + j.ToString(), _gradientArray[i]);
                }
            }
        }

        public void run(BRgait gait)
        {

            _performanceIndexVal = EvalPerformanceIndex(ref gait);


        }
        public double EvalPerformanceIndex(ref BRgait gait)
        {
            double[,] THETA = calculateReducedDynamics.run(gait);
            return 0;
        }
        public void evaluateGradient(Expression[] _gradientArray)
        {
            int len = _gradientArray.Length;
            _gradient = Vector<double>.Build.Dense(len);
            for(int i = 0; i < len; i++)
            {
                _gradient[i] = evaluateFunction(_gradientArray[i]);
            }
        }
        public void evaluateHessian(Expression[,] hessianMatrix)
        {
            int len = _gradientArray.Length;
            _hessian = Matrix<double>.Build.Dense(len, len);
            for (int i = 0; i < len; i++)
            {
                for (int j = 0; j < len; j++)
                {
                    _hessian[i, j] = evaluateFunction(hessianMatrix[i, j]);
                }
            }

        }      
        public double evaluateFunction(Expression exp)
        {
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, exp).RealValue;
        }

    }
}
