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
    public class AGSImpact
    {
        private Expression _performanceIndex;
        private Expression _eqConstraint1;
        private Expression _eqConstraint2;
        private Expression _eqConstraint3;
        private Expression _lagrangian;
        private Expression[] _performanceGradientArray;
        private Expression[,] _hessianMatrix;
        private Expression[] _eqconstraint1GradientArray;
        private Expression[] _eqconstraint2GradientArray;
        private Expression[] _eqconstraint3GradientArray;

        private double[] _eqconstraint1Gradient;
        private double[] _eqconstraint2Gradient;
        private double[] _eqconstraint3Gradient;
        private double[] _gradient;
        private double[,] _hessian;
        private double _performanceIndexVal;
        private double[] _parameterValues;
        private Dictionary<string, FloatingPoint> _parameters;

        public delegate double func(Expression exp);

        public AGSImpact(BRVHC vhc)
        {
            //because of theta entry we take -1. this is used for the gradients and we need the length of all the minimizing parameters
            int len = vhc.phi1Parameters.Count - 1;
            _parameterValues = new double[3*len];

            Expression ddtheta = Expression.Symbol("ddtheta");
            Expression dthetaSquared = Expression.Symbol("dtheta^2");
            Expression theta = Expression.Symbol("theta");

            _parameters = new Dictionary<string, FloatingPoint>();
            _parameters.Add("g", vhc.parameters["g"]);
            _parameters.Add("m1", vhc.parameters["m1"]);
            _parameters.Add("m2", vhc.parameters["m2"]);
            _parameters.Add("m3", vhc.parameters["m3"]);

            _parameters.Add("l1", vhc.parameters["l1"]);
            _parameters.Add("l2", vhc.parameters["l2"]);
            _parameters.Add("l3", vhc.parameters["l3"]);

            _parameters.Add("L1", vhc.parameters["L1"]);
            _parameters.Add("L2", vhc.parameters["L2"]);
            _parameters.Add("L3", vhc.parameters["L3"]);

            _parameters.Add("J1", vhc.parameters["J1"]);
            _parameters.Add("J2", vhc.parameters["J2"]);
            _parameters.Add("J3", vhc.parameters["J3"]);
            _parameters.Add("theta", 0);
            _parameters.Add("dtheta", 0);
            _parameters.Add("ddtheta", 0);

            int k = 0;
            foreach (KeyValuePair<string, FloatingPoint> entry in vhc.phi1Parameters)
            {
                if (entry.Key != "theta")
                {
                    _parameterValues[k] = entry.Value.RealValue;
                    _parameters.Add(entry.Key, entry.Value);
                    k++;
                }
            }
            foreach (KeyValuePair<string, FloatingPoint> entry in vhc.phi2Parameters)
            {
                if (entry.Key != "theta")
                {
                    _parameterValues[k] = entry.Value.RealValue;
                    _parameters.Add(entry.Key, entry.Value);
                    k++;
                }
            }
            foreach (KeyValuePair<string, FloatingPoint> entry in vhc.phi3Parameters)
            {
                if (entry.Key != "theta")
                {
                    _parameterValues[k] = entry.Value.RealValue;
                    _parameters.Add(entry.Key, entry.Value);
                    k++;
                }
            }
            Expression impactPosFirstLine = vhc.impactPosFirstLine;
            impactPosFirstLine = Structure.Substitute("phi1", vhc.phi1, impactPosFirstLine);
            impactPosFirstLine = Structure.Substitute("phi3", vhc.phi3, impactPosFirstLine);
            impactPosFirstLine = Structure.Substitute("dphi3", vhc.dphi3, impactPosFirstLine);
            impactPosFirstLine = Structure.Substitute("theta", "0", impactPosFirstLine);

            Expression impactPosSecondLine = vhc.impactPosSecondLine;
            impactPosSecondLine = Structure.Substitute("phi2", vhc.phi2, impactPosSecondLine);
            impactPosSecondLine = Structure.Substitute("phi3", vhc.phi3, impactPosSecondLine);
            impactPosSecondLine = Structure.Substitute("dphi2", vhc.dphi2, impactPosSecondLine);
            impactPosSecondLine = Structure.Substitute("theta", "0", impactPosSecondLine);

            Expression impactPosThirdLine = vhc.impactPosThirdLine;
            impactPosThirdLine = Structure.Substitute("phi1", vhc.phi1, impactPosThirdLine);
            impactPosThirdLine = Structure.Substitute("phi2", vhc.phi2, impactPosThirdLine);
            impactPosThirdLine = Structure.Substitute("phi3", vhc.phi3, impactPosThirdLine);
            impactPosThirdLine = Structure.Substitute("dphi2", vhc.dphi2, impactPosThirdLine);
            impactPosThirdLine = Structure.Substitute("dphi3", vhc.dphi3, impactPosThirdLine);
            impactPosThirdLine = Structure.Substitute("theta", "0", impactPosThirdLine);

            Expression impactNegFirstLine = vhc.impactNegFirstLine;
            impactNegFirstLine = Structure.Substitute("theta", "1", impactNegFirstLine);

            Expression impactNegSecondLine = vhc.impactNegSecondLine;
            impactNegSecondLine = Structure.Substitute("phi1", vhc.phi1, impactNegSecondLine);
            impactNegSecondLine = Structure.Substitute("phi2", vhc.phi2, impactNegSecondLine);
            impactNegSecondLine = Structure.Substitute("dphi2", vhc.dphi2, impactNegSecondLine);
            impactNegSecondLine = Structure.Substitute("theta", "1", impactNegSecondLine);

            Expression impactNegThirdLine = vhc.impactNegThirdLine;
            impactNegThirdLine = Structure.Substitute("phi1", vhc.phi1, impactNegThirdLine);
            impactNegThirdLine = Structure.Substitute("phi2", vhc.phi2, impactNegThirdLine);
            impactNegThirdLine = Structure.Substitute("phi3", vhc.phi3, impactNegThirdLine);
            impactNegThirdLine = Structure.Substitute("dphi2", vhc.dphi2, impactNegThirdLine);
            impactNegThirdLine = Structure.Substitute("dphi3", vhc.dphi3, impactNegThirdLine);
            impactNegThirdLine = Structure.Substitute("theta", "1", impactNegThirdLine);

            _performanceIndex = Expression.Pow(impactNegFirstLine / impactPosFirstLine - impactNegSecondLine / impactPosSecondLine,"2") - Expression.Pow(impactNegFirstLine / impactPosFirstLine - impactNegThirdLine / impactPosThirdLine,"2");

            Expression phi1Start = vhc.phi1;
            phi1Start = Structure.Substitute("theta", "0", phi1Start);
            Expression phi1End = vhc.phi1;
            phi1End = Structure.Substitute("theta", "1", phi1End);

            Expression phi2Start = vhc.phi2;
            phi2Start = Structure.Substitute("theta", "0", phi2Start);
            Expression phi2End = vhc.phi2;
            phi2End = Structure.Substitute("theta", "1", phi2End);

            Expression phi3Start = vhc.phi3;
            phi3Start = Structure.Substitute("theta", "0", phi3Start);
            Expression phi3End = vhc.phi3;
            phi3End = Structure.Substitute("theta", "1", phi3End);

            _eqConstraint1 = phi1End;
            _eqConstraint2 = phi2End;
            _eqConstraint3 = phi3End;

            _lagrangian = _performanceIndex - "delta1" * _eqConstraint1 - "delta2" * _eqConstraint2 - "delta3" * _eqConstraint3;

            _performanceGradientArray = new Expression[len*3];
            _hessianMatrix = new Expression[len * 3, len * 3];

            _eqconstraint1GradientArray = new Expression[len * 3];
            _eqconstraint2GradientArray = new Expression[len * 3];
            _eqconstraint3GradientArray = new Expression[len * 3];
            for (int i = 0; i < _performanceGradientArray.Length; i++)
            {
                _performanceGradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _performanceIndex);

                _eqconstraint1GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint1);
                _eqconstraint2GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint2);
                _eqconstraint3GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint3);

                Expression expTemp = Calculus.Differentiate("P" + i.ToString(), _lagrangian);
                for (int j = 0; j < _performanceGradientArray.Length; j++)
                {
                    _hessianMatrix[i, j] = Calculus.Differentiate("P" + j.ToString(), expTemp);
                }
            }




        }

        public void run(BRgait gait)
        {
            string a = Infix.Format(_performanceIndex);
            string b = Infix.Format(_eqConstraint1);
            string c = Infix.Format(_eqConstraint2);
            string d = Infix.Format(_eqConstraint3);

            double[] p0 = _parameterValues;
            double[] s = new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            double epsx = 0.00001;
            double radius = 0.1;
            double rho = 50.0;
            int maxits = 0;
            alglib.minnsstate state;
            alglib.minnsreport rep;
            double[] p1;

            alglib.minnscreate(15, p0, out state);
            alglib.minnssetalgoags(state, radius, rho);
            alglib.minnssetcond(state, epsx, maxits);
            alglib.minnssetscale(state, s);

            alglib.minnssetnlc(state, 3, 0);

            alglib.minnsoptimize(state, evaluateObjFuncAndConstraints, null, null);
            alglib.minnsresults(state, out p1, out rep);
            Console.WriteLine("{0}", alglib.ap.format(p1, 15));
            Console.ReadLine();

        }
        
        public void evaluateGradientsAndHessian()
        {
            int len = _performanceGradientArray.Length;
            _hessian = new double[len, len];
            for (int i = 0; i < len; i++)
            {
                for (int j = 0; j < len; j++)
                {
                    _hessian[i, j] = evaluateFunction(_hessianMatrix[i, j]);
                }
            }
        }
        public double evaluateFunction(Expression exp)
        {
            string a = Infix.Format(exp);
            exp = Infix.ParseOrUndefined(a);
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, exp).RealValue;
        }

        public void evaluateObjFuncAndConstraints(double[] p, double[] objFunction, double[,] jacobian, object obj)
        {
            int len = p.Length;
            for (int i = 0; i < len; i++)
            {
                
                jacobian[0, i] = evaluateFunction(_performanceGradientArray[i]);

                jacobian[1, i] = evaluateFunction(_eqconstraint1GradientArray[i]);
                jacobian[2, i] = evaluateFunction(_eqconstraint2GradientArray[i]);
                jacobian[3, i] = evaluateFunction(_eqconstraint3GradientArray[i]);
            }
            objFunction[0] = evaluateFunction(_performanceIndex);
            objFunction[1] = evaluateFunction(_eqConstraint1);
            objFunction[2] = evaluateFunction(_eqConstraint2);
            objFunction[3] = evaluateFunction(_eqConstraint3);
        }

        

    }


    public class SQPPower
    {
        private Expression _performanceIndex;
        private Expression _ineqConstraint1;
        private Expression _ineqConstraint2;
        private Expression _ineqConstraint3;
        private Expression _ineqConstraint4;
        private Expression _lagrangian;
        private Expression[] _performanceGradientArray;
        private Expression[,] _hessianMatrix;
        private Expression[] _constraint1GradientArray;
        private Expression[] _constraint2GradientArray;
        private Expression[] _constraint3GradientArray;
        private Expression[] _constraint4GradientArray;

        private Vector<double> _constraint1Gradient;
        private Vector<double> _constraint2Gradient;
        private Vector<double> _constraint3Gradient;
        private Vector<double> _constraint4Gradient;
        private Vector<double> _gradient;
        private Matrix<double> _hessian;
        private Vector<double> _direction;
        private double _performanceIndexVal;
        private Dictionary<string, FloatingPoint> _parameters;

        public delegate double func(Expression exp);

        public SQPPower(BRVHC vhc)
        {
            Expression ddtheta = Expression.Symbol("ddtheta");
            Expression dthetaSquared = Expression.Symbol("dtheta^2");
            Expression theta = Expression.Symbol("theta");

            _parameters = new Dictionary<string, FloatingPoint>();
            _parameters.Add("g", vhc.parameters["g"]);
            _parameters.Add("m1", vhc.parameters["m1"]);
            _parameters.Add("m2", vhc.parameters["m2"]);
            _parameters.Add("m3", vhc.parameters["m3"]);

            _parameters.Add("l1", vhc.parameters["l1"]);
            _parameters.Add("l2", vhc.parameters["l2"]);
            _parameters.Add("l3", vhc.parameters["l3"]);

            _parameters.Add("L1", vhc.parameters["L1"]);
            _parameters.Add("L2", vhc.parameters["L2"]);
            _parameters.Add("L3", vhc.parameters["L3"]);

            _parameters.Add("J1", vhc.parameters["J1"]);
            _parameters.Add("J2", vhc.parameters["J2"]);
            _parameters.Add("J3", vhc.parameters["J3"]);
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

            int len = vhc.phi1Parameters.Count;
            _performanceGradientArray = new Expression[len * 3];
            _hessianMatrix = new Expression[len * 3, len * 3];
            _constraint1GradientArray = new Expression[len * 3];
            _constraint2GradientArray = new Expression[len * 3];
            _constraint3GradientArray = new Expression[len * 3];
            for (int i = 0; i < _performanceGradientArray.Length; i++)
            {
                _performanceGradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _performanceIndex);
                _constraint1GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _ineqConstraint1);
                _constraint2GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _ineqConstraint2);
                _constraint3GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _ineqConstraint3);
                _constraint4GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _ineqConstraint4);
                Expression expTemp = Calculus.Differentiate("P" + i.ToString(), _lagrangian);
                for (int j = 0; j < _performanceGradientArray.Length; j++)
                {
                    _hessianMatrix[i, j] = Calculus.Differentiate("P" + j.ToString(), expTemp);
                }
            }
        }
    }
}




