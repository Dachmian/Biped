using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Symbolics;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using System.IO;
using System.Globalization;


namespace BipedRobot
{

    //the performanceindex is to be integrated from theta0 to thetaT 
    public class AGS
    {
        private BRgait _gait;
        private Expression _performanceIndex;
        private Expression _impact;
        private List<Expression> _eqConstraints;
        private List<Expression> _ineqConstraints;
        private Expression _lagrangian;
        private Expression[] _performanceGradientArray;
        private Expression[,] _hessianMatrix;
        private List<Expression[]> _eqconstraintsGradientArray;
        private List<Expression[]> _ineqconstraintsGradientArray;

        //private double[] _eqconstraint1Gradient;
        //private double[] _eqconstraint2Gradient;
        //private double[] _eqconstraint3Gradient;
        //private double[] _gradient;
        //private double[,] _hessian;
        private double[] _parameterValues;

        //private double _performanceIndexVal;

        public delegate double func(Expression exp);

        public AGS(BRgait gait, int numOfParams)
        {
            _gait = gait;
            BRVHC vhc = gait.vhc;
            _eqConstraints = new List<Expression>();
            _ineqConstraints = new List<Expression>();
            _eqconstraintsGradientArray = new List<Expression[]>();
            _ineqconstraintsGradientArray = new List<Expression[]>();
            _parameterValues = new double[3*numOfParams];

            string[] pLHS = File.ReadAllLines(@"../../../pLHS.txt");
            string[] pRHS = File.ReadAllLines(@"../../../pRHS.txt");
            string[] posture = File.ReadAllLines(@"../../../posture.txt");
            Matrix<double> Q = Matrix<double>.Build.Dense(3, 3);
            Matrix<double> P = Matrix<double>.Build.Dense(3, 3);
            Matrix<double> pmatrix = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {0,0,1},
                {0,1,0},
                {1,0,0},
            
            });
            Dictionary<string, FloatingPoint> parameters = gait.vhc.parameters;
            parameters.Add("q1", Convert.ToDouble(posture[0], CultureInfo.InvariantCulture));
            parameters.Add("q2", Convert.ToDouble(posture[1], CultureInfo.InvariantCulture));
            parameters.Add("q3", Convert.ToDouble(posture[2], CultureInfo.InvariantCulture));
            for (int i = 0; i < 3; i++)
            {
                string temp = pLHS[i];
                temp = temp.Replace("(double)", "").Replace("0.2e1", "2");
                Expression exp = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
                double val = (double)MathNet.Symbolics.Evaluate.Evaluate(parameters, exp).RealValue;
                Q[0, i] = val;

                temp = pRHS[i];
                temp = temp.Replace("(double)", "").Replace("0.2e1", "2");
                exp = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
                val = (double)MathNet.Symbolics.Evaluate.Evaluate(parameters, exp).RealValue;
                P[0, i] = val;
            }

            for (int i = 0; i < 3; i++)
            {
                string temp = pLHS[i + 3];
                temp = temp.Replace("(double)", "").Replace("0.2e1", "2");
                Expression exp = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
                double val = (double)MathNet.Symbolics.Evaluate.Evaluate(parameters, exp).RealValue;
                Q[1, i] = val;

                temp = pRHS[i + 3];
                temp = temp.Replace("(double)", "").Replace("0.2e1", "2");
                exp = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
                val = (double)MathNet.Symbolics.Evaluate.Evaluate(parameters, exp).RealValue;
                P[1, i] = val;
            }

            for (int i = 0; i < 3; i++)
            {
                string temp = pLHS[i + 6];
                temp = temp.Replace("(double)", "").Replace("0.2e1","2");
                Expression exp = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
                double val = (double)MathNet.Symbolics.Evaluate.Evaluate(parameters, exp).RealValue;
                Q[2, i] = val;

                temp = pRHS[i + 6];
                temp = temp.Replace("(double)", "").Replace("0.2e1", "2");
                exp = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
                val = (double)MathNet.Symbolics.Evaluate.Evaluate(parameters, exp).RealValue;
                P[2, i] = val;
            }


            Matrix<double> A = pmatrix * Q.Solve(P);
    }


        public void setAnalyticalExpressions()
    {
            //int len = _gait.vhc.phi1Parameters.Count;
            //_performanceIndex = 0;
            //_lagrangian = 0;
            //_performanceGradientArray = new Expression[len * 3];
            //_hessianMatrix = new Expression[len * 3, len * 3];

            //_eqconstraint1GradientArray = new Expression[len * 3];
            //_eqconstraint2GradientArray = new Expression[len * 3];
            //_eqconstraint3GradientArray = new Expression[len * 3];
            //_eqconstraint4GradientArray = new Expression[len * 3];
            //_eqconstraint5GradientArray = new Expression[len * 3];

            //_ineqconstraint1GradientArray = new Expression[len * 3];
            //_ineqconstraint2GradientArray = new Expression[len * 3];
            //_ineqconstraint3GradientArray = new Expression[len * 3];
            //for (int i = 0; i < _performanceGradientArray.Length; i++)
            //{
            //    _performanceGradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _performanceIndex);

            //    _eqconstraint1GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint1);
            //    _eqconstraint2GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint2);
            //    _eqconstraint3GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint3);
            //    _eqconstraint4GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint4);
            //    _eqconstraint5GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint5);

            //    _ineqconstraint1GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _ineqConstraint1);
            //    _ineqconstraint2GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _ineqConstraint2);
            //    _ineqconstraint3GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _ineqConstraint3);

            //    Expression expTemp = Calculus.Differentiate("P" + i.ToString(), _lagrangian);
            //    for (int j = 0; j < _performanceGradientArray.Length; j++)
            //    {
            //        _hessianMatrix[i, j] = Calculus.Differentiate("P" + j.ToString(), expTemp);
            //    }
            //}
        }



        public void runAnalytical(BRgait gait)
        {

            double[] p0 = _parameterValues;
            double[] s = new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            double epsx = 0.0001;
            double radius = 0.1;
            double rho = 50.0;
            int maxits = 0;
            alglib.minnsstate state;
            alglib.minnsreport rep;
            double[] p1;

            alglib.minnscreate(p0.Length, p0, out state);
            alglib.minnssetalgoags(state, radius, rho);
            alglib.minnssetcond(state, epsx, maxits);
            alglib.minnssetscale(state, s);
            
            alglib.minnssetnlc(state, 5, 3);

            alglib.minnsoptimize(state, evaluateObjFuncAndConstraintsAnalytical, null, null);
            alglib.minnsresults(state, out p1, out rep);
            Console.WriteLine("{0}", alglib.ap.format(p1, 15));
            
            Console.ReadLine();

        }
        public void runNumerical(BRgait gait)
        {
            double[] p0 = _parameterValues;
            double[] s = new double[] { 1, 1};
            double epsx = 0.001;
            double radius = 0.1;
            double rho = 30.0;
            double diffstep = 0.001;
            int maxits = 0;
            alglib.minnsstate state;
            alglib.minnsreport rep;
            double[] p1;
            alglib.minnscreatef(p0.Length, p0, diffstep, out state);
            alglib.minnssetalgoags(state, radius, rho);
            alglib.minnssetcond(state, epsx, maxits);
            alglib.minnssetscale(state, s);


            alglib.minnssetnlc(state, 0, _ineqConstraints.Count);
            try
            {
                alglib.minnsoptimize(state, evaluateObjFuncAndConstraintsNumerical, null, null);
            }
            catch(Exception ex)
            {
                Console.WriteLine(ex);
            }
            alglib.minnsresults(state, out p1, out rep);
            Console.WriteLine("{0}", alglib.ap.format(p1, 15));
            Console.ReadLine();

            //testImpact(gait, p1);
        }
        
        
        

        public void evaluateObjFuncAndConstraintsAnalytical(double[] p, double[] objFunction, double[,] jacobian, object obj)
        {
            


            //objFunction[0] = evaluateFunction(_performanceIndex);
            //objFunction[1] = evaluateFunction(_eqConstraint1);
            //objFunction[2] = evaluateFunction(_eqConstraint2);
            //objFunction[3] = evaluateFunction(_eqConstraint3);
            //objFunction[4] = evaluateFunction(_eqConstraint4);
            //objFunction[5] = evaluateFunction(_eqConstraint5);
            //objFunction[6] = evaluateFunction(_ineqConstraint1);
            //objFunction[7] = evaluateFunction(_ineqConstraint2);
            //objFunction[8] = evaluateFunction(_ineqConstraint3);
        }
        public void evaluateObjFuncAndConstraintsNumerical(double[] p, double[] objFunction, object obj)
        {
            double dthetaMin = evaluateDthetaConstraint();
            double alphaMax = evaluateAlphaConstraint();
            for (int i = 0; i < p.Length; i++)
            {
                _gait.vhc.parameters["P" + i.ToString()] = p[i];

            }

            
            objFunction[0] = 0;
            objFunction[1] = 0;
            objFunction[2] =0;
            objFunction[3] = 0;
            objFunction[4] = 0;
            objFunction[5] = 0;
            objFunction[6] =0;
            objFunction[7] = 0;

        }



        public double evaluateDthetaConstraint()
        {
            double[][] firstIntegral = RiemannSum.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha);
            double[][] secondIntegral = RiemannSum.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1]);
            int len = firstIntegral[1].Length - 1;
            double dthetaTSquared = (- secondIntegral[1][len])/(1-Math.Exp(-firstIntegral[1][len])*Math.Pow(_gait.impactSecondLine(0,1), 2) );
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(0, 1), 2);
            double val = 0;
            if(double.IsNaN(dthetaTSquared))
            {
                val = (- Math.Pow(10, 300));
                val = (- Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                val = Math.Pow(10, 300);
                val = Math.Pow(10, 300);
            }
            else if(double.IsNaN(dtheta0Squared))
            {
                val = (-Math.Pow(10, 300));
                val = (-Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                val = Math.Pow(10, 300);
                val = Math.Pow(10, 300);
            }
            else
            {
                val = Math.Sqrt(dtheta0Squared);
                for (int i = 1; i < len; i++)
                {
                    val = Math.Sqrt(-secondIntegral[1][i] + Math.Exp(-firstIntegral[1][i]) * dtheta0Squared);
                }
            }
            return val;
        }

        public double evaluateAlphaConstraint()
        {
            

            double dx = 0.01;
            double theta = 0;
            double val = double.MinValue;
            double alphaVal = 0;
            for (int i = 0; i < 1/dx; i++)
            {
                alphaVal = _gait.vhc.evalAlpha(theta);
                if(alphaVal > val)
                {
                    val = alphaVal;
                }
                theta += dx;
            }
            return val;
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

            _ineqConstraint3 = vhc.impactFirstLine - vhc.impactSecondLine;
            _ineqConstraint4 = vhc.impactFirstLine - vhc.impactThirdLine;


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




