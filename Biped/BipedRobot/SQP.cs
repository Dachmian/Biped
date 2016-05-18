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
using NagLibrary;


namespace BipedRobot
{

    //the performanceindex is to be integrated from theta0 to thetaT 
    public class SQPThirdOrder
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

        double _q1End;
        double _q2End;
        double _q3End;

        //private double[] _eqconstraint1Gradient;
        //private double[] _eqconstraint2Gradient;
        //private double[] _eqconstraint3Gradient;
        //private double[] _gradient;
        //private double[,] _hessian;
        private double[] _parameterValues;

        //private double _performanceIndexVal;

        public delegate double func(Expression exp);

        public SQPThirdOrder(BRgait gait, int numOfParams)
        {
            _gait = gait;
            BRVHC vhc = gait.vhc;
            _eqConstraints = new List<Expression>();
            _ineqConstraints = new List<Expression>();
            _eqconstraintsGradientArray = new List<Expression[]>();
            _ineqconstraintsGradientArray = new List<Expression[]>();
            _parameterValues = new double[3*numOfParams - 3];
            string[] parametersFromFile = File.ReadAllLines(@"../../../../parameters.txt");
            for (int i = 0; i < numOfParams - 1; i++)
            {
                _parameterValues[i] = Convert.ToDouble(parametersFromFile[i + 1], CultureInfo.InvariantCulture);
            }
            for (int i = numOfParams - 1; i < 2 * numOfParams- 2; i++)
            {
                _parameterValues[i] = Convert.ToDouble(parametersFromFile[i + 2], CultureInfo.InvariantCulture);

            }
            for (int i = 2 * numOfParams - 2; i < 3 * numOfParams - 3; i++)
            {
                _parameterValues[i] = Convert.ToDouble(parametersFromFile[i + 3], CultureInfo.InvariantCulture);

            }
            string[] pLHS = File.ReadAllLines(@"../../../../pLHS.txt");
            string[] pRHS = File.ReadAllLines(@"../../../../pRHS.txt");
            string[] posture = File.ReadAllLines(@"../../../../posture.txt");
            Matrix<double> Q = Matrix<double>.Build.Dense(3, 3);
            Matrix<double> P = Matrix<double>.Build.Dense(3, 3);
            Matrix<double> pmatrix = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {0,0,1},
                {0,1,0},
                {1,0,0},
            
            });
            Dictionary<string, FloatingPoint> parameters = gait.vhc.parameters;
            _q1End = -Convert.ToDouble(posture[0], CultureInfo.InvariantCulture);
            _q2End = Convert.ToDouble(posture[1], CultureInfo.InvariantCulture);
            _q3End = -Convert.ToDouble(posture[2], CultureInfo.InvariantCulture);
            parameters.Add("q1", _q1End);
            parameters.Add("q2", _q2End);
            parameters.Add("q3", _q3End);
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




        public double[] runNumericalSQP()
        {
            E04.E04WD_CONFUN confunE04WD = new E04.E04WD_CONFUN(confun);
            E04.E04WD_OBJFUN objfunE04WD = new E04.E04WD_OBJFUN(objfun);
            try
            {
                double objf = 0.0;
                int ifail, majits, n, nclin, ncnln;
                n = 9;
                nclin = 3;
                ncnln = 4;
                double[,] a = new double[,]{
                    {1, 1, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 1, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 1, 1},
                };
                double[] bl = new double[] { -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, 2*_q1End, 0, 2*_q3End, 0, 0, 0.2, -1.0E+25 };
                double[] bu = new double[] { 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 2*_q1End, 0, 2*_q3End, 0, 0, 50, -1 };
                double[] ccon = new double[ncnln];
                double[,] cjac = new double[ncnln, n];
                double[] clamda = new double[nclin + ncnln + n];
                double[] grad = new double[n];
                double[,] h = new double[n, n];
                double[] p = _parameterValues;
                int[] istate = new int[nclin + ncnln + n];

                E04.e04wdOptions options = new E04.e04wdOptions();
                options.Set("Print file = 6");
                options.Set("Derivative level = 0");
                //options.Set("Major Step limit = 0.00000000001");
                options.Set("Major Feasibility Tolerance = 1.0E-6");
                //options.Set("Function precision = 1.0E-5");

                E04.e04wd(n, nclin, ncnln, a, bl, bu, confunE04WD, objfunE04WD, out majits, istate,
                    ccon, cjac, clamda, out objf, grad, h, p, options, out ifail);
                Console.WriteLine("");
                Console.WriteLine("  On exit from e04wd, ifail = {0,5}", ifail);
                if (ifail == 0)
                {
                    Console.WriteLine("  Final objective value = {0,11:f3}", objf);
                    Console.WriteLine("");
                    Console.Write("  Output x");
                    for (int i = 1; i <= n; i++)
                    {
                        Console.Write(" " + " {0,11:f3}", p[i - 1]);
                    }
                    Console.WriteLine(" ");
                }
                return p;
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                Console.WriteLine("Exception Raised");
                return null;
            }

        }
        public void objfun(ref int mode, int n, double[] p, ref double objf,
            double[] grad, int nstate)
        {
            int len = p.Length / 3;


            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + (i + 1).ToString()] = p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + (i + 2).ToString()] = p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + (i + 3).ToString()] = p[i];

            }
            double[] values = evalTorques();
            objf = Math.Pow(values[2], 2) + Math.Pow(values[3], 2) +Math.Pow(values[0], 2) + Math.Pow(values[1], 2);
            mode = 0;
           
        }
        public void confun(ref int mode, int ncnln, int n, int[] needc, double[] p,
            double[] ccon, double[,] cjac, int nstate)
        {
            int len = p.Length / 3;
            double val1 = _gait.vhc.phi1Parameters["P0"].RealValue;
            double val2 = _gait.vhc.phi2Parameters["P" + (len + 1).ToString()].RealValue;
            double val3 = _gait.vhc.phi3Parameters["P" + (2 * len + 2).ToString()].RealValue;
            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + (i + 1).ToString()] = p[i];
                val1 += p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + (i + 2).ToString()] = p[i];
                val2 += p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + (i + 3).ToString()] = p[i];
                val3 += p[i];

            }
            double alphaMax = evaluateAlphaConstraint();
            double dthetaMin = evalDtheta();
            ccon[0] = _gait.impactFirstLine(0, 1) - _gait.impactSecondLine(0, 1);
            ccon[1] = _gait.impactSecondLine(0, 1) - _gait.impactThirdLine(0, 1);
            ccon[2] = dthetaMin;
            ccon[3] = alphaMax;
            mode = 0;

        }
        public void evaluateObjFuncAndConstraintsNumerical(double[] p, double[] objFunction, object obj)
        {
            double[] values = new double[5];
            double alphaMax = evaluateAlphaConstraint();
            if (alphaMax >= 0)
            {
                values[0] = 1;
            }
            else
            {
                values = evaluateDthetaAndTorqueConstraint();
            }
            int len = p.Length / 3;

            //add values to calculate endposture constraint
            double val1 = _gait.vhc.phi1Parameters["P0"].RealValue;
            double val2 = _gait.vhc.phi2Parameters["P" + (len + 1).ToString()].RealValue;
            double val3 = _gait.vhc.phi3Parameters["P" + (2*len + 2).ToString()].RealValue;
            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + (i + 1).ToString()] = p[i];
                val1 += p[i];

            }
            for (int i = len; i < 2*len; i++)
            {
                _gait.vhc.phi2Parameters["P" + (i + 2).ToString()] = p[i];
                val2 += p[i];

            }
            for (int i = 2*len; i < 3*len; i++)
            {
                _gait.vhc.phi3Parameters["P" + (i + 3).ToString()] = p[i];
                val3 += p[i];

            }

            if (values[0] == 1)
            {
                objFunction[0] = Math.Pow(10, 300);
                objFunction[1] = 0;
                objFunction[2] = 0;
                objFunction[3] = 0;
                objFunction[4] = 0;
                objFunction[5] = 0;
                objFunction[6] = 0;
                objFunction[7] = 0;
            }
            else
            {
                objFunction[0] = Math.Pow(values[2], 2) + Math.Pow(values[3], 2) + Math.Pow(values[4], 2) + Math.Pow(values[5], 2);
                objFunction[1] = val1 - _q1End;
                objFunction[2] = val2 - _q2End;
                objFunction[3] = val3 - _q3End;
                objFunction[4] = _gait.impactFirstLine(0, 1) - _gait.impactSecondLine(0, 1);
                objFunction[5] = _gait.impactFirstLine(0, 1) - _gait.impactThirdLine(0, 1);
                objFunction[6] = alphaMax;
                objFunction[7] = -values[1];
            }
        }



        public double[] evaluateDthetaAndTorqueConstraint()
        {
            double[][] firstIntegral = RiemannSum.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha);
            double[][] secondIntegral = RiemannSum.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1]);
            int len = firstIntegral[1].Length - 1;
            double dthetaTSquared = (- secondIntegral[1][len])/(1-Math.Exp(-firstIntegral[1][len])*Math.Pow(_gait.impactSecondLine(0,1), 2) );
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(0, 1), 2);
            double dthetaMin = 0;
            double torque1Max = 0;
            double torque2Max = 0;
            double torque1Diff = 0;
            double torque2Diff = 0;
            double infeasible = 0;
            double[] values;
            if (double.IsNaN(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                infeasible = 1;
                values = new double[] { infeasible, dthetaMin, torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                infeasible = 1;
                values = new double[] { infeasible, dthetaMin, torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if(double.IsNaN(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                infeasible = 1;
                values = new double[] { infeasible, dthetaMin, torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                infeasible = 1;
                values = new double[] { infeasible, dthetaMin, torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else
            {
                double[,] THETA = new double[2, len];
                THETA[0, 0] = 0;
                THETA[1, 0] = Math.Sqrt(dtheta0Squared);
                dthetaMin = Math.Sqrt(dtheta0Squared);

                double torque1 = 0;
                double torque2 = 0;
                double ddtheta = 0;
                for (int i = 1; i < len; i++)
                {
                    THETA[0, i] = secondIntegral[0][i];
                    THETA[1, i] = Math.Sqrt(-secondIntegral[1][i] + Math.Exp(-firstIntegral[1][i]) * dtheta0Squared);
                    if (THETA[1, i] < dthetaMin)
                    {
                        dthetaMin = THETA[1, i];
                    }
                }
                for (int i = 0; i < len; i += 5)
                {
                    ddtheta = (-_gait.vhc.evalBeta(THETA[0, i]) * Math.Pow(THETA[1, i], 2) - _gait.vhc.evalGamma(THETA[0, i])) / _gait.vhc.evalAlpha(THETA[0, i]);
                    torque1 = Math.Abs(_gait.vhc.evalAlpha1(THETA[0, i]) * ddtheta + _gait.vhc.evalBeta1(THETA[0, i]) * Math.Pow(THETA[1, i], 2) + _gait.vhc.evalGamma1(THETA[0, i]));
                    torque2 = Math.Abs(_gait.vhc.evalAlpha3(THETA[0, i]) * ddtheta + _gait.vhc.evalBeta3(THETA[0, i]) * Math.Pow(THETA[1, i], 2) + _gait.vhc.evalGamma3(THETA[0, i]));
                    if (torque1 > torque1Max)
                    {
                        torque1Max = torque1;
                    }
                    if (torque2 > torque2Max)
                    {
                        torque2Max = torque2;
                    }
                }
                double ddtheta0 = (-_gait.vhc.evalBeta(THETA[0, 0]) * Math.Pow(THETA[1, 0], 2) - _gait.vhc.evalGamma(THETA[0, 0])) / _gait.vhc.evalAlpha(THETA[0, 0]);

                torque1Diff = _gait.vhc.evalAlpha1(THETA[0, 0]) * ddtheta0 + _gait.vhc.evalBeta1(THETA[0, 0]) * Math.Pow(THETA[1, 0], 2) + _gait.vhc.evalGamma1(THETA[0, 0]) -
                    (_gait.vhc.evalAlpha3(THETA[0, len - 1]) * ddtheta + _gait.vhc.evalBeta3(THETA[0, len - 1]) * Math.Pow(THETA[1, len - 1], 2) + _gait.vhc.evalGamma3(THETA[0, len - 1]));

                torque2Diff = _gait.vhc.evalAlpha3(THETA[0, 0]) * ddtheta0 + _gait.vhc.evalBeta3(THETA[0, 0]) * Math.Pow(THETA[1, 0], 2) + _gait.vhc.evalGamma3(THETA[0, 0]) -
                    (_gait.vhc.evalAlpha1(THETA[0, len - 1]) * ddtheta + _gait.vhc.evalBeta1(THETA[0, len - 1]) * Math.Pow(THETA[1, len - 1], 2) + _gait.vhc.evalGamma1(THETA[0, len - 1]));
                values = new double[] { infeasible, dthetaMin, torque1Max, torque2Max, torque1Diff, torque2Diff };

            }
            return values;
        }

        public double evalDtheta()
        {
            double[][] firstIntegral = RiemannSum.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha);
            double[][] secondIntegral = RiemannSum.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1]);
            int len = firstIntegral[1].Length;
            double dthetaTSquared = (-secondIntegral[1][len - 1]) / (1 - Math.Exp(-firstIntegral[1][len - 1]) * Math.Pow(_gait.impactSecondLine(0, 1), 2));
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(0, 1), 2);
            double dthetaMin = 0;
            if (double.IsNaN(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsNaN(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else
            {
                double[,] THETA = new double[2, len + 1];
                THETA[0, 0] = 0;
                THETA[1, 0] = dtheta0Squared;
                dthetaMin = dtheta0Squared;
                for (int i = 1; i < len; i++)
                {
                    THETA[0, i] = secondIntegral[0][i];
                    THETA[1, i] = -secondIntegral[1][i] + Math.Exp(-firstIntegral[1][i]) * dtheta0Squared;
                    if (THETA[1, i] < dthetaMin)
                    {
                        dthetaMin = THETA[1, i];
                    }
                }
            }
            return dthetaMin;
        }
        public double[] evalTorques()
        {
            double[][] firstIntegral = RiemannSum.calculateFirstIntegral(_gait.firstIntegral);
            double[][] secondIntegral = RiemannSum.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1]);
            int len = firstIntegral[1].Length;
            double dthetaTSquared = (- secondIntegral[1][len - 1])/(1-Math.Exp(-firstIntegral[1][len - 1])*Math.Pow(_gait.impactSecondLine(0,1), 2) );
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(0, 1), 2);
            double torque1Max = 0;
            double torque2Max = 0;
            double torque1Diff = 0;
            double torque2Diff = 0;
            double[] values;
            if (double.IsNaN(dthetaTSquared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] {torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] {torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if(double.IsNaN(dtheta0Squared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] {torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] {torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else
            {
                double[,] THETA = new double[2, len + 1];
                THETA[0, 0] = 0;
                THETA[1, 0] = dtheta0Squared;

                double torque1 = 0;
                double torque2 = 0;
                double theta = 0;
                double dthetaSquared = 0;
                double ddtheta = 0;
                for (int i = 0; i < len; i += 5)
                {
                    theta = secondIntegral[0][i];
                    dthetaSquared = -secondIntegral[1][i] + Math.Exp(-firstIntegral[1][i]) * dtheta0Squared;
                    ddtheta = (-_gait.vhc.evalBeta(theta) * dthetaSquared - _gait.vhc.evalGamma(theta)) / _gait.vhc.evalAlpha(theta);
                    torque1 = Math.Abs(_gait.vhc.evalAlpha1(theta) * ddtheta + _gait.vhc.evalBeta1(theta) * dthetaSquared + _gait.vhc.evalGamma1(theta));
                    torque2 = Math.Abs(_gait.vhc.evalAlpha3(theta) * ddtheta + _gait.vhc.evalBeta3(theta) * dthetaSquared + _gait.vhc.evalGamma3(theta));
                    if (torque1 > torque1Max)
                    {
                        torque1Max = torque1;
                    }
                    if (torque2 > torque2Max)
                    {
                        torque2Max = torque2;
                    }
                }
                double ddtheta0 = (-_gait.vhc.evalBeta(0) * dtheta0Squared - _gait.vhc.evalGamma(0)) / _gait.vhc.evalAlpha(0);

                torque1Diff = _gait.vhc.evalAlpha1(0) * ddtheta0 + _gait.vhc.evalBeta1(0) * dtheta0Squared + _gait.vhc.evalGamma1(0) -
                    (_gait.vhc.evalAlpha3(1) * ddtheta + _gait.vhc.evalBeta3(1) * dthetaTSquared + _gait.vhc.evalGamma3(1));

                torque2Diff = _gait.vhc.evalAlpha3(0) * ddtheta0 + _gait.vhc.evalBeta3(0) * dtheta0Squared + _gait.vhc.evalGamma3(0) -
                    (_gait.vhc.evalAlpha1(1) * ddtheta + _gait.vhc.evalBeta1(1) * dthetaTSquared + _gait.vhc.evalGamma1(1));
                values = new double[] {torque1Max, torque2Max, torque1Diff, torque2Diff };

            }
            return values;
        }

        public double evaluateAlphaConstraint()
        {
            

            double dx = 0.02;
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


    public class SQPFifthOrder
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

        double _q1End;
        double _q2End;
        double _q3End;

        //private double[] _eqconstraint1Gradient;
        //private double[] _eqconstraint2Gradient;
        //private double[] _eqconstraint3Gradient;
        //private double[] _gradient;
        //private double[,] _hessian;
        private double[] _parameterValues;

        //private double _performanceIndexVal;

        public delegate double func(Expression exp);

        public SQPFifthOrder(BRgait gait, int numOfParams)
        {
            _gait = gait;
            BRVHC vhc = gait.vhc;
            _eqConstraints = new List<Expression>();
            _ineqConstraints = new List<Expression>();
            _eqconstraintsGradientArray = new List<Expression[]>();
            _ineqconstraintsGradientArray = new List<Expression[]>();
            _parameterValues = new double[3 * numOfParams - 3];
            string[] parametersFromFile = File.ReadAllLines(@"../../../../parameters.txt");
            for (int i = 0; i < numOfParams - 1; i++)
            {
                _parameterValues[i] = Convert.ToDouble(parametersFromFile[i + 1], CultureInfo.InvariantCulture);
            }
            for (int i = numOfParams - 1; i < 2 * numOfParams - 2; i++)
            {
                _parameterValues[i] = Convert.ToDouble(parametersFromFile[i + 2], CultureInfo.InvariantCulture);

            }
            for (int i = 2 * numOfParams - 2; i < 3 * numOfParams - 3; i++)
            {
                _parameterValues[i] = Convert.ToDouble(parametersFromFile[i + 3], CultureInfo.InvariantCulture);

            }
            string[] pLHS = File.ReadAllLines(@"../../../../pLHS.txt");
            string[] pRHS = File.ReadAllLines(@"../../../../pRHS.txt");
            string[] posture = File.ReadAllLines(@"../../../../posture.txt");
            Matrix<double> Q = Matrix<double>.Build.Dense(3, 3);
            Matrix<double> P = Matrix<double>.Build.Dense(3, 3);
            Matrix<double> pmatrix = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {0,0,1},
                {0,1,0},
                {1,0,0},
            
            });
            Dictionary<string, FloatingPoint> parameters = gait.vhc.parameters;
            _q1End = -Convert.ToDouble(posture[0], CultureInfo.InvariantCulture);
            _q2End = Convert.ToDouble(posture[1], CultureInfo.InvariantCulture);
            _q3End = -Convert.ToDouble(posture[2], CultureInfo.InvariantCulture);
            parameters.Add("q1", _q1End);
            parameters.Add("q2", _q2End);
            parameters.Add("q3", _q3End);
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
                temp = temp.Replace("(double)", "").Replace("0.2e1", "2");
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


        




        public double[] runNumericalSQP()
        {
            E04.E04WD_CONFUN confunE04WD = new E04.E04WD_CONFUN(confun);
            E04.E04WD_OBJFUN objfunE04WD = new E04.E04WD_OBJFUN(objfun);
            try
            {
                double objf = 0.0;
                int ifail, majits, n, nclin, ncnln;
                n = 15;
                nclin = 3;
                ncnln = 4;
                double[,] a = new double[,]{
                    {1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1},
                };
                double[] bl = new double[] { -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, -1.0E+25, 2*_q1End, 0, 2*_q3End, 0, 0, 0.0, -1.0E+25 };
                double[] bu = new double[] { 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 1.0E+25, 2*_q1End, 0, 2*_q3End, 0, 0, 50, 0.0 };
                double[] ccon = new double[ncnln];
                double[,] cjac = new double[ncnln, n];
                double[] clamda = new double[nclin + ncnln + n];
                double[] grad = new double[n];
                double[,] h = new double[n, n];
                double[] p = _parameterValues;
                int[] istate = new int[nclin + ncnln + n];

                E04.e04wdOptions options = new E04.e04wdOptions();
                options.Set("Print file = 6");
                options.Set("Derivative level = 0");
                //options.Set("Major Step limit = 0.00000000001");
                options.Set("Major Feasibility Tolerance = 1.0E-6");
                options.Set("Function precision = 1.0E-5");

                E04.e04wd(n, nclin, ncnln, a, bl, bu, confunE04WD, objfunE04WD, out majits, istate,
                    ccon, cjac, clamda, out objf, grad, h, p, options, out ifail);
                Console.WriteLine("");
                Console.WriteLine("  On exit from e04wd, ifail = {0,5}", ifail);
                if (ifail == 0)
                {
                    Console.WriteLine("  Final objective value = {0,11:f3}", objf);
                    Console.WriteLine("");
                    Console.Write("  Output x");
                    for (int i = 1; i <= n; i++)
                    {
                        Console.Write(" " + " {0,11:f3}", p[i - 1]);
                    }
                    Console.WriteLine(" ");
                }
                return p;
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                Console.WriteLine("Exception Raised");
                return null;
            }

        }
        public void objfun(ref int mode, int n, double[] p, ref double objf,
            double[] grad, int nstate)
        {
            int len = p.Length / 3;


            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + (i + 1).ToString()] = p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + (i + 2).ToString()] = p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + (i + 3).ToString()] = p[i];

            }
            double[] values = evalTorques();
            objf = Math.Pow(values[2], 2) + Math.Pow(values[3], 2);// + Math.Pow(values[2], 2) + Math.Pow(values[3], 2);
            mode = 0;

        }
        public void confun(ref int mode, int ncnln, int n, int[] needc, double[] p,
            double[] ccon, double[,] cjac, int nstate)
        {
            int len = p.Length / 3;
            double val1 = _gait.vhc.phi1Parameters["P0"].RealValue;
            double val2 = _gait.vhc.phi2Parameters["P" + (len + 1).ToString()].RealValue;
            double val3 = _gait.vhc.phi3Parameters["P" + (2 * len + 2).ToString()].RealValue;
            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + (i + 1).ToString()] = p[i];
                val1 += p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + (i + 2).ToString()] = p[i];
                val2 += p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + (i + 3).ToString()] = p[i];
                val3 += p[i];

            }
            double alphaMax = evaluateAlphaConstraint();
            double dthetaMin = evalDtheta();
            ccon[0] = _gait.impactFirstLine(0, 1) - _gait.impactSecondLine(0, 1);
            ccon[1] = _gait.impactSecondLine(0, 1) - _gait.impactThirdLine(0, 1);
            ccon[2] = dthetaMin;
            ccon[3] = alphaMax;
            mode = 0;

        }
        public void evaluateObjFuncAndConstraintsNumerical(double[] p, double[] objFunction, object obj)
        {
            double[] values = new double[5];
            double alphaMax = evaluateAlphaConstraint();
            if (alphaMax >= 0)
            {
                values[0] = 1;
            }
            else
            {
                values = evaluateDthetaAndTorqueConstraint();
            }
            int len = p.Length / 3;

            //add values to calculate endposture constraint
            double val1 = _gait.vhc.phi1Parameters["P0"].RealValue;
            double val2 = _gait.vhc.phi2Parameters["P" + (len + 1).ToString()].RealValue;
            double val3 = _gait.vhc.phi3Parameters["P" + (2 * len + 2).ToString()].RealValue;
            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + (i + 1).ToString()] = p[i];
                val1 += p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + (i + 2).ToString()] = p[i];
                val2 += p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + (i + 3).ToString()] = p[i];
                val3 += p[i];

            }

            if (values[0] == 1)
            {
                objFunction[0] = Math.Pow(10, 300);
                objFunction[1] = 0;
                objFunction[2] = 0;
                objFunction[3] = 0;
                objFunction[4] = 0;
                objFunction[5] = 0;
                objFunction[6] = 0;
                objFunction[7] = 0;
            }
            else
            {
                objFunction[0] = Math.Pow(values[2], 2) + Math.Pow(values[3], 2) + Math.Pow(values[4], 2) + Math.Pow(values[5], 2);
                objFunction[1] = val1 - _q1End;
                objFunction[2] = val2 - _q2End;
                objFunction[3] = val3 - _q3End;
                objFunction[4] = _gait.impactFirstLine(0, 1) - _gait.impactSecondLine(0, 1);
                objFunction[5] = _gait.impactFirstLine(0, 1) - _gait.impactThirdLine(0, 1);
                objFunction[6] = alphaMax;
                objFunction[7] = -values[1];
            }
        }



        public double[] evaluateDthetaAndTorqueConstraint()
        {
            double[][] firstIntegral = RiemannSum.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha);
            double[][] secondIntegral = RiemannSum.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1]);
            int len = firstIntegral[1].Length - 1;
            double dthetaTSquared = (-secondIntegral[1][len]) / (1 - Math.Exp(-firstIntegral[1][len]) * Math.Pow(_gait.impactSecondLine(0, 1), 2));
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(0, 1), 2);
            double dthetaMin = 0;
            double torque1Max = 0;
            double torque2Max = 0;
            double torque1Diff = 0;
            double torque2Diff = 0;
            double infeasible = 0;
            double[] values;
            if (double.IsNaN(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                infeasible = 1;
                values = new double[] { infeasible, dthetaMin, torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                infeasible = 1;
                values = new double[] { infeasible, dthetaMin, torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsNaN(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                infeasible = 1;
                values = new double[] { infeasible, dthetaMin, torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                infeasible = 1;
                values = new double[] { infeasible, dthetaMin, torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else
            {
                double[,] THETA = new double[2, len];
                THETA[0, 0] = 0;
                THETA[1, 0] = Math.Sqrt(dtheta0Squared);
                dthetaMin = Math.Sqrt(dtheta0Squared);

                double torque1 = 0;
                double torque2 = 0;
                double ddtheta = 0;
                for (int i = 1; i < len; i++)
                {
                    THETA[0, i] = secondIntegral[0][i];
                    THETA[1, i] = Math.Sqrt(-secondIntegral[1][i] + Math.Exp(-firstIntegral[1][i]) * dtheta0Squared);
                    if (THETA[1, i] < dthetaMin)
                    {
                        dthetaMin = THETA[1, i];
                    }
                }
                for (int i = 0; i < len; i += 5)
                {
                    ddtheta = (-_gait.vhc.evalBeta(THETA[0, i]) * Math.Pow(THETA[1, i], 2) - _gait.vhc.evalGamma(THETA[0, i])) / _gait.vhc.evalAlpha(THETA[0, i]);
                    torque1 = Math.Abs(_gait.vhc.evalAlpha1(THETA[0, i]) * ddtheta + _gait.vhc.evalBeta1(THETA[0, i]) * Math.Pow(THETA[1, i], 2) + _gait.vhc.evalGamma1(THETA[0, i]));
                    torque2 = Math.Abs(_gait.vhc.evalAlpha3(THETA[0, i]) * ddtheta + _gait.vhc.evalBeta3(THETA[0, i]) * Math.Pow(THETA[1, i], 2) + _gait.vhc.evalGamma3(THETA[0, i]));
                    if (torque1 > torque1Max)
                    {
                        torque1Max = torque1;
                    }
                    if (torque2 > torque2Max)
                    {
                        torque2Max = torque2;
                    }
                }
                double ddtheta0 = (-_gait.vhc.evalBeta(THETA[0, 0]) * Math.Pow(THETA[1, 0], 2) - _gait.vhc.evalGamma(THETA[0, 0])) / _gait.vhc.evalAlpha(THETA[0, 0]);

                torque1Diff = _gait.vhc.evalAlpha1(THETA[0, 0]) * ddtheta0 + _gait.vhc.evalBeta1(THETA[0, 0]) * Math.Pow(THETA[1, 0], 2) + _gait.vhc.evalGamma1(THETA[0, 0]) -
                    (_gait.vhc.evalAlpha3(THETA[0, len - 1]) * ddtheta + _gait.vhc.evalBeta3(THETA[0, len - 1]) * Math.Pow(THETA[1, len - 1], 2) + _gait.vhc.evalGamma3(THETA[0, len - 1]));

                torque2Diff = _gait.vhc.evalAlpha3(THETA[0, 0]) * ddtheta0 + _gait.vhc.evalBeta3(THETA[0, 0]) * Math.Pow(THETA[1, 0], 2) + _gait.vhc.evalGamma3(THETA[0, 0]) -
                    (_gait.vhc.evalAlpha1(THETA[0, len - 1]) * ddtheta + _gait.vhc.evalBeta1(THETA[0, len - 1]) * Math.Pow(THETA[1, len - 1], 2) + _gait.vhc.evalGamma1(THETA[0, len - 1]));
                values = new double[] { infeasible, dthetaMin, torque1Max, torque2Max, torque1Diff, torque2Diff };

            }
            return values;
        }

        public double evalDtheta()
        {
            double[][] firstIntegral = RiemannMiddleSum.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha);
            double[][] secondIntegral = RiemannMiddleSum.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1]);
            int len = firstIntegral[1].Length;
            double dthetaTSquared = (-secondIntegral[1][len - 1]) / (1 - Math.Exp(-firstIntegral[1][len - 1]) * Math.Pow(_gait.impactSecondLine(0, 1), 2));
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(0, 1), 2);
            double dthetaMin = 0;
            if (double.IsNaN(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsNaN(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else
            {
                double[,] THETA = new double[2, len];
                THETA[0, 0] = 0;
                THETA[1, 0] = dtheta0Squared;
                dthetaMin = dtheta0Squared;
                for (int i = 1; i < len; i++)
                {
                    THETA[0, i] = secondIntegral[0][i];
                    THETA[1, i] = -secondIntegral[1][i] + Math.Exp(-firstIntegral[1][i]) * dtheta0Squared;
                    if (THETA[1, i] < dthetaMin)
                    {
                        dthetaMin = THETA[1, i];
                    }
                }
            }
            return dthetaMin;
        }
        public double[] evalTorques()
        {
            double[][] firstIntegral = RiemannMiddleSum.calculateFirstIntegral(_gait.firstIntegral);
            double[][] secondIntegral = RiemannMiddleSum.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1]);
            int len = firstIntegral[1].Length;
            double dthetaTSquared = (-secondIntegral[1][len - 1]) / (1 - Math.Exp(-firstIntegral[1][len - 1]) * Math.Pow(_gait.impactSecondLine(0, 1), 2));
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(0, 1), 2);
            double torque1Max = 0;
            double torque2Max = 0;
            double torque1Diff = 0;
            double torque2Diff = 0;
            double[] values;
            if (double.IsNaN(dthetaTSquared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsNaN(dtheta0Squared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else
            {
                double[,] THETA = new double[2, len];
                THETA[0, 0] = 0;
                THETA[1, 0] = dtheta0Squared;

                double torque1 = 0;
                double torque2 = 0;
                double theta = 0;
                double dthetaSquared = 0;
                double ddtheta = 0;
                for (int i = 0; i < len; i += 5)
                {
                    theta = secondIntegral[0][i];
                    dthetaSquared = -secondIntegral[1][i] + Math.Exp(-firstIntegral[1][i]) * dtheta0Squared;
                    ddtheta = (-_gait.vhc.evalBeta(theta) * dthetaSquared - _gait.vhc.evalGamma(theta)) / _gait.vhc.evalAlpha(theta);
                    torque1 = Math.Abs(_gait.vhc.evalAlpha1(theta) * ddtheta + _gait.vhc.evalBeta1(theta) * dthetaSquared + _gait.vhc.evalGamma1(theta));
                    torque2 = Math.Abs(_gait.vhc.evalAlpha3(theta) * ddtheta + _gait.vhc.evalBeta3(theta) * dthetaSquared + _gait.vhc.evalGamma3(theta));
                    if (torque1 > torque1Max)
                    {
                        torque1Max = torque1;
                    }
                    if (torque2 > torque2Max)
                    {
                        torque2Max = torque2;
                    }
                }
                double ddtheta0 = (-_gait.vhc.evalBeta(0) * dtheta0Squared - _gait.vhc.evalGamma(0)) / _gait.vhc.evalAlpha(0);
                torque1Diff = (_gait.vhc.evalAlpha1(0) * ddtheta0 + _gait.vhc.evalBeta1(0) * dtheta0Squared + _gait.vhc.evalGamma1(0)) -
                    (_gait.vhc.evalAlpha3(1) * ddtheta + _gait.vhc.evalBeta3(1) * dthetaTSquared + _gait.vhc.evalGamma3(1));

                torque2Diff = _gait.vhc.evalAlpha3(0) * ddtheta0 + _gait.vhc.evalBeta3(0) * dtheta0Squared + _gait.vhc.evalGamma3(0) -
                    (_gait.vhc.evalAlpha1(1) * ddtheta + _gait.vhc.evalBeta1(1) * dthetaTSquared + _gait.vhc.evalGamma1(1));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };

            }
            return values;
        }

        public double evaluateAlphaConstraint()
        {


            double dx = 0.02;
            double theta = 0;
            double val = double.MinValue;
            double alphaVal = 0;
            for (int i = 0; i < 1 / dx; i++)
            {
                alphaVal = _gait.vhc.evalAlpha(theta);
                if (alphaVal > val)
                {
                    val = alphaVal;
                }
                theta += dx;
            }
            return val;
        }
    }
}




