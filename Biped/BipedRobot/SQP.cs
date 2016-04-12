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
    public class AGS
    {
        private BRgait _gait;
        private Expression _performanceIndex;
        private Expression _eqConstraint1;
        private Expression _eqConstraint2;
        private Expression _eqConstraint3;
        private Expression _eqConstraint4;
        private Expression _eqConstraint5;
        private Expression _eqConstraint6;
        private Expression _eqConstraint7;
        private Expression _ineqConstraint1;
        private Expression _ineqConstraint2;
        private Expression _ineqConstraint3;
        private Expression _lagrangian;
        private Expression[] _performanceGradientArray;
        private Expression[,] _hessianMatrix;
        private Expression[] _eqconstraint1GradientArray;
        private Expression[] _eqconstraint2GradientArray;
        private Expression[] _eqconstraint3GradientArray;
        private Expression[] _eqconstraint4GradientArray;
        private Expression[] _eqconstraint5GradientArray;
        private Expression[] _ineqconstraint1GradientArray;
        private Expression[] _ineqconstraint2GradientArray;
        private Expression[] _ineqconstraint3GradientArray;

        private double[] _eqconstraint1Gradient;
        private double[] _eqconstraint2Gradient;
        private double[] _eqconstraint3Gradient;
        private double[] _gradient;
        private double[,] _hessian;
        private double[] _parameterValues;

        private double _performanceIndexVal;

        private Dictionary<string, FloatingPoint> _parameters;

        public delegate double func(Expression exp);

        public AGS(BRgait gait)
        {
            _gait = gait;
            BRVHC vhc = gait.vhc;
            //because of theta entry we take -1. this is used for the gradients and we need the length of all the minimizing parameters
            int len = vhc.phi1Parameters.Count;
            _parameterValues = new double[3*len];

            

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
            _parameters.Add("dtheta0Squared", 0);
            _parameters.Add("dthetaTSquared", 0);
            _parameters.Add("ddtheta0", 0);
            _parameters.Add("ddthetaT", 0);

            int k = 0;
            foreach (KeyValuePair<string, FloatingPoint> entry in vhc.phi1Parameters)
            {
                _parameterValues[k] = 1;//entry.Value.RealValue; //1;
                _parameters.Add(entry.Key, entry.Value);
                k++;
            }
            foreach (KeyValuePair<string, FloatingPoint> entry in vhc.phi2Parameters)
            {

                _parameterValues[k] = 1;//entry.Value.RealValue; //1;//
                _parameters.Add(entry.Key, entry.Value);
                k++;

            }
            foreach (KeyValuePair<string, FloatingPoint> entry in vhc.phi3Parameters)
            {
                _parameterValues[k] = 1;//entry.Value.RealValue; //
                _parameters.Add(entry.Key, entry.Value);
                k++;
            }
            StreamReader fs = null;
            fs = new StreamReader(@"../../../alpha1.txt");
            string temp = fs.ReadLine();
            Expression alpha1 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            alpha1 = Structure.Substitute("phi1", vhc.phi1, alpha1);
            alpha1 = Structure.Substitute("phi2", vhc.phi2, alpha1);
            alpha1 = Structure.Substitute("phi3", vhc.phi3, alpha1);
            alpha1 = Structure.Substitute("dphi1", vhc.dphi1, alpha1);
            alpha1 = Structure.Substitute("dphi2", vhc.dphi2, alpha1);
            alpha1 = Structure.Substitute("dphi3", vhc.dphi3, alpha1);
            fs.Close();

            fs = new StreamReader(@"../../../beta1.txt");
            temp = fs.ReadLine();
            Expression beta1 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            beta1 = Structure.Substitute("phi1", vhc.phi1, beta1);
            beta1 = Structure.Substitute("phi2", vhc.phi2, beta1);
            beta1 = Structure.Substitute("phi3", vhc.phi3, beta1);
            beta1 = Structure.Substitute("dphi2", vhc.dphi2, beta1);
            beta1 = Structure.Substitute("dphi3", vhc.dphi3, beta1);
            beta1 = Structure.Substitute("ddphi1", vhc.ddphi1, beta1);
            beta1 = Structure.Substitute("ddphi2", vhc.ddphi2, beta1);
            fs.Close();

            fs = new StreamReader(@"../../../gamma1.txt");
            temp = fs.ReadLine();
            Expression gamma1 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            gamma1 = Structure.Substitute("phi1", vhc.phi1, gamma1);
            fs.Close();

            Expression torque1 = alpha1 * "ddtheta" + beta1 * "dthetaSquared" + gamma1;
            Expression torque1Start = Structure.Substitute("theta", "0", torque1);
            torque1Start = Structure.Substitute("dthetaSquared", "dtheta0Squared", torque1Start);
            torque1Start = Structure.Substitute("ddtheta", "ddtheta0", torque1Start);
            Expression torque1End = Structure.Substitute("theta", "1", torque1);
            torque1End = Structure.Substitute("dthetaSquared", "dthetaTSquared", torque1End);
            torque1End = Structure.Substitute("ddtheta", "ddthetaT", torque1End);

            fs = new StreamReader(@"../../../alpha3.txt");
            temp = fs.ReadLine();
            Expression alpha3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            alpha3 = Structure.Substitute("phi1", vhc.phi1, alpha3);
            alpha3 = Structure.Substitute("phi3", vhc.phi3, alpha3);
            alpha3 = Structure.Substitute("dphi1", vhc.dphi1, alpha3);
            alpha3 = Structure.Substitute("dphi3", vhc.dphi3, alpha3);
            fs.Close();

            fs = new StreamReader(@"../../../beta3.txt");
            temp = fs.ReadLine();
            Expression beta3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            beta3 = Structure.Substitute("phi1", vhc.phi1, beta3);
            beta3 = Structure.Substitute("phi3", vhc.phi3, beta3);
            beta3 = Structure.Substitute("dphi1", vhc.dphi1, beta3);
            beta3 = Structure.Substitute("dphi3", vhc.dphi3, beta3);
            beta3 = Structure.Substitute("ddphi1", vhc.ddphi1, beta3);
            fs.Close();

            fs = new StreamReader(@"../../../gamma3.txt");
            temp = fs.ReadLine();
            Expression gamma3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            gamma3 = Structure.Substitute("phi3", vhc.phi3, gamma3);
            fs.Close();

            Expression torque2 = alpha3 * "ddtheta" + beta3 * "dthetaSquared" + gamma3;
            Expression torque2Start = Structure.Substitute("theta", "0", torque2);
            torque2Start = Structure.Substitute("dthetaSquared", "dtheta0Squared", torque2Start);
            torque2Start = Structure.Substitute("ddtheta", "ddtheta0", torque2Start);
            Expression torque2End = Structure.Substitute("theta", "1", torque2);
            torque2End = Structure.Substitute("dthetaSquared", "dthetaTSquared", torque2End);
            torque2End = Structure.Substitute("ddtheta", "ddthetaT", torque2End);

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

            _performanceIndex = Expression.Pow(impactNegFirstLine / impactPosFirstLine - impactNegSecondLine / impactPosSecondLine,"2") + Expression.Pow(impactNegFirstLine / impactPosFirstLine - impactNegThirdLine / impactPosThirdLine,"2");

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


            //normalisere. gange med 0 gir -1. phi1 starter med å være negativ. phi3 starter med å være positiv
            _eqConstraint1 = (phi1End + phi1Start);
            _eqConstraint2 = (phi2End - phi2Start);
            _eqConstraint3 = (phi3End + phi3Start);
            _eqConstraint4 = (phi1Start + phi3Start);
            _eqConstraint5 = (phi1End + phi3End);

            //torquesa må være like
            _eqConstraint6 = torque1Start - torque2End;
            _eqConstraint7 = torque2Start - torque1End;
            Console.WriteLine(Infix.Format(_eqConstraint6));
            Console.WriteLine(Infix.Format(_eqConstraint7));


            //faktorene må være positive
            _ineqConstraint1 = -(impactNegFirstLine / impactPosFirstLine);
            _ineqConstraint2 = -(impactNegSecondLine / impactPosSecondLine);
            _ineqConstraint3 = -(impactNegThirdLine / impactPosThirdLine);

            


            _lagrangian = _performanceIndex - "delta1" * _eqConstraint1 - "delta2" * _eqConstraint2 - "delta3" * _eqConstraint3;

            _performanceGradientArray = new Expression[len*3];
            _hessianMatrix = new Expression[len * 3, len * 3];

            _eqconstraint1GradientArray = new Expression[len * 3];
            _eqconstraint2GradientArray = new Expression[len * 3];
            _eqconstraint3GradientArray = new Expression[len * 3];
            _eqconstraint4GradientArray = new Expression[len * 3];
            _eqconstraint5GradientArray = new Expression[len * 3];

            _ineqconstraint1GradientArray = new Expression[len * 3];
            _ineqconstraint2GradientArray = new Expression[len * 3];
            _ineqconstraint3GradientArray = new Expression[len * 3];
            for (int i = 0; i < _performanceGradientArray.Length; i++)
            {
                _performanceGradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _performanceIndex);

                _eqconstraint1GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint1);
                _eqconstraint2GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint2);
                _eqconstraint3GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint3);
                _eqconstraint4GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint4);
                _eqconstraint5GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _eqConstraint5);

                _ineqconstraint1GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _ineqConstraint1);
                _ineqconstraint2GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _ineqConstraint2);
                _ineqconstraint3GradientArray[i] = Calculus.Differentiate("P" + i.ToString(), _ineqConstraint3);

                Expression expTemp = Calculus.Differentiate("P" + i.ToString(), _lagrangian);
                for (int j = 0; j < _performanceGradientArray.Length; j++)
                {
                    _hessianMatrix[i, j] = Calculus.Differentiate("P" + j.ToString(), expTemp);
                }
            }




        }

        public void runAnalytical(BRgait gait)
        {
            string a = Infix.Format(_performanceIndex);
            string b = Infix.Format(_eqConstraint1);
            string c = Infix.Format(_eqConstraint2);
            string d = Infix.Format(_eqConstraint3);

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
            testImpact(gait, p1);

        }
        public void runNumerical(BRgait gait)
        {
            double[] p0 = _parameterValues;
            double[] s = new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            double epsx = 0.001;
            double radius = 0.1;
            double rho = 50.0;
            double diffstep = 0.001;
            int maxits = 0;
            alglib.minnsstate state;
            alglib.minnsreport rep;
            double[] p1;
            alglib.minnscreatef(p0.Length, p0, diffstep, out state);
            alglib.minnssetalgoags(state, radius, rho);
            alglib.minnssetcond(state, epsx, maxits);
            alglib.minnssetscale(state, s);

            alglib.minnssetnlc(state, 5, 4);
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
        
        
        public double evaluateFunction(Expression exp)
        {
            string a = Infix.Format(exp);
            exp = Infix.ParseOrUndefined(a);
            try
            {
                return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, exp).RealValue;
            }
            catch (Exception)
            {
                return Math.Pow(10, 300);
            }
        }

        public void evaluateObjFuncAndConstraintsAnalytical(double[] p, double[] objFunction, double[,] jacobian, object obj)
        {
            for (int i = 0; i < p.Length; i++)
            {
                _parameters["P"+i.ToString()] = p[i];
                
            }
            int len = p.Length;
            for (int i = 0; i < len; i++)
            {
                
                jacobian[0, i] = evaluateFunction(_performanceGradientArray[i]);

                jacobian[1, i] = evaluateFunction(_eqconstraint1GradientArray[i]);
                jacobian[2, i] = evaluateFunction(_eqconstraint2GradientArray[i]);
                jacobian[3, i] = evaluateFunction(_eqconstraint3GradientArray[i]);
                jacobian[4, i] = evaluateFunction(_eqconstraint4GradientArray[i]);
                jacobian[5, i] = evaluateFunction(_eqconstraint5GradientArray[i]);

                jacobian[6, i] = evaluateFunction(_ineqconstraint1GradientArray[i]);
                jacobian[7, i] = evaluateFunction(_ineqconstraint2GradientArray[i]);
                jacobian[8, i] = evaluateFunction(_ineqconstraint3GradientArray[i]);
            }


            objFunction[0] = evaluateFunction(_performanceIndex);
            objFunction[1] = evaluateFunction(_eqConstraint1);
            objFunction[2] = evaluateFunction(_eqConstraint2);
            objFunction[3] = evaluateFunction(_eqConstraint3);
            objFunction[4] = evaluateFunction(_eqConstraint4);
            objFunction[5] = evaluateFunction(_eqConstraint5);
            objFunction[6] = evaluateFunction(_ineqConstraint1);
            objFunction[7] = evaluateFunction(_ineqConstraint2);
            objFunction[8] = evaluateFunction(_ineqConstraint3);
        }
        public void evaluateObjFuncAndConstraintsNumerical(double[] p, double[] objFunction, object obj)
        {
            double dthetaMin = evaluateDthetaConstraint();
            double alphaMax = evaluateAlphaConstraint();
            for (int i = 0; i < p.Length; i++)
            {
                _parameters["P" + i.ToString()] = p[i];

            }


            objFunction[0] = evaluateFunction(_performanceIndex);
            objFunction[1] = evaluateFunction(_eqConstraint1);
            objFunction[2] = evaluateFunction(_eqConstraint2);
            objFunction[3] = evaluateFunction(_eqConstraint3);
            objFunction[4] = evaluateFunction(_eqConstraint4);
            objFunction[5] = evaluateFunction(_eqConstraint5);
            //objFunction[6] = evaluateFunction(_eqConstraint6);
            //objFunction[7] = evaluateFunction(_eqConstraint7);
            objFunction[6] = evaluateFunction(_ineqConstraint1);
            objFunction[7] = evaluateFunction(_ineqConstraint2);
            objFunction[8] = evaluateFunction(_ineqConstraint3);
            //objFunction[11] = -dthetaMin;
            objFunction[9] = -alphaMax;

        }



        public double evaluateDthetaConstraint()
        {
            double[,] firstIntegral = RiemannSum.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha);
            double[,] secondIntegral = RiemannSum.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral);

            int len = secondIntegral.Length/secondIntegral.Rank;
            double dthetaTSquared = -secondIntegral[len - 1, 1] / (1 - Math.Exp(-firstIntegral[len - 1, 1]) * Math.Pow(-evaluateFunction(_ineqConstraint1), 2));
            double dtheta0Squared = Math.Pow(-evaluateFunction(_ineqConstraint1), 2) * dthetaTSquared;


            //set parameters for use later
            
            if(double.IsNaN(dthetaTSquared))
            {
                _parameters["dthetaTSquared"] = (FloatingPoint)(- Math.Pow(10, 300));
                _parameters["ddthetaT"] = (FloatingPoint)(- Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                _parameters["dthetaTSquared"] = (FloatingPoint)Math.Pow(10, 300);
                _parameters["ddthetaT"] = (FloatingPoint)Math.Pow(10, 300);
            }
            else
            {
                _parameters["dthetaTSquared"] = dthetaTSquared;
                _parameters["ddthetaT"] = (FloatingPoint)BRReducedDynamics.rhs1D(Vector<double>.Build.Dense(new double[] { 1, dthetaTSquared }), _gait.vhc.evalAlpha, _gait.vhc.evalBeta, _gait.vhc.evalGamma);
            }
            if(double.IsNaN(dtheta0Squared))
            {
                _parameters["dtheta0Squared"] = (FloatingPoint)(-Math.Pow(10, 300));
                _parameters["ddtheta0"] = (FloatingPoint)(-Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                _parameters["dtheta0Squared"] = (FloatingPoint)Math.Pow(10, 300);
                _parameters["ddtheta0"] = (FloatingPoint)Math.Pow(10, 300);
            }
            else
            {
                _parameters["dtheta0Squared"] = (FloatingPoint)dtheta0Squared;
                _parameters["ddtheta0"] = (FloatingPoint)BRReducedDynamics.rhs1D(Vector<double>.Build.Dense(new double[] { 0, dtheta0Squared }), _gait.vhc.evalAlpha, _gait.vhc.evalBeta, _gait.vhc.evalGamma);
            }
            double[] dthetaSquare = new double[len];
            double min = 0;
            double dthetatSquare = 0;
            for (int i = 0; i < len; i++)
            {
                dthetatSquare = -secondIntegral[i, 1] + Math.Exp(firstIntegral[i, 1]) * dtheta0Squared;
                if(dthetatSquare < min)
                {
                    min = dthetatSquare;
                }
            }
            if (double.IsNegativeInfinity(min)){
                return (-Math.Pow(10, 300));
            }
            if (double.IsPositiveInfinity(min))
            {
                return (Math.Pow(10, 300));
            }
            if (double.IsNaN(min))
            {
                return (-Math.Pow(10, 300));
            }
            return min;
        }

        public double evaluateAlphaConstraint()
        {
            _gait.vhc.phi1Parameters = _parameters;
            _gait.vhc.phi2Parameters = _parameters;
            _gait.vhc.phi3Parameters = _parameters;

            double dx = 0.002;
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

        public void testImpact(BRgait gait, double[] p1)
        {
            int numberOfPoints = 6;
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);
            tuple.Item1["P0"] = p1[0];
            tuple.Item1["P1"] = p1[1];
            tuple.Item1["P2"] = p1[2];
            tuple.Item1["P3"] = p1[3];
            tuple.Item1["P4"] = p1[4];
            tuple.Item1["P5"] = p1[5];
            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            tuple.Item1["P6"] = p1[6];
            tuple.Item1["P7"] = p1[7];
            tuple.Item1["P8"] = p1[8];
            tuple.Item1["P9"] = p1[9];
            tuple.Item1["P10"] = p1[10];
            tuple.Item1["P11"] = p1[11];
            gait.vhc.phi2Parameters = tuple.Item1;

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            tuple.Item1["P12"] = p1[12];
            tuple.Item1["P13"] = p1[13];
            tuple.Item1["P14"] = p1[14];
            tuple.Item1["P15"] = p1[15];
            tuple.Item1["P16"] = p1[16];
            tuple.Item1["P17"] = p1[17];
            gait.vhc.phi3Parameters = tuple.Item1;

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());
            Console.WriteLine(gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactSecondLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactThirdLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
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




