using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.Linq;
using System.IO;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;
using MathNet.Symbolics;

namespace BipedRobot
{
    public class BRGaitParameters
    {
        private Tuple<Vector<double>, Vector<double>, Vector<double>> _gaitParameters;
        private double _objFunVal;
        private double _intervalStart;
        private double _intervalEnd;
        private Vector<double> _initialControlPointsq1;
        private Vector<double> _initialControlPointsq2;
        private Vector<double> _initialControlPointsq3;
        

        public BRGaitParameters(int numOfPoints)
        {
            /*
            StreamReader fs = null;
            fs = new StreamReader(@"../../../parameters.txt");
            string temp = fs.ReadLine();
            temp = temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1);
            int tempint = Int32.Parse(temp);
            fs.Close();
            */
            Vector<double> q1 = Vector<double>.Build.Dense(numOfPoints);

            Vector<double> q2 = Vector<double>.Build.Dense(numOfPoints);

            Vector<double> q3 = Vector<double>.Build.Dense(numOfPoints);
            
            _gaitParameters = new Tuple<Vector<double>, Vector<double>, Vector<double>>(
                q1, q2, q3);

            _initialControlPointsq1 = Vector<double>.Build.Dense(numOfPoints);
            _initialControlPointsq2 = Vector<double>.Build.Dense(numOfPoints);
            _initialControlPointsq3 = Vector<double>.Build.Dense(numOfPoints);
            setPosture();
            setInitialControlPoints(numOfPoints);
        }

        public void setInitialControlPoints(int numOfPoints)
        {
            _initialControlPointsq1[0] = _intervalStart;
            _initialControlPointsq1[1] = _intervalStart/2;
            _initialControlPointsq1[2] = 0;
            _initialControlPointsq1[3] = _intervalEnd/2;
            _initialControlPointsq1[4] = _intervalEnd;


            _initialControlPointsq2[0] = 0;
            _initialControlPointsq2[1] = Math.PI/18;
            _initialControlPointsq2[2] = 0;
            _initialControlPointsq2[3] = -Math.PI/18;
            _initialControlPointsq2[4] = 0;


            _initialControlPointsq3[0] = _intervalEnd;
            _initialControlPointsq3[1] = _intervalEnd / 2;
            _initialControlPointsq3[2] = 0;
            _initialControlPointsq3[3] = _intervalStart/2;
            _initialControlPointsq3[4] = _intervalStart;
        }
        public void setPosture()
        {
            _intervalStart = -Math.PI / 8;
            _intervalEnd = Math.PI / 8;
        }

        public Tuple<Vector<double>, Vector<double>, Vector<double>> gaitparameters
        {
            get
            {
                return _gaitParameters;
            }
            set
            {
                _gaitParameters = value;
            }
        }
        public double objFunVal
        {
            get
            {
                return _objFunVal;
            }
            set
            {
                _objFunVal = value;
            }
        }
        public double intervalStart
        {
            get
            {
                return _intervalStart;
            }
            set
            {
                _intervalStart = value;
            }
        }
        public double intervalEnd
        {
            get
            {
                return _intervalEnd;
            }
            set
            {
                _intervalEnd = value;
            }
        }
        public Vector<double> initialControlPointsq1
        {
            get
            {
                return _initialControlPointsq1;
            }
            set
            {
                _initialControlPointsq1 = value;
            }
        }
        public Vector<double> initialControlPointsq2
        {
            get
            {
                return _initialControlPointsq2;
            }
            set
            {
                _initialControlPointsq2 = value;
            }
        }
        public Vector<double> initialControlPointsq3
        {
            get
            {
                return _initialControlPointsq3;
            }
            set
            {
                _initialControlPointsq3 = value;
            }
        }
    }

    public static class gaitSearch
    {
        public static void run(ref Biped biped)
        {
            //first time running use rand to find a valid gait
            int numberOfPoints = 5;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            for (int i = 0; i < 200000; i++)
            {
                setParametersRandom(ref gait);
                setVHC(ref gait, numberOfPoints);
                if (verifyParameters(gait) && verifyImpact(gait))
                {
                    Console.WriteLine("found gait");
                    Console.WriteLine(i);
                    Infix.Format(gait.vhc.phi1);
                    Infix.Format(gait.vhc.phi2);
                    Infix.Format(gait.vhc.phi3);
                    biped.gaits.Add(gait);
                    break;
                }
            }
            //begin search
            SQP sqp = new SQP(gait.vhc);





            //write to file
            StreamWriter fs = null;
            fs = new StreamWriter(@"../../../foundGaits.txt");
            foreach (BRgait g in biped.gaits)
            {
                fs.WriteLine(Infix.Format(g.vhc.phi1));
                fs.WriteLine(Infix.Format(g.vhc.phi1));
                fs.WriteLine(Infix.Format(g.vhc.phi1));
                fs.WriteLine(g.impactFirstLine(gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd));
                fs.WriteLine(g.impactSecondLine(gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd));
                fs.WriteLine(g.impactThirdLine(gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd));
                fs.WriteLine();
            }
            fs.Close();
            
        }
        public static void setParametersRandom(ref BRgait gait)
        {
            Random rndm = new Random();
            double legFactor = gait.gaitParam.initialControlPointsq1[0];

            gait.gaitParam.gaitparameters.Item1[0] = gait.gaitParam.initialControlPointsq1[0];
            for (int i = 1; i < gait.gaitParam.gaitparameters.Item1.Count; i++)
            {
                gait.gaitParam.gaitparameters.Item1[i] = gait.gaitParam.initialControlPointsq1[i] + 2 * legFactor * (rndm.NextDouble() - 0.5);
                
                    
            }
            
            gait.gaitParam.gaitparameters.Item2[0] = gait.gaitParam.initialControlPointsq2[0];
            for (int i = 1; i < gait.gaitParam.gaitparameters.Item2.Count; i++)
            {
                gait.gaitParam.gaitparameters.Item2[i] = gait.gaitParam.initialControlPointsq2[i] + 2 * (rndm.NextDouble() - 0.5);


            }
            gait.gaitParam.gaitparameters.Item3[0] = gait.gaitParam.initialControlPointsq3[0];
            for (int i = 1; i < gait.gaitParam.gaitparameters.Item3.Count; i++)
            {
                gait.gaitParam.gaitparameters.Item3[i] = gait.gaitParam.initialControlPointsq3[i] + 2 * legFactor * (rndm.NextDouble() - 0.5);


            }
        }
        public static void setVHC(ref BRgait gait, int numberOfPoints)
        {
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);
            Tuple<Dictionary<string, double>, string> tuple = new Tuple<Dictionary<string, double>, string>();
            gait.vhc.phi1 = Infix.ParseOrUndefined(brCrv.phi1ToString());
            gait.vhc.phi2 = Infix.ParseOrUndefined(brCrv.phi2ToString());
            gait.vhc.phi3 = Infix.ParseOrUndefined(brCrv.phi3ToString());

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());
        }
        public static bool verifyParameters(BRgait gait)
        {
            double thetaDotAtTSquare = (-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.secondIntegral, gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd, 32)) /
                (1 - Math.Exp(-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.firstIntegral, gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd, 32)) 
                * Math.Pow(gait.impactFirstLine(gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd), 2));
            //Console.WriteLine(thetaDotAtTSquare);
            if (thetaDotAtTSquare < 0.1 || thetaDotAtTSquare > 10000 || thetaDotAtTSquare == double.NaN)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        public static bool verifyImpact(BRgait gait)
        {
            double firstAndSecond = gait.impactFirstLine(gait.gaitParam.intervalStart,gait.gaitParam.intervalEnd) - 
                gait.impactSecondLine(gait.gaitParam.intervalStart,gait.gaitParam.intervalEnd); 
            double firstAndThird = gait.impactFirstLine(gait.gaitParam.intervalStart,gait.gaitParam.intervalEnd) - 
                gait.impactThirdLine(gait.gaitParam.intervalStart,gait.gaitParam.intervalEnd);
            if ((firstAndSecond <= 0.001 && firstAndSecond >= -0.001) && (firstAndThird <= 0.001 && firstAndThird >= -0.001)) 
            {
                Console.WriteLine(gait.impactFirstLine(gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd));
                Console.WriteLine(gait.impactSecondLine(gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd));
                Console.WriteLine(gait.impactThirdLine(gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd));
                return true;
            }
            return false;
        }
    }

    public class BRVHC
    {
        public BRVHC(BRParameters param, int numberOfPoints)
        {
            _parameters = new Dictionary<string, FloatingPoint>();

            StreamReader fs = null;
            fs = new StreamReader(@"../../../q3.txt");
            string temp = fs.ReadLine();
            //temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _q3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../ddq3.txt");
            temp = fs.ReadLine();
            //temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _ddq3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../alpha.txt");
            temp = fs.ReadLine();
            //temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _alpha = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../beta.txt");
            temp = fs.ReadLine();
            //temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _beta = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../gamma.txt");
            temp = fs.ReadLine();
            //temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _gamma = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactPosFirstLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("(double)", "").Replace("0.2e1","2");
            _impactPosFirstLine = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactPosSecondLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("(double)", "").Replace("0.2e1", "2");
            _impactPosSecondLine = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactPosThirdLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("(double)", "").Replace("0.2e1", "2");
            _impactPosThirdLine = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactNegFirstLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("(double)", "");
            _impactNegFirstLine = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactNegSecondLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("(double)", "");
            _impactNegSecondLine = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactNegThirdLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("(double)", "");
            _impactNegThirdLine = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();
            setPhysicalParameters(param, numberOfPoints);
        }
        #region parameters
        private Dictionary<string, FloatingPoint> _parameters;
        private Dictionary<string, FloatingPoint> _phi1Parameters;
        private Dictionary<string, FloatingPoint> _phi2Parameters;
        private Dictionary<string, FloatingPoint> _phi3Parameters;

        public Dictionary<string, FloatingPoint> parameters
        {
            get
            {
                return _parameters;
            }
            set
            {
                _parameters = value;
            }
        }
        public Dictionary<string, FloatingPoint> phi1Parameters
        {
            get
            {
                return _phi1Parameters;
            }
            set
            {
                _phi1Parameters = value;
            }
        }
        public Dictionary<string, FloatingPoint> phi2Parameters
        {
            get
            {
                return _phi2Parameters;
            }
            set
            {
                _phi2Parameters = value;
            }
        }
        public Dictionary<string, FloatingPoint> phi3Parameters
        {
            get
            {
                return _phi3Parameters;
            }
            set
            {
                _phi3Parameters = value;
            }
        }
        #endregion
        #region general coordinates
        private Expression _q1;
        private Expression _q2;
        private Expression _q3;

        private Expression _dq1;
        private Expression _dq2;
        private Expression _dq3;

        private Expression _ddq1;
        private Expression _ddq2;
        private Expression _ddq3;

        public double evalq1(double theta)
        {
            return 0;
        }
        public double evalq2(double theta)
        {
            return Math.Cos(theta);
        }
        public double evalq3(double theta)
        {
            _parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _q3).RealValue;
        }

        public double evaldq1(double theta, double dtheta)
        {
            return Math.Cos(theta);
        }
        public double evaldq2(double theta, double dtheta)
        {
            return Math.Cos(theta);
        }
        public double evaldq3(double theta, double dtheta)
        {
            return 0;
        }

        public double evalddq1(double theta, double dtheta, double ddtheta)
        {
            return Math.Cos(theta);
        }
        public double evalddq2(double theta, double dtheta, double ddtheta)
        {
            return Math.Cos(theta);
        }
        public double evalddq3(double theta, double dtheta, double ddtheta)
        {
            _parameters["theta"] = theta;
            _parameters["dtheta"] = dtheta;
            _parameters["ddtheta"] = ddtheta;

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _ddq3).RealValue;
        }

        public Expression q1
        {
            get
            {
                return _q1;
            }
            set
            {
                _q1 = value;
            }
        }
        public Expression q2
        {
            get
            {
                return _q2;
            }
            set
            {
                _q2 = value;
            }
        }
        public Expression q3
        {
            get
            {
                return _q3;
            }
            set
            {
                _q3 = value;
            }
        }

        public Expression dq1
        {
            get
            {
                return _dq1;
            }
            set
            {
                _dq1 = value;
            }
        }
        public Expression dq2
        {
            get
            {
                return _dq2;
            }
            set
            {
                _dq2 = value;
            }
        }
        public Expression dq3
        {
            get
            {
                return _dq3;
            }
            set
            {
                _dq3 = value;
            }
        }

        public Expression ddq1
        {
            get
            {
                return _ddq1;
            }
            set
            {
                _ddq1 = value;
            }
        }
        public Expression ddq2
        {
            get
            {
                return _ddq2;
            }
            set
            {
                _ddq2 = value;
            }
        }
        public Expression ddq3
        {
            get
            {
                return _ddq3;
            }
            set
            {
                _ddq3 = value;
            }
        }
        #endregion
        #region constraints
        private Expression _phi1;
        private Expression _phi2;
        private Expression _phi3;

        private Expression _dphi1;
        private Expression _dphi2;
        private Expression _dphi3;

        private Expression _ddphi1;
        private Expression _ddphi2;
        private Expression _ddphi3;

        public double evalPhi1(double theta)
        {
            _parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _phi1).RealValue;
        }
        public double evalPhi2(double theta)
        {
            _parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _phi2).RealValue;
        }
        public double evalPhi3(double theta)
        {
            _parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _phi3).RealValue;
        }

        public double evalDphi1(double theta)
        {
            _parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _dphi1).RealValue;
        }
        public double evalDphi2(double theta)
        {
            _parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _dphi2).RealValue;
        }
        public double evalDphi3(double theta)
        {
            _parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _dphi3).RealValue;
        }

        public double evalDdphi1(double theta)
        {
            _parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _ddphi1).RealValue;
        }
        public double evalDdphi2(double theta)
        {
            _parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _ddphi2).RealValue;
        }
        public double evalDdphi3(double theta)
        {
            _parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _ddphi3).RealValue;
        }

        public Expression phi1
        {
            get
            {
                return _phi1;
            }
            set
            {
                _phi1 = value;
            }
        }
        public Expression phi2
        {
            get
            {
                return _phi2;
            }
            set
            {
                _phi2 = value;
            }
        }
        public Expression phi3
        {
            get
            {
                return _phi3;
            }
            set
            {
                _phi3 = value;
            }
        }

        public Expression dphi1
        {
            get
            {
                return _dphi1;
            }
            set
            {
                _dphi1 = value;
            }
        }
        public Expression dphi2
        {
            get
            {
                return _dphi2;
            }
            set
            {
                _dphi2 = value;
            }
        }
        public Expression dphi3
        {
            get
            {
                return _dphi3;
            }
            set
            {
                _dphi3 = value;
            }
        }

        public Expression ddphi1
        {
            get
            {
                return _ddphi1;
            }
            set
            {
                _ddphi1 = value;
            }
        }
        public Expression ddphi2
        {
            get
            {
                return _ddphi2;
            }
            set
            {
                _ddphi2 = value;
            }
        }
        public Expression ddphi3
        {
            get
            {
                return _ddphi3;
            }
            set
            {
                _ddphi3 = value;
            }
        }
        #endregion
        #region zerodynamics
        private Expression _alpha;
        private Expression _beta;
        private Expression _gamma;

        public double evalAlpha(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _alpha).RealValue;
        }
        public double evalBeta(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _beta).RealValue;
        }
        public double evalGamma(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _gamma).RealValue;
        }

        public double evalTwoTimesBetaDividedByAlpha(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);

            return (2 * (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _beta).RealValue / (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _alpha).RealValue);
        }
        public double evalTwoTimesGammaDividedByAlpha(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);

            return (2 * (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _gamma).RealValue / (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _alpha).RealValue);
        }
        #endregion
        #region impacts
        private Expression _impactPosFirstLine;
        private Expression _impactPosSecondLine;
        private Expression _impactPosThirdLine;

        private Expression _impactNegFirstLine;
        private Expression _impactNegSecondLine;
        private Expression _impactNegThirdLine;


        public double evalImpactPosFirstLine(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _impactPosFirstLine).RealValue;
        }
        public double evalImpactPosSecondLine(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _impactPosSecondLine).RealValue;
        }
        public double evalImpactPosThirdLine(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _impactPosThirdLine).RealValue;
        }

        public double evalImpactNegFirstLine(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _impactNegFirstLine).RealValue;
        }
        public double evalImpactNegSecondLine(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _impactNegSecondLine).RealValue;
        }
        public double evalImpactNegThirdLine(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDphi1(theta);
            _parameters["ddphi2"] = evalDphi2(theta);
            _parameters["ddphi3"] = evalDphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _impactNegThirdLine).RealValue;
        }

        public Expression impactPosFirstLine
        {
            get
            {
                return _impactPosFirstLine;
            }
            set
            {
                _impactPosFirstLine = value;
            }
        }
        public Expression impactPosSecondLine
        {
            get
            {
                return _impactPosSecondLine;
            }
            set
            {
                _impactPosSecondLine = value;
            }
        }
        public Expression impactPosThirdLine
        {
            get
            {
                return _impactPosThirdLine;
            }
            set
            {
                _impactPosThirdLine = value;
            }
        }

        public Expression impactNegFirstLine
        {
            get
            {
                return _impactNegFirstLine;
            }
            set
            {
                _impactNegFirstLine = value;
            }
        }
        public Expression impactNegSecondLine
        {
            get
            {
                return _impactNegSecondLine;
            }
            set
            {
                _impactNegSecondLine = value;
            }
        }
        public Expression impactNegThirdLine
        {
            get
            {
                return _impactNegThirdLine;
            }
            set
            {
                _impactNegThirdLine = value;
            }
        }
        #endregion
        public void setPhysicalParameters(BRParameters param, int numberOfPoints)
        {
            _parameters.Add("g", param.g);
            _parameters.Add("m1", param.m1);
            _parameters.Add("m2", param.m2);
            _parameters.Add("m3", param.m3);

            _parameters.Add("l1", param.l1);
            _parameters.Add("l2", param.l2);
            _parameters.Add("l3", param.l3);

            _parameters.Add("L1", param.L1);
            _parameters.Add("L2", param.L2);
            _parameters.Add("L3", param.L3);

            _parameters.Add("J1", param.J1);
            _parameters.Add("J2", param.J2);
            _parameters.Add("J3", param.J3);

            _parameters.Add("phi1", 0);
            _parameters.Add("phi2", 0);
            _parameters.Add("phi3", 0);

            _parameters.Add("dphi1", 0);
            _parameters.Add("dphi2", 0);
            _parameters.Add("dphi3", 0);

            _parameters.Add("ddphi1", 0);
            _parameters.Add("ddphi2", 0);
            _parameters.Add("ddphi3", 0);

            _parameters.Add("theta", 0);
            _parameters.Add("dtheta", 0);
            _parameters.Add("ddtheta", 0);

            for (int i = 1; i <= numberOfPoints; i++)
            {
                _parameters.Add("P"+ i.ToString(), 0);
            }
        }


    }

  
    public class BRgait
    {
        private BRVHC _vhc;
        private BRGaitParameters _gaitParam;

        public BRgait(BRParameters param, int numberOfPoints)
        {
            _vhc = new BRVHC(param, numberOfPoints);
            _gaitParam = new BRGaitParameters(numberOfPoints);
        }
        public BRVHC vhc
        {
            get
            {
                return _vhc;
            }
            set
            {
                _vhc = value;
            }
        }
        public BRGaitParameters gaitParam
        {
            get
            {
                return _gaitParam;
            }
            set
            {
                _gaitParam = value;
            }
        }

        public double firstIntegral(double theta)
        {

            return _vhc.evalTwoTimesBetaDividedByAlpha(theta);
        }

        public double secondIntegral(double theta)
        {
            double a = MathNet.Numerics.Integration.GaussLegendreRule.Integrate(firstIntegral, _gaitParam.intervalStart, theta, 32);
            return Math.Exp(a) * _vhc.evalTwoTimesGammaDividedByAlpha(theta);
        }

        public double impactFirstLine(double thetaStart, double thetaEnd)
        {
            return _vhc.evalImpactNegFirstLine(thetaEnd)/_vhc.evalImpactPosFirstLine(thetaStart);
        }
        public double impactSecondLine(double thetaStart, double thetaEnd)
        {
            return _vhc.evalImpactNegSecondLine(thetaEnd)/ _vhc.evalImpactPosSecondLine(thetaStart);
        }
        public double impactThirdLine(double thetaStart, double thetaEnd)
        {
            return _vhc.evalImpactNegThirdLine(thetaEnd)/ _vhc.evalImpactPosThirdLine(thetaStart);
        }

    }


    //the performanceindex is to be integrated from theta0 to thetaT 
    public class SQP
    {
        private Expression _performanceIndex;
        private Expression _ineqConstraint1;
        private Expression _ineqConstraint2;
        private Expression _eqConstraint1;
        private Expression _eqConstraint2;

        public SQP(BRVHC vhc)
        {
            
            Expression ddtheta = Expression.Symbol("ddtheta");
            Expression dthetaSquared = Expression.Symbol("dtheta^2");
            Expression theta = Expression.Symbol("theta");

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

            _ineqConstraint2 = alpha3 * ddtheta + beta3 * dthetaSquared + gamma3 -150;

            _performanceIndex = (alpha1 * ddtheta + beta1 * dthetaSquared + gamma1) * vhc.dphi1 + (alpha3 * ddtheta + beta3 * dthetaSquared + gamma3) * vhc.dphi3;

            _eqConstraint1 = vhc.impactNegFirstLine / vhc.impactPosFirstLine - vhc.impactNegSecondLine/vhc.impactPosSecondLine;
            _eqConstraint2 = vhc.impactNegFirstLine / vhc.impactPosFirstLine - vhc.impactNegThirdLine / vhc.impactPosThirdLine;
        }


        public double EvalPerformanceIndex()
        {
            return 0;
        }

        public void run()
        {
            
        }


    }
}
