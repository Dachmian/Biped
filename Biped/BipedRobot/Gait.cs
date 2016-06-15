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
        private double _startAngle;
        private double _endAngle;
        private double _theta0;
        private double _thetaT;
        private double _dtheta0;
        private double _dthetaT;

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
            
            

            _initialControlPointsq1 = Vector<double>.Build.Dense(numOfPoints);
            _initialControlPointsq2 = Vector<double>.Build.Dense(numOfPoints);
            _initialControlPointsq3 = Vector<double>.Build.Dense(numOfPoints);
            setPosture();
            setInitialControlPoints(numOfPoints);
            setInitialTheta();

            _gaitParameters = new Tuple<Vector<double>, Vector<double>, Vector<double>>(
                _initialControlPointsq1, _initialControlPointsq2, _initialControlPointsq3);
        }

        public void setInitialControlPoints(int numOfPoints)
        {
            double factor = (endAngle - startAngle) / numOfPoints;
            for (int i = 0; i < numOfPoints; i++)
            {
                _initialControlPointsq1[i] = _startAngle + i*factor;
            }


            
            for (int i = 0; i < numOfPoints; i++)
            {
                _initialControlPointsq2[i] = 0;
            }


            for (int i = 0; i < numOfPoints; i++)
            {
                _initialControlPointsq1[i] = _endAngle - i * factor;
            }
        }
        public void setPosture()
        {
            _startAngle = -Math.PI / 10;
            _endAngle = Math.PI / 10;
        }
        public void setInitialTheta()
        {
            _theta0 = 0;
            _thetaT = 1;
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
        public double startAngle
        {
            get
            {
                return _startAngle;
            }
            set
            {
                _startAngle = value;
            }
        }
        public double endAngle
        {
            get
            {
                return _endAngle;
            }
            set
            {
                _endAngle = value;
            }
        }
        public double theta0
        {
            get
            {
                return _theta0;
            }
            set
            {
                _theta0 = value;
            }
        }
        public double thetaT
        {
            get
            {
                return _thetaT;
            }
            set
            {
                _thetaT = value;
            }
        }
        public double dtheta0
        {
            get
            {
                return _dtheta0;
            }
            set
            {
                _dtheta0 = value;
            }
        }
        public double dthetaT
        {
            get
            {
                return _dthetaT;
            }
            set
            {
                _dthetaT = value;
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
            int numberOfPoints = 10;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            /*while (testImpact(gait))
            {
                setParametersRandom(ref gait);
            }*/
            Console.WriteLine("ute av impact");
            setVHC(ref gait, numberOfPoints);
            //begin search

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

            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);
            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            gait.vhc.phi2Parameters = tuple.Item1;

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            gait.vhc.phi3Parameters = tuple.Item1;

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());
        }
        public static bool setAndVerifyParameters(ref BRgait gait)
        {
            double dthetaTSquare = (-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.secondIntegral, gait.gaitParam.theta0, gait.gaitParam.thetaT, 32)) /
                (1 - Math.Exp(-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.firstIntegral, gait.gaitParam.theta0, gait.gaitParam.thetaT, 32))
                * Math.Pow(gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT), 2));
            Console.WriteLine(dthetaTSquare);
            if (dthetaTSquare < 0.1 || dthetaTSquare > 10000 || dthetaTSquare == double.NaN)
            {
                return false;
            }
            else
            {
                gait.gaitParam.dthetaT = Math.Sqrt(dthetaTSquare);
                gait.gaitParam.dtheta0 = gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT) * Math.Sqrt(dthetaTSquare);
                return true;
            }
        }
        public static bool verifyImpact(BRgait gait)
        {
            double firstAndSecond = gait.impactFirstLine(0,1) - 
                gait.impactSecondLine(0,1); 
            double firstAndThird = gait.impactFirstLine(0,1) - 
                gait.impactThirdLine(0,1);
            if ((firstAndSecond <= 0.001 && firstAndSecond >= -0.001) && (firstAndThird <= 0.001 && firstAndThird >= -0.001)) 
            {
                Console.WriteLine(gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
                Console.WriteLine(gait.impactSecondLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
                Console.WriteLine(gait.impactThirdLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
                return true;
            }
            return false;
        }

        public static bool testImpact(BRgait gait)
        {
            int numberOfPoints = 6;
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);
            Vector<double> first = gait.gaitParam.gaitparameters.Item1;
            Vector<double> second = gait.gaitParam.gaitparameters.Item2;
            Vector<double> third = gait.gaitParam.gaitparameters.Item3;
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);
            tuple.Item1["P0"] =first[0];
            tuple.Item1["P1"] = first[1];
            tuple.Item1["P2"] = first[2];
            tuple.Item1["P3"] = first[3];
            tuple.Item1["P4"] = first[4];
            tuple.Item1["P5"] = first[5];
            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            tuple.Item1["P6"] = second[0];
            tuple.Item1["P7"] = second[1];
            tuple.Item1["P8"] = second[2];
            tuple.Item1["P9"] = second[3];
            tuple.Item1["P10"] = second[4];
            tuple.Item1["P11"] = second[5];
            gait.vhc.phi2Parameters = tuple.Item1;

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            tuple.Item1["P12"] = third[0];
            tuple.Item1["P13"] = third[1];
            tuple.Item1["P14"] = third[2];
            tuple.Item1["P15"] = third[3];
            tuple.Item1["P16"] = third[4];
            tuple.Item1["P17"] = third[5];
            gait.vhc.phi3Parameters = tuple.Item1;

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());
            if( (gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT)>0) && (gait.impactSecondLine(gait.gaitParam.theta0, gait.gaitParam.thetaT)>0) && (gait.impactThirdLine(gait.gaitParam.theta0, gait.gaitParam.thetaT) > 0))
            {
                return false;
            }
            return true;
        }
    }

    public class BRVHC
    {
        public BRVHC(BRParameters param, int numberOfPoints)
        {
            _parameters = new Dictionary<string, FloatingPoint>();
            _phi1Parameters = new Dictionary<string, FloatingPoint>();
            _phi2Parameters = new Dictionary<string, FloatingPoint>();
            _phi3Parameters = new Dictionary<string, FloatingPoint>();

            StreamReader fs = null;
            string temp;

            fs = new StreamReader(@"../../../../alpha.txt");
            temp = fs.ReadLine();
            //temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _alpha = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();


            fs = new StreamReader(@"../../../../gamma.txt");
            temp = fs.ReadLine();
            //temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _gamma = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../../beta.txt");
            temp = fs.ReadLine();
            _beta = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../../alpha1.txt");
            temp = fs.ReadLine();
            _alpha1 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../../beta1.txt");
            temp = fs.ReadLine();
            _beta1 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../../gamma1.txt");
            temp = fs.ReadLine();
            _gamma1 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../../alpha3.txt");
            temp = fs.ReadLine();
            _alpha3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../../beta3.txt");
            temp = fs.ReadLine();
            _beta3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../../gamma3.txt");
            temp = fs.ReadLine();
            _gamma3 = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../../impactFirstLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("0.2e1","2").Replace("pow","");
            temp = temp.Replace(", 2)", ")^2");
            _impactFirstLine = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../../impactSecondLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("0.2e1", "2").Replace("pow", "");
            temp = temp.Replace(", 2)", ")^2");
            _impactSecondLine = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../../impactThirdLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("0.2e1", "2").Replace("pow", "");
            temp = temp.Replace(", 2)", ")^2");
            _impactThirdLine = Infix.ParseOrUndefined(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
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
            _phi1Parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_phi1Parameters, _phi1).RealValue;
        }
        public double evalPhi2(double theta)
        {
            _phi2Parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_phi2Parameters, _phi2).RealValue;
        }
        public double evalPhi3(double theta)
        {
            _phi3Parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_phi3Parameters, _phi3).RealValue;
        }

        public double evalDphi1(double theta)
        {
            _phi1Parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_phi1Parameters, _dphi1).RealValue;
        }
        public double evalDphi2(double theta)
        {
            _phi2Parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_phi2Parameters, _dphi2).RealValue;
        }
        public double evalDphi3(double theta)
        {
            _phi3Parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_phi3Parameters, _dphi3).RealValue;
        }

        public double evalDdphi1(double theta)
        {
            _phi1Parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_phi1Parameters, _ddphi1).RealValue;
        }
        public double evalDdphi2(double theta)
        {
            _phi2Parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_phi2Parameters, _ddphi2).RealValue;
        }
        public double evalDdphi3(double theta)
        {
            _phi3Parameters["theta"] = theta;
            return (double)MathNet.Symbolics.Evaluate.Evaluate(_phi3Parameters, _ddphi3).RealValue;
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

        private Expression _alpha1;
        private Expression _beta1;
        private Expression _gamma1;

        private Expression _alpha3;
        private Expression _beta3;
        private Expression _gamma3;

        public double evalAlpha(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

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
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

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
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _gamma).RealValue;
        }

        public double evalAlpha1(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _alpha1).RealValue;
        }
        public double evalBeta1(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _beta1).RealValue;
        }
        public double evalGamma1(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _gamma1).RealValue;
        }

        public double evalAlpha3(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _alpha3).RealValue;
        }
        public double evalBeta3(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _beta3).RealValue;
        }
        public double evalGamma3(double theta)
        {
            _parameters["theta"] = theta;
            _parameters["phi1"] = evalPhi1(theta);
            _parameters["phi2"] = evalPhi2(theta);
            _parameters["phi3"] = evalPhi3(theta);
            _parameters["dphi1"] = evalDphi1(theta);
            _parameters["dphi2"] = evalDphi2(theta);
            _parameters["dphi3"] = evalDphi3(theta);
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _gamma3).RealValue;
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
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

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
            _parameters["ddphi1"] = evalDdphi1(theta);
            _parameters["ddphi2"] = evalDdphi2(theta);
            _parameters["ddphi3"] = evalDdphi3(theta);

            return (2 * (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _gamma).RealValue / (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _alpha).RealValue);
        }

        public Expression alpha
        {
            get
            {
                return _alpha;
            }
            set
            {
                _alpha = value;
            }
        }
        public Expression beta
        {
            get
            {
                return _beta;
            }
            set
            {
                _beta = value;
            }
        }
        public Expression gamma
        {
            get
            {
                return _gamma;
            }
            set
            {
                _gamma = value;
            }
        }

        public Expression alpha1
        {
            get
            {
                return _alpha1;
            }
            set
            {
                _alpha1 = value;
            }
        }
        public Expression beta1
        {
            get
            {
                return _beta1;
            }
            set
            {
                _beta1 = value;
            }
        }
        public Expression gamma1
        {
            get
            {
                return _gamma1;
            }
            set
            {
                _gamma1 = value;
            }
        }

        public Expression alpha3
        {
            get
            {
                return _alpha3;
            }
            set
            {
                _alpha3 = value;
            }
        }
        public Expression beta3
        {
            get
            {
                return _beta3;
            }
            set
            {
                _beta3 = value;
            }
        }
        public Expression gamma3
        {
            get
            {
                return _gamma3;
            }
            set
            {
                _gamma3 = value;
            }
        }
        #endregion
        #region impacts
        private Expression _impactFirstLine;
        private Expression _impactSecondLine;
        private Expression _impactThirdLine;


        public double evalImpactFirstLine(double theta)
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

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _impactFirstLine).RealValue;
        }
        public double evalImpactSecondLine(double theta)
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

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _impactSecondLine).RealValue;
        }
        public double evalImpactThirdLine(double theta)
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

            return (double)MathNet.Symbolics.Evaluate.Evaluate(_parameters, _impactThirdLine).RealValue;
        }

        public Expression impactFirstLine
        {
            get
            {
                return _impactFirstLine;
            }
            set
            {
                _impactFirstLine = value;
            }
        }
        public Expression impactSecondLine
        {
            get
            {
                return _impactSecondLine;
            }
            set
            {
                _impactSecondLine = value;
            }
        }
        public Expression impactThirdLine
        {
            get
            {
                return _impactThirdLine;
            }
            set
            {
                _impactThirdLine = value;
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
        private BRReducedSimulationData _data;

        public BRgait(BRParameters param, int numberOfPoints)
        {
            _gaitParam = new BRGaitParameters(numberOfPoints);
            _vhc = new BRVHC(param, numberOfPoints);
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
        public BRReducedSimulationData data
        {
            get
            {
                return _data;
            }
            set
            {
                _data = value;
            }
        }

        public double firstIntegral(double theta)
        {

            return _vhc.evalTwoTimesBetaDividedByAlpha(theta);
        }

        public double secondIntegral(double theta)
        {
            double a = MathNet.Numerics.Integration.GaussLegendreRule.Integrate(firstIntegral, 1, theta, 100);
            return Math.Exp(a) * _vhc.evalTwoTimesGammaDividedByAlpha(theta);
        }
        public double secondIntegralPartial(double theta)
        {
            return _vhc.evalTwoTimesGammaDividedByAlpha(theta);
        }

        public double impactFirstLine(double thetaStart, double thetaEnd)
        {
            return _vhc.evalImpactFirstLine(thetaEnd)/_vhc.evalDphi1(thetaStart);
        }
        public double impactSecondLine(double thetaStart, double thetaEnd)
        {
            return _vhc.evalImpactSecondLine(thetaEnd)/_vhc.evalDphi2(thetaStart);
        }
        public double impactThirdLine(double thetaStart, double thetaEnd)
        {
            return _vhc.evalImpactThirdLine(thetaEnd)/_vhc.evalDphi3(thetaStart);
        }

    }

    

}
