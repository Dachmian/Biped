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
using NCalc;

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
            _initialControlPointsq1[1] = -Math.PI/8;
            _initialControlPointsq1[2] = 0;
            _initialControlPointsq1[3] = Math.PI/8;
            _initialControlPointsq1[4] = _intervalEnd;


            _initialControlPointsq2[0] = 0;
            _initialControlPointsq2[1] = Math.PI/9;
            _initialControlPointsq2[2] = 0;
            _initialControlPointsq2[3] = -Math.PI/9;
            _initialControlPointsq2[4] = 0;


            _initialControlPointsq3[0] = _intervalEnd;
            _initialControlPointsq3[1] = Math.PI/8;
            _initialControlPointsq3[2] = 0;
            _initialControlPointsq3[3] = -Math.PI/8;
            _initialControlPointsq3[4] = _intervalStart;
        }
        public void setPosture()
        {
            _intervalStart = -Math.PI / 4;
            _intervalEnd = Math.PI / 4;
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
            for (int i = 0; i < 10; i++)
            {
                setParametersRandom(ref gait);
                setVHC(ref gait, numberOfPoints);
                verifyParameters(gait);
            }
        }
        public static void setPosture(ref BRgait gait)
        {
            
        }
        public static void setParametersRandom(ref BRgait gait)
        {
            Random rndm = new Random();
            gait.gaitParam.gaitparameters.Item1[0] = gait.gaitParam.initialControlPointsq1[0];
            for (int i = 1; i < gait.gaitParam.gaitparameters.Item1.Count; i++)
            {
                gait.gaitParam.gaitparameters.Item1[i] = gait.gaitParam.initialControlPointsq1[i] + 2 * (rndm.NextDouble() - 0.5);
                
                    
            }
            gait.gaitParam.gaitparameters.Item2[0] = gait.gaitParam.initialControlPointsq2[0];
            for (int i = 1; i < gait.gaitParam.gaitparameters.Item2.Count; i++)
            {
                gait.gaitParam.gaitparameters.Item2[i] = gait.gaitParam.initialControlPointsq2[i] + 2 * (rndm.NextDouble() - 0.5);


            }
            gait.gaitParam.gaitparameters.Item3[0] = gait.gaitParam.initialControlPointsq3[0];
            for (int i = 1; i < gait.gaitParam.gaitparameters.Item3.Count; i++)
            {
                gait.gaitParam.gaitparameters.Item3[i] = gait.gaitParam.initialControlPointsq3[i] + 2 * (rndm.NextDouble() - 0.5);


            }
        }
        public static void setVHC(ref BRgait gait, int numberOfPoints)
        {
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);
            gait.vhc.phi1 = new Expression(brCrv.phi1ToString());
            gait.vhc.phi2 = new Expression(brCrv.phi2ToString());
            gait.vhc.phi3 = new Expression(brCrv.phi3ToString());

            gait.vhc.dphi1 = new Expression(brCrv.dphi1ToString());
            gait.vhc.dphi2 = new Expression(brCrv.dphi2ToString());
            gait.vhc.dphi3 = new Expression(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = new Expression(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = new Expression(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = new Expression(brCrv.ddphi3ToString());
        }
        public static bool verifyParameters(BRgait gait)
        {
            double thetaDotAtTSquare = (-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.secondIntegral, gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd, 20)) /
                (1 - Math.Exp(-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.firstIntegral, gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd, 20)) 
                * Math.Pow(gait.impactFirstLine(gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd), 2));
            Console.WriteLine(thetaDotAtTSquare);
            if (thetaDotAtTSquare < 0)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }

    public class BRVHC
    {
        public BRVHC(BRParameters param)
        {
            StreamReader fs = null;
            fs = new StreamReader(@"../../../q3.txt");
            string temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _q3 = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _q3);
            fs.Close();

            fs = new StreamReader(@"../../../ddq3.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _ddq3 = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _ddq3);
            fs.Close();

            fs = new StreamReader(@"../../../alpha.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _alpha = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _alpha);
            fs.Close();

            fs = new StreamReader(@"../../../beta.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _beta = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _beta);
            fs.Close();

            fs = new StreamReader(@"../../../gamma.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _gamma = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _gamma);
            fs.Close();

            fs = new StreamReader(@"../../../impactPosFirstLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactPosFirstLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _impactPosFirstLine);
            fs.Close();

            fs = new StreamReader(@"../../../impactPosSecondLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactPosSecondLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _impactPosSecondLine);
            fs.Close();

            fs = new StreamReader(@"../../../impactPosThirdLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactPosThirdLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _impactPosThirdLine);
            fs.Close();

            fs = new StreamReader(@"../../../impactNegFirstLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactNegFirstLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _impactNegFirstLine);
            fs.Close();

            fs = new StreamReader(@"../../../impactNegSecondLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactNegSecondLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _impactNegSecondLine);
            fs.Close();

            fs = new StreamReader(@"../../../impactNegThirdLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactNegThirdLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            setPhysicalParameters(param, ref _impactNegThirdLine);
            fs.Close();
        }
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
            _q3.Parameters["theta"] = theta;
            return (double)_q3.Evaluate();
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
            _ddq3.Parameters["theta"] = theta;
            _ddq3.Parameters["dtheta"] = dtheta;
            _ddq3.Parameters["ddtheta"] = ddtheta;

            return (double)_ddq3.Evaluate();
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
            _phi1.Parameters["theta"] = theta;
            return (double)_phi1.Evaluate();
        }
        public double evalPhi2(double theta)
        {
            _phi2.Parameters["theta"] = theta;
            return (double)_phi2.Evaluate();
        }
        public double evalPhi3(double theta)
        {
            _phi3.Parameters["theta"] = theta;
            return (double)_phi3.Evaluate();
        }

        public double evalDphi1(double theta)
        {
            _dphi1.Parameters["theta"] = theta;
            return (double)_dphi1.Evaluate();
        }
        public double evalDphi2(double theta)
        {
            _dphi2.Parameters["theta"] = theta;
            return (double)_dphi2.Evaluate();
        }
        public double evalDphi3(double theta)
        {
            _dphi3.Parameters["theta"] = theta;
            return (double)_dphi3.Evaluate();
        }

        public double evalDdphi1(double theta)
        {
            _ddphi1.Parameters["theta"] = theta;
            return (double)_ddphi1.Evaluate();
        }
        public double evalDdphi2(double theta)
        {
            _ddphi2.Parameters["theta"] = theta;
            return (double)_ddphi2.Evaluate();
        }
        public double evalDdphi3(double theta)
        {
            _ddphi3.Parameters["theta"] = theta;
            return (double)_ddphi3.Evaluate();
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
            _alpha.Parameters["theta"] = theta;
            _alpha.Parameters["phi1"] = evalPhi1(theta);
            _alpha.Parameters["phi2"] = evalPhi2(theta);
            _alpha.Parameters["phi3"] = evalPhi3(theta);
            _alpha.Parameters["dphi1"] = evalDphi1(theta);
            _alpha.Parameters["dphi2"] = evalDphi2(theta);
            _alpha.Parameters["dphi3"] = evalDphi3(theta);
            _alpha.Parameters["ddphi1"] = evalDphi1(theta);
            _alpha.Parameters["ddphi2"] = evalDphi2(theta);
            _alpha.Parameters["ddphi3"] = evalDphi3(theta);

            return (double)_alpha.Evaluate();
        }
        public double evalBeta(double theta)
        {
            _beta.Parameters["theta"] = theta;
            _beta.Parameters["phi1"] = evalPhi1(theta);
            _beta.Parameters["phi2"] = evalPhi2(theta);
            _beta.Parameters["phi3"] = evalPhi3(theta);
            _beta.Parameters["dphi1"] = evalDphi1(theta);
            _beta.Parameters["dphi2"] = evalDphi2(theta);
            _beta.Parameters["dphi3"] = evalDphi3(theta);
            _beta.Parameters["ddphi1"] = evalDphi1(theta);
            _beta.Parameters["ddphi2"] = evalDphi2(theta);
            _beta.Parameters["ddphi3"] = evalDphi3(theta);

            return (double)_beta.Evaluate();
        }
        public double evalGamma(double theta)
        {
            _gamma.Parameters["theta"] = theta;
            _gamma.Parameters["phi1"] = evalPhi1(theta);
            _gamma.Parameters["phi2"] = evalPhi2(theta);
            _gamma.Parameters["phi3"] = evalPhi3(theta);
            _gamma.Parameters["dphi1"] = evalDphi1(theta);
            _gamma.Parameters["dphi2"] = evalDphi2(theta);
            _gamma.Parameters["dphi3"] = evalDphi3(theta);
            _gamma.Parameters["ddphi1"] = evalDphi1(theta);
            _gamma.Parameters["ddphi2"] = evalDphi2(theta);
            _gamma.Parameters["ddphi3"] = evalDphi3(theta);

            return (double)_gamma.Evaluate();
        }

        public double evalTwoTimesBetaDividedByAlpha(double theta)
        {
            _alpha.Parameters["theta"] = theta;
            _alpha.Parameters["phi1"] = evalPhi1(theta);
            _alpha.Parameters["phi2"] = evalPhi2(theta);
            _alpha.Parameters["phi3"] = evalPhi3(theta);
            _alpha.Parameters["dphi1"] = evalDphi1(theta);
            _alpha.Parameters["dphi2"] = evalDphi2(theta);
            _alpha.Parameters["dphi3"] = evalDphi3(theta);
            _alpha.Parameters["ddphi1"] = evalDphi1(theta);
            _alpha.Parameters["ddphi2"] = evalDphi2(theta);
            _alpha.Parameters["ddphi3"] = evalDphi3(theta);

            _beta.Parameters["theta"] = theta;
            _beta.Parameters["phi1"] = evalPhi1(theta);
            _beta.Parameters["phi2"] = evalPhi2(theta);
            _beta.Parameters["phi3"] = evalPhi3(theta);
            _beta.Parameters["dphi1"] = evalDphi1(theta);
            _beta.Parameters["dphi2"] = evalDphi2(theta);
            _beta.Parameters["dphi3"] = evalDphi3(theta);
            _beta.Parameters["ddphi1"] = evalDphi1(theta);
            _beta.Parameters["ddphi2"] = evalDphi2(theta);
            _beta.Parameters["ddphi3"] = evalDphi3(theta);

            return (2 * (double)_beta.Evaluate() / (double)_alpha.Evaluate());
        }
        public double evalTwoTimesGammaDividedByAlpha(double theta)
        {
            _gamma.Parameters["theta"] = theta;
            _gamma.Parameters["phi1"] = evalPhi1(theta);
            _gamma.Parameters["phi2"] = evalPhi2(theta);
            _gamma.Parameters["phi3"] = evalPhi3(theta);
            _gamma.Parameters["dphi1"] = evalDphi1(theta);
            _gamma.Parameters["dphi2"] = evalDphi2(theta);
            _gamma.Parameters["dphi3"] = evalDphi3(theta);
            _gamma.Parameters["ddphi1"] = evalDphi1(theta);
            _gamma.Parameters["ddphi2"] = evalDphi2(theta);
            _gamma.Parameters["ddphi3"] = evalDphi3(theta);

            _alpha.Parameters["theta"] = theta;
            _alpha.Parameters["phi1"] = evalPhi1(theta);
            _alpha.Parameters["phi2"] = evalPhi2(theta);
            _alpha.Parameters["phi3"] = evalPhi3(theta);
            _alpha.Parameters["dphi1"] = evalDphi1(theta);
            _alpha.Parameters["dphi2"] = evalDphi2(theta);
            _alpha.Parameters["dphi3"] = evalDphi3(theta);
            _alpha.Parameters["ddphi1"] = evalDphi1(theta);
            _alpha.Parameters["ddphi2"] = evalDphi2(theta);
            _alpha.Parameters["ddphi3"] = evalDphi3(theta);



            return (2 * (double)_gamma.Evaluate() / (double)_alpha.Evaluate());
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
            _impactPosFirstLine.Parameters["theta"] = theta;
            _impactPosFirstLine.Parameters["phi1"] = evalPhi1(theta);
            _impactPosFirstLine.Parameters["phi2"] = evalPhi2(theta);
            _impactPosFirstLine.Parameters["phi3"] = evalPhi3(theta);
            _impactPosFirstLine.Parameters["dphi1"] = evalDphi1(theta);
            _impactPosFirstLine.Parameters["dphi2"] = evalDphi2(theta);
            _impactPosFirstLine.Parameters["dphi3"] = evalDphi3(theta);
            _impactPosFirstLine.Parameters["ddphi1"] = evalDphi1(theta);
            _impactPosFirstLine.Parameters["ddphi2"] = evalDphi2(theta);
            _impactPosFirstLine.Parameters["ddphi3"] = evalDphi3(theta);

            return (double)_impactPosFirstLine.Evaluate();
        }
        public double evalImpactPosSecondLine(double theta)
        {
            _impactPosSecondLine.Parameters["theta"] = theta;
            _impactPosSecondLine.Parameters["phi1"] = evalPhi1(theta);
            _impactPosSecondLine.Parameters["phi2"] = evalPhi2(theta);
            _impactPosSecondLine.Parameters["phi3"] = evalPhi3(theta);
            _impactPosSecondLine.Parameters["dphi1"] = evalDphi1(theta);
            _impactPosSecondLine.Parameters["dphi2"] = evalDphi2(theta);
            _impactPosSecondLine.Parameters["dphi3"] = evalDphi3(theta);
            _impactPosSecondLine.Parameters["ddphi1"] = evalDphi1(theta);
            _impactPosSecondLine.Parameters["ddphi2"] = evalDphi2(theta);
            _impactPosSecondLine.Parameters["ddphi3"] = evalDphi3(theta);
            
            return (double)_impactPosSecondLine.Evaluate();
        }
        public double evalImpactPosThirdLine(double theta)
        {
            _impactPosThirdLine.Parameters["theta"] = theta;
            _impactPosThirdLine.Parameters["phi1"] = evalPhi1(theta);
            _impactPosThirdLine.Parameters["phi2"] = evalPhi2(theta);
            _impactPosThirdLine.Parameters["phi3"] = evalPhi3(theta);
            _impactPosThirdLine.Parameters["dphi1"] = evalDphi1(theta);
            _impactPosThirdLine.Parameters["dphi2"] = evalDphi2(theta);
            _impactPosThirdLine.Parameters["dphi3"] = evalDphi3(theta);
            _impactPosThirdLine.Parameters["ddphi1"] = evalDphi1(theta);
            _impactPosThirdLine.Parameters["ddphi2"] = evalDphi2(theta);
            _impactPosThirdLine.Parameters["ddphi3"] = evalDphi3(theta);
            return (double)_impactPosThirdLine.Evaluate();
        }

        public double evalImpactNegFirstLine(double theta)
        {
            _impactNegFirstLine.Parameters["theta"] = theta;
            _impactNegFirstLine.Parameters["phi1"] = evalPhi1(theta);
            _impactNegFirstLine.Parameters["phi2"] = evalPhi2(theta);
            _impactNegFirstLine.Parameters["phi3"] = evalPhi3(theta);
            _impactNegFirstLine.Parameters["dphi1"] = evalDphi1(theta);
            _impactNegFirstLine.Parameters["dphi2"] = evalDphi2(theta);
            _impactNegFirstLine.Parameters["dphi3"] = evalDphi3(theta);
            _impactNegFirstLine.Parameters["ddphi1"] = evalDphi1(theta);
            _impactNegFirstLine.Parameters["ddphi2"] = evalDphi2(theta);
            _impactNegFirstLine.Parameters["ddphi3"] = evalDphi3(theta);
            
            return (double)_impactNegFirstLine.Evaluate();
        }
        public double evalImpactNegSecondLine(double theta)
        {
            _impactNegSecondLine.Parameters["theta"] = theta;
            _impactNegSecondLine.Parameters["phi1"] = evalPhi1(theta);
            _impactNegSecondLine.Parameters["phi2"] = evalPhi2(theta);
            _impactNegSecondLine.Parameters["phi3"] = evalPhi3(theta);
            _impactNegSecondLine.Parameters["dphi1"] = evalDphi1(theta);
            _impactNegSecondLine.Parameters["dphi2"] = evalDphi2(theta);
            _impactNegSecondLine.Parameters["dphi3"] = evalDphi3(theta);
            _impactNegSecondLine.Parameters["ddphi1"] = evalDphi1(theta);
            _impactNegSecondLine.Parameters["ddphi2"] = evalDphi2(theta);
            _impactNegSecondLine.Parameters["ddphi3"] = evalDphi3(theta);
            return (double)_impactNegSecondLine.Evaluate();
        }
        public double evalImpactNegThirdLine(double theta)
        {
            _impactNegThirdLine.Parameters["theta"] = theta;
            _impactNegThirdLine.Parameters["phi1"] = evalPhi1(theta);
            _impactNegThirdLine.Parameters["phi2"] = evalPhi2(theta);
            _impactNegThirdLine.Parameters["phi3"] = evalPhi3(theta);
            _impactNegThirdLine.Parameters["dphi1"] = evalDphi1(theta);
            _impactNegThirdLine.Parameters["dphi2"] = evalDphi2(theta);
            _impactNegThirdLine.Parameters["dphi3"] = evalDphi3(theta);
            _impactNegThirdLine.Parameters["ddphi1"] = evalDphi1(theta);
            _impactNegThirdLine.Parameters["ddphi2"] = evalDphi2(theta);
            _impactNegThirdLine.Parameters["ddphi3"] = evalDphi3(theta);
            return (double)_impactNegThirdLine.Evaluate();
        }
        #endregion
        public void setPhysicalParameters(BRParameters param, ref Expression exp)
        {
            exp.Parameters["g"] = param.g;
            exp.Parameters["m1"] = param.m1;
            exp.Parameters["m2"] = param.m2;
            exp.Parameters["m3"] = param.m3;

            exp.Parameters["l1"] = param.l1;
            exp.Parameters["l2"] = param.l2;
            exp.Parameters["l3"] = param.l3;

            exp.Parameters["L1"] = param.L1;
            exp.Parameters["L2"] = param.L2;
            exp.Parameters["L3"] = param.L3;

            exp.Parameters["J1"] = param.J1;
            exp.Parameters["J2"] = param.J2;
            exp.Parameters["J3"] = param.J3;
        }


    }

  
    public class BRgait
    {
        private BRVHC _vhc;
        private BRGaitParameters _gaitParam;

        public BRgait(BRParameters param, int numberOfPoints)
        {
            _vhc = new BRVHC(param);
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
            double a = MathNet.Numerics.Integration.GaussLegendreRule.Integrate(firstIntegral, _gaitParam.intervalStart, theta, 2);
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
}
