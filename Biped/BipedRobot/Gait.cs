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
        Dictionary<string, double> _gaitParameters;
        private double _objFunVal;
        private double _intervalStart;
        private double _intervalEnd;

        public BRGaitParameters()
        {
            StreamReader fs = null;
            fs = new StreamReader(@"../../../parameters.txt");
            string temp = fs.ReadLine();
            temp = temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1);
            int tempint = Int32.Parse(temp);
            _gaitParameters = new Dictionary<string, double>();
            for (int i = 0; i < tempint; i++)
            {
                _gaitParameters.Add("P" + i.ToString(), 0.0);
            }

            fs.Close();
        }

        public Dictionary<string, double> gaitparameters
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
    }

    public static class gaitSearch
    {
        public static void run(ref Biped biped)
        {
            //first time running use rand to find a valid gait
            BRgait gait = new BRgait(biped.param);
            Console.WriteLine(gait.impactFirstLine(-0.2,0.2));
            Console.WriteLine(gait.impactSecondLine(-0.2, 0.2));
            Console.WriteLine(gait.impactThirdLine(-0.2, 0.2));
            setPosture(ref gait);
            verifyParameters(gait);
        }
        public static void setPosture(ref BRgait gait)
        {
            gait.gaitParam.intervalStart = -Math.PI / 4;
            gait.gaitParam.intervalEnd = Math.PI / 4;
        }
        public static void setParametersRandom(ref BRgait gait)
        {
            foreach (var pair in gait.gaitParam.gaitparameters)
            {

            }
        }
        public static bool verifyParameters(BRgait gait)
        {
            double thetaDotAtTSquare = (-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.secondIntegral, gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd, 20)) /
                (1 - Math.Exp(-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.firstIntegral, gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd, 20)) 
                * Math.Pow(gait.impactFirstLine(gait.gaitParam.intervalStart, gait.gaitParam.intervalEnd), 20));
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
        private Expression _q3;
        private Expression _ddq3;

        private Expression _alpha;
        private Expression _beta;
        private Expression _gamma;

        private Expression _impactPosFirstLine;
        private Expression _impactPosSecondLine;
        private Expression _impactPosThirdLine;

        private Expression _impactNegFirstLine;
        private Expression _impactNegSecondLine;
        private Expression _impactNegThirdLine;
        public BRVHC()
        {
            StreamReader fs = null;
            fs = new StreamReader(@"../../../q3.txt");
            string temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _q3 = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../ddq3.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _ddq3 = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../alpha.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _alpha = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../beta.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _beta = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../gamma.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _gamma = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactPosFirstLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactPosFirstLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactPosSecondLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactPosSecondLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactPosThirdLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactPosThirdLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactNegFirstLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactNegFirstLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactNegSecondLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactNegSecondLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();

            fs = new StreamReader(@"../../../impactNegThirdLine.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(double)", "");
            _impactNegThirdLine = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            fs.Close();
        }
        public double q1(double theta)
        {
            return 0;
        }
        public double q2(double theta)
        {
            return Math.Cos(theta);
        }
        public double q3(double theta)
        {
            _q3.Parameters["theta"] = theta;
            return (double)_q3.Evaluate();
        }

        public double dq1(double theta)
        {
            return Math.Cos(theta);
        }
        public double dq2(double theta)
        {
            return Math.Cos(theta);
        }
        public double dq3(double theta)
        {
            return 0;
        }

        public double ddq1(double theta)
        {
            return Math.Cos(theta);
        }
        public double ddq2(double theta)
        {
            return Math.Cos(theta);
        }
        public double ddq3(double theta, double dtheta, double ddtheta, Dictionary<string, double> gaitParam, BRParameters physicalParam)
        {
            _ddq3.Parameters["theta"] = theta;
            _ddq3.Parameters["dtheta"] = dtheta;
            _ddq3.Parameters["ddtheta"] = ddtheta;
            foreach (var pair in gaitParam)
            {
                _ddq3.Parameters[pair.Key] = pair.Value;
            }
            return (double)_ddq3.Evaluate();
        }

        public double alpha(double theta, Dictionary<string, double> gaitParam)
        {
            _alpha.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _alpha.Parameters[pair.Key] = pair.Value;
            }
            return (double)_alpha.Evaluate();
        }
        public double beta(double theta, Dictionary<string, double> gaitParam)
        {
            _beta.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _beta.Parameters[pair.Key] = pair.Value;
            }
            return (double)_beta.Evaluate();
        }
        public double gamma(double theta, Dictionary<string, double> gaitParam)
        {
            _gamma.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _gamma.Parameters[pair.Key] = pair.Value;
            }
            return (double)_gamma.Evaluate();
        }

        public double twoTimesBetaDividedByAlpha(double theta, Dictionary<string, double> gaitParam, BRParameters physicalParam)
        {
            _alpha.Parameters["theta"] = theta;
            _beta.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _alpha.Parameters[pair.Key] = pair.Value;
                _beta.Parameters[pair.Key] = pair.Value;
            }
            _alpha.Parameters["g"] = physicalParam.g;
            _alpha.Parameters["m1"] = physicalParam.m1;
            _alpha.Parameters["m2"] = physicalParam.m2;
            _alpha.Parameters["m3"] = physicalParam.m3;

            _alpha.Parameters["l1"] = physicalParam.l1;
            _alpha.Parameters["l2"] = physicalParam.l2;
            _alpha.Parameters["l3"] = physicalParam.l3;

            _alpha.Parameters["L1"] = physicalParam.L1;
            _alpha.Parameters["L2"] = physicalParam.L2;
            _alpha.Parameters["L3"] = physicalParam.L3;

            _alpha.Parameters["J1"] = physicalParam.J1;
            _alpha.Parameters["J2"] = physicalParam.J2;
            _alpha.Parameters["J3"] = physicalParam.J3;

            _beta.Parameters["g"] = physicalParam.g;
            _beta.Parameters["m1"] = physicalParam.m1;
            _beta.Parameters["m2"] = physicalParam.m2;
            _beta.Parameters["m3"] = physicalParam.m3;

            _beta.Parameters["l1"] = physicalParam.l1;
            _beta.Parameters["l2"] = physicalParam.l2;
            _beta.Parameters["l3"] = physicalParam.l3;

            _beta.Parameters["L1"] = physicalParam.L1;
            _beta.Parameters["L2"] = physicalParam.L2;
            _beta.Parameters["L3"] = physicalParam.L3;

            _beta.Parameters["J1"] = physicalParam.J1;
            _beta.Parameters["J2"] = physicalParam.J2;
            _beta.Parameters["J3"] = physicalParam.J3;
            return (2 * (double)_beta.Evaluate() / (double)_alpha.Evaluate());
        }
        public double twoTimesGammaDividedByAlpha(double theta, Dictionary<string, double> gaitParam, BRParameters physicalParam)
        {
            _gamma.Parameters["theta"] = theta;
            _alpha.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _alpha.Parameters[pair.Key] = pair.Value;
                _gamma.Parameters[pair.Key] = pair.Value;
            }
            _alpha.Parameters["g"] = physicalParam.g;
            _alpha.Parameters["m1"] = physicalParam.m1;
            _alpha.Parameters["m2"] = physicalParam.m2;
            _alpha.Parameters["m3"] = physicalParam.m3;

            _alpha.Parameters["l1"] = physicalParam.l1;
            _alpha.Parameters["l2"] = physicalParam.l2;
            _alpha.Parameters["l3"] = physicalParam.l3;

            _alpha.Parameters["L1"] = physicalParam.L1;
            _alpha.Parameters["L2"] = physicalParam.L2;
            _alpha.Parameters["L3"] = physicalParam.L3;

            _alpha.Parameters["J1"] = physicalParam.J1;
            _alpha.Parameters["J2"] = physicalParam.J2;
            _alpha.Parameters["J3"] = physicalParam.J3;

            _gamma.Parameters["g"] = physicalParam.g;
            _gamma.Parameters["m1"] = physicalParam.m1;
            _gamma.Parameters["m2"] = physicalParam.m2;
            _gamma.Parameters["m3"] = physicalParam.m3;

            _gamma.Parameters["l1"] = physicalParam.l1;
            _gamma.Parameters["l2"] = physicalParam.l2;
            _gamma.Parameters["l3"] = physicalParam.l3;

            _gamma.Parameters["L1"] = physicalParam.L1;
            _gamma.Parameters["L2"] = physicalParam.L2;
            _gamma.Parameters["L3"] = physicalParam.L3;

            _gamma.Parameters["J1"] = physicalParam.J1;
            _gamma.Parameters["J2"] = physicalParam.J2;
            _gamma.Parameters["J3"] = physicalParam.J3;
            return (2 * (double)_gamma.Evaluate() / (double)_alpha.Evaluate());
        }

        public double impactPosFirstLine(double theta, Dictionary<string, double> gaitParam, BRParameters physicalParam)
        {
            _impactPosFirstLine.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _impactPosFirstLine.Parameters[pair.Key] = pair.Value;
            }

            _impactPosFirstLine.Parameters["g"] = physicalParam.g;
            _impactPosFirstLine.Parameters["m1"] = physicalParam.m1;
            _impactPosFirstLine.Parameters["m2"] = physicalParam.m2;
            _impactPosFirstLine.Parameters["m3"] = physicalParam.m3;

            _impactPosFirstLine.Parameters["l1"] = physicalParam.l1;
            _impactPosFirstLine.Parameters["l2"] = physicalParam.l2;
            _impactPosFirstLine.Parameters["l3"] = physicalParam.l3;

            _impactPosFirstLine.Parameters["L1"] = physicalParam.L1;
            _impactPosFirstLine.Parameters["L2"] = physicalParam.L2;
            _impactPosFirstLine.Parameters["L3"] = physicalParam.L3;

            _impactPosFirstLine.Parameters["J1"] = physicalParam.J1;
            _impactPosFirstLine.Parameters["J2"] = physicalParam.J2;
            _impactPosFirstLine.Parameters["J3"] = physicalParam.J3;
            return (double)_impactPosFirstLine.Evaluate();
        }
        public double impactPosSecondLine(double theta, Dictionary<string, double> gaitParam, BRParameters physicalParam)
        {
            _impactPosSecondLine.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _impactPosSecondLine.Parameters[pair.Key] = pair.Value;
            }

            _impactPosSecondLine.Parameters["g"] = physicalParam.g;
            _impactPosSecondLine.Parameters["m1"] = physicalParam.m1;
            _impactPosSecondLine.Parameters["m2"] = physicalParam.m2;
            _impactPosSecondLine.Parameters["m3"] = physicalParam.m3;

            _impactPosSecondLine.Parameters["l1"] = physicalParam.l1;
            _impactPosSecondLine.Parameters["l2"] = physicalParam.l2;
            _impactPosSecondLine.Parameters["l3"] = physicalParam.l3;

            _impactPosSecondLine.Parameters["L1"] = physicalParam.L1;
            _impactPosSecondLine.Parameters["L2"] = physicalParam.L2;
            _impactPosSecondLine.Parameters["L3"] = physicalParam.L3;

            _impactPosSecondLine.Parameters["J1"] = physicalParam.J1;
            _impactPosSecondLine.Parameters["J2"] = physicalParam.J2;
            _impactPosSecondLine.Parameters["J3"] = physicalParam.J3;
            return (double)_impactPosSecondLine.Evaluate();
        }
        public double impactPosThirdLine(double theta, Dictionary<string, double> gaitParam, BRParameters physicalParam)
        {
            _impactPosThirdLine.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _impactPosThirdLine.Parameters[pair.Key] = pair.Value;
            }

            _impactPosThirdLine.Parameters["g"] = physicalParam.g;
            _impactPosThirdLine.Parameters["m1"] = physicalParam.m1;
            _impactPosThirdLine.Parameters["m2"] = physicalParam.m2;
            _impactPosThirdLine.Parameters["m3"] = physicalParam.m3;

            _impactPosThirdLine.Parameters["l1"] = physicalParam.l1;
            _impactPosThirdLine.Parameters["l2"] = physicalParam.l2;
            _impactPosThirdLine.Parameters["l3"] = physicalParam.l3;

            _impactPosThirdLine.Parameters["L1"] = physicalParam.L1;
            _impactPosThirdLine.Parameters["L2"] = physicalParam.L2;
            _impactPosThirdLine.Parameters["L3"] = physicalParam.L3;

            _impactPosThirdLine.Parameters["J1"] = physicalParam.J1;
            _impactPosThirdLine.Parameters["J2"] = physicalParam.J2;
            _impactPosThirdLine.Parameters["J3"] = physicalParam.J3;
            return (double)_impactPosThirdLine.Evaluate();
        }

        public double impactNegFirstLine(double theta, Dictionary<string, double> gaitParam, BRParameters physicalParam)
        {
            _impactNegFirstLine.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _impactNegFirstLine.Parameters[pair.Key] = pair.Value;
            }

            _impactNegFirstLine.Parameters["g"] = physicalParam.g;
            _impactNegFirstLine.Parameters["m1"] = physicalParam.m1;
            _impactNegFirstLine.Parameters["m2"] = physicalParam.m2;
            _impactNegFirstLine.Parameters["m3"] = physicalParam.m3;

            _impactNegFirstLine.Parameters["l1"] = physicalParam.l1;
            _impactNegFirstLine.Parameters["l2"] = physicalParam.l2;
            _impactNegFirstLine.Parameters["l3"] = physicalParam.l3;

            _impactNegFirstLine.Parameters["L1"] = physicalParam.L1;
            _impactNegFirstLine.Parameters["L2"] = physicalParam.L2;
            _impactNegFirstLine.Parameters["L3"] = physicalParam.L3;

            _impactNegFirstLine.Parameters["J1"] = physicalParam.J1;
            _impactNegFirstLine.Parameters["J2"] = physicalParam.J2;
            _impactNegFirstLine.Parameters["J3"] = physicalParam.J3;
            return (double)_impactNegFirstLine.Evaluate();
        }
        public double impactNegSecondLine(double theta, Dictionary<string, double> gaitParam, BRParameters physicalParam)
        {
            _impactNegSecondLine.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _impactNegSecondLine.Parameters[pair.Key] = pair.Value;
            }

            _impactNegSecondLine.Parameters["g"] = physicalParam.g;
            _impactNegSecondLine.Parameters["m1"] = physicalParam.m1;
            _impactNegSecondLine.Parameters["m2"] = physicalParam.m2;
            _impactNegSecondLine.Parameters["m3"] = physicalParam.m3;

            _impactNegSecondLine.Parameters["l1"] = physicalParam.l1;
            _impactNegSecondLine.Parameters["l2"] = physicalParam.l2;
            _impactNegSecondLine.Parameters["l3"] = physicalParam.l3;

            _impactNegSecondLine.Parameters["L1"] = physicalParam.L1;
            _impactNegSecondLine.Parameters["L2"] = physicalParam.L2;
            _impactNegSecondLine.Parameters["L3"] = physicalParam.L3;

            _impactNegSecondLine.Parameters["J1"] = physicalParam.J1;
            _impactNegSecondLine.Parameters["J2"] = physicalParam.J2;
            _impactNegSecondLine.Parameters["J3"] = physicalParam.J3;
            return (double)_impactNegSecondLine.Evaluate();
        }
        public double impactNegThirdLine(double theta, Dictionary<string, double> gaitParam, BRParameters physicalParam)
        {
            _impactNegThirdLine.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _impactNegThirdLine.Parameters[pair.Key] = pair.Value;
            }

            _impactNegThirdLine.Parameters["g"] = physicalParam.g;
            _impactNegThirdLine.Parameters["m1"] = physicalParam.m1;
            _impactNegThirdLine.Parameters["m2"] = physicalParam.m2;
            _impactNegThirdLine.Parameters["m3"] = physicalParam.m3;

            _impactNegThirdLine.Parameters["l1"] = physicalParam.l1;
            _impactNegThirdLine.Parameters["l2"] = physicalParam.l2;
            _impactNegThirdLine.Parameters["l3"] = physicalParam.l3;

            _impactNegThirdLine.Parameters["L1"] = physicalParam.L1;
            _impactNegThirdLine.Parameters["L2"] = physicalParam.L2;
            _impactNegThirdLine.Parameters["L3"] = physicalParam.L3;

            _impactNegThirdLine.Parameters["J1"] = physicalParam.J1;
            _impactNegThirdLine.Parameters["J2"] = physicalParam.J2;
            _impactNegThirdLine.Parameters["J3"] = physicalParam.J3;
            return (double)_impactNegThirdLine.Evaluate();
        }
    }

    public class BRgait
    {
        private BRVHC _vhc;
        private BRGaitParameters _gaitParam;
        private BRParameters _physicalParam;

        public BRgait(BRParameters param)
        {
            _vhc = new BRVHC();
            _gaitParam = new BRGaitParameters();
            _physicalParam = param;
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

            return _vhc.twoTimesBetaDividedByAlpha(theta, _gaitParam.gaitparameters, _physicalParam);
        }

        public double secondIntegral(double theta)
        {
            double a = MathNet.Numerics.Integration.GaussLegendreRule.Integrate(firstIntegral, _gaitParam.intervalStart, theta, 2);
            return Math.Exp(a) * _vhc.twoTimesGammaDividedByAlpha(theta, _gaitParam.gaitparameters, _physicalParam);
        }

        public double impactFirstLine(double thetaStart, double thetaEnd)
        {
            return _vhc.impactNegFirstLine(thetaEnd, _gaitParam.gaitparameters, _physicalParam)/_vhc.impactPosFirstLine(thetaStart, _gaitParam.gaitparameters, _physicalParam);
        }
        public double impactSecondLine(double thetaStart, double thetaEnd)
        {
            return _vhc.impactNegSecondLine(thetaEnd, _gaitParam.gaitparameters, _physicalParam)/ _vhc.impactPosSecondLine(thetaStart, _gaitParam.gaitparameters, _physicalParam);
        }
        public double impactThirdLine(double thetaStart, double thetaEnd)
        {
            return _vhc.impactNegThirdLine(thetaEnd, _gaitParam.gaitparameters, _physicalParam)/ _vhc.impactPosThirdLine(thetaStart, _gaitParam.gaitparameters, _physicalParam);
        }

    }
}
