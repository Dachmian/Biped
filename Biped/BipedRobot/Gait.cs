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
            gait.impact(0);
        }
        public static void setPosture(ref BRgait gait)
        {
            gait.gaitParam.intervalStart = -Math.PI / 4;
            gait.gaitParam.intervalStart = Math.PI / 4;
        }
        public static void setParametersRandom(ref BRgait gait)
        {
            foreach (var pair in gait.gaitParam.gaitparameters)
            {

            }
        }
        public static bool verifyParameters(BRgait gait)
        {
            double thetaDotAtTSquare = (-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.secondIntegral, gait.param.intervalStart, gait.param.intervalEnd, 2)) /
                (1 - Math.Exp(-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.firstIntegral, gait.param.intervalStart, gait.param.intervalEnd, 2)) * Math.Pow(gait.impact(gait.param.intervalEnd), 2));
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

        private Expression _impact;

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

            fs = new StreamReader(@"../../../impact.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow").Replace("(System.Double)","");
            _impact = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
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

        public double impact(double theta, Dictionary<string, double> gaitParam, BRParameters physicalParam)
        {
            _impact.Parameters["theta"] = theta;
            foreach (var pair in gaitParam)
            {
                _impact.Parameters[pair.Key] = pair.Value;
            }

            _impact.Parameters["g"] = physicalParam.g;
            _impact.Parameters["m1"] = physicalParam.m1;
            _impact.Parameters["m2"] = physicalParam.m2;
            _impact.Parameters["m3"] = physicalParam.m3;

            _impact.Parameters["l1"] = physicalParam.l1;
            _impact.Parameters["l2"] = physicalParam.l2;
            _impact.Parameters["l3"] = physicalParam.l3;

            _impact.Parameters["L1"] = physicalParam.L1;
            _impact.Parameters["L2"] = physicalParam.L2;
            _impact.Parameters["L3"] = physicalParam.L3;

            _impact.Parameters["J1"] = physicalParam.J1;
            _impact.Parameters["J2"] = physicalParam.J2;
            _impact.Parameters["J3"] = physicalParam.J3;
            return (double)_impact.Evaluate();
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

        public double impact(double theta)
        {
            return _vhc.impact(theta, _gaitParam.gaitparameters, _physicalParam);
        }

    }
}
