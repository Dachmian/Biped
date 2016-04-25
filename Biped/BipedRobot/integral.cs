using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BipedRobot
{
    public class GaussianQuadrature
    {
        private double[] _weigths;
        private double[] _abscissae;
        private double[,] _THETA;
        public delegate double function(double theta);

        public double[] weigths
        {
            get
            {
                return _weigths;
            }
            set
            {
                _weigths = value;
            }
        }
        public double[] abscissae
        {
            get
            {
                return _abscissae;
            }
            set
            {
                _abscissae = value;
            }
        }
        public double[,] THETA
        {
            get
            {
                return _THETA;
            }
            set
            {
                _THETA = value;
            }
        }

        public GaussianQuadrature()
        {
            _weigths = new double[64];
            _abscissae = new double[64];

            _weigths[0] = 0.0486909570091397;
            _weigths[1] = 0.0486909570091397;
            _weigths[2] = 0.0485754674415034;
            _weigths[3] = 0.0485754674415034;
            _weigths[4] = 0.0483447622348030;
            _weigths[5] = 0.0483447622348030;
            _weigths[6] = 0.0479993885964583;
            _weigths[7] = 0.0479993885964583;
            _weigths[8] = 0.0475401657148303;
            _weigths[9] = 0.0475401657148303;

            _weigths[10] = 0.0469681828162100;
            _weigths[11] = 0.0469681828162100;
            _weigths[12] = 0.0462847965813144;
            _weigths[13] = 0.0462847965813144;
            _weigths[14] = 0.0454916279274181;
            _weigths[15] = 0.0454916279274181;
            _weigths[16] = 0.0445905581637566;
            _weigths[17] = 0.0445905581637566;
            _weigths[18] = 0.0435837245293235;
            _weigths[19] = 0.0435837245293235;

            _weigths[20] = 0.0424735151236536;
            _weigths[21] = 0.0424735151236536;
            _weigths[22] = 0.0412625632426235;
            _weigths[23] = 0.0412625632426235;
            _weigths[24] = 0.0399537411327203;
            _weigths[25] = 0.0399537411327203;
            _weigths[26] = 0.0385501531786156;
            _weigths[27] = 0.0385501531786156;
            _weigths[28] = 0.0370551285402400;
            _weigths[29] = 0.0370551285402400;

            _weigths[30] = 0.0354722132568824;
            _weigths[31] = 0.0354722132568824;
            _weigths[32] = 0.0338051618371416;
            _weigths[33] = 0.0338051618371416;
            _weigths[34] = 0.0320579283548516;
            _weigths[35] = 0.0320579283548516;
            _weigths[36] = 0.0302346570724025;
            _weigths[37] = 0.0302346570724025;
            _weigths[38] = 0.0283396726142595;
            _weigths[39] = 0.0283396726142595;

            _weigths[40] = 0.0263774697150547;
            _weigths[41] = 0.0263774697150547;
            _weigths[42] = 0.0243527025687109;
            _weigths[43] = 0.0243527025687109;
            _weigths[44] = 0.0222701738083833;
            _weigths[45] = 0.0222701738083833;
            _weigths[46] = 0.0201348231535302;
            _weigths[47] = 0.0201348231535302;
            _weigths[48] = 0.0179517157756973;
            _weigths[49] = 0.0179517157756973;

            _weigths[50] = 0.0157260304760247;
            _weigths[51] = 0.0157260304760247;
            _weigths[52] = 0.0134630478967186;
            _weigths[53] = 0.0134630478967186;
            _weigths[54] = 0.0111681394601311;
            _weigths[55] = 0.0111681394601311;
            _weigths[56] = 0.0088467598263639;
            _weigths[57] = 0.0088467598263639;
            _weigths[58] = 0.0065044579689784;
            _weigths[59] = 0.0065044579689784;

            _weigths[60] = 0.0041470332605625;
            _weigths[61] = 0.0041470332605625;
            _weigths[62] = 0.0017832807216964;
            _weigths[63] = 0.0017832807216964;

            _abscissae[0] = -0.0243502926634244;
            _abscissae[1] = 0.0243502926634244;
            _abscissae[2] = -0.0729931217877990;
            _abscissae[3] = 0.0729931217877990;
            _abscissae[4] = -0.1214628192961206;
            _abscissae[5] = 0.1214628192961206;
            _abscissae[6] = -0.1696444204239928;
            _abscissae[7] = 0.1696444204239928;
            _abscissae[8] = -0.2174236437400071;
            _abscissae[9] = 0.2174236437400071;

            _abscissae[10] = -0.2646871622087674;
            _abscissae[11] = 0.2646871622087674;
            _abscissae[12] = -0.3113228719902110;
            _abscissae[13] = 0.3113228719902110;
            _abscissae[14] = -0.3572201583376681;
            _abscissae[15] = 0.3572201583376681;
            _abscissae[16] = -0.4022701579639916;
            _abscissae[17] = 0.4022701579639916;
            _abscissae[18] = -0.4463660172534641;
            _abscissae[19] = 0.4463660172534641;

            _abscissae[20] = -0.4894031457070530;
            _abscissae[21] = 0.4894031457070530;
            _abscissae[22] = -0.5312794640198946;
            _abscissae[23] = 0.5312794640198946;
            _abscissae[24] = -0.5718956462026340;
            _abscissae[25] = 0.5718956462026340;
            _abscissae[26] = -0.6111553551723933;
            _abscissae[27] = 0.6111553551723933;
            _abscissae[28] = -0.6489654712546573;
            _abscissae[29] = 0.6489654712546573;

            _abscissae[30] = -0.6852363130542333;
            _abscissae[31] = 0.6852363130542333;
            _abscissae[32] = -0.7198818501716109;
            _abscissae[33] = 0.7198818501716109;
            _abscissae[34] = -0.7528199072605319;
            _abscissae[35] = 0.7528199072605319;
            _abscissae[36] = -0.7839723589433414;
            _abscissae[37] = 0.7839723589433414;
            _abscissae[38] = -0.8132653151227975;
            _abscissae[39] = 0.8132653151227975;

            _abscissae[40] = -0.8406292962525803;
            _abscissae[41] = 0.8406292962525803;
            _abscissae[42] = -0.8659993981540928;
            _abscissae[43] = 0.8659993981540928;
            _abscissae[44] = -0.8893154459951141;
            _abscissae[45] = 0.8893154459951141;
            _abscissae[46] = -0.9105221370785028;
            _abscissae[47] = 0.9105221370785028;
            _abscissae[48] = -0.9295691721319396;
            _abscissae[49] = 0.9295691721319396;

            _abscissae[50] = -0.9464113748584028;
            _abscissae[51] = 0.9464113748584028;
            _abscissae[52] = -0.9610087996520538;
            _abscissae[53] = 0.9610087996520538;
            _abscissae[54] = -0.9733268277899110;
            _abscissae[55] = 0.9733268277899110;
            _abscissae[56] = -0.9833362538846260;
            _abscissae[57] = 0.9833362538846260;
            _abscissae[58] = -0.9910133714767443;
            _abscissae[59] = 0.9910133714767443;

            _abscissae[60] = -0.9963401167719553;
            _abscissae[61] = 0.9963401167719553;
            _abscissae[62] = -0.9993050417357722;
            _abscissae[63] = -0.9993050417357722;
        }
        //the values in the array are theta,integralvalue
        public double[,] run(BRgait gait, function f)
        {
            _THETA = new double[64,2];
            double a = gait.gaitParam.theta0;
            double b = gait.gaitParam.thetaT;
            double dtheta0 = gait.gaitParam.dtheta0;
            double val = 0;
            for (int i = 0; i < 64; i += 2)
            {
                val += ((a - b) / 2) * _weigths[i] * f(((a - b) / 2) * _abscissae[i] + ((a + b) / 2));
                _THETA[64 / 2 - 1 - i / 2, 0] = ((a - b) / 2) * _abscissae[i] + ((a + b) / 2);
                _THETA[64 / 2 - 1 - i / 2, 1] = val;

                val += ((a - b) / 2) * _weigths[i+1] * f(((a - b) / 2) * _abscissae[i+1] + ((a + b) / 2));
                _THETA[64 / 2 + i / 2, 0] = ((a - b) / 2) * _abscissae[i] + ((a + b) / 2);
                _THETA[64 / 2 + i / 2, 1] = val;
            }
            return _THETA;
        }

    }
    public static class RiemannSum
    {
        public delegate double integralFunc(double theta);

        public static double[,] calculateFirstIntegral(integralFunc f)
        {
            double dx = 0.002;
            double b = 1;
            double val = 0;
            double theta = 0;
            double[,] RES = new double[(int)(b / dx), 2];
            for(int i = 0;i < b/ dx; i++)
            {
                val += f(theta) * dx;
                RES[i, 0] = theta;
                RES[i, 1] = val;

                theta += dx;
            }
            return RES;
        }
        public static double[,] calculateSecondIntegral(integralFunc f, double[,] firstIntegralVal)
        {
            double dx = 0.002;
            double b = 1;
            double val = 0;
            double theta = 0;
            double[,] RES = new double[(int)(b / dx), 2];
            for (int i = 0; i < b / dx; i++)
            {

                val += (Math.Exp(firstIntegralVal[i, 1]) * f(theta))*dx;

                RES[i, 0] = theta;
                RES[i, 1] = val;

                theta += dx;
            }
            return RES;
        }
    }
}
