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
    public class Biped
    {
        private BRParameters _param;
        private BRData _data;
        private BRgait[] _gaits;

        public Biped(string fileName)
        {
            _param = XMLBRParser.getXMLData(fileName);
            _data = new BRData();
        }

        public BRParameters param
        {
            get
            {
                return _param;
            }
            set
            {
                _param = value;
            }
        }

        public BRData data
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

        public BRgait[] gaits
        {
            get
            {
                return _gaits;
            }
            set
            {
                _gaits = value;
            }
        }

    }

    public class BRParameters
    {
        private string _description;
        private double _g;
        private double _m1;
        private double _L1;
        private double _l1;
        private double _J1;
        private double _m2;
        private double _L2;
        private double _l2;
        private double _J2;
        private double _m3;
        private double _L3;
        private double _l3;
        private double _J3;
        private double _mh;

        private double _JLM1;
        private double _JLM2;
        private double _JLM3;


        public BRParameters(string description, double g, double m1, double L1, double l1, double J1,
            double m2, double L2, double l2, double J2, double m3, double L3, double l3, double J3, double mh)
        {
            _description = description;
            _g = g;
            _mh = mh;
            _m1 = m1;
            _L1 = L1;
            _l1 = l1;
            _J1 = J1;

            _m2 = m2;
            _L2 = L2;
            _l2 = l2;
            _J2 = J2;

            _m3 = m3;
            _L3 = L3;
            _l3 = l3;
            _J3 = J3;

            _JLM1 = Math.Pow(L1, 2) * m2 + Math.Pow(L1, 2) * m3 + Math.Pow(l1, 2) * m1 + Math.Pow(l2, 2) * m2 + Math.Pow(l3, 2) * m3 + J1 + J2 + J3;
            _JLM2 = Math.Pow(l2, 2) * m2 + Math.Pow(l3, 2) * m3 + J2 + J3;
            _JLM3 = Math.Pow(l3, 2) * m3 + J3;


        }

        public string description
        {
            get
            {
                return _description;
            }
            set
            {
                _description = value;
            }
        }
        public double g
        {
            get
            {
                return _g;
            }
            set
            {
                _g = value;
            }
        }
        public double mh
        {
            get
            {
                return _mh;
            }
            set
            {
                _mh = value;
            }
        }
        public double m1
        {
            get
            {
                return _m1;
            }
            set
            {
                _m1 = value;
            }
        }
        public double L1
        {
            get
            {
                return _L1;
            }
            set
            {
                _L1 = value;
            }
        }
        public double l1
        {
            get
            {
                return _l1;
            }
            set
            {
                _l1 = value;
            }
        }
        public double J1
        {
            get
            {
                return _J1;
            }
            set
            {
                _J1 = value;
            }
        }

        public double m2
        {
            get
            {
                return _m2;
            }
            set
            {
                _m2 = value;
            }
        }
        public double L2
        {
            get
            {
                return _L2;
            }
            set
            {
                _L2 = value;
            }
        }
        public double l2
        {
            get
            {
                return _l2;
            }
            set
            {
                _l2 = value;
            }
        }
        public double J2
        {
            get
            {
                return _J2;
            }
            set
            {
                _J2 = value;
            }
        }

        public double m3
        {
            get
            {
                return _m3;
            }
            set
            {
                _m3 = value;
            }
        }
        public double L3
        {
            get
            {
                return _L3;
            }
            set
            {
                _L3 = value;
            }
        }
        public double l3
        {
            get
            {
                return _l3;
            }
            set
            {
                _l3 = value;
            }
        }
        public double J3
        {
            get
            {
                return _J3;
            }
            set
            {
                _J3 = value;
            }
        }

        public double JLM1
        {
            get
            {
                return _JLM1;
            }
            set
            {
                _JLM1 = value;
            }
        }
        public double JLM2
        {
            get
            {
                return _JLM2;
            }
            set
            {
                _JLM2 = value;
            }
        }
        public double JLM3
        {
            get
            {
                return _JLM3;
            }
            set
            {
                _JLM3 = value;
            }
        }
    }

    public class BRData
    {

        private Vector<double> _currentQ;
        private Vector<double> _currentDQ;
        private Vector<double> _currentDDQ;

        private double[] _time;
        private double _timestep;
        private Tuple<Vector<double>, double>[] _RES;

        public Vector<double> currentQ
        {
            get
            {
                return _currentQ;
            }
            set
            {
                _currentQ = value;
            }
        }
        public Vector<double> currentDQ
        {
            get
            {
                return _currentDQ;
            }
            set
            {
                _currentDQ = value;
            }
        }
        public Vector<double> currentDDQ
        {
            get
            {
                return _currentDDQ;
            }
            set
            {
                _currentDDQ = value;
            }
        }

        public double[] time
        {
            get
            {
                return _time;
            }
            set
            {
                _time = value;
            }
        }

        public double timestep
        {
            get
            {
                return _timestep;
            }
            set
            {
                _timestep = value;
            }
        }

        public Tuple<Vector<double>, double>[] RES
        {
            get
            {
                return _RES;
            }
            set
            {
                _RES = value;
            }
        }



    }

    public static class BRDynamics
    {
        public static Vector<double> rhs3D(Biped biped, Vector<double> dx)
        {
            BRData data = biped.data;
            BRParameters param = biped.param;
            double q1 = data.currentQ[0] + dx[0];
            double q2 = data.currentQ[1] + dx[1];
            double q3 = data.currentQ[2] + dx[2];

            double dq1 = data.currentDQ[0] + dx[3];
            double dq2 = data.currentDQ[1] + dx[4];
            double dq3 = data.currentDQ[2] + dx[5];

            double m11 = Math.Pow(param.L1, 2) * param.m2 + Math.Pow(param.L1, 2) * param.m3 + param.m1 * Math.Pow(param.l1, 2) + param.J1;
            double m12 = param.L1 * Math.Cos(-q2 + q1) * param.l2 * param.m2;
            double m13 = -param.L1 * Math.Cos(-q3 + q1) * param.l3 * param.m3;

            double m21 = param.L1 * Math.Cos(-q2 + q1) * param.l2 * param.m2;
            double m22 = Math.Pow(param.l2, 2) * param.m2 + param.J2;
            double m23 = 0;

            double m31 = -param.L1 * Math.Cos(-q3 + q1) * param.l3 * param.m3;
            double m32 = 0;
            double m33 = Math.Pow(param.l3, 2) * param.m3 + param.J3;

            var M = Matrix<double>.Build.DenseOfArray(new double[,] {
                {m11,m12,m13},
                {m21,m22,m23},
                {m31,m32,m33}
            });

            double c11 = 0;
            double c12 = param.L1 * Math.Sin(-q2 + q1) * param.l2 * param.m2 * dq2;
            double c13 = -param.L1 * Math.Sin(-q3 + q1) * param.l3 * param.m3 * dq3;

            double c21 = -param.L1 * Math.Sin(-q2 + q1) * param.l2 * param.m2 * dq1;
            double c22 = 0;
            double c23 = 0;

            double c31 = param.L1 * Math.Sin(-q3 + q1) * param.l3 * param.m3 * dq1;
            double c32 = 0;
            double c33 = 0;

            var C = Matrix<double>.Build.DenseOfArray(new double[,] {
                {c11,c12,c13},
                {c21,c22,c23},
                {c31,c32,c33},
            });

            double g1 = -(param.L1 * param.g * param.m2 + param.L1 * param.g * param.m3 + param.g * param.l1 * param.m1) * Math.Sin(q1);
            double g2 = -Math.Sin(q2) * param.g * param.l2 * param.m2;
            double g3 = Math.Sin(q3) * param.g * param.l3 * param.m3;

            var G = Vector<double>.Build.Dense(new double[] { g1, g2, g3 });

            //Vector<double> DDQ = M.Solve(-C.Multiply(Vector<double>.Build.Dense(new double[] { dq1, dq2, dq3 })) - G);
            Vector<double> DDQ = M.Inverse() * (-C.Multiply(Vector<double>.Build.Dense(new double[] { dq1, dq2, dq3 })) - G);
            return DDQ;
        }
        public static Vector<double> rhs6D(Biped biped, Vector<double> dx)
        {
            BRData data = biped.data;
            BRParameters param = biped.param;
            Vector<double> rhs3d = rhs3D(biped, dx);
            Vector<double> rhs6d = Vector.Build.Dense(new double[] { data.currentDQ[0] + dx[3], data.currentDQ[1] + dx[4], data.currentDQ[2] + dx[5], rhs3d[0], rhs3d[1], rhs3d[2] });
            return rhs6d;
        }
        public static double potentialEnergyAbsoluteAngles(Biped biped)
        {
            BRData data = biped.data;
            BRParameters param = biped.param;
            double q1 = data.currentQ[0];
            double q2 = data.currentQ[1];
            double q3 = data.currentQ[2];

            double P1 = param.m1 * param.g * param.l1 * Math.Cos(q1);
            double P2 = param.m2 * param.g * (param.L1 * Math.Cos(q1) + param.l2 * Math.Cos(q2));
            double P3 = param.m3 * param.g * (param.L1 * Math.Cos(q1) - Math.Cos(q3) * param.l3);
            double P = P1 + P2 + P3;
            return P;
        }
        public static double kineticEnergyAbsoluteAngles(Biped biped)
        {
            BRData data = biped.data;
            BRParameters param = biped.param;
            double K1 = 0.5 * param.m1 * Math.Pow(param.l1, 2) * Math.Pow(data.currentDQ[0], 2) + 0.5 * param.J1 * Math.Pow(data.currentDQ[0], 2);

            double K2 = 0.5 * Math.Pow(param.L1, 2) * Math.Pow(data.currentDQ[0], 2) * param.m2 + param.L1 * Math.Cos(-data.currentQ[1] + data.currentQ[0]) *
                data.currentDQ[1] * data.currentDQ[0] * param.l2 * param.m2 + (0.5 * Math.Pow(param.l2, 2) * param.m2 + 0.5 * param.J2) * Math.Pow(data.currentDQ[1], 2);

            double K3 = 0.5 * Math.Pow(param.L1, 2) * Math.Pow(data.currentDQ[0], 2) * param.m3 - param.L1 * data.currentDQ[0] * Math.Cos(-data.currentQ[2] + data.currentQ[0]) *
                data.currentDQ[2] * param.l3 * param.m3 + (0.5 * Math.Pow(param.l3, 2) * param.m3 + 0.5 * param.J3) * Math.Pow(data.currentDQ[2], 2);

            double K = K1 + K2 + K3;
            return K;
        }

        public static Vector<double> impactMapAngularMomentum(Biped biped)
        {
            BRData data = biped.data;
            BRParameters param = biped.param;
            double q1 = data.currentQ[0];
            double q2 = data.currentQ[1];
            double q3 = data.currentQ[2];

            double dq1 = data.currentDQ[0];
            double dq2 = data.currentDQ[1];
            double dq3 = data.currentDQ[2];

            double Qpos11 = -param.L1 * param.L3 * param.m1 * Math.Cos(q1 - q3) + param.L3 * param.l1 * param.m1 * Math.Cos(q1 - q3);
            double Qpos12 = 0;
            double Qpos13 = Math.Pow(param.L1, 2) * param.m1 - 2 * param.L1 * param.l1 * param.m1 + Math.Pow(param.l1, 2) * param.m1 + param.J1;

            double Qpos21 = param.L3 * param.l2 * param.m1 * Math.Cos(q2 - q3);
            double Qpos22 = Math.Pow(param.l2, 2) * param.m1 + param.J2;
            double Qpos23 = 0;

            double Qpos31 = Math.Pow(param.l3, 2) * param.m3 + Math.Pow(param.L3, 2) * param.m2 + Math.Pow(param.L3, 2) * param.m3 +
                Math.Pow(param.L3, 2) * param.m1 - 2 * param.L3 * param.l3 * param.m3 + param.L3 * param.l2 * param.m2 * Math.Cos(q2 - q3) +
                param.L3 * param.l1 * param.m1 * Math.Cos(q1 - q3) - param.L1 * param.L3 * param.m1 * Math.Cos(q1 - q3) + param.J3;
            double Qpos32 = param.L3 * param.l2 * param.m2 * Math.Cos(q2 - q3) + Math.Pow(param.l2, 2) * param.m2 + param.J2;
            double Qpos33 = Math.Pow(param.L1, 2) * param.m1 - param.L1 * param.L3 * param.m1 * Math.Cos(q1 - q3) -
                2 * param.L1 * param.l1 * param.m1 + param.L3 * param.l1 * param.m1 * Math.Cos(q1 - q3) + Math.Pow(param.l1, 2) * param.m1 + param.J1;


            double Qneg11 = -param.L1 * param.l1 * param.m1 + param.m1 * param.l1 * param.l1 + param.J1;
            double Qneg12 = 0;
            double Qneg13 = 0;
            double Qneg21 = param.L1 * param.l2 * param.m1 * Math.Cos(q1 - q2);
            double Qneg22 = param.l2 * param.l2 * param.m1 + param.J2;
            double Qneg23 = 0;
            double Qneg31 = param.J1 - param.L1 * param.l1 * param.m1 + param.L1 * param.L3 * param.m2 * Math.Cos(q1 - q3) + param.L1 * param.L3 * param.m3 * Math.Cos(q1 - q3) -
                param.L1 * Math.Cos(q1 - q3) * param.l3 * param.m3 + param.L1 * Math.Cos(q1 - q2) * param.l2 * param.m2 + param.L3 * param.l1 * param.m1 * Math.Cos(q1 - q3) + 
                param.m1 * param.l1 * param.l1;
            double Qneg32 = param.L3 * param.l2 * param.m2 * Math.Cos(-q3 + q2) + param.l2 * param.l2 * param.m2 + param.J2;
            double Qneg33 = -param.L3 * param.l3 * param.m3 + param.l3 * param.l3 * param.m3 + param.J3;

            Matrix <double> Qpos = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {Qpos11,Qpos12,Qpos13 },
                {Qpos21,Qpos22,Qpos23 },
                {Qpos31,Qpos32,Qpos33 }
            });

            Matrix<double> Qneg = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {Qneg11,Qneg12,Qneg13 },
                {Qneg21,Qneg22,Qneg23 },
                {Qneg31,Qneg32,Qneg33 }
            });

            Vector<double> qpos = Qpos.Solve(Qneg*data.currentDQ);

            return qpos;
        }

        public static Vector<double> stanceSwitch(Biped biped)
        {
            Matrix<double> ss = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {0,0,1 },
                {0,1,0 },
                {1,0,0 }
            });

            return ss.Multiply(biped.data.currentQ);
        }
    }
}
