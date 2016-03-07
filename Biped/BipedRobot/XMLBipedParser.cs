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

namespace BipedRobot{

    public static class XMLBRParser
    {

        public static BRParameters getXMLData(string XMLFileName)
        {
            FileStream fs = null;
            try
            {
                fs = new FileStream(XMLFileName, FileMode.Open, FileAccess.Read);
                XElement model = XElement.Load(fs);
                XElement physicalParameters = model.Element("physical_parameters");

                double g = (double)physicalParameters.Element("gravitational_constant");
                string description = (string)physicalParameters.Element("description");

                XElement stanceLeg = physicalParameters.Element("stance_leg_parameters");
                double m1 = (double)stanceLeg.Element("mass");
                double L1 = (double)stanceLeg.Element("length");
                double l1 = (double)stanceLeg.Element("CoM_position");
                double J1 = (double)stanceLeg.Element("inertia");

                XElement body = physicalParameters.Element("body_parameters");
                double m2 = (double)body.Element("mass");
				double L2 = (double)body.Element("length");
				double l2 = (double)body.Element("CoM_position");
				double J2 = (double)body.Element("inertia");

                XElement swingLeg = physicalParameters.Element("swing_leg_parameters");
                double m3 = (double)swingLeg.Element("mass");
				double L3 = (double)swingLeg.Element("length");
				double l3 = (double)swingLeg.Element("CoM_position");
				double J3 = (double)swingLeg.Element("inertia");

                XElement hip = physicalParameters.Element("hips_mass_parameter");
                double mh = (double)hip.Element("mass");

                //even if mh == 0.0
                double lH2 = l2 * m2 / (m2 + mh);
                m2 = m2 + mh;
                l2 = lH2;
                J2 = J2 + mh * Math.Pow(lH2,2) + m2 * Math.Pow((l2 - lH2),2);

                BRParameters param = new BRParameters(description, g, m1, L1, l1, J1, m2, L2, l2, J2, m3, L3, l3, J3, mh);
                return param;
            }
            catch (SystemException e)
            {
                Console.WriteLine(e.Message);
                Console.WriteLine(e.GetType());
                return null;
            }
            finally
            {
                if (fs != null)
                {
                    fs.Close();
                }
            }

        }


    }

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

        public BRgait gaits { get; set; }

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
            _JLM2 = Math.Pow(l2,2) * m2 + Math.Pow(l3,2) * m3 + J2 + J3;
            _JLM3 = Math.Pow(l3, 2) * m3 + J3;


        }

        public string description
        {
            get{
                return _description;
            }
            set{
                _description = value;
            }
        }
        public double g{
            get{
                return _g;
            }
            set{
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
            get{
                return _m1;
            }
            set{
                _m1 = value;
            }
        }
        public double L1
        {
            get{
                return _L1;
            }
            set{
                _L1 = value;
            }
        }
        public double l1
        {
            get{
                return _l1;
            }
            set{
                _l1 = value;
            }
        }
        public double J1
        {
            get{
                return _J1;
            }
            set{
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
		private Tuple<Vector<double>,double>[] _RES;

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

		public double[] time{
			get{
				return _time;
			}
			set{
				_time = value;
			}
		}

		public double timestep{
			get{
				return _timestep;
			}
			set{
				_timestep = value;
			}
		}

		public Tuple<Vector<double>,double>[] RES {
			get {
				return _RES;
			}
			set {
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
            Vector<double> DDQ = M.Inverse()*(-C.Multiply(Vector<double>.Build.Dense(new double[] { dq1, dq2, dq3 })) - G);
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


            double Qneg1 = -dq1 * param.L1 * param.l1 * param.m1 + dq1 * Math.Pow(param.l1, 2) * param.m1 + param.J1 * dq1;
            double Qneg2 = dq1 * param.L1 * param.l2 * param.m1 * Math.Cos(-q2 + q1) + param.l2 * param.l2 * param.m1 * dq2 + param.J2 * dq2;
            double Qneg3 = dq3 * param.l3 * param.l3 * param.m3 + dq1 * param.l1 * param.l1 * param.m1 + dq2 * param.l2 * param.l2 * param.m2 + 
                dq1 * param.L3 * param.l1 * param.m1 * Math.Cos(q1 - q3)+ dq2 * param.L3 * param.l2 * param.m2 * Math.Cos(q2 - q3) + 
                dq1 * param.L1 * param.L3 * param.m3 * Math.Cos(q1 - q3) + dq1 * param.L1 * param.l2 * param.m2 * Math.Cos(-q2 + q1) - 
                dq1 * param.L1 * param.l3 * param.m3 * Math.Cos(q1 - q3) - dq3 * param.L3 * param.l3 * param.m3 + 
                dq1 * param.L1 * param.L3 * param.m2 * Math.Cos(q1 - q3) - dq1 * param.L1 * param.l1 * param.m1 +
                param.J3 * dq3 + param.J2 * dq2 + param.J1 * dq1;

            Matrix<double> Qpos = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {Qpos11,Qpos12,Qpos13 },
                {Qpos21,Qpos22,Qpos23 },
                {Qpos31,Qpos32,Qpos33 }
            });

            Vector<double> Qneg = Vector<double>.Build.Dense(new double[] { Qneg1, Qneg2, Qneg3 });

            Vector<double> qpos = Qpos.Solve(Qneg);

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



    public static class integration
    {

		public static void setInitialConditions(ref Biped biped)
		{
            double q1 = 0;
            double q2 = 0;
            double q3 = 0.4;
            double dq1 = -0.2;
            double dq2 = 0;
            double dq3 = 0.1;
            Vector<double> initialConditions = Vector<double>.Build.Dense(new double[] {
				q1,
				q2,
				q3,
				dq1,
				dq2,
				dq3,
			});
			double[] time = Generate.LinearSpaced (1001, 0.0, 2.0);
			Tuple<Vector<double>,double>[] RES = new Tuple<Vector<double>, double>[time.Length];
			RES [0] = new Tuple<Vector<double>, double>(initialConditions, time [0]);

			biped.data.time = time;
			biped.data.timestep = (time.Last()-time.First()) / time.Length;
			biped.data.RES = RES;
			biped.data.currentQ = Vector<double>.Build.Dense (new double[] {initialConditions[0],initialConditions[1],initialConditions[2] });
			biped.data.currentDQ = Vector<double>.Build.Dense (new double[] {initialConditions[3],initialConditions[4],initialConditions[5] });
		}

		public static double stoppingConditionREL(Vector<double> currentState)
		{
			double q1 = currentState [0];
			double q2 = currentState [1];
			double q3 = currentState [2];
			double val = (Math.Sin (q1) + Math.Sin (q1 + q2 + q3));
			return val;
		}

        public static double stoppingConditionABS(Vector<double> currentState)
        {
            double q1 = currentState[0];
            double q2 = currentState[1];
            double q3 = currentState[2];
            double val = Math.Cos(q1) - Math.Cos(q3);
            Console.WriteLine(val);
            return val;
        }


        public static void integrateUntilConditionGEO(ref Biped biped)
		{
			if (stoppingConditionABS (biped.data.currentQ) < 0) {
				Console.WriteLine ("ugyldig start");
			} else {
				integrateWhileConditionGEO (ref biped);
				}
			}

		public static void integrateWhileConditionGEO(ref Biped biped)
		{
			Vector<double> oneStep;
			for (int i = 0; i<biped.data.time.Length-1; i++) {
				//oneStep = rk4 (biped);
                oneStep = constantStep(biped);
				
				biped.data.RES [i+1] = new Tuple<Vector<double>,double> (oneStep, biped.data.time [i+1]);
				biped.data.currentQ = Vector<double>.Build.Dense (new double[] {oneStep[0], oneStep[1], oneStep[2] });
				biped.data.currentDQ = Vector<double>.Build.Dense (new double[] {oneStep[3], oneStep[4], oneStep[5] });
				if (stoppingConditionABS (biped.data.currentQ) < 0.0) {
					break;
				}
			}


		}

		public static void run(ref Biped biped){
			setInitialConditions (ref biped);

			integrateUntilConditionGEO(ref biped);

            updateAfterImpact(ref biped);
		}





		public static Vector<double> rk4(Biped biped)
		{
			double dx = biped.data.timestep;
			double halfdx = 0.5 * dx;
			double sixth = 1.0 / 6.0;

			Vector<double> k0 = Vector<double>.Build.Dense (new double[] { 0, 0, 0, 0, 0, 0 });
			Vector<double> k1 = dx * BRDynamics.rhs6D(biped, k0);
			Vector<double> k2 = dx * BRDynamics.rhs6D(biped, k1*halfdx);
			Vector<double> k3 = dx * BRDynamics.rhs6D(biped, k2*halfdx);
			Vector<double> k4 = dx * BRDynamics.rhs6D(biped, k3*dx);

			return (Vector<double>.Build.Dense (new double[] {biped.data.currentQ[0], biped.data.currentQ[1], biped.data.currentQ[2],
				biped.data.currentDQ[0], biped.data.currentDQ[1], biped.data.currentDQ[2]
			}) +
				sixth * (k1 + 2 * k2 + 2 * k3 + k4));
		}




        public static Vector<double> constantStep(Biped biped)
        {
            double dx = biped.data.timestep;
            Vector<double> k0 = Vector<double>.Build.Dense(new double[] { 0, 0, 0, 0, 0, 0 });
            Vector<double> dq =biped.data.currentDQ + dx * BRDynamics.rhs3D(biped, k0);
            Vector<double> q = biped.data.currentQ + dx * 0.5 * (dq + biped.data.currentDQ);

            return Vector<double>.Build.Dense(new double[] { q[0], q[1], q[2], dq[0], dq[1], dq[2] });
        }


        public static void updateAfterImpact(ref Biped biped)
        {
            Vector<double> DQ = BRDynamics.impactMapAngularMomentum(biped);
            Vector<double> Q = BRDynamics.stanceSwitch(biped);

            biped.data.currentQ = Q;
            biped.data.currentDQ = DQ;
        }


    }


    public static class plotting
    {
        public static void plotStates(Biped biped)
        {
            graph stateGraph = new graph(biped);
        }
    }


    public class BRGaitParameters
    {
        private Vector<double> _gaitParameters;
        private double _objFunVal;
        private double _intervalStart;
        private double _intervalEnd;

        public Vector<double> gaitparameters { get; set; }
        public double objFunVal { get; set; }
        public double intervalStart { get; set; }
        public double intervalEnd { get; set; }
    }

    public static class gaitSearch
    {
        public static void run(ref Biped biped)
        {
            //first time running use rand to find a valid gait
            BRgait gait = new BRgait();
        }

        public static bool verifyParameters(BRgait gait)
        {
            double thetaDotAtTSquare = (-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.secondIntegral, gait.param.intervalStart, gait.param.intervalEnd, 2)) /
                (1 - Math.Exp(-MathNet.Numerics.Integration.GaussLegendreRule.Integrate(gait.firstIntegral, gait.param.intervalStart, gait.param.intervalEnd, 2)) * Math.Pow(gait.impact(gait.param.intervalEnd), 2));
            if (thetaDotAtTSquare < 0){
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
        private Expression _q1;
        private Expression _q3;
        private Expression _dq3;
        private Expression _ddq3;

        private Expression _alpha;
        private Expression _beta;
        private Expression _gamma;

        private Expression _impact;
        public delegate double function(double theta);

        public BRVHC()
        {
            StreamReader fs = null;
            fs = new StreamReader(@"../../../q3.txt");
            string temp = fs.ReadLine();
            temp =temp.Replace("tan", "Tan").Replace("cos","Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _q3 = new Expression(temp.Substring(temp.IndexOf('=')+1,temp.LastIndexOf(';')- temp.IndexOf('=')-1));
            Console.WriteLine(q3(0));
            fs.Close();

            fs = new StreamReader(@"../../../ddq3.txt");
            temp = fs.ReadLine();
            temp = temp.Replace("tan", "Tan").Replace("cos", "Cos").Replace("sin", "Sin").Replace("pow", "Pow");
            _ddq3 = new Expression(temp.Substring(temp.IndexOf('=') + 1, temp.LastIndexOf(';') - temp.IndexOf('=') - 1));
            Console.WriteLine(ddq3(0,0,1));
            fs.Close();
        }
        public double q1(double theta)
        {
            _q1.Parameters["theta"] = theta;
            return (double)_q1.Evaluate();
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
            _dq3.Parameters["theta"] = theta;
            return (double)_q3.Evaluate();
        }

        public double ddq1(double theta)
        {
            return Math.Cos(theta);
        }
        public double ddq2(double theta)
        {
            return Math.Cos(theta);
        }
        public double ddq3(double theta, double dtheta, double ddtheta)
        {
            _ddq3.Parameters["theta"] = theta;
            _ddq3.Parameters["dtheta"] = dtheta;
            _ddq3.Parameters["ddtheta"] = ddtheta;
            return (double)_ddq3.Evaluate();
        }

        public double alpha(double theta)
        {
            _alpha.Parameters["theta"] = theta;
            return (double)_alpha.Evaluate();
        }
        public double beta(double theta)
        {
            _beta.Parameters["theta"] = theta;
            return (double)_beta.Evaluate();
        }
        public double gamma(double theta)
        {
            _gamma.Parameters["theta"] = theta;
            return (double)_gamma.Evaluate();
        }

        public double twoTimesBetaDividedByAlpha(double theta)
        {
            _alpha.Parameters["theta"] = theta;
            _beta.Parameters["theta"] = theta;
            return (2*(double)_beta.Evaluate()/(double)_alpha.Evaluate());
        }
        public double twoTimeGammaDividedByAlpha(double theta)
        {
            _gamma.Parameters["theta"] = theta;
            _alpha.Parameters["theta"] = theta;
            return (2 * (double)_gamma.Evaluate() / (double)_alpha.Evaluate());
        }

        public double impact(double theta)
        {
            _impact.Parameters["theta"] = theta;
            return (double)_impact.Evaluate();
        }
    }

    public class BRgait
    {
        private BRVHC _vhc;
        private BRGaitParameters _param;

        public BRgait()
        {
            _vhc = new BRVHC();
            _param = new BRGaitParameters();
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
        public BRGaitParameters param { get; set; }

        public double firstIntegral(double theta)
        {

            return _vhc.twoTimeGammaDividedByAlpha(theta);
        }

        public double secondIntegral(double theta)
        {
            double a = MathNet.Numerics.Integration.GaussLegendreRule.Integrate(_vhc.twoTimesBetaDividedByAlpha, _param.intervalStart, theta, 2);
            return Math.Exp(a) * _vhc.twoTimeGammaDividedByAlpha(theta);
        }

        public double impact(double theta)
        {
            return _vhc.impact(theta);
        }

    }
}









