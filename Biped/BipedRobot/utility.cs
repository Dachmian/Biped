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


    public static class IntegrationFullDynamics
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
				oneStep = rk4 (biped);
                //oneStep = constantStep(biped);
				
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

    public static class integrationReducedDynamics
    {
        

        public static void run(BRVHC vhc)
        {
            
        }

        public static Vector<double> rk4(Vector<double> THETA, BRVHC vhc)
        {
            double dx = 0.01;
            double halfdx = 0.5 * dx;
            double sixth = 1.0 / 6.0;

            Vector<double> k0 = Vector<double>.Build.Dense(new double[] { 0, 0, 0, 0, 0, 0 });
            Vector<double> k1 = dx * BRReducedDynamics.rhs2D(THETA + k0, vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma);
            Vector<double> k2 = dx * BRReducedDynamics.rhs2D(THETA + k1 * halfdx, vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma);
            Vector<double> k3 = dx * BRReducedDynamics.rhs2D(THETA + k2 * halfdx, vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma);
            Vector<double> k4 = dx * BRReducedDynamics.rhs2D(THETA + k3 * dx, vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma);

            return + sixth * (k1 + 2 * k2 + 2 * k3 + k4);
        }
    }

    public static class plotting
    {
        public static void plotStates(Biped biped)
        {
            graph stateGraph = new graph(biped);
        }
    }



}








