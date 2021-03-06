﻿using System;
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
using System.Globalization;

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
                J2 = J2 + mh * Math.Pow(lH2, 2) + (m2 - mh) * Math.Pow((lH2 - l2), 2);
                l2 = lH2;

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
        public static void writeXMLDataToFile(string txtFileName, BRParameters parameters)
        {
            string m = String.Format(CultureInfo.InvariantCulture, "{0} {1} {2}", parameters.m1, parameters.m2, parameters.m3);
            string l = String.Format(CultureInfo.InvariantCulture, "{0} {1} {2}", parameters.l1, parameters.l2, parameters.l3);
            string L = String.Format(CultureInfo.InvariantCulture, "{0} {1} {2}", parameters.L1, parameters.L2, parameters.L3);
            string J = String.Format(CultureInfo.InvariantCulture, "{0} {1} {2}", parameters.J1, parameters.J2, parameters.J3);
            string etc = String.Format(CultureInfo.InvariantCulture, "{0} {1} {2}", parameters.g, parameters.mh, 0);

            string[] data = { m, l, L, J, etc };
            File.WriteAllLines(@"../../../../"+txtFileName+".txt", data);

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

			biped.simulationData.time = time;
            biped.simulationData.timestep = (time.Last() - time.First()) / time.Length;
            biped.simulationData.RES = RES;
            biped.simulationData.currentQ = Vector<double>.Build.Dense(new double[] { initialConditions[0], initialConditions[1], initialConditions[2] });
            biped.simulationData.currentDQ = Vector<double>.Build.Dense(new double[] { initialConditions[3], initialConditions[4], initialConditions[5] });
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
            if (stoppingConditionABS(biped.simulationData.currentQ) < 0)
            {
				Console.WriteLine ("ugyldig start");
			} else {
				integrateWhileConditionGEO (ref biped);
				}
			}

		public static void integrateWhileConditionGEO(ref Biped biped)
		{
			Vector<double> oneStep;
            for (int i = 0; i < biped.simulationData.time.Length - 1; i++)
            {
				oneStep = rk4 (biped);
                //oneStep = constantStep(biped);

                biped.simulationData.RES[i + 1] = new Tuple<Vector<double>, double>(oneStep, biped.simulationData.time[i + 1]);
                biped.simulationData.currentQ = Vector<double>.Build.Dense(new double[] { oneStep[0], oneStep[1], oneStep[2] });
                biped.simulationData.currentDQ = Vector<double>.Build.Dense(new double[] { oneStep[3], oneStep[4], oneStep[5] });
                //check if both legs hits the ground and a heel strike is occuring (dq3 < 0)
                if ((stoppingConditionABS(biped.simulationData.currentQ) < 0.0) && oneStep[5] < 0.0)
                {
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
            double dx = biped.simulationData.timestep;
			double halfdx = 0.5 * dx;
			double sixth = 1.0 / 6.0;

			Vector<double> k0 = Vector<double>.Build.Dense (new double[] { 0, 0, 0, 0, 0, 0 });
			Vector<double> k1 = dx * BRDynamics.rhs6D(biped, k0);
			Vector<double> k2 = dx * BRDynamics.rhs6D(biped, k1*halfdx);
			Vector<double> k3 = dx * BRDynamics.rhs6D(biped, k2*halfdx);
			Vector<double> k4 = dx * BRDynamics.rhs6D(biped, k3*dx);

            return (Vector<double>.Build.Dense(new double[] {biped.simulationData.currentQ[0], biped.simulationData.currentQ[1], biped.simulationData.currentQ[2],
				biped.simulationData.currentDQ[0], biped.simulationData.currentDQ[1], biped.simulationData.currentDQ[2]
			}) +
				sixth * (k1 + 2 * k2 + 2 * k3 + k4));
		}

        public static Vector<double> constantStep(Biped biped)
        {
            double dx = biped.simulationData.timestep;
            Vector<double> k0 = Vector<double>.Build.Dense(new double[] { 0, 0, 0, 0, 0, 0 });
            Vector<double> dq = biped.simulationData.currentDQ + dx * BRDynamics.rhs3D(biped, k0);
            Vector<double> q = biped.simulationData.currentQ + dx * 0.5 * (dq + biped.simulationData.currentDQ);

            return Vector<double>.Build.Dense(new double[] { q[0], q[1], q[2], dq[0], dq[1], dq[2] });
        }

        public static void updateAfterImpact(ref Biped biped)
        {
            Vector<double> DQ = BRDynamics.impactMapAngularMomentum(biped);
            Vector<double> Q = BRDynamics.stanceSwitch(biped);

            biped.simulationData.currentQ = Q;
            biped.simulationData.currentDQ = DQ;
        }
    }

    public static class integrationReducedDynamics
    {
        public static BRReducedSimulationData run(BRVHC vhc, Vector<double> integrationStart, Vector<double> integrationEnd)
        {
            BRReducedSimulationData data = new BRReducedSimulationData(0.001, Vector<double>.Build.Dense(new double[] { integrationStart[0], integrationStart[1],
                BRReducedDynamics.rhs1D(integrationStart, vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma)}));
            int i = 1;
            Vector<double> oneStep = integrationStart;
            while (true){
				oneStep = rk4(oneStep, data.timestep, vhc);
                double ddtheta = BRReducedDynamics.rhs1D(oneStep, vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma);
				data.RES.Add(new Tuple<Vector<double>,double> (Vector<double>.Build.Dense(new double[]{oneStep[0], oneStep[1], ddtheta}), i*data.timestep));
				if (oneStep[0] > integrationEnd[0]) {
					break;
				}
                i++;
                if(i > 20000)
                {
                    break;
                }
            }
            return data;
        }

        public static Vector<double> rk4(Vector<double> THETA, double dx, BRVHC vhc)
        {
            double halfdx = 0.5 * dx;
            double sixthdx = 1.0 / 6.0 * dx;

            Vector<double> k0 = Vector<double>.Build.Dense(new double[] { 0, 0});
            Vector<double> k1 = BRReducedDynamics.rhs2D(THETA + k0, vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma);
            Vector<double> k2 = BRReducedDynamics.rhs2D(THETA + k1 * halfdx, vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma);
            Vector<double> k3 = BRReducedDynamics.rhs2D(THETA + k2 * halfdx, vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma);
            Vector<double> k4 = BRReducedDynamics.rhs2D(THETA + k3 * dx, vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma);

            return THETA + sixthdx * (k1 + 2 * k2 + 2 * k3 + k4);
        }
    }
    

    public static class calculateReducedDynamics
    {
        public static double[,] run(BRgait gait)
        {
            GaussianQuadrature quad = new GaussianQuadrature();
            double[,] firstIntegral = quad.run(gait, gait.firstIntegral);
            double[,] secondIntegral = quad.run(gait, gait.secondIntegral);
            double[,] THETA = new double[firstIntegral.Length, 3];
            THETA[0, 0] = 0;
            THETA[0, 1] = 1;
            double sumFirstIntegral = 0;
            double sumSecondIntegral = 0;
            double dtheta0 = gait.gaitParam.dtheta0;
            double dtheta;
            for (int i = 1; i < firstIntegral.Length; i++)
            {
                sumFirstIntegral += firstIntegral[i,1];
                sumSecondIntegral += secondIntegral[i,1];
                dtheta = Math.Sqrt(Math.Exp(sumFirstIntegral)) * Math.Pow(dtheta0, 2) + sumSecondIntegral;
                THETA[i, 0] = firstIntegral[i, 0];
                THETA[i, 1] = dtheta;
                THETA[i, 2] = BRReducedDynamics.rhs1D(Vector<double>.Build.Dense(new[] { THETA[i, 0], dtheta }), gait.vhc.evalAlpha, gait.vhc.evalBeta, gait.vhc.evalGamma);
            }
            return THETA;
        }
    }

    public static class calculateTorques
    {
        public static BRTorques run(BRVHC vhc, BRReducedSimulationData data)
        {
            BRTorques torques = new BRTorques(data.RES.Count);

            for (int i = 0; i < data.RES.Count; i++)
            {
                torques.torque1[i] = vhc.evalAlpha1(data.RES[i].Item1[0]) * data.RES[i].Item1[2] + vhc.evalBeta1(data.RES[i].Item1[0]) * Math.Pow(data.RES[i].Item1[1],2) + vhc.evalGamma1(data.RES[i].Item1[0]);
                torques.torque2[i] = vhc.evalAlpha3(data.RES[i].Item1[0]) * data.RES[i].Item1[2] + vhc.evalBeta3(data.RES[i].Item1[0]) * Math.Pow(data.RES[i].Item1[1], 2) + vhc.evalGamma3(data.RES[i].Item1[0]);
                torques.theta[i] = data.RES[i].Item1[0];
                torques.dtheta[i] = data.RES[i].Item1[1];
                torques.ddtheta[i] = data.RES[i].Item1[2];

            }
            return torques;
        }

    }

    public static class calculateTorques2
    {
        public static BRTorques run(BRVHC vhc, double[,] data)
        {
            BRTorques torques = new BRTorques(data.Length / data.Rank);
            double ddtheta;
            for (int i = 0; i < data.Length / data.Rank; i++)
            {
                ddtheta = (-vhc.evalBeta(data[0, i]) * Math.Pow(data[1, i], 2) - vhc.evalGamma(data[0, i])) / vhc.evalAlpha(data[0, i]);
                torques.torque1[i] = vhc.evalAlpha1(data[0, i]) * ddtheta + vhc.evalBeta1(data[0, i]) * Math.Pow(data[1, i], 2) + vhc.evalGamma1(data[0, i]);
                torques.torque2[i] = vhc.evalAlpha3(data[0, i]) * ddtheta + vhc.evalBeta3(data[0, i]) * Math.Pow(data[1, i], 2) + vhc.evalGamma3(data[0, i]);
                torques.theta[i] = data[0, i];
                torques.dtheta[i] = data[1, i];
                torques.ddtheta[i] = ddtheta;

            }
            return torques;
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








