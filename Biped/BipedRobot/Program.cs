using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Symbolics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Integration;
using System.IO;
using System.Globalization;

namespace BipedRobot{
    static class Program{
        [STAThread]
        static void Main(string[] args)
        {
            //MainForm mainForm = new MainForm();
            //basic_randomized_set
            //CAD_params_ILeg
            Biped biped = new Biped(@"../../../CAD_params_ILeg.xml");
            //IntegrationFullDynamics.run(ref biped);
            //plotting.plotStates(biped);
            //findGaitsManually(biped);
            test(biped);
            //XMLBRParser.writeXMLDataToFile("physicalParameters", biped.param);
            //findCoMGait(biped);
            //findq1Gait(biped);
            //testq1(biped);
            
        }

        static void test(Biped biped)
        {
            int numberOfPoints = 6;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);


            string[] parameters = File.ReadAllLines(@"../../../../parameters.txt");
            //double dtheta0 = Convert.ToDouble(File.ReadAllLines(@"../../../../dtheta0.txt")[0], CultureInfo.InvariantCulture);
            //double dthetaT = Convert.ToDouble(File.ReadAllLines(@"../../../../dthetaT.txt")[0], CultureInfo.InvariantCulture);
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);

            for (int i = 0; i < 6; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 6; i < 12; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi2Parameters = tuple.Item1;

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 12; i < 18; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi3Parameters = tuple.Item1;

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());


            double dtheta0Squared = 0;// Math.Pow(dtheta0, 2);
            double dthetaTSquared = 0;// Math.Pow(dthetaT, 2);

            double[,] THETA = evalDthetaConstraint(gait, ref dtheta0Squared, ref dthetaTSquared);
            //writeThetaToFile(THETA, @"../../../../../../../../MATLAB/theta.csv", gait.vhc);
            Phaseportrait plot = new Phaseportrait(THETA);

            BRReducedSimulationData data = integrationReducedDynamics.run(gait.vhc, Vector<double>.Build.Dense(new double[] { 0, Math.Sqrt(dtheta0Squared) }),
                Vector<double>.Build.Dense(new double[] { 1, Math.Sqrt(dthetaTSquared) }));
            biped.reducedSimulationData = data;
            //writeThetaEvolutionToFile(data, @"../../../../../../../../MATLAB/thetaTime.csv");
            graph graph = new graph(biped);

            ErrorGraph errorGraph = new ErrorGraph(gait, THETA, data);

            BRTorques torques = calculateTorques2.run(gait.vhc, THETA);
            TorquesGraph torquesGraph = new TorquesGraph(torques);


            TrajectoryGraph trajectory = new TrajectoryGraph(gait, THETA);

            TrajectoryPortrait trajectoryPortrait = new TrajectoryPortrait(gait, THETA);

            Console.WriteLine(gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactSecondLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactThirdLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));

            
            
        }

        static void findGaitsManually(Biped biped)
        {
            int numberOfPoints = 6;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);


            string[] parameters = File.ReadAllLines(@"../../../../parameters.txt");
            //double dtheta0 = Convert.ToDouble(File.ReadAllLines(@"../../../../dtheta0.txt")[0], CultureInfo.InvariantCulture);
            //double dthetaT = Convert.ToDouble(File.ReadAllLines(@"../../../../dthetaT.txt")[0], CultureInfo.InvariantCulture);
            for (int i = 0; i < parameters.Length; i++)
            {
                Console.WriteLine(parameters[i]);
            }
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);

            for (int i = 0; i < 6; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 6; i < 12; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi2Parameters = tuple.Item1;

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 12; i < 18; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi3Parameters = tuple.Item1;

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());




            Console.WriteLine(gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactSecondLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactThirdLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));


            double dtheta0Squared = 0;// Math.Pow(dtheta0, 2);
            double dthetaTSquared = 0;// Math.Pow(dthetaT, 2);


            //SLSQPFifthOrderTorqueEquality slSQPTE = new SLSQPFifthOrderTorqueEquality(gait, numberOfPoints);
            //slSQPTE.runNumericalSLSQP();

            //SLSQPFifthOrderTorque slSQPT = new SLSQPFifthOrderTorque(gait, numberOfPoints);
            //slSQPT.runNumericalSLSQP();

            //SLSQPFifthOrder slSQP = new SLSQPFifthOrder(gait, numberOfPoints);
            //slSQP.runNumericalSLSQP();

            //SLSQPFifthOrderTorqueDtheta slSQPTD = new SLSQPFifthOrderTorqueDtheta(gait, numberOfPoints);
            //slSQPTD.runNumericalSLSQP();

            //SLSQPFifthOrderDtheta slSQPTD = new SLSQPFifthOrderDtheta(gait, numberOfPoints);
            //slSQPTD.runNumericalSLSQP();

            SLSQPFifthOrderTorqueConstraint slSQPTD = new SLSQPFifthOrderTorqueConstraint(gait, numberOfPoints);
            slSQPTD.runNumericalSLSQP();
            double[,] THETA = evalDthetaConstraint(gait, ref dtheta0Squared, ref dthetaTSquared);

            Phaseportrait plot = new Phaseportrait(THETA);

            BRReducedSimulationData data = integrationReducedDynamics.run(gait.vhc, Vector<double>.Build.Dense(new double[] { 0, Math.Sqrt(dtheta0Squared) }),
                Vector<double>.Build.Dense(new double[] { 1, Math.Sqrt(dthetaTSquared) }));
            biped.reducedSimulationData = data;

            graph graph = new graph(biped);


            using (StreamWriter file =
                                        new System.IO.StreamWriter(@"../../../../foundgaits.txt", true))
            {
                file.WriteLine(gait.vhc.phi1Parameters["P0"].RealValue);
                file.WriteLine(gait.vhc.phi1Parameters["P1"].RealValue);
                file.WriteLine(gait.vhc.phi1Parameters["P2"].RealValue);
                file.WriteLine(gait.vhc.phi1Parameters["P3"].RealValue);
                file.WriteLine(gait.vhc.phi1Parameters["P4"].RealValue);
                file.WriteLine(gait.vhc.phi1Parameters["P5"].RealValue);

                file.WriteLine(gait.vhc.phi2Parameters["P6"].RealValue);
                file.WriteLine(gait.vhc.phi2Parameters["P7"].RealValue);
                file.WriteLine(gait.vhc.phi2Parameters["P8"].RealValue);
                file.WriteLine(gait.vhc.phi2Parameters["P9"].RealValue);
                file.WriteLine(gait.vhc.phi2Parameters["P10"].RealValue);
                file.WriteLine(gait.vhc.phi2Parameters["P11"].RealValue);

                file.WriteLine(gait.vhc.phi3Parameters["P12"].RealValue);
                file.WriteLine(gait.vhc.phi3Parameters["P13"].RealValue);
                file.WriteLine(gait.vhc.phi3Parameters["P14"].RealValue);
                file.WriteLine(gait.vhc.phi3Parameters["P15"].RealValue);
                file.WriteLine(gait.vhc.phi3Parameters["P16"].RealValue);
                file.WriteLine(gait.vhc.phi3Parameters["P17"].RealValue);
                file.WriteLine("NEW GAIT");
                file.WriteLine("");
            }

            BRTorques torques = calculateTorques2.run(gait.vhc, THETA);
            TorquesGraph torquesGraph = new TorquesGraph(torques);
        }

        public static double[,] evalDthetaConstraint(BRgait gait, ref double dtheta0Squared, ref double dthetaTSquared)
        {
            double[][] firstIntegral = TrapezoidalSum2.calculateFirstIntegral(gait.vhc.evalTwoTimesBetaDividedByAlpha, 0, 1, 200);
            double firstIntegralVALUE = GaussLegendreRule.Integrate(gait.firstIntegral, 0, 1, 100);
            double secondIntegralVALUE = GaussLegendreRule.Integrate(gait.secondIntegral, 0, 1, 100);
            double[][] secondIntegral = TrapezoidalSum2.calculateSecondIntegral(gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1], 0, 1, 200);            
            int len = firstIntegral[1].Length;
            dthetaTSquared = (-secondIntegral[1][len - 1] * Math.Exp(firstIntegral[1][len - 1])) / (1 - Math.Exp(firstIntegral[1][len - 1]) * Math.Pow(gait.impactSecondLine(0, 1), 2));
            dtheta0Squared = dthetaTSquared * Math.Pow(gait.impactSecondLine(0, 1), 2);
            double dthetaTSquaredTest = (-secondIntegralVALUE) / (1 - Math.Exp(-firstIntegralVALUE) * Math.Pow(gait.impactSecondLine(0, 1), 2));
            double dtheta0squaredTest = dthetaTSquaredTest * Math.Pow(gait.impactSecondLine(0, 1), 2);
            int numberOfCycles = 1;
            double[,] THETA = new double[2, numberOfCycles * (len)];
            double temp = 0;
            //write firstintegral and secondintegral values to matlab
            writeRhoToFile(firstIntegral, secondIntegral, @"../../../../../../../../MATLAB/rho.csv", gait.vhc);
            for (int j = 0; j < numberOfCycles; j++)
            {
                THETA[0, j*len] = 0;
                THETA[1, j*len] = Math.Sqrt(dtheta0Squared);
                for (int i = 1; i < len; i++)
                {
                    THETA[0, i + j*len] = secondIntegral[0][i];
                    THETA[1, i + j*len] = Math.Sqrt(-secondIntegral[1][i] * Math.Exp(firstIntegral[1][i]) + Math.Exp(firstIntegral[1][i]) * dtheta0Squared);
                    temp = THETA[1, i];
                }
                dtheta0Squared = Math.Pow(gait.impactSecondLine(0, 1), 2) * Math.Pow(temp, 2);
            }
            using (StreamWriter file =
                        new System.IO.StreamWriter(@"../../../../VALS.txt", true))
            {
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(THETA[0, i]);
                }
                file.WriteLine("");
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(THETA[1, i]);
                }
                file.WriteLine("");
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(Math.Exp(firstIntegral[1][i]));
                }
                file.WriteLine("");
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(secondIntegral[1][i]);
                }
                file.WriteLine("");
            }
            return THETA;
        }
        public static void writeThetaToFile(double[,] THETA, string fileName, BRVHC vhc)
        {
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName))
            {
                for (int i = 0; i < THETA.Rank; i++)
                {
                    for (int j = 0; j < THETA.Length / (THETA.Rank * 8); j += 2)
                    {
                        file.Write(THETA[i, j].ToString(CultureInfo.InvariantCulture) + ",");
                    }
                    file.Write("0"+Environment.NewLine);
                }
                for (int i = 0; i < THETA.Length / (THETA.Rank * 8); i += 2)
                {
                    file.Write(BRReducedDynamics.rhs1D(Vector<double>.Build.Dense(new double[] { THETA[0, i], THETA[1, i] }), vhc.evalAlpha, vhc.evalBeta, vhc.evalGamma).ToString(CultureInfo.InvariantCulture) + ",");
                }
                file.Write("0" + Environment.NewLine);
            }
        }
        public static void writeRhoToFile(double[][] rho1, double[][] rho2, string fileName, BRVHC vhc)
        {
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName))
            {
                for (int i = 0; i < rho1[1].Length; i += 2)
                {
                    file.Write(rho1[1][i].ToString(CultureInfo.InvariantCulture) + ",");
                }
                file.Write("0" + Environment.NewLine);
                for (int i = 0; i < rho2[1].Length; i += 2)
                {
                    file.Write(rho2[1][i].ToString(CultureInfo.InvariantCulture) + ",");
                }
                file.Write("0" + Environment.NewLine);
            }
        }
        public static void writeThetaEvolutionToFile(BRReducedSimulationData data, string fileName)
        {
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(fileName))
            {
                for (int j = 0; j < data.RES.Count; j ++)
                {
                    file.Write(data.RES[j].Item1[0].ToString(CultureInfo.InvariantCulture) + ",");
                }
                file.Write("0" + Environment.NewLine);
                for (int j = 0; j < data.RES.Count; j++)
                {
                    file.Write(data.RES[j].Item1[1].ToString(CultureInfo.InvariantCulture) + ",");
                }
                file.Write("0" + Environment.NewLine);
                for (int j = 0; j < data.RES.Count; j++)
                {
                    file.Write(data.RES[j].Item1[2].ToString(CultureInfo.InvariantCulture) + ",");
                }
                file.Write("0" + Environment.NewLine);
                for (int j = 0; j < data.RES.Count; j++)
                {
                    file.Write(data.RES[j].Item2.ToString(CultureInfo.InvariantCulture) + ",");
                }
                file.Write("0" + Environment.NewLine);
            }
        }

        static void findCoMGait(Biped biped)
        {
            int numberOfPoints = 6;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);


            string[] parameters = File.ReadAllLines(@"../../../../parameters.txt");
            //double dtheta0 = Convert.ToDouble(File.ReadAllLines(@"../../../../dtheta0.txt")[0], CultureInfo.InvariantCulture);
            //double dthetaT = Convert.ToDouble(File.ReadAllLines(@"../../../../dthetaT.txt")[0], CultureInfo.InvariantCulture);
            for (int i = 0; i < parameters.Length; i++)
            {
                Console.WriteLine(parameters[i]);
            }
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);

            for (int i = 0; i < 6; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 6; i < 12; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi2Parameters = tuple.Item1;

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 12; i < 18; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi3Parameters = tuple.Item1;

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());

            biped.gait = gait;
            double dtheta0Squared = 0;// Math.Pow(dtheta0, 2);
            double dthetaTSquared = 0;// Math.Pow(dthetaT, 2);
            while (true)
            {
                SLSQPCoM slSQP = new SLSQPCoM(biped, numberOfPoints);
                double[] theta = thetaRange(slSQP.q1, slSQP.q2, slSQP.q3, biped.param);
                double[,] THETA = evalDthetaCoM(gait, ref dtheta0Squared, ref dthetaTSquared, theta[0], theta[1]);
                Phaseportrait plot = new Phaseportrait(THETA);
                TrajectoryGraph trajectory = new TrajectoryGraph(gait, THETA);
                slSQP.runAlpha();
                theta = thetaRange(slSQP.q1, slSQP.q2, slSQP.q3, biped.param);
                double test2 = evaluateAlphaConstraint(theta[0], theta[1], gait);
                THETA = evalDthetaCoM(gait, ref dtheta0Squared, ref dthetaTSquared, theta[0], theta[1]);
                plot = new Phaseportrait(THETA);
                trajectory = new TrajectoryGraph(gait, THETA);
                if (evaluateAlphaConstraint(theta[0], theta[1], gait) > 0)
                {
                    continue;
                }
                slSQP.runDtheta();
                theta = thetaRange(slSQP.q1, slSQP.q2, slSQP.q3, biped.param);
                test2 = evaluateAlphaConstraint(theta[0], theta[1], gait);
                double testere = evaluateDthetaConstraintCoM(theta[0], theta[1], gait);
                THETA = evalDthetaCoM(gait, ref dtheta0Squared, ref dthetaTSquared, theta[0], theta[1]);
                plot = new Phaseportrait(THETA);
                trajectory = new TrajectoryGraph(gait, THETA);
                if ((evaluateAlphaConstraint(theta[0], theta[1], gait) > 0) || (evaluateDthetaConstraintCoM(theta[0], theta[1], gait) < 0))
                {
                    continue;
                }
                slSQP.runImpact();
                theta = thetaRange(slSQP.q1, slSQP.q2, slSQP.q3, biped.param);
                Console.WriteLine(gait.impactFirstLine(theta[0], theta[1]));
                Console.WriteLine(gait.impactSecondLine(theta[0], theta[1]));
                Console.WriteLine(gait.impactThirdLine(theta[0], theta[1]));

                THETA = evalDthetaCoM(gait, ref dtheta0Squared, ref dthetaTSquared, theta[0], theta[1]);
                plot = new Phaseportrait(THETA);
                trajectory = new TrajectoryGraph(gait, THETA);
                using (StreamWriter file =
                                        new System.IO.StreamWriter(@"../../../../ImpactSatisfied.txt", true))
                {
                    file.WriteLine(gait.vhc.phi1Parameters["P0"].RealValue);
                    file.WriteLine(gait.vhc.phi1Parameters["P1"].RealValue);
                    file.WriteLine(gait.vhc.phi1Parameters["P2"].RealValue);
                    file.WriteLine(gait.vhc.phi1Parameters["P3"].RealValue);
                    file.WriteLine(gait.vhc.phi1Parameters["P4"].RealValue);
                    file.WriteLine(gait.vhc.phi1Parameters["P5"].RealValue);

                    file.WriteLine(gait.vhc.phi2Parameters["P6"].RealValue);
                    file.WriteLine(gait.vhc.phi2Parameters["P7"].RealValue);
                    file.WriteLine(gait.vhc.phi2Parameters["P8"].RealValue);
                    file.WriteLine(gait.vhc.phi2Parameters["P9"].RealValue);
                    file.WriteLine(gait.vhc.phi2Parameters["P10"].RealValue);
                    file.WriteLine(gait.vhc.phi2Parameters["P11"].RealValue);

                    file.WriteLine(gait.vhc.phi3Parameters["P12"].RealValue);
                    file.WriteLine(gait.vhc.phi3Parameters["P13"].RealValue);
                    file.WriteLine(gait.vhc.phi3Parameters["P14"].RealValue);
                    file.WriteLine(gait.vhc.phi3Parameters["P15"].RealValue);
                    file.WriteLine(gait.vhc.phi3Parameters["P16"].RealValue);
                    file.WriteLine(gait.vhc.phi3Parameters["P17"].RealValue);
                    file.WriteLine("NEW GAIT");
                    file.WriteLine("");
                }
                slSQP.run();

                using (StreamWriter file =
                                        new System.IO.StreamWriter(@"../../../../foundgaits.txt", true))
                {
                    file.WriteLine(gait.vhc.phi1Parameters["P0"].RealValue);
                    file.WriteLine(gait.vhc.phi1Parameters["P1"].RealValue);
                    file.WriteLine(gait.vhc.phi1Parameters["P2"].RealValue);
                    file.WriteLine(gait.vhc.phi1Parameters["P3"].RealValue);
                    file.WriteLine(gait.vhc.phi1Parameters["P4"].RealValue);
                    file.WriteLine(gait.vhc.phi1Parameters["P5"].RealValue);

                    file.WriteLine(gait.vhc.phi2Parameters["P6"].RealValue);
                    file.WriteLine(gait.vhc.phi2Parameters["P7"].RealValue);
                    file.WriteLine(gait.vhc.phi2Parameters["P8"].RealValue);
                    file.WriteLine(gait.vhc.phi2Parameters["P9"].RealValue);
                    file.WriteLine(gait.vhc.phi2Parameters["P10"].RealValue);
                    file.WriteLine(gait.vhc.phi2Parameters["P11"].RealValue);

                    file.WriteLine(gait.vhc.phi3Parameters["P12"].RealValue);
                    file.WriteLine(gait.vhc.phi3Parameters["P13"].RealValue);
                    file.WriteLine(gait.vhc.phi3Parameters["P14"].RealValue);
                    file.WriteLine(gait.vhc.phi3Parameters["P15"].RealValue);
                    file.WriteLine(gait.vhc.phi3Parameters["P16"].RealValue);
                    file.WriteLine(gait.vhc.phi3Parameters["P17"].RealValue);
                    file.WriteLine("NEW GAIT");
                    file.WriteLine("");
                }
            }
            //BRReducedSimulationData data = integrationReducedDynamics.run(gait.vhc, Vector<double>.Build.Dense(new double[] { 0, Math.Sqrt(dtheta0Squared) }),
            //    Vector<double>.Build.Dense(new double[] { 1, Math.Sqrt(dthetaTSquared) }));
            //biped.reducedSimulationData = data;

            //graph graph = new graph(biped);

            //BRTorques torques = calculateTorques2.run(gait.vhc, THETA);
            //TorquesGraph torquesGraph = new TorquesGraph(torques);
        }

        public static double[,] evalDthetaCoM(BRgait gait, ref double dtheta0Squared, ref double dthetaTSquared, double thetaMin, double thetaMax)
        {
            double[][] firstIntegral = TrapezoidalSum2.calculateFirstIntegral(gait.vhc.evalTwoTimesBetaDividedByAlpha, thetaMin, thetaMax, 200);
            double firstIntegralVALUE = GaussLegendreRule.Integrate(gait.firstIntegral, thetaMin, thetaMax, 100);
            double secondIntegralVALUE = GaussLegendreRule.Integrate(gait.secondIntegral, thetaMin, thetaMax, 100);
            double[][] secondIntegral = TrapezoidalSum2.calculateSecondIntegral(gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1], thetaMin, thetaMax, 200);
            int len = firstIntegral[1].Length;
            dthetaTSquared = (-secondIntegral[1][len - 1] * Math.Exp(firstIntegral[1][len - 1])) / (1 - Math.Exp(firstIntegral[1][len - 1]) * Math.Pow(gait.impactSecondLine(thetaMin, thetaMax), 2));
            dtheta0Squared = dthetaTSquared * Math.Pow(gait.impactSecondLine(thetaMin, thetaMax), 2);
            double dthetaTSquaredTest = (-secondIntegralVALUE) / (1 - Math.Exp(-firstIntegralVALUE) * Math.Pow(gait.impactSecondLine(thetaMin, thetaMax), 2));
            double dtheta0squaredTest = dthetaTSquaredTest * Math.Pow(gait.impactSecondLine(thetaMin, thetaMax), 2);
            double[,] THETA = new double[2, 8 * (len)];
            double temp = 0;
            for (int j = 0; j < 8; j++)
            {
                THETA[0, j * len] = thetaMin;
                THETA[1, j * len] = Math.Sqrt(dtheta0Squared);
                for (int i = 1; i < len; i++)
                {
                    THETA[0, i + j * len] = secondIntegral[0][i];
                    THETA[1, i + j * len] = Math.Sqrt(-secondIntegral[1][i] * Math.Exp(firstIntegral[1][i]) + Math.Exp(firstIntegral[1][i]) * dtheta0Squared);
                    temp = THETA[1, i];
                }
                dtheta0Squared = Math.Pow(gait.impactSecondLine(0, 1), 2) * Math.Pow(temp, 2);
            }
            using (StreamWriter file =
                        new System.IO.StreamWriter(@"../../../../VALS.txt", true))
            {
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(THETA[0, i]);
                }
                file.WriteLine("");
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(THETA[1, i]);
                }
                file.WriteLine("");
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(Math.Exp(firstIntegral[1][i]));
                }
                file.WriteLine("");
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(secondIntegral[1][i]);
                }
                file.WriteLine("");
            }
            return THETA;
        }
        public static double[] thetaRange(double q1, double q2, double q3, BRParameters _param)
        {
            double a = -(_param.m1 * _param.l1 + _param.m2 * _param.L1 + _param.m3 * _param.L1) / (_param.m1 + _param.m2 + _param.m3);
            double b = -(_param.m2 * _param.l2) / (_param.m1 + _param.m2 + _param.m3);
            double c = (_param.m3 * _param.l3) / (_param.m1 + _param.m2 + _param.m3);
            double thetaMin = a * q1 + b * q2 + c * q3;
            double thetaMax = a * (-q1) + b * (q2) + c * (-q3);
            return new double[] { thetaMin, thetaMax };
        }
        public static double evaluateAlphaConstraint(double thetaMin, double thetaMax, BRgait _gait)
        {


            double dx = 0.02;
            double theta = thetaMin;
            double val = double.MinValue;
            double alphaVal = 0;
            for (int i = 0; i < thetaMax / dx; i++)
            {
                alphaVal = _gait.vhc.evalAlpha(theta);
                if (alphaVal > val)
                {
                    val = alphaVal;
                }
                theta += dx;
            }
            return val;
        }
        public static double evaluateDthetaConstraintCoM(double thetaMin, double thetaMax, BRgait _gait)
        {
            double[][] firstIntegral = TrapezoidalSum2.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha, thetaMin, thetaMax, 200);
            double[][] secondIntegral = TrapezoidalSum2.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1], thetaMin, thetaMax, 200);
            int len = firstIntegral[1].Length;
            double dthetaTSquared = (-secondIntegral[1][len - 1] * Math.Exp(firstIntegral[1][len - 1])) / (1 - Math.Exp(firstIntegral[1][len - 1]) * Math.Pow(_gait.impactSecondLine(thetaMin, thetaMax), 2));
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(thetaMin, thetaMax), 2);
            double dthetaMin = 0;
            if (double.IsNaN(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsNaN(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else
            {
                double[,] THETA = new double[2, len];
                THETA[0, 0] = thetaMin;
                THETA[1, 0] = dtheta0Squared;
                dthetaMin = dtheta0Squared;
                for (int i = 1; i < len; i++)
                {
                    THETA[0, i] = secondIntegral[0][i];
                    THETA[1, i] = -secondIntegral[1][i] * Math.Exp(firstIntegral[1][i]) + Math.Exp(firstIntegral[1][i]) * dtheta0Squared;
                    if (THETA[1, i] < dthetaMin)
                    {
                        dthetaMin = THETA[1, i];
                    }
                }
            }
            return dthetaMin;
        }


        public static void donothing()
        {
            int a = 1;
            //int len = THETA.Length / (THETA.Rank * 8);
            //using (StreamWriter file =
            //            new System.IO.StreamWriter(@"../../../../VALS.txt", true))
            //{
            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(gait.vhc.evalAlpha(THETA[0, i]));
            //    }
            //    file.WriteLine("");

            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(gait.vhc.evalBeta(THETA[0, i]));
            //    }
            //    file.WriteLine("");
            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(gait.vhc.evalGamma(THETA[0, i]));
            //    }
            //    file.WriteLine("");

            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(gait.vhc.evalAlpha1(THETA[0, i]));
            //    }
            //    file.WriteLine("");

            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(gait.vhc.evalBeta1(THETA[0, i]));
            //    }
            //    file.WriteLine("");
            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(gait.vhc.evalGamma1(THETA[0, i]));
            //    }
            //    file.WriteLine("");

            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(gait.vhc.evalAlpha3(THETA[0, i]));
            //    }
            //    file.WriteLine("");

            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(gait.vhc.evalBeta3(THETA[0, i]));
            //    }
            //    file.WriteLine("");
            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(gait.vhc.evalGamma3(THETA[0, i]));
            //    }
            //    file.WriteLine("");

            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(torques.torque1[i]);
            //    }
            //    file.WriteLine("");

            //    for (int i = 0; i < len; i += 20)
            //    {
            //        file.WriteLine(torques.torque2[i]);
            //    }
        }

        static void findq1Gait(Biped biped)
        {
            int numberOfPoints = 6;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);


            string[] parameters = File.ReadAllLines(@"../../../../parameters.txt");
            //double dtheta0 = Convert.ToDouble(File.ReadAllLines(@"../../../../dtheta0.txt")[0], CultureInfo.InvariantCulture);
            //double dthetaT = Convert.ToDouble(File.ReadAllLines(@"../../../../dthetaT.txt")[0], CultureInfo.InvariantCulture);
            for (int i = 0; i < parameters.Length; i++)
            {
                Console.WriteLine(parameters[i]);
            }
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);

            tuple.Item1["P0"] = Convert.ToDouble(parameters[0], CultureInfo.InvariantCulture);
            tuple.Item1["P1"] = -2*Convert.ToDouble(parameters[0], CultureInfo.InvariantCulture);
            for (int i = 2; i < 6; i++)
            {
                tuple.Item1["P" + i.ToString()] = 0;

            }

            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 6; i < 12; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi2Parameters = tuple.Item1;

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 12; i < 18; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi3Parameters = tuple.Item1;

            gait.vhc.phi2Parameters.Add("q1end", -Convert.ToDouble(parameters[0], CultureInfo.InvariantCulture));
            gait.vhc.phi3Parameters.Add("q1end", -Convert.ToDouble(parameters[0], CultureInfo.InvariantCulture));

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined("(P7 + 2 * P8 * theta + 3 * P9 * theta ^ 2 + 4 * P10 * theta ^ 3 + 5 * P11 * theta ^ 4) / (2 * q1end)");
            gait.vhc.dphi3 = Infix.ParseOrUndefined("(P13 + 2 * P14 * theta + 3 * P15 * theta ^ 2 + 4 * P16 * theta ^ 3 + 5 * P17 * theta ^ 4) / (2 * q1end)");

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined("(2 * P8 + 6 * P9 * theta + 12 * P10 * theta ^ 2 + 20 * P11 * theta ^ 3) / (2 * q1end)^2");
            gait.vhc.ddphi3 = Infix.ParseOrUndefined("(2 * P14 + 6 * P15 * theta + 12 * P16 * theta ^ 2 + 20 * P17 * theta ^ 3) / (2 * q1end)^2");

            Console.WriteLine(gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactSecondLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactThirdLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));


            double dtheta0Squared = 0;// Math.Pow(dtheta0, 2);
            double dthetaTSquared = 0;// Math.Pow(dthetaT, 2);


            //SLSQPFifthOrderTorqueEquality slSQPTE = new SLSQPFifthOrderTorqueEquality(gait, numberOfPoints);
            //slSQPTE.runNumericalSLSQP();

            //SLSQPFifthOrderTorque slSQPT = new SLSQPFifthOrderTorque(gait, numberOfPoints);
            //slSQPT.runNumericalSLSQP();

            //SLSQPFifthOrder slSQP = new SLSQPFifthOrder(gait, numberOfPoints);
            //slSQP.runNumericalSLSQP();

            //SLSQPFifthOrderTorqueDtheta slSQPTD = new SLSQPFifthOrderTorqueDtheta(gait, numberOfPoints);
            //slSQPTD.runNumericalSLSQP();

            //SLSQPFifthOrderDtheta slSQPTD = new SLSQPFifthOrderDtheta(gait, numberOfPoints);
            //slSQPTD.runNumericalSLSQP();

            SLSQPFifthOrderq1 slSQPq1 = new SLSQPFifthOrderq1(gait, numberOfPoints, gait.vhc.phi1Parameters["P0"].RealValue);
            //slSQPq1.runNumericalSLSQP();
            double[,] THETA = evalDthetaConstraint(gait, ref dtheta0Squared, ref dthetaTSquared);

            Phaseportrait plot = new Phaseportrait(THETA);

            BRReducedSimulationData data = integrationReducedDynamics.run(gait.vhc, Vector<double>.Build.Dense(new double[] { 0, Math.Sqrt(dtheta0Squared) }),
                Vector<double>.Build.Dense(new double[] { 1, Math.Sqrt(dthetaTSquared) }));
            biped.reducedSimulationData = data;

            graph graph = new graph(biped);


            using (StreamWriter file =
                                        new System.IO.StreamWriter(@"../../../../foundgaits.txt", true))
            {
                file.WriteLine(gait.vhc.phi1Parameters["P0"].RealValue);
                file.WriteLine(gait.vhc.phi1Parameters["P1"].RealValue);
                file.WriteLine(gait.vhc.phi1Parameters["P2"].RealValue);
                file.WriteLine(gait.vhc.phi1Parameters["P3"].RealValue);
                file.WriteLine(gait.vhc.phi1Parameters["P4"].RealValue);
                file.WriteLine(gait.vhc.phi1Parameters["P5"].RealValue);

                file.WriteLine(gait.vhc.phi2Parameters["P6"].RealValue);
                file.WriteLine(gait.vhc.phi2Parameters["P7"].RealValue);
                file.WriteLine(gait.vhc.phi2Parameters["P8"].RealValue);
                file.WriteLine(gait.vhc.phi2Parameters["P9"].RealValue);
                file.WriteLine(gait.vhc.phi2Parameters["P10"].RealValue);
                file.WriteLine(gait.vhc.phi2Parameters["P11"].RealValue);

                file.WriteLine(gait.vhc.phi3Parameters["P12"].RealValue);
                file.WriteLine(gait.vhc.phi3Parameters["P13"].RealValue);
                file.WriteLine(gait.vhc.phi3Parameters["P14"].RealValue);
                file.WriteLine(gait.vhc.phi3Parameters["P15"].RealValue);
                file.WriteLine(gait.vhc.phi3Parameters["P16"].RealValue);
                file.WriteLine(gait.vhc.phi3Parameters["P17"].RealValue);
                file.WriteLine("NEW GAIT");
                file.WriteLine("");
            }

            BRTorques torques = calculateTorques2.run(gait.vhc, THETA);
            TorquesGraph torquesGraph = new TorquesGraph(torques);
        }
        static void testq1(Biped biped)
        {
            int numberOfPoints = 8;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);


            string[] parameters = File.ReadAllLines(@"../../../../parameters.txt");
            //double dtheta0 = Convert.ToDouble(File.ReadAllLines(@"../../../../dtheta0.txt")[0], CultureInfo.InvariantCulture);
            //double dthetaT = Convert.ToDouble(File.ReadAllLines(@"../../../../dthetaT.txt")[0], CultureInfo.InvariantCulture);
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);

            for (int i = 0; i < 8; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 8; i < 16; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi2Parameters = tuple.Item1;

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 15; i < 24; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi3Parameters = tuple.Item1;

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());

            Infix.Format(gait.vhc.phi1);
            Infix.Format(gait.vhc.phi2);
            Infix.Format(gait.vhc.phi3);
            double dtheta0Squared = 0;// Math.Pow(dtheta0, 2);
            double dthetaTSquared = 0;// Math.Pow(dthetaT, 2);

            double[,] THETA = evalDthetaConstraintq1(gait, ref dtheta0Squared, ref dthetaTSquared);
            writeThetaToFile(THETA, @"../../../../../../../../MATLAB/theta.csv", gait.vhc);
            Phaseportrait plot = new Phaseportrait(THETA);
            double thetaStart = 0.083430086;
            double thetaEnd = -0.083430086;
            double dtheta0 = -0.1943;
            BRReducedSimulationData data = integrationReducedDynamics.run(gait.vhc, Vector<double>.Build.Dense(new double[] { thetaStart, dtheta0 }),
                Vector<double>.Build.Dense(new double[] { thetaEnd, Math.Sqrt(dthetaTSquared) }));
            biped.reducedSimulationData = data;
            writeThetaEvolutionToFile(data, @"../../../../../../../../MATLAB/thetaTime.csv");
            graph graph = new graph(biped);

            BRTorques torques = calculateTorques2.run(gait.vhc, THETA);
            TorquesGraph torquesGraph = new TorquesGraph(torques);


            TrajectoryGraph trajectory = new TrajectoryGraph(gait, THETA);
            
            Console.WriteLine(gait.impactFirstLine(thetaStart, thetaEnd));
            Console.WriteLine(gait.impactSecondLine(thetaStart, thetaEnd));
            Console.WriteLine(gait.impactThirdLine(thetaStart, thetaEnd));



        }
        public static double[,] evalDthetaConstraintq1(BRgait gait, ref double dtheta0Squared, ref double dthetaTSquared)
        {
            double thetaStart = 0.083430086;
            double thetaEnd = -0.083430086;
            double[][] firstIntegral = TrapezoidalSum2.calculateFirstIntegral(gait.vhc.evalTwoTimesBetaDividedByAlpha, thetaStart, thetaEnd, 200);
            double[][] secondIntegral = TrapezoidalSum2.calculateSecondIntegral(gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1], thetaStart, thetaEnd, 200);
            int len = firstIntegral[1].Length;
            dthetaTSquared = (-secondIntegral[1][len - 1] * Math.Exp(firstIntegral[1][len - 1])) / (1 - Math.Exp(firstIntegral[1][len - 1]) * Math.Pow(gait.impactSecondLine(thetaStart, thetaEnd), 2));
            dtheta0Squared = dthetaTSquared * Math.Pow(gait.impactSecondLine(thetaStart, thetaEnd), 2);

            int numberOfCycles = 1;
            double[,] THETA = new double[2, numberOfCycles * (len)];
            double temp = 0;
            //write firstintegral and secondintegral values to matlab
            writeRhoToFile(firstIntegral, secondIntegral, @"../../../../../../../../MATLAB/rho.csv", gait.vhc);
            for (int j = 0; j < numberOfCycles; j++)
            {
                THETA[0, j * len] = 0;
                THETA[1, j * len] = Math.Sqrt(dtheta0Squared);
                for (int i = 1; i < len; i++)
                {
                    THETA[0, i + j * len] = secondIntegral[0][i];
                    THETA[1, i + j * len] = Math.Sqrt(-secondIntegral[1][i] * Math.Exp(firstIntegral[1][i]) + Math.Exp(firstIntegral[1][i]) * dtheta0Squared);
                    temp = THETA[1, i];
                }
                dtheta0Squared = Math.Pow(gait.impactSecondLine(thetaStart, thetaEnd), 2) * Math.Pow(temp, 2);
            }
            using (StreamWriter file =
                        new System.IO.StreamWriter(@"../../../../VALS.txt", true))
            {
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(THETA[0, i]);
                }
                file.WriteLine("");
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(THETA[1, i]);
                }
                file.WriteLine("");
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(Math.Exp(firstIntegral[1][i]));
                }
                file.WriteLine("");
                for (int i = 0; i < len; i += 20)
                {
                    file.WriteLine(secondIntegral[1][i]);
                }
                file.WriteLine("");
            }
            return THETA;
        }



    }

    


}
