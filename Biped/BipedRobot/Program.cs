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
        static void Main(string[] args)
        {
            Biped biped = new Biped(@"../../../CAD_params_ILeg.xml");
            //IntegrationFullDynamics.run(ref biped);
            //plotting.plotStates(biped);
            findGaitsManually(biped);
            //test(biped);
            
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

            Phaseportrait plot = new Phaseportrait(THETA);

            BRReducedData data = integrationReducedDynamics.run(gait.vhc, Vector<double>.Build.Dense(new double[] { 0, Math.Sqrt(dtheta0Squared) }),
                Vector<double>.Build.Dense(new double[] { 1, Math.Sqrt(dthetaTSquared) }));
            biped.reducedData = data;

            graph graph = new graph(biped);

            BRTorques torques = calculateTorques2.run(gait.vhc, THETA);
            TorquesGraph torquesGraph = new TorquesGraph(torques);

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

            BRReducedData data = integrationReducedDynamics.run(gait.vhc, Vector<double>.Build.Dense(new double[] { 0, Math.Sqrt(dtheta0Squared) }),
                Vector<double>.Build.Dense(new double[] { 1, Math.Sqrt(dthetaTSquared) }));
            biped.reducedData = data;

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
            double[,] THETA = new double[2, 8*(len)];
            double temp = 0;
            for (int j = 0; j < 8; j++)
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
            return THETA;
        }


        
       }

    
}
