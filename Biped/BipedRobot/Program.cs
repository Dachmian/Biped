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
            Biped biped = new Biped (@"../../basic_randomized_set.xml");
            //IntegrationFullDynamics.run(ref biped);
            //plotting.plotStates(biped);
            //gaitSearch.run(ref biped);
            findGaitsManually(biped);
            //findGaitsHere(biped);
            
        }

        static void test()
        {
            Dictionary<string, FloatingPoint>  _parameters = new Dictionary<string, FloatingPoint>();
            _parameters.Add("g", 1);
            _parameters.Add("m1", 1);
            _parameters.Add("m2", 1);
            _parameters.Add("m3", 1);

            _parameters.Add("l1", 1);
            _parameters.Add("l2", 1);
            _parameters.Add("l3", 10);

            _parameters.Add("L1", 10);
            _parameters.Add("L2", 1);
            _parameters.Add("L3", 1);

            _parameters.Add("J1", 10);
            _parameters.Add("J2", 1);
            _parameters.Add("J3", 10);
            _parameters.Add("theta", 0);
            _parameters.Add("dtheta", 0);
            _parameters.Add("ddtheta", 0);
            _parameters.Add("P0", 1);
            _parameters.Add("P1", 1);
            _parameters.Add("P2", 0);
            _parameters.Add("P3", 0);
            _parameters.Add("P4", 0);
            _parameters.Add("P5", 0);
            _parameters.Add("P6", 0);
            _parameters.Add("P7", 0);
            _parameters.Add("P8", 0);
            _parameters.Add("P9", 0);
            _parameters.Add("P10", 0);
            _parameters.Add("P11", 0);
            _parameters.Add("P12", 0);
            _parameters.Add("P13", 0);
            _parameters.Add("P14", 1);

            Expression test = Infix.ParseOrUndefined("2*(-((J1 - L1*l1*m1 + l1^2*m1)*(0^4*L1*L3*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) - 0^4*L3*l1*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)))/((4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1) - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4))^2 + (1^4*L1*l2*m1*sin(1^4*P0 + 1^3*P1 + 1^2*P2 + 1*P3 + P4 - (1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)))/((4*0^3*P5 + 3*0^2*P6 + 2*0*P7 + P8)*(J2 + l2^2*m1) + L3*l2*m1*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9)))*((J1 - L1*l1*m1 + l1^2*m1)/((4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1) - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)) - ((4*1^3*P5 + 3*1^2*P6 + 2*1*P7 + P8)*(J2 + l2^2*m1) + L1*l2*m1*cos(1^4*P0 + 1^3*P1 + 1^2*P2 + 1*P3 + P4 - (1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)))/((4*0^3*P5 + 3*0^2*P6 + 2*0*P7 + P8)*(J2 + l2^2*m1) + L3*l2*m1*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9)))^(-1 + 2) - 2*((J1 - L1*l1*m1 + l1^2*m1)/((4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1) - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)) - (J1 - L1*l1*m1 + l1^2*m1 + (4*1^3*P10 + 3*1^2*P11 + 2*1*P12 + P13)*(J3 - L3*l3*m3 + l3^2*m3) + L3*l1*m1*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + L1*L3*m2*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + L1*L3*m3*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) - L1*l3*m3*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + (4*1^3*P5 + 3*1^2*P6 + 2*1*P7 + P8)*(J2 + l2^2*m2 + L3*l2*m2*cos(-(1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)) + L1*l2*m2*cos(1^4*P0 + 1^3*P1 + 1^2*P2 + 1*P3 + P4 - (1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)))/(J3 + L3^2*m1 + L3^2*m2 + L3^2*m3 - 2*L3*l3*m3 + l3^2*m3 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + (4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)) + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9) + (4*0^3*P5 + 3*0^2*P6 + 2*0*P7 + P8)*(J2 + l2^2*m2 + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9))))^(-1 + 2)*(-((J1 - L1*l1*m1 + l1^2*m1)*(0^4*L1*L3*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) - 0^4*L3*l1*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)))/((4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1) - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4))^2 - ((-1^4*L3*l1*m1*sin(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) - 1^4*L1*L3*m2*sin(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) - 1^4*L1*L3*m3*sin(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + 1^4*L1*l3*m3*sin(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) - 1^4*L1*l2*m2*sin(1^4*P0 + 1^3*P1 + 1^2*P2 + 1*P3 + P4 - (1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)))/(J3 + L3^2*m1 + L3^2*m2 + L3^2*m3 - 2*L3*l3*m3 + l3^2*m3 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + (4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)) + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9) + (4*0^3*P5 + 3*0^2*P6 + 2*0*P7 + P8)*(J2 + l2^2*m2 + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9))) - ((0^4*L1*L3*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) - 0^4*L3*l1*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + (4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(0^4*L1*L3*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) - 0^4*L3*l1*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)))*(J1 - L1*l1*m1 + l1^2*m1 + (4*1^3*P10 + 3*1^2*P11 + 2*1*P12 + P13)*(J3 - L3*l3*m3 + l3^2*m3) + L3*l1*m1*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + L1*L3*m2*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + L1*L3*m3*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) - L1*l3*m3*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + (4*1^3*P5 + 3*1^2*P6 + 2*1*P7 + P8)*(J2 + l2^2*m2 + L3*l2*m2*cos(-(1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)) + L1*l2*m2*cos(1^4*P0 + 1^3*P1 + 1^2*P2 + 1*P3 + P4 - (1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9))))/(J3 + L3^2*m1 + L3^2*m2 + L3^2*m3 - 2*L3*l3*m3 + l3^2*m3 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + (4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)) + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9) + (4*0^3*P5 + 3*0^2*P6 + 2*0*P7 + P8)*(J2 + l2^2*m2 + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9)))^2))");
            Console.WriteLine(MathNet.Symbolics.Evaluate.Evaluate(_parameters, test).RealValue);
            
        }

        static void findGaitsManually(Biped biped)
        {
            int numberOfPoints = 4;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);


            string[] parameters = File.ReadAllLines(@"../../../parameters.txt");
            double dtheta0 = Convert.ToDouble(File.ReadAllLines(@"../../../dtheta0.txt")[0], CultureInfo.InvariantCulture);
            double dthetaT = Convert.ToDouble(File.ReadAllLines(@"../../../dthetaT.txt")[0], CultureInfo.InvariantCulture);
            for (int i = 0; i < parameters.Length; i++)
            {
                Console.WriteLine(parameters[i]);
            }
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);

            for (int i = 0; i < 4; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 4; i < 8; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi2Parameters = tuple.Item1;

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 8; i < 12; i++)
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


            double dtheta0Squared = Math.Pow(dtheta0, 2);
            double dthetaTSquared = Math.Pow(dthetaT, 2);

            double[,] THETA = evalDthetaConstraint(gait, ref dtheta0Squared, ref dthetaTSquared);

            Phaseportrait plot = new Phaseportrait(THETA);

            BRReducedData data = integrationReducedDynamics.run(gait.vhc, Vector<double>.Build.Dense(new double[] { 0, Math.Sqrt(dtheta0Squared) }),
                Vector<double>.Build.Dense(new double[] { 1, Math.Sqrt(dthetaTSquared) }));
            biped.reducedData = data;

            graph graph = new graph(biped);
            //using (StreamWriter file =
            //                            new System.IO.StreamWriter(@"../../../foundgaits.txt", true))
            //{
            //    file.WriteLine(gait.vhc.phi1Parameters["P0"].RealValue);
            //    file.WriteLine(gait.vhc.phi1Parameters["P1"].RealValue);
            //    file.WriteLine(gait.vhc.phi1Parameters["P2"].RealValue);
            //    file.WriteLine(gait.vhc.phi1Parameters["P3"].RealValue);
            //    file.WriteLine(gait.vhc.phi2Parameters["P4"].RealValue);
            //    file.WriteLine(gait.vhc.phi2Parameters["P5"].RealValue);
            //    file.WriteLine(gait.vhc.phi2Parameters["P6"].RealValue);
            //    file.WriteLine(gait.vhc.phi2Parameters["P7"].RealValue);
            //    file.WriteLine(gait.vhc.phi3Parameters["P8"].RealValue);
            //    file.WriteLine(gait.vhc.phi3Parameters["P9"].RealValue);
            //    file.WriteLine(gait.vhc.phi3Parameters["P10"].RealValue);
            //    file.WriteLine(gait.vhc.phi3Parameters["P11"].RealValue);
            //    file.WriteLine("NEW GAIT");
            //    file.WriteLine("");
            //}
            BRTorques torques = calculateTorques.run(gait.vhc, data);
            TorquesGraph torquesGraph = new TorquesGraph(torques);

            AGS ags = new AGS(gait, numberOfPoints);
        }

        public static double[,] evalDthetaConstraint(BRgait gait, ref double dtheta0Squared, ref double dthetaTSquared)
        {
            double[][] firstIntegral = RiemannSum.calculateFirstIntegral(gait.firstIntegral);
            //double firstIntegralVALUE = GaussLegendreRule.Integrate(gait.firstIntegral, 0, 1, 1000);
            //double secondIntegralVALUE = GaussLegendreRule.Integrate(gait.secondIntegral, 0, 1, 100);
            double[][] secondIntegral = RiemannSum.calculateSecondIntegral(gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1]);            
            int len = firstIntegral[1].Length - 1;
            //dthetaTSquared = (- secondIntegral[1][len])/(1-Math.Exp(-firstIntegral[1][len])*Math.Pow(gait.impactSecondLine(0,1), 2) );
            //dtheta0Squared = dthetaTSquared * Math.Pow(gait.impactSecondLine(0, 1), 2);
            double[,] THETA = new double[2, len];
            THETA[0, 0] = 0;
            THETA[1, 0] = Math.Sqrt(dtheta0Squared);
            for (int i = 1; i < len; i++)
            {
                THETA[0, i] = secondIntegral[0][i];
                THETA[1, i] = Math.Sqrt(-secondIntegral[1][i] + Math.Exp(-firstIntegral[1][i]) * dtheta0Squared);
            }
            return THETA;
        }


        
       }

    
}
